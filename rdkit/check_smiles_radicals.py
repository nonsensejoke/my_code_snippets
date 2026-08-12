#!/usr/bin/env python3
"""
check_smiles_radicals.py

检测 SMILES 中「冗余方括号导致隐式氢丢失、产生非预期自由基」的问题。

背景
----
SMILES 有两套原子写法：
  * 有机子集原子 (B C N O P S F Cl Br I) 裸写时，按标准价态自动补充隐式氢；
  * 写进方括号后变成 bracket atom，氢数必须显式声明，**默认为 0**。

因此 `[N]` 不等于 `N`：前者是 0 个氢的氮，剩余价电子被 RDKit 记为
自由基电子 (radical electrons)。这类 SMILES 语法完全合法，
sanitization 不报错，会静默通过常规校验 —— 所以需要专门检测。

判定逻辑（结构层面，不依赖正则匹配字符串）
----------------------------------------
对每个原子，若同时满足：
  1. GetNoImplicit() 为 True        → 该原子来自方括号
  2. GetNumExplicitHs() == 0        → 作者没有写氢数
  3. GetNumRadicalElectrons() > 0   → 实际产生了未配对电子
  4. 元素属于有机子集              → 去掉方括号后语法仍合法
则判为 SUSPECT（高置信度书写错误）。

若该原子另外带有形式电荷或同位素标记，说明作者是有意在方括号里写东西，
降级为 AMBIGUOUS（需人工确认，可能是真实的自由基/离子物种）。

裸写原子产生的自由基（如 `C[CH2]` 写成 `CC[CH2]`… 或显式 `[CH3]`）
不会被判为 SUSPECT —— 显式写了氢数就是作者的主动选择。

用法
----
【只想知道会不会产生自由基 —— 用这个】

  python check_smiles_radicals.py library/ --strict

  --strict 的含义就是「任何自由基都算失败」，适用于生物活性分子/药物库
  这类本不应出现自由基的场景。开启后有三种失败状态，它们**都代表存在
  自由基**，区别只在成因（用于定位上游 bug 是哪一类）：

      SUSPECT    方括号原子没写氢数，也没有电荷/同位素   → 几乎必然是书写错误
      AMBIGUOUS  方括号原子带电荷或同位素标记            → 需人工确认
      RADICAL    氢数是显式写出来的（如 [CH2]CC）        → 另一类上游 bug

  另有 PARSE_ERROR / EMPTY 两种状态，属于「SMILES 根本读不出来」，
  跟自由基无关，但同样计入失败（读不出来就无法断言它没有自由基）。

  只要统计数字、不要逐条明细：
      python check_smiles_radicals.py library/ --strict --quiet

  留档以便逐条排查：
      python check_smiles_radicals.py library/ --strict -o report.csv

  ⚠️ 不加 --strict 时，显式声明氢数的自由基只会记为 OTHER_RADICAL 而
  **不计入失败、退出码仍为 0**。所以纯粹的「查自由基」用途必须加 --strict。

【在自己的代码里调用】

  from check_smiles_radicals import has_radical
  has_radical('[N][N][C]1C(=O)C=CC1=O')   # -> True
  has_radical('NNC1C(=O)C=CC1=O')         # -> False
  has_radical('C1CC')                     # -> None（解析失败，需自行处理）

【其他】

  python check_smiles_radicals.py mols/*.sdf                # 通配符
  python check_smiles_radicals.py mols/ --crosscheck        # 额外核对 mol block
  python check_smiles_radicals.py mols/ --strict --all      # 连通过的也逐条打印
  echo 'CCO' | python check_smiles_radicals.py - --strict   # 从 stdin 读 SMILES

退出码：0 = 全部通过；1 = 发现问题（方便接进 CI）；2 = 用法错误。

已知边界
--------
一氧化氮 (`[N]=O`)、氮氧自由基 (TEMPO 类探针/造影剂) 是**真实存在的
自由基类药物**，会被 --strict 标记。数量极少，人工过一遍报告的
bad_atoms 列即可，不建议为此放宽规则。
"""

import argparse
import csv
import glob
import os
import sys

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors
except ImportError:
    sys.exit("需要 RDKit：pip install rdkit")

RDLogger.DisableLog("rdApp.*")  # 抑制 RDKit 自带的 stderr 噪音

# 可以裸写（自动补氢）的有机子集元素
ORGANIC_SUBSET = {"B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"}

STATUS_OK = "OK"
STATUS_SUSPECT = "SUSPECT"
STATUS_AMBIGUOUS = "AMBIGUOUS"
STATUS_OTHER_RADICAL = "OTHER_RADICAL"
STATUS_PARSE_ERROR = "PARSE_ERROR"
STATUS_EMPTY = "EMPTY"
STATUS_RADICAL = "RADICAL"             # 仅 --strict 时可能出现
STATUS_MOLBLOCK_BAD = "MOLBLOCK_BAD"   # 仅 --crosscheck 时可能出现
STATUS_MISMATCH = "SMILES_MOLBLOCK_MISMATCH"

PROBLEM_STATUSES = {
    STATUS_SUSPECT, STATUS_AMBIGUOUS, STATUS_PARSE_ERROR,
    STATUS_RADICAL, STATUS_MOLBLOCK_BAD, STATUS_MISMATCH,
}


# --------------------------------------------------------------------------
# 最小接口：只判断「会不会产生自由基」
# --------------------------------------------------------------------------
def has_radical(smi):
    """
    单一职责的布尔判定，供其他脚本 import 使用。

    返回值：
      True  -> 该 SMILES 会产生自由基电子（不区分成因）
      False -> 不会
      None  -> SMILES 空或无法解析，无从判断（调用方必须自行处理这种情况，
               不要把 None 当成 False，否则坏数据会被静默放过）

    等价于命令行的 --strict 模式，只是不做成因分类、不给修正建议。
    需要知道是哪个原子出问题、或需要修正后的 SMILES 时，改用 analyse_smiles()。
    """
    if not smi or not smi.strip():
        return None
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return Descriptors.NumRadicalElectrons(mol) > 0


# --------------------------------------------------------------------------
# 核心检测
# --------------------------------------------------------------------------
def analyse_smiles(smi, strict=False):
    """
    分析单个 SMILES，返回结果 dict。

    strict=True 适用于「库内不应存在任何自由基」的场景（如生物活性分子/药物库）：
    此时显式声明了氢数的自由基（如 [CH2]CC）也判为失败，并一并尝试修正。
    原子级别的成因分类仍然保留在 bad_atoms 里，用于定位上游 bug 的类型。
    """
    result = {
        "smiles": smi,
        "status": STATUS_OK,
        "n_radical_electrons": 0,
        "bad_atoms": "",
        "formula": "",
        "fixed_smiles": "",
        "note": "",
    }

    if not smi or not smi.strip():
        result["status"] = STATUS_EMPTY
        result["note"] = "空的 SMILES 字段"
        return result

    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        result["status"] = STATUS_PARSE_ERROR
        result["note"] = "RDKit 无法解析"
        return result

    from rdkit.Chem import rdMolDescriptors

    result["formula"] = rdMolDescriptors.CalcMolFormula(mol)
    result["n_radical_electrons"] = Descriptors.NumRadicalElectrons(mol)

    suspects, ambiguous, others = [], [], []
    for atom in mol.GetAtoms():
        n_rad = atom.GetNumRadicalElectrons()
        if n_rad == 0:
            continue

        desc = (
            f"idx{atom.GetIdx()}:{atom.GetSymbol()}"
            f"(rad={n_rad},H={atom.GetTotalNumHs()},deg={atom.GetDegree()})"
        )

        from_bracket = atom.GetNoImplicit()
        no_h_written = atom.GetNumExplicitHs() == 0
        in_subset = atom.GetSymbol() in ORGANIC_SUBSET
        annotated = atom.GetFormalCharge() != 0 or atom.GetIsotope() != 0

        if from_bracket and no_h_written and in_subset:
            (ambiguous if annotated else suspects).append(desc)
        else:
            others.append(desc)

    if suspects:
        result["status"] = STATUS_SUSPECT
        result["bad_atoms"] = " ".join(suspects + ambiguous + others)
        result["note"] = (
            f"{len(suspects)} 个方括号原子未声明氢数，隐式氢被抑制为 0"
        )
        result["fixed_smiles"] = repair(mol, all_radicals=strict)
    elif ambiguous:
        result["status"] = STATUS_AMBIGUOUS
        result["bad_atoms"] = " ".join(ambiguous + others)
        result["note"] = "方括号原子带电荷/同位素标记，可能是有意的自由基物种，请人工确认"
        result["fixed_smiles"] = repair(mol, all_radicals=strict)
    elif others:
        result["bad_atoms"] = " ".join(others)
        if strict:
            result["status"] = STATUS_RADICAL
            result["note"] = "存在自由基（氢数为显式声明）；strict 模式下一律判为失败"
            result["fixed_smiles"] = repair(mol, all_radicals=True)
        else:
            result["status"] = STATUS_OTHER_RADICAL
            result["note"] = "存在自由基，但氢数是显式声明的，推测为有意书写"

    return result


def repair(mol, all_radicals=False):
    """
    把自由基原子恢复为「按标准价态自动补氢」，返回修正后的 canonical SMILES。

    默认只处理「未声明氢数的方括号原子」。all_radicals=True 时处理所有自由基原子
    （包括 [CH2] 这类已声明氢数的），此时 RDKit 会在已有显式氢之上补足隐式氢。

    在分子对象上操作而不是做字符串替换，因此重原子骨架和键序保证不变。
    修正失败时返回空字符串（不猜）。

    注意：补足到标准价态只是「一种」可能的原意。上游也可能本该写成
    电荷分离形式或多一个双键。所以修正结果应当复核，最可靠的参照物是
    SDF 内部的 mol block（见 --crosscheck）。
    """
    rw = Chem.RWMol(mol)
    changed = False
    for atom in rw.GetAtoms():
        if atom.GetNumRadicalElectrons() == 0:
            continue
        if atom.GetSymbol() not in ORGANIC_SUBSET:
            continue  # 金属等留给人工处理
        if not all_radicals and not (atom.GetNoImplicit() and atom.GetNumExplicitHs() == 0):
            continue
        atom.SetNumRadicalElectrons(0)
        atom.SetNoImplicit(False)
        changed = True
    if not changed:
        return ""
    try:
        fixed = rw.GetMol()
        Chem.SanitizeMol(fixed)
    except Exception:
        return ""
    # 双重保险：骨架必须没变，且自由基必须清零
    if fixed.GetNumAtoms() != mol.GetNumAtoms() or fixed.GetNumBonds() != mol.GetNumBonds():
        return ""
    if Descriptors.NumRadicalElectrons(fixed) > 0:
        return ""
    return Chem.MolToSmiles(fixed)


# --------------------------------------------------------------------------
# 文件读取
# --------------------------------------------------------------------------
def read_titles(path):
    """
    读取 SDF，按 $$$$ 分记录，取每条记录的第一行（title line）。
    返回 [(record_index, title_line), ...]。单分子文件即只有一条。
    """
    with open(path, "r", errors="replace") as fh:
        text = fh.read()

    records = []
    for chunk in text.split("$$$$"):
        lines = chunk.lstrip("\n").splitlines()
        if not lines:
            continue
        if not any(line.strip() for line in lines):
            continue
        records.append(lines[0].strip())
    return list(enumerate(records, start=1))


def crosscheck_molblock(path, record_index):
    """
    可选：把 title 行的 SMILES 与 SDF 内部 mol block 的实际结构做比对。
    这能抓出所有类型的不一致，而不只是本脚本针对的氢数问题。

    返回 (smiles, error)：
      * 正常     -> (canonical_smiles, "")
      * 价键非法 -> (不做 sanitize 得到的 smiles 或 "", 错误说明)
    mol block 本身无法 sanitize 也是必须报出来的问题，不能静默跳过。
    """
    supplier = Chem.SDMolSupplier(path, sanitize=True, removeHs=False)
    mols = list(supplier)
    if record_index > len(mols):
        return "", "mol block 记录数不足"

    mol = mols[record_index - 1]
    if mol is not None:
        return Chem.MolToSmiles(mol), ""

    # sanitize 失败 —— 退回不校验读取，至少把连接表捞出来给人看
    raw = Chem.SDMolSupplier(path, sanitize=False, removeHs=False)
    raw_mols = list(raw)
    raw_mol = raw_mols[record_index - 1] if record_index <= len(raw_mols) else None
    if raw_mol is None:
        return "", "mol block 无法读取"

    problems = Chem.DetectChemistryProblems(raw_mol)
    detail = "; ".join(p.Message() for p in problems) or "未知 sanitization 失败"
    try:
        rough = Chem.MolToSmiles(raw_mol)
    except Exception:
        rough = ""
    return rough, f"mol block 本身无法 sanitize：{detail}"


def expand_inputs(args_paths):
    """把目录 / 通配符 / 文件名混合的输入展开成文件列表。"""
    files = []
    for p in args_paths:
        if p == "-":
            files.append("-")
        elif os.path.isdir(p):
            files.extend(sorted(glob.glob(os.path.join(p, "**", "*.sdf"), recursive=True)))
        elif any(ch in p for ch in "*?["):
            files.extend(sorted(glob.glob(p, recursive=True)))
        else:
            files.append(p)
    return files


# --------------------------------------------------------------------------
# 主流程
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="检测 SMILES 中冗余方括号导致的非预期自由基（SMILES 取自 SDF 第一行）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("paths", nargs="+", help="SDF 文件、目录、通配符，或 - 表示从 stdin 逐行读 SMILES")
    ap.add_argument("-o", "--output", help="把完整报告写成 CSV")
    ap.add_argument("--crosscheck", action="store_true",
                    help="额外把 title SMILES 与 SDF 内部 mol block 结构比对")
    ap.add_argument("--strict", action="store_true",
                    help="【查自由基就用这个】任何自由基一律判为失败，\n"
                         "包括 [CH2] 这类显式声明了氢数的。生物活性分子/药物库建议常开")
    ap.add_argument("--quiet", action="store_true", help="只打印统计摘要")
    ap.add_argument("--all", action="store_true", help="逐条打印（含通过的分子）")
    args = ap.parse_args()

    files = expand_inputs(args.paths)
    if not files:
        sys.exit(2)

    rows = []
    for path in files:
        if path == "-":
            for i, line in enumerate(sys.stdin, start=1):
                r = analyse_smiles(line.strip(), strict=args.strict)
                r.update(file="<stdin>", record=i)
                rows.append(r)
            continue

        if not os.path.isfile(path):
            rows.append({
                "file": path, "record": 0, "smiles": "", "status": STATUS_PARSE_ERROR,
                "n_radical_electrons": 0, "bad_atoms": "", "formula": "",
                "fixed_smiles": "", "note": "文件不存在",
            })
            continue

        for idx, title in read_titles(path):
            r = analyse_smiles(title, strict=args.strict)
            r.update(file=path, record=idx)
            if args.crosscheck:
                mb, mb_err = crosscheck_molblock(path, idx)
                parsed = Chem.MolFromSmiles(r["smiles"]) if r["smiles"] else None
                ref = r["fixed_smiles"] or (Chem.MolToSmiles(parsed) if parsed else "")

                if mb_err:
                    msg = f"{mb_err}" + (f" (连接表约为 {mb})" if mb else "")
                    r["status"] = STATUS_MOLBLOCK_BAD if r["status"] == STATUS_OK else r["status"]
                elif mb and ref and mb != ref:
                    msg = f"title SMILES 与 mol block 不一致 (molblock={mb})"
                    r["status"] = STATUS_MISMATCH if r["status"] == STATUS_OK else r["status"]
                elif mb and ref:
                    msg = "与 mol block 一致"
                else:
                    msg = "无法比对"
                r["note"] = (r["note"] + " | " if r["note"] else "") + msg
            rows.append(r)

    # ---- 输出 ----
    counts = {}
    for r in rows:
        counts[r["status"]] = counts.get(r["status"], 0) + 1

    if not args.quiet:
        shown = rows if args.all else [r for r in rows if r["status"] != STATUS_OK]
        for r in shown:
            mark = "OK " if r["status"] == STATUS_OK else "!! "
            print(f"{mark}[{r['status']}] {r['file']}#{r['record']}")
            print(f"    smiles : {r['smiles']}")
            if r["formula"]:
                print(f"    formula: {r['formula']}  radicals={r['n_radical_electrons']}")
            if r["bad_atoms"]:
                print(f"    atoms  : {r['bad_atoms']}")
            if r["note"]:
                print(f"    note   : {r['note']}")
            if r["fixed_smiles"]:
                print(f"    fix -> : {r['fixed_smiles']}")
            print()

    total = len(rows)
    n_bad = sum(counts.get(s, 0) for s in PROBLEM_STATUSES)
    print(f"共 {total} 条 | " + " | ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    print(f"需要处理：{n_bad} 条")

    if args.output:
        fields = ["file", "record", "status", "smiles", "fixed_smiles",
                  "formula", "n_radical_electrons", "bad_atoms", "note"]
        with open(args.output, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        print(f"报告已写入 {args.output}")

    sys.exit(1 if n_bad else 0)


if __name__ == "__main__":
    main()
