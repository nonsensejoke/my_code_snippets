#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
topo_check.py — 检查 SDF 中 title 行的 SMILES 与 connection table 的重原子拓扑是否一致

问题本质：这是一个「带节点标签(元素)、边无标签」的图同构判定问题。
必须忽略的东西：bond order / 芳香性 / 形式电荷 / 氢原子 / 原子编号顺序 / 3D 坐标。

流程（逐级加严，越早失败越便宜）：
  Stage 0  解析：SMILES 与 molblock 都用 sanitize=False 读，避免因价键异常直接报错
  Stage 1  廉价不变量：分子式、原子数、键数、连通分支数、环数(E-V+C)、逐元素度分布
  Stage 2  WL(Weisfeiler-Lehman)哈希：必要非充分条件，但能挡掉几乎所有不一致
  Stage 3  精确同构判定：VF2(networkx) 或内置回溯算法，给出最终判决
  Stage 4  不一致时的定位诊断：原子环境(元素+邻居元素多重集)的差集

用法:
    python topo_check.py mol_127311.sdf mol_347022.sdf
    python topo_check.py /path/to/sdf_dir --csv result.csv
    python topo_check.py *.sdf --quiet          # 只打印不一致的
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from collections import Counter, deque
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

try:
    import networkx as nx
except ImportError:  # 可选依赖
    nx = None


# ────────────────────────────────────────────────────────────────────────────
# 图结构
# ────────────────────────────────────────────────────────────────────────────
@dataclass
class Graph:
    """重原子图：节点带元素标签，边无标签(平行边已折叠)。"""
    symbols: list
    adj: list  # list[set[int]]

    @property
    def n(self) -> int:
        return len(self.symbols)

    @property
    def n_bonds(self) -> int:
        return sum(len(s) for s in self.adj) // 2

    def degrees(self) -> list:
        return [len(s) for s in self.adj]

    def n_components(self) -> int:
        seen, comp = set(), 0
        for s in range(self.n):
            if s in seen:
                continue
            comp += 1
            q = deque([s])
            seen.add(s)
            while q:
                u = q.popleft()
                for v in self.adj[u]:
                    if v not in seen:
                        seen.add(v)
                        q.append(v)
        return comp

    def n_rings(self) -> int:
        """SSSR 环数 = E - V + C（圈秩），比调用 RDKit 环感知更鲁棒。"""
        return self.n_bonds - self.n + self.n_components()

    def formula(self) -> str:
        c = Counter(self.symbols)
        head = [f"{e}{c.pop(e)}" for e in ("C", "N", "O") if e in c]
        return "".join(head + [f"{e}{c[e]}" for e in sorted(c)])

    def env_counter(self) -> Counter:
        """原子环境签名：元素(排序后的邻居元素)，用于定位差异。"""
        return Counter(
            f"{self.symbols[i]}({','.join(sorted(self.symbols[j] for j in self.adj[i]))})"
            for i in range(self.n)
        )

    def bfs_order(self) -> list:
        order, seen = [], set()
        for s in range(self.n):
            if s in seen:
                continue
            seen.add(s)
            q = deque([s])
            while q:
                u = q.popleft()
                order.append(u)
                for v in sorted(self.adj[u]):
                    if v not in seen:
                        seen.add(v)
                        q.append(v)
        return order


def mol_to_graph(mol) -> Graph:
    """RDKit Mol -> 重原子图。丢弃显式 H、折叠重复键、忽略键级/电荷/立体。"""
    keep, idx_map, symbols = [], {}, []
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 1 and a.GetDegree() > 0:  # 显式氢，丢弃
            continue
        idx_map[a.GetIdx()] = len(symbols)
        symbols.append(a.GetSymbol() if a.GetAtomicNum() else "*")
        keep.append(a.GetIdx())
    adj = [set() for _ in symbols]
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i not in idx_map or j not in idx_map:
            continue
        u, v = idx_map[i], idx_map[j]
        if u == v:  # 自环，忽略
            continue
        adj[u].add(v)
        adj[v].add(u)
    return Graph(symbols, adj)


# ────────────────────────────────────────────────────────────────────────────
# WL 哈希
# ────────────────────────────────────────────────────────────────────────────
def _h(s: str) -> str:
    return hashlib.blake2b(s.encode(), digest_size=8).hexdigest()


def wl_labels(g: Graph, iterations: int = 4) -> list:
    labels = [f"{s}:{len(g.adj[i])}" for i, s in enumerate(g.symbols)]
    for _ in range(iterations):
        labels = [
            _h(labels[i] + "|" + ",".join(sorted(labels[j] for j in g.adj[i])))
            for i in range(g.n)
        ]
    return labels


def wl_hash(g: Graph, iterations: int = 4) -> str:
    """整图哈希：累计每一轮的标签直方图（比只取最后一轮更强）。"""
    labels = [f"{s}:{len(g.adj[i])}" for i, s in enumerate(g.symbols)]
    acc = ["|".join(sorted(labels))]
    for _ in range(iterations):
        labels = [
            _h(labels[i] + "|" + ",".join(sorted(labels[j] for j in g.adj[i])))
            for i in range(g.n)
        ]
        acc.append("|".join(sorted(labels)))
    return _h("#".join(acc))


# ────────────────────────────────────────────────────────────────────────────
# 精确同构判定
# ────────────────────────────────────────────────────────────────────────────
def isomorphic(g1: Graph, g2: Graph) -> bool:
    if g1.n != g2.n or g1.n_bonds != g2.n_bonds:
        return False
    if Counter(g1.symbols) != Counter(g2.symbols):
        return False
    if nx is not None:
        return nx.is_isomorphic(
            _to_nx(g1), _to_nx(g2), node_match=lambda a, b: a["el"] == b["el"]
        )
    return _backtrack_iso(g1, g2)


def _to_nx(g: Graph):
    G = nx.Graph()
    for i, s in enumerate(g.symbols):
        G.add_node(i, el=s)
    for i in range(g.n):
        for j in g.adj[i]:
            if i < j:
                G.add_edge(i, j)
    return G


def _backtrack_iso(g1: Graph, g2: Graph) -> bool:
    """无 networkx 时的兜底：WL 标签剪枝 + 回溯，检查边与非边都保持。"""
    l1, l2 = wl_labels(g1), wl_labels(g2)
    if Counter(l1) != Counter(l2):
        return False
    cand = {i: [j for j in range(g2.n) if l2[j] == l1[i]] for i in range(g1.n)}
    order = g1.bfs_order()
    mapping, used = {}, set()

    def bt(pos: int) -> bool:
        if pos == g1.n:
            return True
        i = order[pos]
        for j in cand[i]:
            if j in used:
                continue
            if all((mapping[k] in g2.adj[j]) == (k in g1.adj[i]) for k in mapping):
                mapping[i], _ = j, used.add(j)
                if bt(pos + 1):
                    return True
                del mapping[i]
                used.discard(j)
        return False

    return bt(0)


# ────────────────────────────────────────────────────────────────────────────
# SDF 解析
# ────────────────────────────────────────────────────────────────────────────
def iter_records(path: Path):
    """产出 (record_index, title, molblock)。支持多分子 SDF。"""
    text = path.read_text(errors="replace")
    for k, chunk in enumerate(text.split("$$$$")):
        lines = chunk.splitlines()
        while lines and not lines[0].strip():
            lines.pop(0)
        if not lines:
            continue
        title = lines[0].strip()
        end = next((i for i, l in enumerate(lines) if l.strip() == "M  END"), None)
        block = "\n".join(lines[: end + 1] if end is not None else lines)
        yield k, title, block


def parse_smiles(smi: str):
    m = Chem.MolFromSmiles(smi)          # 先按标准解析
    if m is None:
        m = Chem.MolFromSmiles(smi, sanitize=False)  # 退化：只要连接性
    return m


def parse_block(block: str):
    return Chem.MolFromMolBlock(
        block, sanitize=False, removeHs=False, strictParsing=False
    )


# ────────────────────────────────────────────────────────────────────────────
# 单条记录检查
# ────────────────────────────────────────────────────────────────────────────
def check_record(title: str, block: str) -> dict:
    r = {
        "smiles": title,
        "status": "",
        "consistent": None,
        "reasons": [],
        "info": {},
    }
    if not title:
        r["status"], r["reasons"] = "NO_SMILES", ["title 行为空"]
        return r

    ms, mb = parse_smiles(title), parse_block(block)
    if ms is None:
        r["status"], r["reasons"] = "PARSE_ERROR", ["SMILES 无法解析"]
        return r
    if mb is None:
        r["status"], r["reasons"] = "PARSE_ERROR", ["molblock 无法解析"]
        return r

    g1, g2 = mol_to_graph(ms), mol_to_graph(mb)
    r["info"] = {
        "n_atoms": f"{g1.n}/{g2.n}",
        "n_bonds": f"{g1.n_bonds}/{g2.n_bonds}",
        "n_rings": f"{g1.n_rings()}/{g2.n_rings()}",
        "n_frags": f"{g1.n_components()}/{g2.n_components()}",
        "formula": f"{g1.formula()}/{g2.formula()}",
    }

    # Stage 1: 廉价不变量（必要条件）
    reasons = []
    if g1.n != g2.n:
        reasons.append(f"重原子数不同: SMILES {g1.n} vs SDF {g2.n}")
    if Counter(g1.symbols) != Counter(g2.symbols):
        reasons.append(f"分子式不同: {g1.formula()} vs {g2.formula()}")
    if g1.n_bonds != g2.n_bonds:
        reasons.append(f"键数不同: SMILES {g1.n_bonds} vs SDF {g2.n_bonds}")
    if g1.n_components() != g2.n_components():
        reasons.append(f"连通分支数不同: {g1.n_components()} vs {g2.n_components()}")
    if g1.n_rings() != g2.n_rings():
        reasons.append(f"环数(E-V+C)不同: SMILES {g1.n_rings()} vs SDF {g2.n_rings()}")

    if not reasons:
        # Stage 2: WL 哈希
        if wl_hash(g1) != wl_hash(g2):
            reasons.append("WL 哈希不同：原子数/键数相同但连接方式不同")
        # Stage 3: 精确同构（最终判决）
        elif not isomorphic(g1, g2):
            reasons.append("精确同构判定失败（WL 哈希碰撞）")

    if reasons:
        # Stage 4: 定位诊断
        d1, d2 = g1.env_counter(), g2.env_counter()
        only1, only2 = d1 - d2, d2 - d1
        if only1:
            reasons.append("仅出现在 SMILES 的原子环境: " + _fmt(only1))
        if only2:
            reasons.append("仅出现在 SDF 的原子环境: " + _fmt(only2))
        r["status"], r["consistent"] = "MISMATCH", False
    else:
        r["status"], r["consistent"] = "OK", True
        r["info"]["atom_order_preserved"] = str(g1.symbols == g2.symbols)

    r["reasons"] = reasons
    return r


def _fmt(c: Counter) -> str:
    return ", ".join(f"{k}×{v}" if v > 1 else k for k, v in sorted(c.items()))


# ────────────────────────────────────────────────────────────────────────────
# CLI
# ────────────────────────────────────────────────────────────────────────────
def collect(paths):
    out = []
    for p in paths:
        p = Path(p)
        out.extend(sorted(p.rglob("*.sdf"))) if p.is_dir() else out.append(p)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="检查 SDF title 行 SMILES 与 CT 块的拓扑一致性")
    ap.add_argument("paths", nargs="+", help="sdf 文件或目录")
    ap.add_argument("--csv", help="结果写入 CSV")
    ap.add_argument("-q", "--quiet", action="store_true", help="只打印不一致/出错的记录")
    args = ap.parse_args()

    rows, n_bad = [], 0
    for f in collect(args.paths):
        for k, title, block in iter_records(f):
            res = check_record(title, block)
            bad = res["status"] != "OK"
            n_bad += bad
            rows.append(
                {
                    "file": f.name,
                    "record": k,
                    "smiles": res["smiles"],
                    "status": res["status"],
                    "reason": " | ".join(res["reasons"]),
                    **res["info"],
                }
            )
            if args.quiet and not bad:
                continue
            mark = {"OK": "✓", "MISMATCH": "✗"}.get(res["status"], "!")
            print(f"{mark} {f.name}[{k}]  {res['status']}  {res['smiles']}")
            for key, val in res["info"].items():
                print(f"      {key:<20} {val}   (SMILES/SDF)" if "/" in val else f"      {key:<20} {val}")
            for reason in res["reasons"]:
                print(f"      - {reason}")

    print(f"\n共 {len(rows)} 条记录，{n_bad} 条不一致/异常。")
    if args.csv and rows:
        keys = sorted({k for r in rows for k in r})
        head = ["file", "record", "smiles", "status", "reason"]
        keys = head + [k for k in keys if k not in head]
        with open(args.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)
        print(f"已写入 {args.csv}")
    return 1 if n_bad else 0


if __name__ == "__main__":
    sys.exit(main())
