#!/usr/bin/env python3
"""
patch_similarity.py
用法: python patch_similarity.py <pdb1> <pdb2> [--max_dist 16.0] [--bin_width 1.0] [--sigma 1.0]
"""

import argparse
import sys
import numpy as np
from itertools import combinations_with_replacement
from pathlib import Path

# ─────────────────────────────────────────────
# 1. 氨基酸映射（含 AMBER 力场变体）
# ─────────────────────────────────────────────

THREE_TO_ONE = {
    # 标准 20 种
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    # HIS 变体
    'HSD': 'H', 'HSE': 'H', 'HSP': 'H',
    'HIE': 'H', 'HID': 'H', 'HIP': 'H',
    # 其他常见变体
    'MSE': 'M', 'CSE': 'C', 'SEC': 'C',
    # AMBER 力场质子化态变体
    'LYN': 'K',   # 中性赖氨酸
    'GLH': 'E',   # 质子化谷氨酸
    'ASH': 'D',   # 质子化天冬氨酸
    'CYX': 'C',   # 参与二硫键的半胱氨酸
}

AA_GROUPS = {
    'hydrophobic': set('AVLIM'),
    'aromatic':    set('FWY'),
    'polar':       set('STNQ'),
    'positive':    set('KRH'),
    'negative':    set('DE'),
    'sulfur':      set('C'),
    'special':     set('GP'),
    'unknown':     set('X'),
}
GROUP_ID   = {aa: i for i, (_, s) in enumerate(AA_GROUPS.items()) for aa in s}
GROUP_NAME = {i: name for i, (name, _) in enumerate(AA_GROUPS.items())}
K          = len(AA_GROUPS)
PAIRS      = list(combinations_with_replacement(range(K), 2))
PAIR_IDX   = {p: i for i, p in enumerate(PAIRS)}

# ─────────────────────────────────────────────
# 2. PDB 解析
# ─────────────────────────────────────────────

def parse_pdb(path):
    """
    规则：
    - 有 CA 原子：
        - resname 在 THREE_TO_ONE 中 → 使用对应单字母
        - resname 不在字典中         → 使用 'X'
    - 无 CA 原子：跳过该残基

    返回 list of (one_letter: str, xyz: np.ndarray(3,))
    每个 res_id 只取第一个 CA（多构象取 A）
    """
    # 第一遍：收集每个 res_id 的 resname 和 CA 坐标
    res_info   = {}   # res_id → resname（第一次见到时记录）
    res_ca     = {}   # res_id → xyz（第一个 CA）

    with open(path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            resname  = line[17:20].strip()
            chain    = line[21]
            res_id   = (chain, line[22:26].strip(), line[26].strip())
            atom     = line[12:16].strip()

            if res_id not in res_info:
                res_info[res_id] = resname

            if atom == 'CA' and res_id not in res_ca:
                xyz = np.array([float(line[30:38]),
                                float(line[38:46]),
                                float(line[46:54])])
                res_ca[res_id] = xyz

    # 第二遍：只保留有 CA 的残基，resname 未知则 'X'
    residues = []
    unknown_resnames = set()
    for res_id in res_info:
        if res_id not in res_ca:
            continue   # 无 CA → 跳过
        resname    = res_info[res_id]
        one_letter = THREE_TO_ONE.get(resname)
        if one_letter is None:
            one_letter = 'X'
            unknown_resnames.add(resname)
        residues.append((one_letter, res_ca[res_id]))

    return residues, unknown_resnames

# ─────────────────────────────────────────────
# 3. 指纹计算
# ─────────────────────────────────────────────

def interface_fingerprint(residues, max_dist=16.0, bin_width=1.0, sigma=None):
    """
    sigma=None → 硬 bin；否则高斯展宽（单位 Å）
    返回归一化 1-D numpy 数组，长度 = n_pairs × n_bins
    """
    bins        = np.arange(0, max_dist + bin_width, bin_width)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    n_bins      = len(bin_centers)

    labels = np.array([GROUP_ID.get(aa, GROUP_ID['X']) for aa, _ in residues])
    coords = np.array([xyz for _, xyz in residues], dtype=float)
    D      = np.linalg.norm(coords[:, None] - coords[None, :], axis=-1)

    fp = np.zeros((len(PAIRS), n_bins))
    n  = len(residues)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = sorted((labels[i], labels[j]))
            d    = D[i, j]
            if d >= bins[-1]:
                continue
            if sigma is None:
                fp[PAIR_IDX[(a, b)], np.searchsorted(bins, d) - 1] += 1
            else:
                fp[PAIR_IDX[(a, b)]] += np.exp(-0.5 * ((bin_centers - d) / sigma) ** 2)

    fp /= (fp.sum() + 1e-8)
    return fp.flatten()

def cosine_similarity(a, b):
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-8))

# ─────────────────────────────────────────────
# 4. CLI
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Compute interface fingerprint similarity between two protein patches.')
    parser.add_argument('pdb1',                  help='第一个 PDB 文件路径')
    parser.add_argument('pdb2',                  help='第二个 PDB 文件路径')
    parser.add_argument('--max_dist',  type=float, default=16.0, help='距离直方图上限 (Å)，默认 16.0')
    parser.add_argument('--bin_width', type=float, default=1.0,  help='bin 宽度 (Å)，默认 1.0')
    parser.add_argument('--sigma',     type=float, default=1.0,  help='高斯展宽 σ (Å)，默认 1.0；设为 0 则使用硬 bin')
    parser.add_argument('--output',    type=str,   default=None, help='输出文本文件路径（默认写到 stdout）')
    args = parser.parse_args()

    sigma = args.sigma if args.sigma > 0 else None

    # 解析
    res1, unk1 = parse_pdb(args.pdb1)
    res2, unk2 = parse_pdb(args.pdb2)

    lines = []
    lines.append("=" * 60)
    lines.append("Patch Fingerprint Similarity")
    lines.append("=" * 60)

    n_bins  = int(args.max_dist / args.bin_width)
    n_pairs = len(PAIRS)
    sigma_str = f"{args.sigma} Å" if sigma is not None else "None (hard bin)"
    lines.append(f"参数：max_dist={args.max_dist} Å, bin_width={args.bin_width} Å, σ={sigma_str}")
    lines.append(f"指纹维度：{n_pairs} pairs × {n_bins} bins = {n_pairs * n_bins} 维")
    lines.append("")

    for label, path, res, unk in [
        (Path(args.pdb1).name, args.pdb1, res1, unk1),
        (Path(args.pdb2).name, args.pdb2, res2, unk2),
    ]:
        aa_str = ' '.join(f"{r[0]}({GROUP_NAME[GROUP_ID.get(r[0], GROUP_ID['X'])]})"
                          for r in res)
        lines.append(f"{label}: {len(res)} 残基 → {aa_str}")
        if unk:
            lines.append(f"  [警告] 未识别 resname（已映射为 X）: {', '.join(sorted(unk))}")

    lines.append("")

    # 计算指纹与相似度
    fp1  = interface_fingerprint(res1, args.max_dist, args.bin_width, sigma)
    fp2  = interface_fingerprint(res2, args.max_dist, args.bin_width, sigma)
    sim  = cosine_similarity(fp1, fp2)

    lines.append(f"Cosine Similarity: {sim:.4f}")
    lines.append("=" * 60)

    output = "\n".join(lines)

    if args.output:
        with open(args.output, 'w') as f:
            f.write(output + "\n")
        print(f"结果已写入 {args.output}")
    else:
        print(output)

if __name__ == '__main__':
    main()