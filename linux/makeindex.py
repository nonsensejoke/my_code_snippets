#!/usr/bin/env python3
"""
文件索引生成脚本
用法:
    python make_index.py <相对路径> [选项]

示例:
    python make_index.py a/b/c
    python make_index.py a/b/c --root /data/molecules --output /data/indexes
    python make_index.py a/b/c --ext .pdb .mol2 .sdf
"""

import os
import sys
import argparse
from datetime import datetime
from pathlib import Path

# ──────────────────────────────────────────────
# 默认配置（可直接修改这里）
# ──────────────────────────────────────────────
DEFAULT_ROOT        = "/path/to/source"       # 默认 source root 目录
DEFAULT_OUTPUT      = "/path/to/output"       # 默认输出目录
DEFAULT_EXTENSIONS  = {".pdb", ".pdbqt"}      # 默认支持的文件后缀
DEFAULT_INDEX_NAME  = "index.txt"             # 默认索引文件名
# ──────────────────────────────────────────────


def human_size(nbytes: int) -> str:
    """将字节数转换为人类可读格式（如 1.23 MB）"""
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if nbytes < 1024:
            return f"{nbytes:.2f} {unit}" if unit != "B" else f"{nbytes} B"
        nbytes /= 1024
    return f"{nbytes:.2f} PB"


def build_index(
    rel_path: str,
    root: str,
    output: str,
    extensions: set[str],
    index_name: str,
    dry_run: bool = False,
) -> int:
    """
    在 root/rel_path/ 下扫描文件并生成索引。

    返回值: 写入的文件行数（不含表头）
    """
    src_dir = Path(root) / rel_path
    out_dir = Path(output) / rel_path

    if not src_dir.is_dir():
        print(f"[错误] 源目录不存在: {src_dir}", file=sys.stderr)
        sys.exit(1)

    # 扫描符合后缀的文件（只取直接子文件，不递归）
    entries = []
    for item in sorted(src_dir.iterdir()):
        if item.is_file() and item.suffix.lower() in extensions:
            stat = item.stat()
            size_str  = human_size(stat.st_size)
            mtime_str = datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d %H:%M:%S")
            entries.append((item.name, size_str, mtime_str))

    if not entries:
        print(f"[警告] 未找到符合条件的文件: {src_dir}")
        print(f"       支持的后缀: {', '.join(sorted(extensions))}")
        return 0

    index_path = out_dir / index_name

    if dry_run:
        print(f"[预览] 将写入 {index_path}  ({len(entries)} 条记录)")
        print(f"{'文件名':<40}\t{'大小':>12}\t最后修改时间")
        print("-" * 72)
        for name, size, mtime in entries:
            print(f"{name:<40}\t{size:>12}\t{mtime}")
        return len(entries)

    # 创建输出目录
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(index_path, "w", encoding="utf-8", newline="") as f:
        # TSV 表头
        f.write("filename\tsize\tlast_modified\n")
        for name, size, mtime in entries:
            f.write(f"{name}\t{size}\t{mtime}\n")

    print(f"[完成] 索引已写入: {index_path}  ({len(entries)} 个文件)")
    return len(entries)


def main():
    parser = argparse.ArgumentParser(
        description="为指定子目录生成文件索引（TSV 格式）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "rel_path",
        help="相对于 root 的子目录路径，例如: a/b/c",
    )
    parser.add_argument(
        "--root", "-r",
        default=DEFAULT_ROOT,
        metavar="DIR",
        help=f"source root 目录（默认: {DEFAULT_ROOT}）",
    )
    parser.add_argument(
        "--output", "-o",
        default=DEFAULT_OUTPUT,
        metavar="DIR",
        help=f"索引输出根目录（默认: {DEFAULT_OUTPUT}）",
    )
    parser.add_argument(
        "--ext", "-e",
        nargs="+",
        default=sorted(DEFAULT_EXTENSIONS),
        metavar="EXT",
        help=f"要索引的文件后缀（默认: {' '.join(sorted(DEFAULT_EXTENSIONS))}）",
    )
    parser.add_argument(
        "--index-name", "-n",
        default=DEFAULT_INDEX_NAME,
        metavar="NAME",
        help=f"索引文件名（默认: {DEFAULT_INDEX_NAME}）",
    )
    parser.add_argument(
        "--dry-run", "-d",
        action="store_true",
        help="预览模式：只打印结果，不写入文件",
    )

    args = parser.parse_args()

    # 统一后缀格式（确保以 . 开头，小写）
    extensions = set()
    for e in args.ext:
        e = e.strip()
        if not e.startswith("."):
            e = "." + e
        extensions.add(e.lower())

    build_index(
        rel_path   = args.rel_path,
        root       = args.root,
        output     = args.output,
        extensions = extensions,
        index_name = args.index_name,
        dry_run    = args.dry_run,
    )


if __name__ == "__main__":
    main()
  
