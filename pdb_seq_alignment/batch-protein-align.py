#!/usr/bin/env python3
"""
批量蛋白质序列比对工具
读取 TargetInfo.dat 文件，批量运行 protein_alignment.py 脚本，
将 L1~L10 位置的 PDB ID 替换为相似度百分比，输出到 TargetInfo_seq_similarity.dat
"""

import os
import sys
import subprocess
import re
from typing import List, Tuple, Optional


def parse_target_info(file_path: str) -> List[List[str]]:
    """
    解析 TargetInfo.dat 文件
    
    Args:
        file_path: TargetInfo.dat 文件路径
        
    Returns:
        每行数据的列表，每个元素是一行的所有列
    """
    data = []
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                # 跳过注释行和空行
                if line.startswith('#') or not line:
                    continue
                
                # 分割列数据
                columns = line.split()
                if columns:
                    data.append(columns)
                    
    except FileNotFoundError:
        print(f"错误：找不到文件 {file_path}")
        return []
    except Exception as e:
        print(f"解析文件时出错：{e}")
        return []
        
    return data


def run_protein_alignment(pdb_id1: str, pdb_id2: str) -> Optional[float]:
    """
    运行 protein_alignment.py 脚本并提取相似度百分比
    
    Args:
        pdb_id1: 第一个PDB ID
        pdb_id2: 第二个PDB ID
        
    Returns:
        相似度百分比，如果失败则返回None
    """
    try:
        # 运行 protein_alignment.py 脚本
        cmd = ['python3', 'protein_alignment.py', pdb_id1, pdb_id2]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        if result.returncode != 0:
            print(f"警告：比对 {pdb_id1} vs {pdb_id2} 失败: {result.stderr.strip()}")
            return None
            
        # 从输出中提取相似度百分比
        output = result.stdout
        
        # 查找 "相似度: XX.X%" 的模式
        similarity_pattern = r'相似度:\s*(\d+\.?\d*)%'
        match = re.search(similarity_pattern, output)
        
        if match:
            similarity = float(match.group(1))
            return similarity
        else:
            print(f"警告：无法从输出中提取相似度信息 ({pdb_id1} vs {pdb_id2})")
            return None
            
    except subprocess.TimeoutExpired:
        print(f"警告：比对 {pdb_id1} vs {pdb_id2} 超时")
        return None
    except Exception as e:
        print(f"警告：运行比对时出错 ({pdb_id1} vs {pdb_id2}): {e}")
        return None


def format_similarity_percentage(similarity: Optional[float]) -> str:
    """
    格式化相似度百分比为右对齐的3位整数+%
    
    Args:
        similarity: 相似度百分比
        
    Returns:
        格式化后的字符串
    """
    if similarity is None:
        return "  -%"  # 失败时显示为 "  -%"
    
    # 四舍五入到整数
    rounded_similarity = round(similarity)
    
    # 格式化为右对齐的3位数字+%
    return f"{rounded_similarity:3d}%"


def format_pdb_with_similarity(pdb_id: str, similarity: Optional[float]) -> str:
    """
    格式化PDB ID和相似度百分比为 pdb_id(similarity%) 的格式
    
    Args:
        pdb_id: PDB ID
        similarity: 相似度百分比
        
    Returns:
        格式化后的字符串，如 "3ejr(100%)" 或 "3dx2( 90%)"
    """
    if similarity is None:
        return f"{pdb_id}(  -%)"  # 失败时显示为 "pdb_id(  -%)"
    
    # 四舍五入到整数
    rounded_similarity = round(similarity)
    
    # 格式化为 pdb_id(similarity%)，相似度右对齐3位
    return f"{pdb_id}({rounded_similarity:3d}%)"


def process_batch_alignment(input_file: str, output_file1: str, output_file2: str):
    """
    批量处理蛋白质序列比对
    
    Args:
        input_file: 输入文件路径 (TargetInfo.dat)
        output_file1: 输出文件路径1 (TargetInfo_seq_similarity.dat)
        output_file2: 输出文件路径2 (TargetInfo_pdbid_seq_similarity.dat)
    """
    print("批量蛋白质序列比对工具")
    print("=" * 50)
    
    # 解析输入文件
    print(f"正在读取 {input_file}...")
    data = parse_target_info(input_file)
    
    if not data:
        print("错误：无法读取输入文件或文件为空")
        return
        
    print(f"共找到 {len(data)} 个目标蛋白")
    
    # 处理每一行数据
    processed_data1 = []  # 用于 TargetInfo_seq_similarity.dat
    processed_data2 = []  # 用于 TargetInfo_pdbid_seq_similarity.dat
    
    for i, row in enumerate(data, 1):
        if len(row) < 2:
            print(f"警告：第 {i} 行数据不完整，跳过")
            continue
            
        target_pdb = row[0]  # 第一列是目标蛋白
        print(f"\n处理第 {i}/{len(data)} 个目标蛋白: {target_pdb}")
        
        # 处理后的行数据
        new_row1 = [target_pdb]  # 第一列保持不变 (仅相似度)
        new_row2 = [target_pdb]  # 第一列保持不变 (PDB ID + 相似度)
        
        # 处理 L1~L10 列（从第2列开始）
        for j, ligand_pdb in enumerate(row[1:], 1):
            print(f"  正在比对 {target_pdb} vs {ligand_pdb} (L{j})...", end=" ")
            
            # 运行比对
            similarity = run_protein_alignment(target_pdb, ligand_pdb)
            
            # 格式化结果
            formatted_similarity = format_similarity_percentage(similarity)
            formatted_pdb_similarity = format_pdb_with_similarity(ligand_pdb, similarity)
            
            new_row1.append(formatted_similarity)
            new_row2.append(formatted_pdb_similarity)
            
            if similarity is not None:
                print(f"相似度: {similarity:.1f}% -> {formatted_similarity}")
            else:
                print("失败")
        
        processed_data1.append(new_row1)
        processed_data2.append(new_row2)
    
    # 写入输出文件1 (仅相似度百分比)
    print(f"\n正在写入结果到 {output_file1}...")
    
    try:
        with open(output_file1, 'w') as f:
            # 写入头部注释
            f.write("#" + "=" * 79 + "\n")
            f.write("# Target proteins and their sequence similarity percentages\n")
            f.write("# Generated by batch-protein-align.py\n")
            f.write("# Format: Target_PDB  L1_similarity%  L2_similarity%  ...\n")
            f.write("# Similarity percentages are right-aligned 3-digit integers + %\n")
            f.write("#" + "=" * 79 + "\n")
            f.write("#T    L1    L2    L3    L4    L5    L6    L7    L8    L9    L10\n")
            
            # 写入数据
            for row in processed_data1:
                # 使用制表符分隔，确保对齐
                line = "  ".join(f"{col:>4}" for col in row)
                f.write(line + "\n")
                
        print(f"结果已成功写入 {output_file1}")
        
    except Exception as e:
        print(f"错误：写入输出文件1时出错: {e}")
    
    # 写入输出文件2 (PDB ID + 相似度百分比)
    print(f"正在写入结果到 {output_file2}...")
    
    try:
        with open(output_file2, 'w') as f:
            # 写入头部注释
            f.write("#" + "=" * 79 + "\n")
            f.write("# Target proteins and their PDB IDs with sequence similarity percentages\n")
            f.write("# Generated by batch-protein-align.py\n")
            f.write("# Format: Target_PDB  L1_pdb(similarity%)  L2_pdb(similarity%)  ...\n")
            f.write("# Example: 3ejr  3ejr(100%)  3dx2( 90%)  3d4z(100%)  ...\n")
            f.write("#" + "=" * 79 + "\n")
            f.write("#T    L1    L2    L3    L4    L5    L6    L7    L8    L9    L10\n")
            
            # 写入数据
            for row in processed_data2:
                # 使用空格分隔，第一列左对齐，其他列适当间隔
                line = row[0]  # 第一列 (target PDB)
                for col in row[1:]:
                    line += f"  {col}"
                f.write(line + "\n")
                
        print(f"结果已成功写入 {output_file2}")
        print(f"共处理 {len(processed_data1)} 个目标蛋白")
        
    except Exception as e:
        print(f"错误：写入输出文件2时出错: {e}")


def main():
    """主函数"""
    # 检查 protein_alignment.py 是否存在
    if not os.path.exists('protein_alignment.py'):
        print("错误：找不到 protein_alignment.py 脚本")
        print("请确保 protein_alignment.py 在当前目录中")
        sys.exit(1)
    
    # 设置输入输出文件
    input_file = 'TargetInfo.dat'
    output_file1 = 'TargetInfo_seq_similarity.dat'
    output_file2 = 'TargetInfo_pdbid_seq_similarity.dat'
    
    # 检查输入文件是否存在
    if not os.path.exists(input_file):
        print(f"错误：找不到输入文件 {input_file}")
        sys.exit(1)
    
    # 开始批量处理
    process_batch_alignment(input_file, output_file1, output_file2)
    
    print("\n批量处理完成！")


if __name__ == '__main__':
    main()