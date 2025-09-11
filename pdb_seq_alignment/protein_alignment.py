#!/usr/bin/env python3
"""
蛋白质序列比对工具
基于PDBbind数据库的蛋白质结构比较
"""

import os
import sys
from collections import defaultdict
from typing import List, Tuple, Dict


def parse_pdb_file(pdb_file_path: str) -> List[Dict]:
    """
    解析PDB文件，提取ATOM行信息

    Args:
        pdb_file_path: PDB文件路径

    Returns:
        包含原子信息的字典列表
    """
    atoms = []

    try:
        with open(pdb_file_path, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    # PDB格式说明：http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html
                    # 解析ATOM行
                    atom_info = {
                        'record_type': line[0:6].strip(),  # 记录类型
                        'serial': int(line[6:11].strip()),  # 原子序号
                        'atom_name': line[12:16].strip(),  # 原子名称
                        'alt_loc': line[16],  # 交替位置指示符
                        'res_name': line[17:20].strip(),  # 残基名称
                        'chain_id': line[21],  # 链标识符
                        'res_seq': int(line[22:26].strip()),  # 残基序号
                        'insertion': line[26],  # 插入代码
                        'x': float(line[30:38].strip()),  # X坐标
                        'y': float(line[38:46].strip()),  # Y坐标
                        'z': float(line[46:54].strip()),  # Z坐标
                        'occupancy': float(line[54:60].strip() or '1.0'),  # 占有率
                        'temp_factor': float(line[60:66].strip() or '0.0'),  # 温度因子
                        'element': line[76:78].strip(),  # 元素符号
                        'charge': line[78:80].strip()  # 电荷
                    }
                    atoms.append(atom_info)

    except FileNotFoundError:
        print(f"错误：找不到文件 {pdb_file_path}")
        return []
    except Exception as e:
        print(f"解析PDB文件时出错：{e}")
        return []

    return atoms


def get_chain_residue_counts(atoms: List[Dict]) -> Dict[str, int]:
    """
    获取每个链的残基数量统计

    Args:
        atoms: 原子信息列表

    Returns:
        链ID到残基数量的映射
    """
    chain_residues = defaultdict(set)

    for atom in atoms:
        chain_id = atom['chain_id']
        res_seq = atom['res_seq']
        res_name = atom['res_name']
        # 使用(res_seq, res_name)作为唯一标识，避免重复计数
        chain_residues[chain_id].add((res_seq, res_name))

    return {chain: len(residues) for chain, residues in chain_residues.items()}


def select_main_chain(atoms: List[Dict]) -> str:
    """
    选择残基数量最多的链作为主要链

    Args:
        atoms: 原子信息列表

    Returns:
        主要链的ID
    """
    chain_counts = get_chain_residue_counts(atoms)

    if not chain_counts:
        return ''

    # 选择残基数量最多的链
    main_chain = max(chain_counts.items(), key=lambda x: x[1])
    return main_chain[0]


def extract_residue_sequence(atoms: List[Dict], chain_id: str) -> List[Dict]:
    """
    从指定链中提取残基序列

    Args:
        atoms: 原子信息列表
        chain_id: 链标识符

    Returns:
        残基信息列表（包含res_seq和res_name）
    """
    # 收集指定链的所有残基
    residues = {}
    for atom in atoms:
        if atom['chain_id'] == chain_id:
            res_seq = atom['res_seq']
            res_name = atom['res_name']
            if res_seq not in residues:
                residues[res_seq] = res_name

    # 按残基序号排序并返回残基信息列表
    sorted_residues = sorted(residues.items())
    return [{'res_seq': res_seq, 'res_name': res_name} for res_seq, res_name in sorted_residues]

def find_longest_fuzzy_common_subsequence(seq1: List[Dict], seq2: List[Dict]) -> Tuple[List[Dict], List[Dict], int, int, int]:
    """
    寻找两个列表的最长公共子序列（模糊匹配版本）

    Args:
        seq1: 第一个残基序列 (list of dicts with 'res_name')
        seq2: 第二个残基序列 (list of dicts with 'res_name')

    Returns:
        (seq1中的公共子序列, seq2中的公共子序列, 最长公共子序列长度, seq1中的起始索引, seq2中的起始索引)
    """
    # 将残基序列转换为残基名称列表
    list1 = [res['res_name'] for res in seq1]
    list2 = [res['res_name'] for res in seq2]
    
    m = len(list1)
    n = len(list2)
    max_length = 0
    common_sub1 = []
    common_sub2 = []
    index1 = -1
    index2 = -1
    
    # 从最大可能长度递减到 2
    for k in range(min(m, n), 1, -1):
        # 收集 list1 中所有长度 k 的 (start, end) 对
        pairs1 = set()
        for i in range(m - k + 1):
            start = list1[i]
            end = list1[i + k - 1]
            pairs1.add((start, end))
        
        # 收集 list2 中所有长度 k 的 (start, end) 对
        pairs2 = set()
        for j in range(n - k + 1):
            start = list2[j]
            end = list2[j + k - 1]
            pairs2.add((start, end))
        
        # 检查交集
        intersection = pairs1 & pairs2
        if intersection:
            # 找到一个共同的 (start, end)
            start_val, end_val = next(iter(intersection))
            
            # 在 list1 中找到对应的 i
            for i in range(m - k + 1):
                if list1[i] == start_val and list1[i + k - 1] == end_val:
                    common_sub1 = seq1[i:i + k]
                    index1 = i
                    break
            
            # 在 list2 中找到对应的 j
            for j in range(n - k + 1):
                if list2[j] == start_val and list2[j + k - 1] == end_val:
                    common_sub2 = seq2[j:j + k]
                    index2 = j
                    break
            
            max_length = k
            break  # 既然是从大到小，已是最大
    
    return common_sub1, common_sub2, max_length, index1, index2


def calculate_sequence_similarity(seq1: List[Dict], seq2: List[Dict]) -> Dict:
    """
    计算两个序列的相似度信息

    Args:
        seq1: 第一个残基序列 (list of dicts)
        seq2: 第二个残基序列 (list of dicts)

    Returns:
        包含相似度统计信息的字典
    """
    if not seq1 or not seq2:
        return {
            'similarity_percentage': 0.0,
            'common_subsequence1': [],
            'common_subsequence2': [],
            'common_length': 0,
            'seq1_length': len(seq1),
            'seq2_length': len(seq2)
        }

    # 使用模糊匹配版本
    common_sub1, common_sub2, common_length, start_index1, start_index2 = find_longest_fuzzy_common_subsequence(seq1, seq2)

    # 计算相似度百分比（基于较短序列的长度）
    min_length = min(len(seq1), len(seq2))
    similarity_percentage = (common_length / min_length * 100) if min_length > 0 else 0.0

    result = {
        'similarity_percentage': similarity_percentage,
        'common_length': common_length,
        'seq1_length': len(seq1),
        'seq2_length': len(seq2),
        'common_subsequence1': [],
        'common_subsequence2': [],
        'residue_ids': {}
    }

    if common_length > 0:
        result['common_subsequence1'] = [res['res_name'] for res in common_sub1]
        result['common_subsequence2'] = [res['res_name'] for res in common_sub2]
        
        result['residue_ids'] = {
            'seq1': {
                'start': seq1[start_index1]['res_seq'],
                'end': seq1[start_index1 + common_length - 1]['res_seq']
            },
            'seq2': {
                'start': seq2[start_index2]['res_seq'],
                'end': seq2[start_index2 + common_length - 1]['res_seq']
            }
        }

    return result


def compare_protein_structures(pdb_id1: str, pdb_id2: str, base_dir: str = '.') -> Dict:
    """
    比较两个蛋白质结构的序列相似度

    Args:
        pdb_id1: 第一个PDB ID
        pdb_id2: 第二个PDB ID
        base_dir: 基础目录

    Returns:
        比较结果字典
    """
    result = {
        'pdb_id1': pdb_id1,
        'pdb_id2': pdb_id2,
        'success': False,
        'error_message': '',
        'comparison': {}
    }

    try:
        # 构建文件路径
        pocket1_path = os.path.join(base_dir, pdb_id1, f'{pdb_id1}_pocket.pdb')
        protein1_path = os.path.join(base_dir, pdb_id1, f'{pdb_id1}_protein.pdb')
        pocket2_path = os.path.join(base_dir, pdb_id2, f'{pdb_id2}_pocket.pdb')
        protein2_path = os.path.join(base_dir, pdb_id2, f'{pdb_id2}_protein.pdb')

        # 检查文件是否存在
        for path, name in [(pocket1_path, f'{pdb_id1}_pocket.pdb'),
                          (protein1_path, f'{pdb_id1}_protein.pdb'),
                          (pocket2_path, f'{pdb_id2}_pocket.pdb'),
                          (protein2_path, f'{pdb_id2}_protein.pdb')]:
            if not os.path.exists(path):
                result['error_message'] = f'文件不存在：{name}'
                return result

        # 解析PDB文件
        print(f"正在解析 {pdb_id1} 的PDB文件...")
        pocket1_atoms = parse_pdb_file(pocket1_path)
        protein1_atoms = parse_pdb_file(protein1_path)

        print(f"正在解析 {pdb_id2} 的PDB文件...")
        pocket2_atoms = parse_pdb_file(pocket2_path)
        protein2_atoms = parse_pdb_file(protein2_path)

        # 选择主要链
        chain1 = select_main_chain(pocket1_atoms)
        chain2 = select_main_chain(pocket2_atoms)

        print(f"{pdb_id1} 选择链：{chain1}")
        print(f"{pdb_id2} 选择链：{chain2}")

        # 提取残基序列
        seq1 = extract_residue_sequence(protein1_atoms, chain1)
        seq2 = extract_residue_sequence(protein2_atoms, chain2)

        print(f"{pdb_id1} 残基序列原始长度：{len(seq1)}")
        #print(seq1)
        print(f"{pdb_id2} 残基序列原始长度：{len(seq2)}")

        # 检查 seq1/seq2 中 res_seq 是否连续, 假如出现了res_seq 不连续的情况，那么就在对应的地方添加 一个 'res_name' = '-' 的字典 (最后seq1/seq2还是需要排序)
        # 另外，遇到这种情况的话，还需要额外打印 seq1/seq2的最终长度
        # 处理seq1的不连续性
        if seq1:
            min_res_seq = seq1[0]['res_seq']
            max_res_seq = seq1[-1]['res_seq']
            continuous_seq1 = []
            
            for i in range(min_res_seq, max_res_seq + 1):
                found = False
                for res in seq1:
                    if res['res_seq'] == i:
                        continuous_seq1.append(res)
                        found = True
                        break
                if not found:
                    continuous_seq1.append({'res_seq': i, 'res_name': '-'})
            
            if len(continuous_seq1) != len(seq1):
                print(f"{pdb_id1} 残基序列补充缺失后的长度：{len(continuous_seq1)}")
            seq1 = continuous_seq1

        # 处理seq2的不连续性
        if seq2:
            min_res_seq = seq2[0]['res_seq']
            max_res_seq = seq2[-1]['res_seq']
            continuous_seq2 = []
            
            for i in range(min_res_seq, max_res_seq + 1):
                found = False
                for res in seq2:
                    if res['res_seq'] == i:
                        continuous_seq2.append(res)
                        found = True
                        break
                if not found:
                    continuous_seq2.append({'res_seq': i, 'res_name': '-'})
            
            if len(continuous_seq2) != len(seq2):
                print(f"{pdb_id2} 残基序列补充缺失后的长度：{len(continuous_seq2)}")
            seq2 = continuous_seq2




        # 计算相似度
        comparison = calculate_sequence_similarity(seq1, seq2)

        result.update({
            'success': True,
            'chain1': chain1,
            'chain2': chain2,
            'seq1': [res['res_name'] for res in seq1], # Keep original behavior for this key
            'seq2': [res['res_name'] for res in seq2], # Keep original behavior for this key
            'comparison': comparison
        })

    except Exception as e:
        result['error_message'] = f'比较过程中出错：{str(e)}'

    return result


def print_comparison_result(result: Dict):
    """
    打印比较结果

    Args:
        result: 比较结果字典
    """
    if not result['success']:
        print(f"比较失败：{result['error_message']}")
        return

    print("\n" + "="*60)
    print("蛋白质结构序列比对结果")
    print("="*60)
    print(f"PDB ID 1: {result['pdb_id1']} (链 {result['chain1']})")
    print(f"PDB ID 2: {result['pdb_id2']} (链 {result['chain2']})")
    print()

    comp = result['comparison']
    print(f"序列 1 长度: {comp['seq1_length']}")
    print(f"序列 2 长度: {comp['seq2_length']}")
    print(f"最长公共子序列长度: {comp['common_length']}")
    print(f"子序列覆盖度: {comp['similarity_percentage']:.1f}%")

    if comp['common_subsequence1'] and comp['common_subsequence2']:

        if True: 
            print("\n最长公共子序列:") # too long ...
            print(f"序列 1: {' -> '.join(comp['common_subsequence1'])}")
            print(f"序列 2: {' -> '.join(comp['common_subsequence2'])}")



            # 统计 comp['common_subsequence1'] 和 comp['common_subsequence2'] 有多少个元素在相同位置是相同的
            
            # Calculate exact matches within the common subsequence
            exact_matches = sum(1 for a, b in zip(comp['common_subsequence1'], comp['common_subsequence2']) if a == b)
            print(f"  - 精确匹配数量: {exact_matches} / {comp['common_length']}")
            # 统计 comp['common_subsequence1'] 和 comp['common_subsequence2'] 有多少个元素在相同位置是相同的

            # 计算 exact_matches 除以 result['seq1'] 的长度，得到一个百分比
            exact_match_percentage_seq1 = (exact_matches / len(result['seq1']) * 100) if len(result['seq1']) > 0 else 0.0
            exact_match_percentage_seq2 = (exact_matches / len(result['seq2']) * 100) if len(result['seq2']) > 0 else 0.0
            print(f"相似度: {exact_match_percentage_seq1:.1f}%")
            #print(f"  - 精确匹配占序列2的比例: {exact_match_percentage_seq2:.1f}%")
            

        print("\n最长公共子序列(开头3个和末尾3个):")
        # 获取序列1的开头和结尾
        seq1_start = ' -> '.join(comp['common_subsequence1'][:3])
        seq1_end = ' -> '.join(comp['common_subsequence1'][-3:])
        print(f"序列 1: {seq1_start} ... {seq1_end}")

        # 获取序列2的开头和结尾
        seq2_start = ' -> '.join(comp['common_subsequence2'][:3])
        seq2_end = ' -> '.join(comp['common_subsequence2'][-3:])
        print(f"序列 2: {seq2_start} ... {seq2_end}")

        residue_ids = comp.get('residue_ids', {})
        if residue_ids:
            print(f"  - 在 {result['pdb_id1']} (链 {result['chain1']}) 中对应的残基ID范围: {residue_ids['seq1']['start']} -> {residue_ids['seq1']['end']}")
            print(f"  - 在 {result['pdb_id2']} (链 {result['chain2']}) 中对应的残基ID范围: {residue_ids['seq2']['start']} -> {residue_ids['seq2']['end']}")

    else:
        print("\n未找到公共子序列")

    print("="*60)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("用法: python protein_alignment.py <pdb_id1> <pdb_id2>")
        sys.exit(1)

    pdb_id1 = sys.argv[1]
    pdb_id2 = sys.argv[2]

    # The base directory is the current directory
    base_dir = './data_dir' # fixed and no need to change!

    print("蛋白质序列比对工具")
    print("=" * 40)

    print(f"\n正在比较 {pdb_id1} 和 {pdb_id2}...")
    result = compare_protein_structures(pdb_id1, pdb_id2, base_dir)

    print_comparison_result(result)

    print("\n感谢使用蛋白质序列比对工具！")