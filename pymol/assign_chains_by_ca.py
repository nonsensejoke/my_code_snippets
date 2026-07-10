# -*- coding: utf-8 -*-
"""
PyMOL：根据相邻 CA 原子距离自动补全蛋白质 Chain ID
====================================================

用途
----
本脚本用于处理蛋白质结构中 chain ID 缺失或错误的情况。
它依次检查相邻蛋白质残基 CA 原子之间的距离：

- 距离小于或等于 cutoff：视为属于同一条链；
- 距离大于 cutoff：视为出现链边界，后一个残基开始新链。

检测完成后，各连续片段依次被写为 chain A、B、C、D……
即使原结构的第一条链没有 chain ID，第一段也会被明确写回 chain A，
不会继续保留空的 chain 字段。

基本用法
--------
    run /完整路径/assign_chains_by_ca.py
    assign_chains_by_ca system, cutoff=5.0

参数
----
object_name：需要处理的 PyMOL object 名称。
cutoff：相邻 CA 原子间的距离阈值（Å），默认 5.0。
state：用于计算距离的结构状态，默认 1。
selection：需要处理的选择，默认 polymer.protein。
quiet：quiet=0 显示详细结果；quiet=1 减少输出。

检查结果
--------
    print(cmd.get_chains("system"))
    select empty_chain, system and polymer.protein and chain ""
    count_atoms empty_chain

注意事项
--------
1. 几何断点不一定完全等同于真实生物学链。
2. 链内缺失残基可能被误判为新链。
3. 两条链首尾若很接近，可能无法仅凭 CA 距离识别。
4. 脚本按照 cmd.get_model() 返回的 CA 原子顺序判断。
5. 缺少 CA 的残基不会直接参与边界判断。
6. 多状态结构使用指定 state 判断，但 chain 会写回整个 object。
7. 单字符 chain ID 最多支持 62 条链：A-Z、a-z、0-9。
8. 修改 chain 后会自动调用 cmd.sort()。
"""

from math import sqrt
from pymol import cmd


def _distance(coord1, coord2):
    """计算两个三维坐标之间的欧氏距离。"""
    return sqrt(
        (coord1[0] - coord2[0]) ** 2
        + (coord1[1] - coord2[1]) ** 2
        + (coord1[2] - coord2[2]) ** 2
    )


def _make_chain_ids(number):
    """生成适合传统 PDB 使用的单字符 chain ID。"""
    alphabet = (
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz"
        "0123456789"
    )

    if number > len(alphabet):
        raise ValueError(
            "More than 62 chains were detected. "
            "Single-character PDB chain IDs are insufficient."
        )

    return list(alphabet[:number])


def assign_chains_by_ca(
    object_name,
    cutoff=5.0,
    state=1,
    selection="polymer.protein",
    quiet=0,
):
    """
    根据相邻蛋白质残基 CA 原子的距离检测链边界并补全 chain ID。

    第一段始终被写为 chain A，即使原始 chain 字段为空；
    后续连续片段依次写为 B、C、D……
    """

    cutoff = float(cutoff)
    state = int(state)
    quiet = int(quiet)

    if object_name not in cmd.get_object_list():
        raise ValueError(
            "Object {!r} was not found.".format(object_name)
        )

    base_selection = "({}) and ({})".format(
        object_name,
        selection,
    )

    ca_selection = (
        "({}) and name CA and alt ''+A".format(base_selection)
    )

    ca_model = cmd.get_model(ca_selection, state)

    if not ca_model.atom:
        raise ValueError(
            "No CA atoms were found in selection: {}".format(
                ca_selection
            )
        )

    ca_atoms = list(ca_model.atom)

    ca_index_to_chain_number = {}
    chain_number = 0
    previous_atom = None
    breaks = []

    for atom in ca_atoms:
        if previous_atom is not None:
            distance = _distance(
                previous_atom.coord,
                atom.coord,
            )

            if distance > cutoff:
                chain_number += 1
                breaks.append(
                    {
                        "previous_index": previous_atom.index,
                        "previous_resi": previous_atom.resi,
                        "current_index": atom.index,
                        "current_resi": atom.resi,
                        "distance": distance,
                    }
                )

        ca_index_to_chain_number[atom.index] = chain_number
        previous_atom = atom

    number_of_chains = chain_number + 1
    chain_ids = _make_chain_ids(number_of_chains)

    residue_chain_map = {}

    for atom in ca_atoms:
        residue_key = (
            atom.segi,
            atom.resi,
        )

        residue_chain_map[residue_key] = chain_ids[
            ca_index_to_chain_number[atom.index]
        ]

    all_atom_model = cmd.get_model(base_selection, state)
    atom_index_to_chain = {}

    for atom in all_atom_model.atom:
        residue_key = (
            atom.segi,
            atom.resi,
        )

        if residue_key in residue_chain_map:
            atom_index_to_chain[atom.index] = residue_chain_map[
                residue_key
            ]

    if not atom_index_to_chain:
        raise RuntimeError(
            "No protein atoms could be mapped to detected residues."
        )

    cmd.alter(
        base_selection,
        "chain = atom_index_to_chain[index]",
        space={
            "atom_index_to_chain": atom_index_to_chain,
        },
    )

    cmd.sort(object_name)

    if not quiet:
        print("")
        print("Chain assignment completed")
        print("--------------------------")
        print("Object:          {}".format(object_name))
        print("Cutoff:          {:.3f} A".format(cutoff))
        print("State:           {}".format(state))
        print("CA atoms:        {}".format(len(ca_atoms)))
        print("Chains assigned: {}".format(number_of_chains))
        print("Chain IDs:       {}".format(", ".join(chain_ids)))

        if breaks:
            print("")
            print("Detected chain boundaries:")

            for number, item in enumerate(breaks, 1):
                print(
                    "  {}. resi {} -> resi {}: {:.3f} A".format(
                        number,
                        item["previous_resi"],
                        item["current_resi"],
                        item["distance"],
                    )
                )
        else:
            print("")
            print(
                "No distance break was detected; "
                "all selected residues were assigned chain A."
            )

    return {
        "object": object_name,
        "cutoff": cutoff,
        "state": state,
        "chain_ids": chain_ids,
        "breaks": breaks,
        "ca_count": len(ca_atoms),
    }


cmd.extend(
    "assign_chains_by_ca",
    assign_chains_by_ca,
)
