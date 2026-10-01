"""
build_contact_graph.py — 构建蛋白质残基接触图
================================================
输入: PDB 文件（或含 PDB 文件路径的列表）
输出: 每个蛋白一个 .npz 文件:
    nodes      (N, d)  残基特征矩阵
    edge_index (2, E)  接触边（无向，Cβ-Cβ < 8Å）
    edge_attr  (E, c)  边特征（接触类型）
    metadata   dict    wt_seq, chain, 残基编号等

接触定义: 侧链质心距离 < 8 Å（甘氨酸无侧链，用 Cα）。
这是关键修正：原用 Cβ-Cβ，但 FoldX 侧链重排不改骨架，Cβ 几乎不动，
导致突变体接触图与野生型完全相同（边编辑=0）。侧链质心对侧链重排敏感。
（与 ThermoMPNN 等结构方法的残基接触定义一致）

依赖: biopython, numpy
安装: pip install biopython numpy
"""
import argparse
import json
import os

import numpy as np

try:
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1
except ImportError:
    raise SystemExit("请先安装 biopython: pip install biopython")

AA3TO1 = protein_letters_3to1  # 三字母→单字母字典（如 'MET'->'M'）

# 20 种标准氨基酸的理化属性: [疏水性(Kyte-Doolittle), 体积(Å³), 电荷(pH7)]
AA_PROPERTIES = {
    'A': [1.8,  88.6,  0], 'R': [-4.5, 173.4, 1], 'N': [-3.5, 114.1, 0],
    'D': [-3.5, 111.1, -1], 'C': [2.5, 108.5, 0], 'Q': [-3.5, 143.8, 0],
    'E': [-3.5, 138.4, -1], 'G': [-0.4, 60.1,  0], 'H': [-3.2, 153.2, 0],
    'I': [4.5, 166.7,  0], 'L': [3.8, 166.7,  0], 'K': [-3.9, 168.6, 1],
    'M': [1.9, 162.9,  0], 'F': [2.8, 189.9,  0], 'P': [-1.6, 112.7, 0],
    'S': [-0.8, 89.0,  0], 'T': [-0.7, 116.1,  0], 'W': [-0.9, 227.8, 0],
    'Y': [-1.3, 193.6,  0], 'V': [4.2, 140.0,  0],
}
AA_ORDER = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
AA_INDEX = {aa: i for i, aa in enumerate(AA_ORDER)}

CONTACT_CUTOFF = 8.0  # Å


def residue_node_feature(resname):
    """节点特征 = 20 维 one-hot + 3 维理化属性（归一化）"""
    onehot = np.zeros(20, dtype=np.float32)
    if resname in AA_INDEX:
        onehot[AA_INDEX[resname]] = 1.0
    props = np.array(AA_PROPERTIES.get(resname, [0.0, 100.0, 0.0]), dtype=np.float32)
    # 简单归一化（体积除以 250，疏水/电荷保持原尺度）
    props[1] = props[1] / 250.0
    return np.concatenate([onehot, props]).astype(np.float32)  # 23 维


def get_cb_coord(residue):
    """取侧链质心坐标（side-chain centroid）。

    关键修正（2026-09-29）：原用 Cβ，但 FoldX 侧链重排不改骨架，
    Cβ 位置几乎不变，导致 WT/MT 接触图完全一样（边编辑=0）。
    改用侧链质心（除主链外的重原子质心），对侧链重排敏感。

    ★ 关键修正 2（本次审计）：必须排除氢原子。
    原实现只按原子名排除主链，未排除 H。实验结构常含氢（如 1BFM 有 582 个 H），
    而 FoldX 生成的突变体结构不含氢。若不排除 H，WT 质心（含氢）与 MT 质心（无氢）
    不是同一个量，会产生大量**伪断边**（1BFM M35W：伪断边 26 条，真实为 0 条）。

    甘氨酸（无侧链）回退到 Cα。
    """
    backbone = {'N', 'CA', 'C', 'O', 'OXT'}
    atoms = [a for a in residue.get_atoms()
             if a.get_name() not in backbone and a.element != 'H']
    if atoms:
        return np.mean([a.get_coord() for a in atoms], axis=0)
    if residue.has_id('CA'):
        return residue['CA'].get_coord()
    return None


def classify_contact(aa_i, aa_j):
    """边特征：接触类型 one-hot。
    [氢键供体-受体, 疏水-疏水, 静电(异号电荷), 主链-主链, 其他]
    简化规则：基于理化属性。可用 DSSP 派生更精确版本。
    """
    p_i, p_j = AA_PROPERTIES.get(aa_i, [0,0,0]), AA_PROPERTIES.get(aa_j, [0,0,0])
    hydrophobic = 1.0 if (p_i[0] > 1.5 and p_j[0] > 1.5) else 0.0
    electrostatic = 1.0 if (p_i[2] * p_j[2] < 0) else 0.0  # 异号电荷
    return np.array([hydrophobic, electrostatic], dtype=np.float32)


def build_graph(pdb_path, chain_id=None, offset=0, length=None, return_coords=False):
    """从 PDB 构建接触图。返回 (nodes, edge_index, edge_attr, meta)。

    处理：
      ① 只取第一个模型（NMR 多模型结构取 model 0）；
      ② 残基名三字母→单字母；
      ③ 若给定 offset/length，只切出该片段（对应 wt_seq 在链上的位置）。
         —— 这是序列-结构对齐的关键：wt_seq 常是 PDB 链的一段，不是整条链。
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('prot', pdb_path)
    model = structure[0]  # 只取第一个模型

    # 选链
    chains = [c for c in model.get_chains()]
    if chain_id is not None:
        chains = [c for c in chains if c.id == chain_id]
    if not chains:
        return None
    chain = chains[0]

    residues = []
    aa1_list = []
    for res in chain:
        resname3 = res.get_resname().strip()
        # 特殊处理：MSE（硒代甲硫氨酸）按甲硫氨酸 M 处理（X 射线结构的常见衍生物）
        if resname3 == 'MSE':
            aa1 = 'M'
        else:
            # 用 BioPython 的 is_aa 判断标准氨基酸，protein_letters_3to1 转单字母
            if not is_aa(res, standard=True):
                continue
            aa1 = AA3TO1.get(resname3)
            if aa1 is None or aa1 not in AA_INDEX:
                continue
        if get_cb_coord(res) is None:
            continue
        residues.append(res)
        aa1_list.append(aa1)

    # 切片段（offset 是链上 0-based 起始，length 是残基数）
    if offset > 0 or length is not None:
        end = (offset + length) if length is not None else None
        residues = residues[offset:end]
        aa1_list = aa1_list[offset:end]

    n = len(residues)
    if n < 2:
        return None

    nodes = np.stack([residue_node_feature(aa1) for aa1 in aa1_list])
    coords = np.array([get_cb_coord(res) for res in residues], dtype=np.float32)

    # 距离矩阵 → 接触边
    edges_i, edges_j, edges_attr = [], [], []
    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(coords[i] - coords[j])
            if d < CONTACT_CUTOFF:
                edges_i += [i, j]
                edges_j += [j, i]
                attr = classify_contact(aa1_list[i], aa1_list[j])
                edges_attr += [attr, attr]

    edge_index = np.array([edges_i, edges_j], dtype=np.int64) if edges_i else np.zeros((2, 0), dtype=np.int64)
    edge_attr = np.array(edges_attr, dtype=np.float32) if edges_attr else np.zeros((0, 2), dtype=np.float32)

    meta = {
        'pdb': os.path.basename(pdb_path),
        'chain': chain.id,
        'seq': ''.join(aa1_list),
        'n_residues': n,
        'n_contacts': len(edges_i) // 2,
        'offset': offset,
    }
    if return_coords:
        # 供坐标感知模型（EGNN）使用；coords 与 nodes 行一一对应
        return nodes, edge_index, edge_attr, meta, coords
    return nodes, edge_index, edge_attr, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--pdb_dir', required=True, help='PDB 文件目录')
    ap.add_argument('--out_dir', required=True, help='输出 .npz 目录')
    ap.add_argument('--chain', default=None, help='指定链（默认取第一条链）')
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    for fn in sorted(os.listdir(args.pdb_dir)):
        if not fn.lower().endswith(('.pdb', '.ent')):
            continue
        pdb_path = os.path.join(args.pdb_dir, fn)
        try:
            g = build_graph(pdb_path, args.chain)
        except Exception as e:
            print(f'[跳过] {fn}: {e}')
            continue
        if g is None:
            print(f'[跳过] {fn}: 残基不足')
            continue
        nodes, edge_index, edge_attr, meta = g
        out_path = os.path.join(args.out_dir, os.path.splitext(fn)[0] + '.npz')
        np.savez(out_path, nodes=nodes, edge_index=edge_index, edge_attr=edge_attr,
                 metadata=np.array([json.dumps(meta)]))
        print(f'[完成] {fn}: {meta["n_residues"]} 残基, {meta["n_contacts"]} 接触')


if __name__ == '__main__':
    main()
