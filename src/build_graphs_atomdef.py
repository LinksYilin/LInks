"""
build_graphs_atomdef.py — 按指定原子定义构建接触图（用于 Cα/Cβ 对照）
=====================================================================
把图写到 `contact_graphs_{set}_{def}` 目录，格式与既有 _sc 图一致
（nodes / edge_index / edge_attr / metadata），供 GNN 训练与评估复用。
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(__file__))
from contact_graph_defs import AA_INDEX, AA_PROPERTIES, extract_residues, repr_coords
from paths import DATA, FIGURES, PUB_FIGURES, ROOT, SRC, TOOLS  # 路径集中管理

SRC_PATH = str(SRC)
DATA_PATH = str(DATA)
ROOT_PATH = str(ROOT)
FIG_PATH = str(FIGURES)
PUB_PATH = str(PUB_FIGURES)
TOOLS_PATH = str(TOOLS)


AA_ORDER = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


def node_feature(aa):
    onehot = np.zeros(20, dtype=np.float32)
    if aa in AA_INDEX:
        onehot[AA_INDEX[aa]] = 1.0
    props = np.array(AA_PROPERTIES.get(aa, [0.0, 100.0, 0.0]), dtype=np.float32)
    props[1] = props[1] / 250.0
    return np.concatenate([onehot, props]).astype(np.float32)


def build(pdb_path, chain, offset, length, atom_def, cutoff=8.0):
    residues, seq = extract_residues(pdb_path, chain, offset, length)
    if residues is None or len(residues) < 2:
        return None
    if atom_def == 'allatom':
        from contact_graph_defs import allatom_min_pairs

        pairs_dict = allatom_min_pairs(residues, cutoff)
        pairs = [(i, j) for (i, j), d in pairs_dict.items() if d < cutoff]
        # 节点位置：allatom 定义下用侧链质心作为该残基的代表点（供 EGNN 使用）
        coords = repr_coords(residues, 'centroid')
    else:
        coords = repr_coords(residues, atom_def)
        tree = cKDTree(coords)
        pairs = tree.query_pairs(r=cutoff, output_type='ndarray')
        pairs = [(int(a), int(b)) for a, b in pairs]
    ei, ea = [], []
    for i, j in pairs:
        ei += [[i, j], [j, i]]
        prop_i = AA_PROPERTIES.get(seq[i], [0, 0, 0])
        prop_j = AA_PROPERTIES.get(seq[j], [0, 0, 0])
        hydro = 1.0 if (prop_i[0] > 1.5 and prop_j[0] > 1.5) else 0.0
        elec = 1.0 if (prop_i[2] * prop_j[2] < 0) else 0.0
        ea += [[hydro, elec], [hydro, elec]]
    nodes = np.stack([node_feature(aa) for aa in seq])
    edge_index = np.array(ei, dtype=np.int64).T if ei else np.zeros((2, 0), dtype=np.int64)
    edge_attr = np.array(ea, dtype=np.float32) if ea else np.zeros((0, 2), dtype=np.float32)
    # coords: (N,3) float32 —— 供坐标感知模型（如 EGNN）使用
    return nodes, edge_index, edge_attr, seq, np.asarray(coords, dtype=np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--atom_def', required=True, choices=['ca', 'cb', 'centroid', 'allatom'])
    ap.add_argument('--cutoff', type=float, default=8.0)
    args = ap.parse_args()
    # 阈值非 8 Å 时把阈值写进目录名，避免覆盖默认图
    sfx_extra = '' if abs(args.cutoff - 8.0) < 1e-6 else f'{args.cutoff:g}A'

    # ---- S669 与 ssym（从基准 CSV + 对齐）----
    for tag, label_csv, align_csv, out_dir in [
        ('s669', 'benchmarks_s669_clean.csv', 'alignment_map_s669.csv',
         os.path.join(DATA, f'contact_graphs_s669_{args.atom_def}{sfx_extra}')),
        ('ssym', 'benchmarks_ssym_clean.csv', 'alignment_map_ssym.csv',
         os.path.join(DATA, f'contact_graphs_ssym_{args.atom_def}{sfx_extra}')),
    ]:
        df = pd.read_csv(os.path.join(DATA, label_csv))
        align = pd.read_csv(os.path.join(DATA, align_csv))
        amap = {r['pdb_id']: r for _, r in align.iterrows()}
        wtseq = {}
        for _, r in df.iterrows():
            wtseq.setdefault(r['pdb_id'], r['wt_seq'])
        os.makedirs(out_dir, exist_ok=True)
        ok = 0
        for pid, seq in wtseq.items():
            a = amap.get(pid)
            if a is None:
                continue
            fp = os.path.join(DATA, 'structures', f'pdb{pid.lower()}.ent')
            if not os.path.exists(fp):
                continue
            g = build(fp, a['chain'], int(a['offset']), len(seq), args.atom_def, args.cutoff)
            if g is None:
                continue
            nodes, ei, ea, s, xy = g
            if s != seq:
                continue
            np.savez(os.path.join(out_dir, f'{pid}.npz'), nodes=nodes, edge_index=ei, edge_attr=ea,
                     coords=xy,
                     metadata=np.array([json.dumps({'pdb': pid, 'chain': a['chain'],
                                                    'atom_def': args.atom_def, 'cutoff': args.cutoff})]))
            ok += 1
        print(f'{tag} ({args.atom_def}): {ok} 图 -> {out_dir}')

    # ---- ThermoMutDB（训练集）----
    tm = pd.read_csv(os.path.join(DATA, 'thermomutdb_aligned_sc.csv'))
    out_tm = os.path.join(DATA, f'contact_graphs_thermomutdb_{args.atom_def}{sfx_extra}')
    os.makedirs(out_tm, exist_ok=True)
    # 需要对齐信息：复用已有的 thermomutdb 对齐结果（_pdb_res_idx 在 aligned 文件里）
    ok = 0
    for pid, grp in tm.groupby('pdb_id'):
        fp = os.path.join(DATA, 'structures', f'pdb{pid.lower()}.ent')
        if not os.path.exists(fp):
            continue
        # 用行内 _chain/offset 信息（若无则尝试从已有质心图目录读取 metadata）
        ref = os.path.join(DATA, 'contact_graphs_thermomutdb_sc', f'{pid}.npz')
        if not os.path.exists(ref):
            continue
        md = json.loads(np.load(ref, allow_pickle=True)['metadata'][0])
        chain = md.get('chain')
        n_res = md.get('n_residues')
        offset = md.get('offset', 0)
        g = build(fp, chain, offset, n_res, args.atom_def, args.cutoff)
        if g is None:
            continue
        nodes, ei, ea, s, xy = g
        np.savez(os.path.join(out_tm, f'{pid}.npz'), nodes=nodes, edge_index=ei, edge_attr=ea,
                 coords=xy,
                 metadata=np.array([json.dumps({'pdb': pid, 'chain': chain,
                                                'atom_def': args.atom_def, 'cutoff': args.cutoff})]))
        ok += 1
    print(f'thermomutdb ({args.atom_def}): {ok} 图 -> {out_tm}')


if __name__ == '__main__':
    main()
