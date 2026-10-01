"""
build_megascale_graphs_atomdef.py — 按原子定义构建 MegaScale 图（Cβ 对照用）
============================================================================
复用 build_megascale_graphs.py 的对齐逻辑，但边定义换成指定 atom_def。
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from build_megascale_graphs import align_ref_to_built
from contact_graph_defs import AA_INDEX, AA_PROPERTIES, extract_residues, repr_coords

TMPDIR = os.path.join(os.environ.get('TEMP', '.'), 'ms_graphs_atomdef')

AA_ORDER = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']


def node_feature(aa):
    onehot = np.zeros(20, dtype=np.float32)
    if aa in AA_INDEX:
        onehot[AA_INDEX[aa]] = 1.0
    props = np.array(AA_PROPERTIES.get(aa, [0.0, 100.0, 0.0]), dtype=np.float32)
    props[1] = props[1] / 250.0
    return np.concatenate([onehot, props]).astype(np.float32)


def build_from_text(pdb_text, name, atom_def, cutoff=8.0):
    os.makedirs(TMPDIR, exist_ok=True)
    tmp = os.path.join(TMPDIR, f'{name.replace("/", "_")}.pdb')
    with open(tmp, 'w') as f:
        f.write(pdb_text)
    try:
        residues, seq = extract_residues(tmp, None, 0, None)
        if residues is None or len(residues) < 2:
            return None
        from scipy.spatial import cKDTree
        if atom_def == 'allatom':
            from contact_graph_defs import allatom_min_pairs
            pd_ = allatom_min_pairs(residues, cutoff)
            pairs = [(i, j) for (i, j), d in pd_.items() if d < cutoff]
            coords = repr_coords(residues, 'centroid')  # 节点位置供 EGNN 使用
        else:
            coords = repr_coords(residues, atom_def)
            tree = cKDTree(coords)
            pairs = [(int(a), int(b)) for a, b in tree.query_pairs(r=cutoff, output_type='ndarray')]
        ei, ea = [], []
        for i, j in pairs:
            ei += [[i, j], [j, i]]
            pi = AA_PROPERTIES.get(seq[i], [0, 0, 0]); pj = AA_PROPERTIES.get(seq[j], [0, 0, 0])
            hy = 1.0 if (pi[0] > 1.5 and pj[0] > 1.5) else 0.0
            el = 1.0 if (pi[2] * pj[2] < 0) else 0.0
            ea += [[hy, el], [hy, el]]
        nodes = np.stack([node_feature(x) for x in seq])
        eidx = np.array(ei, dtype=np.int64).T if ei else np.zeros((2, 0), dtype=np.int64)
        eattr = np.array(ea, dtype=np.float32) if ea else np.zeros((0, 2), dtype=np.float32)
        xy = np.asarray(coords, dtype=np.float32)
        return nodes, eidx, eattr, seq, xy
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--struct_parquet', required=True)
    ap.add_argument('--train_parquet', required=True)
    ap.add_argument('--out_dir', required=True)
    ap.add_argument('--atom_def', required=True)
    ap.add_argument('--cutoff', type=float, default=8.0)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    struct_df = pd.read_parquet(args.struct_parquet)
    struct_map = {r['name']: r['pdb'] for _, r in struct_df.iterrows()}
    train_df = pd.read_parquet(args.train_parquet, columns=['WT_name', 'aa_seq'])
    wt_ref = {}
    for _, r in train_df.iterrows():
        wt_ref.setdefault(r['WT_name'], r['aa_seq'])

    ok = 0
    for wt_name, ref_seq in sorted(wt_ref.items()):
        if wt_name not in struct_map:
            continue
        g = build_from_text(struct_map[wt_name], wt_name, args.atom_def, args.cutoff)
        if g is None:
            continue
        nodes, ei, ea, built_seq, xy = g
        off = align_ref_to_built(ref_seq, built_seq)
        if off is None:
            continue
        meta = {'pdb': wt_name, 'seq': built_seq, 'n_residues': len(built_seq),
                'n_contacts': int(ei.shape[1] // 2), 'atom_def': args.atom_def,
                'ref_len': len(ref_seq), 'align_offset': off}
        np.savez(os.path.join(args.out_dir, wt_name.replace('/', '_').replace('.pdb', '') + '.npz'),
                 nodes=nodes, edge_index=ei, edge_attr=ea, coords=xy,
                 metadata=np.array([json.dumps(meta)]))
        ok += 1
    print(f'MegaScale ({args.atom_def}): {ok} 图 -> {args.out_dir}')


if __name__ == '__main__':
    main()
