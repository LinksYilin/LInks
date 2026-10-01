"""
build_mutant_contact_graphs.py — 构建突变体接触图 + 计算真实编辑
==================================================================
批量建模完成后，对每个突变体结构构建接触图，
对比 WT vs MT 得到真实断边/成边，供结构编辑恢复评估。

流程：
  1. 读 S669 干净集（含 pdb_id, mut_info, _chain, _pdb_res_idx）
  2. 读对齐 map（含 offset）
  3. 对每个已生成的突变体结构（mutant_structures_s669/{pid}/{mut}/out/{base}_1.pdb）
     用与 WT 相同的 (chain, offset, length) 构建 MT 接触图
  4. 对比 WT vs MT 接触图 → 真实断边/成边
  5. 输出真实编辑清单

用法：
  python build_mutant_contact_graphs.py \
    --label_csv data/benchmarks_s669_clean.csv \
    --alignment data/alignment_map_s669.csv \
    --wt_graph_dir data/contact_graphs_s669 \
    --mutant_dir data/mutant_structures_s669 \
    --out data/true_edits_s669.csv
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from build_contact_graph import build_graph


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--label_csv', required=True)
    ap.add_argument('--alignment', required=True)
    ap.add_argument('--wt_graph_dir', required=True)
    ap.add_argument('--mutant_dir', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.label_csv)
    align = pd.read_csv(args.alignment)
    align_map = {r['pdb_id']: r for _, r in align.iterrows()}

    rows = []
    for _, row in df.iterrows():
        pid = row['pdb_id']
        mut = row['mut_info']
        a = align_map.get(pid)
        if a is None:
            continue
        # 突变体结构路径
        base = f'pdb{pid.lower()}'
        mt_pdb = os.path.join(args.mutant_dir, pid, mut, 'work', 'out', f'{base}_1.pdb')
        if not os.path.exists(mt_pdb):
            continue
        # 构建 MT 接触图（用与 WT 相同的 chain + offset + length）
        wt_seq_len = len(row['wt_seq'])
        mt_graph = build_graph(mt_pdb, chain_id=a['chain'], offset=int(a['offset']), length=wt_seq_len)
        if mt_graph is None:
            continue
        mt_nodes, mt_ei, mt_ea, mt_meta = mt_graph

        # 读 WT 接触图（侧链质心版）
        wt_path = os.path.join(args.wt_graph_dir, f'{pid}.npz')
        if not os.path.exists(wt_path):
            continue
        wt_data = np.load(wt_path, allow_pickle=True)
        wt_ei = wt_data['edge_index']

        # 对比接触边（无向，归一化为 (min, max) 集合）
        def edge_set(ei):
            s = set()
            for a, b in ei.T:
                s.add((min(int(a), int(b)), max(int(a), int(b))))
            return s

        wt_edges = edge_set(wt_ei)
        mt_edges = edge_set(mt_ei)
        broken = wt_edges - mt_edges
        formed = mt_edges - wt_edges

        rows.append({
            'pdb_id': pid,
            'mut_info': mut,
            'mut_idx': int(row['_pdb_res_idx']),
            'n_broken': len(broken),
            'n_formed': len(formed),
            'broken_edges': sorted(broken),
            'formed_edges': sorted(formed),
        })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(args.out, index=False)
    print(f'处理 {len(out_df)} 个突变（有突变体结构的）')
    if len(out_df) > 0:
        print(f'断边数: 均值 {out_df["n_broken"].mean():.1f}, 中位数 {out_df["n_broken"].median():.0f}')
        print(f'成边数: 均值 {out_df["n_formed"].mean():.1f}, 中位数 {out_df["n_formed"].median():.0f}')
    print(f'已写入 {args.out}')


if __name__ == '__main__':
    main()
