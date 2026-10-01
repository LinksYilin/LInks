"""
merge_training.py — 整合 MegaScale + ThermoMutDB 训练集
========================================================
把两个训练集统一成 (接触图 + 突变位点 + ddg) 格式，供 GNN/GEDMut 训练。

统一字段：
  - pdb_id / wt_name：蛋白标识
  - mut_idx：突变残基在接触图中的节点索引（0-based）
  - ddg：标准化后的 ΔΔG（正值=失稳）

用法：
  python merge_training.py \
    --megascale_label data/megascale_clean.csv \
    --megascale_graph_dir data/contact_graphs_megascale \
    --thermomutdb_label data/thermomutdb_aligned.csv \
    --thermomutdb_graph_dir data/contact_graphs_thermomutdb \
    --out data/training_merged.csv
"""
import argparse
import os

import pandas as pd


def merge_megascale(label_csv, graph_dir):
    """MegaScale: pos_1based → mut_idx = pos-1。保留 mt_aa。"""
    df = pd.read_csv(label_csv)
    rows = []
    for _, r in df.iterrows():
        wt = r['WT_name']
        key = wt.replace('/', '_').replace('.pdb', '')
        if not os.path.exists(os.path.join(graph_dir, f'{key}.npz')):
            continue
        rows.append({
            'protein': key,
            'mut_idx': int(r['pos_1based']) - 1,
            'mt_aa': r['mt_aa'],
            'ddg': float(r['ddg']),
            'source': 'megascale',
        })
    return pd.DataFrame(rows)


def merge_thermomutdb(label_csv, graph_dir):
    """ThermoMutDB: _pdb_res_idx 已是接触图中的节点索引。mt_aa 从 mutation_code 提取。"""
    df = pd.read_csv(label_csv)
    rows = []
    for _, r in df.iterrows():
        pid = r['pdb_id']
        if not os.path.exists(os.path.join(graph_dir, f'{pid}.npz')):
            continue
        # mutation_code 如 E49M，mt_aa 是最后一个字符
        mt_aa = str(r['mutation_code'])[-1]
        rows.append({
            'protein': pid,
            'mut_idx': int(r['_pdb_res_idx']),
            'mt_aa': mt_aa,
            'ddg': float(r['ddg']),
            'source': 'thermomutdb',
        })
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--megascale_label', required=True)
    ap.add_argument('--megascale_graph_dir', required=True)
    ap.add_argument('--thermomutdb_label', required=True)
    ap.add_argument('--thermomutdb_graph_dir', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    ms = merge_megascale(args.megascale_label, args.megascale_graph_dir)
    tm = merge_thermomutdb(args.thermomutdb_label, args.thermomutdb_graph_dir)
    merged = pd.concat([ms, tm], ignore_index=True)
    merged.to_csv(args.out, index=False)

    print(f'MegaScale: {len(ms)} 条, ThermoMutDB: {len(tm)} 条')
    print(f'合并: {len(merged)} 条, 蛋白数: {merged["protein"].nunique()}')
    print(f'来源分布: {merged["source"].value_counts().to_dict()}')


if __name__ == '__main__':
    main()
