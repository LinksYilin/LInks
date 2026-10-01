"""
build_contact_graphs.py — 批量构建接触图（结合对齐结果）
========================================================
对干净评估集里的每个蛋白，用对齐的 (chain, offset) 切出 wt_seq 对应片段，
构建接触图并校验：构建出的序列必须与标签 wt_seq 完全一致。

用法：
  python build_contact_graphs.py \
    --clean_csv data/benchmarks_s669_clean.csv \
    --alignment data/alignment_map_s669.csv \
    --struct_dir data/structures \
    --out_dir data/contact_graphs_s669
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

from build_contact_graph import build_graph


def pdb_from_filename(fn):
    return fn.split('.')[0].replace('pdb', '').upper()


def find_pdb_file(struct_dir, pid):
    for fn in os.listdir(struct_dir):
        if (fn.lower().endswith(('.ent', '.pdb'))) and (pdb_from_filename(fn) == pid.upper()):
                return os.path.join(struct_dir, fn)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--clean_csv', required=True)
    ap.add_argument('--alignment', required=True)
    ap.add_argument('--struct_dir', required=True)
    ap.add_argument('--out_dir', required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    df = pd.read_csv(args.clean_csv)
    align = pd.read_csv(args.alignment)
    align_map = {r['pdb_id']: r for _, r in align.iterrows()}

    # 每个蛋白一条 wt_seq（取第一个）
    wtseq_by_pdb = {}
    for _, row in df.iterrows():
        wtseq_by_pdb.setdefault(row['pdb_id'], row['wt_seq'])

    ok, fail = [], []
    for pid, wt_seq in sorted(wtseq_by_pdb.items()):
        a = align_map.get(pid)
        if a is None or a['status'] != 'OK':
            fail.append((pid, '无对齐'))
            continue
        pdb_path = find_pdb_file(args.struct_dir, pid)
        if pdb_path is None:
            fail.append((pid, '无结构文件'))
            continue
        g = build_graph(pdb_path, chain_id=a['chain'], offset=int(a['offset']), length=len(wt_seq))
        if g is None:
            fail.append((pid, '构建失败'))
            continue
        nodes, edge_index, edge_attr, meta = g
        if meta['seq'] != wt_seq:
            fail.append((pid, f"序列不一致 built={len(meta['seq'])} vs wt={len(wt_seq)}"))
            continue
        # 保存
        out_path = os.path.join(args.out_dir, f"{pid}.npz")
        np.savez(out_path, nodes=nodes, edge_index=edge_index, edge_attr=edge_attr,
                 metadata=np.array([json.dumps(meta)]))
        ok.append((pid, meta['n_residues'], meta['n_contacts']))

    print(f'成功构建接触图: {len(ok)}/{len(wtseq_by_pdb)}')
    for pid, n, c in ok[:5]:
        print(f'  ✅ {pid}: {n} 残基, {c} 接触')
    if fail:
        print(f'失败: {len(fail)}')
        for pid, msg in fail:
            print(f'  ❌ {pid}: {msg}')


if __name__ == '__main__':
    main()
