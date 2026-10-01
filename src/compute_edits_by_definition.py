"""
compute_edits_by_definition.py — 4 种原子定义 × 5 个阈值的接触编辑计算
======================================================================
一致定义：所有原子定义均排除氢原子（修正既有实现的 WT含氢/MT无氢 不一致）。

输出 data/edits_by_definition.csv：
  pdb_id, mut_info, mut_idx, wt_has_h, atom_def, threshold,
  wt_edges, mt_edges, n_broken, n_formed
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from contact_graph_defs import compute_all_defs, edges_at, extract_residues

DATA = r'D:\GED_mutation\data'
ATOM_DEFS = ['ca', 'cb', 'centroid', 'allatom']
THRESHOLDS = [6.0, 7.0, 8.0, 9.0, 10.0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', default=os.path.join(DATA, 'edits_by_definition.csv'))
    args = ap.parse_args()

    df = pd.read_csv(os.path.join(DATA, 'benchmarks_s669_clean.csv'))
    align = pd.read_csv(os.path.join(DATA, 'alignment_map_s669.csv'))
    amap = {r['pdb_id']: r for _, r in align.iterrows()}
    mt_dir = os.path.join(DATA, 'mutant_structures_s669')

    wt_cache = {}
    rows = []
    n_done = 0
    for _, row in df.iterrows():
        pid, mut = row['pdb_id'], row['mut_info']
        a = amap.get(pid)
        if a is None:
            continue
        base = f'pdb{pid.lower()}'
        mt_path = os.path.join(mt_dir, pid, mut, 'work', 'out', f'{base}_1.pdb')
        wt_path = os.path.join(DATA, 'structures', f'{base}.ent')
        if not (os.path.exists(mt_path) and os.path.exists(wt_path)):
            continue
        L = len(row['wt_seq'])
        chain, off = a['chain'], int(a['offset'])

        if pid not in wt_cache:
            dw, sw = compute_all_defs(wt_path, chain, off, L, max_cut=max(THRESHOLDS))
            wt_cache[pid] = (dw, sw)
        dw, sw = wt_cache[pid]
        if dw is None:
            continue
        dm, sm = compute_all_defs(mt_path, chain, off, L, max_cut=max(THRESHOLDS))
        if dm is None or sm is None:
            continue
        if len(sw) != len(sm):
            pass  # 残基数不一致时仍比较（索引对齐）

        # WT 是否含氢
        res, _ = extract_residues(wt_path, chain, off, L)
        wt_has_h = 1 if res and any(a_.element == 'H' for r in res for a_ in r.get_atoms()) else 0

        for ad in ATOM_DEFS:
            for th in THRESHOLDS:
                W = edges_at(dw[ad], th)
                M = edges_at(dm[ad], th)
                rows.append({
                    'pdb_id': pid, 'mut_info': mut, 'mut_idx': int(row['_pdb_res_idx']),
                    'wt_has_h': wt_has_h, 'atom_def': ad, 'threshold': th,
                    'wt_edges': len(W), 'mt_edges': len(M),
                    'n_broken': len(W - M), 'n_formed': len(M - W),
                })
        n_done += 1
        if n_done % 50 == 0:
            print(f'  已处理 {n_done} 个突变')

    out = pd.DataFrame(rows)
    out.to_csv(args.out, index=False)
    print(f'\n已写入 {args.out}（{len(out)} 行，{n_done} 个突变 × {len(ATOM_DEFS)} 定义 × {len(THRESHOLDS)} 阈值）')

    # 摘要
    print('\n=== 8 Å 下各定义编辑统计 ===')
    sub = out[out['threshold'] == 8.0]
    print(f'{"atom_def":<10} {"断边均值":>9} {"成边均值":>9} {"≥1变化%":>9} {"WT边均值":>9} {"MT边均值":>9}')
    for ad in ATOM_DEFS:
        s = sub[sub['atom_def'] == ad]
        tot = s['n_broken'] + s['n_formed']
        print(f'{ad:<10} {s["n_broken"].mean():>9.2f} {s["n_formed"].mean():>9.2f} '
              f'{(tot > 0).mean()*100:>8.1f}% {s["wt_edges"].mean():>9.1f} {s["mt_edges"].mean():>9.1f}')


if __name__ == '__main__':
    main()
