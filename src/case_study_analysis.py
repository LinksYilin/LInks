"""case_study_analysis.py — 候选案例的详细分析（替换 1BFM M35W）"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from contact_graph_defs import extract_residues, pairwise_within, repr_coords

from pathlib import Path as _P
DATA = str(_P(__file__).resolve().parent.parent / 'data')
TH, MAXC = 8.0, 10.0
AA_VOL = {'A':88.6,'R':173.4,'N':114.1,'D':111.1,'C':108.5,'Q':143.8,'E':138.4,'G':60.1,'H':153.2,
          'I':166.7,'L':166.7,'K':168.6,'M':162.9,'F':189.9,'P':112.7,'S':89.0,'T':116.1,'W':227.8,'Y':193.6,'V':140.0}

CASES = [('1R2Y', 'R244E'), ('3O39', 'L32P'), ('1XZO', 'W36A'), ('1BA3', 'H461D'), ('1PFL', 'M113T')]

align = pd.read_csv(os.path.join(DATA, 'alignment_map_s669.csv'))
amap = {r['pdb_id']: r for _, r in align.iterrows()}
ss = pd.read_csv(os.path.join(DATA, 'benchmarks_s669_clean.csv'))
ed = pd.read_csv(os.path.join(DATA, 'edits_corrected.csv'))


def edges(pairs, th=TH):
    return {k for k, d in pairs.items() if d < th}


for pid, mut in CASES:
    a = amap.get(pid)
    src = ss[(ss['pdb_id'] == pid) & (ss['mut_info'] == mut)]
    if a is None or len(src) == 0:
        print(f'{pid} {mut}: 缺数据')
        continue
    src = src.iloc[0]
    base = f'pdb{pid.lower()}'
    wt_path = os.path.join(DATA, 'structures', f'{base}.ent')
    mt_path = os.path.join(DATA, 'mutant_structures_s669', pid, mut, 'work', 'out', f'{base}_1.pdb')
    L = len(src['wt_seq'])
    ch, off = a['chain'], int(a['offset'])
    wr, sw = extract_residues(wt_path, ch, off, L)
    mr, sm = extract_residues(mt_path, ch, off, L)
    if wr is None or mr is None or len(sw) != len(sm):
        print(f'{pid} {mut}: 残基数不匹配')
        continue
    if '_node_idx' not in src.index or pd.isna(src['_node_idx']):
        raise ValueError('缺少 _node_idx')
    mi = int(src['_node_idx'])  # 修正：切片后索引
    wc = repr_coords(wr, 'centroid')
    mc = repr_coords(mr, 'centroid')

    print(f'=== {pid} {mut} (ddg = {src["ddg"]:.2f} kcal/mol) ===')
    print(f'  链 {ch}, 残基 {len(sw)}, 突变位点 idx {mi} ({sw[mi]}->{sm[mi]})')
    print(f'  体积: {AA_VOL.get(sw[mi],0):.1f} -> {AA_VOL.get(sm[mi],0):.1f} Å³ '
          f'(Δ={AA_VOL.get(sm[mi],0)-AA_VOL.get(sw[mi],0):+.1f})')

    res = {}
    for ad in ['ca', 'cb', 'centroid', 'allatom']:
        if ad == 'allatom':
            from contact_graph_defs import allatom_min_pairs
            W = edges(allatom_min_pairs(wr, MAXC))
            M = edges(allatom_min_pairs(mr, MAXC))
        else:
            W = edges(pairwise_within(repr_coords(wr, ad), MAXC))
            M = edges(pairwise_within(repr_coords(mr, ad), MAXC))
        res[ad] = (len(W), len(M), len(W - M), len(M - W))
        print(f'  {ad:<9} WT={res[ad][0]:>5} MT={res[ad][1]:>5} 断={res[ad][2]:>3} 成={res[ad][3]:>3}')

    # 断边距离分布（质心定义）
    W = edges(pairwise_within(wc, MAXC))
    M = edges(pairwise_within(mc, MAXC))
    broken = sorted(W - M)
    kept = W & M
    mc0 = wc[mi]
    db = [min(np.linalg.norm(wc[i] - mc0), np.linalg.norm(wc[j] - mc0)) for i, j in broken]
    dk = [min(np.linalg.norm(wc[i] - mc0), np.linalg.norm(wc[j] - mc0)) for i, j in kept]
    if db:
        print(f'  断边到突变位点: 均值 {np.mean(db):.2f} Å, 中位 {np.median(db):.2f} Å, '
              f'<5Å {np.mean(np.array(db) < 5)*100:.0f}%')
    if dk:
        print(f'  保持边到突变位点: 均值 {np.mean(dk):.2f} Å')
    # 最近的断边
    if db:
        order = np.argsort(db)[:3]
        for k in order:
            i, j = broken[k]
            print(f'    最近断边: {i}-{j} 距离 {db[k]:.2f} Å  (涉及突变位点: {mi in (i, j)})')
    print()
