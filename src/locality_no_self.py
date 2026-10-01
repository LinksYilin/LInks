# -*- coding: utf-8 -*-
"""locality_no_self.py — 局域性统计：排除与突变残基自身相连的边

原 locality_corrected.py 用 min(||wcoord[i]-mc||, ||wcoord[j]-mc||)，
对与突变残基相连的边恒为 0，导致 358 个突变中 124 个 mean_d_broken = 0，
把"断裂接触平均距离"从 ~7.7 Å 拉到 4.03 Å。

本脚本分别报告：
  A) 全部断裂边（原口径，含自身边）
  B) 排除自身边后（本应报告的"其他接触"距离）
"""
import os

import numpy as np
import pandas as pd

from contact_graph_defs import extract_residues, pairwise_within, repr_coords

D = r'D:\GED_mutation\data'
TH, MAXC = 8.0, 10.0
bench = pd.read_csv(os.path.join(D, 'benchmarks_s669_clean.csv'))
align = pd.read_csv(os.path.join(D, 'alignment_map_s669.csv'))
amap = {r['pdb_id']: r for _, r in align.iterrows()}
mt_dir = os.path.join(D, 'mutant_structures_s669')
cache = {}


def boot_ci(vals, groups, B=2000, seed=0):
    rng = np.random.default_rng(seed)
    g = np.asarray(groups)
    uniq = np.unique(g)
    obs = float(np.mean(vals))
    out = []
    for _ in range(B):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([np.where(g == u)[0] for u in pick])
        if len(idx):
            out.append(float(np.mean(np.asarray(vals)[idx])))
    return obs, float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


rows = []
for _, row in bench.iterrows():
    if '_node_idx_ok' in bench.columns and not bool(row['_node_idx_ok']):
        continue
    pid, mut = row['pdb_id'], str(row['mut_info'])
    a = amap.get(pid)
    if a is None:
        continue
    base = f'pdb{pid.lower()}'
    wt_p = os.path.join(D, 'structures', f'{base}.ent')
    mt_p = os.path.join(mt_dir, pid, mut, 'work', 'out', f'{base}_1.pdb')
    if not (os.path.exists(wt_p) and os.path.exists(mt_p)):
        continue
    L = len(row['wt_seq'])
    ch, off = a['chain'], int(a['offset'])
    if pid not in cache:
        cache[pid] = extract_residues(wt_p, ch, off, L)
    wr, sw = cache[pid]
    if wr is None:
        continue
    mr, sm = extract_residues(mt_p, ch, off, L)
    if mr is None or len(sw) != len(sm):
        continue
    wcoord = repr_coords(wr, 'centroid')
    mcoord = repr_coords(mr, 'centroid')
    W = {k for k, d in pairwise_within(wcoord, MAXC).items() if d < TH}
    M = {k for k, d in pairwise_within(mcoord, MAXC).items() if d < TH}
    broken, kept = W - M, W & M
    mi = int(row['_node_idx'])
    if mi >= len(wcoord):
        continue
    mc = wcoord[mi]

    def d_all(edges):
        return [min(np.linalg.norm(wcoord[i] - mc), np.linalg.norm(wcoord[j] - mc))
                for (i, j) in edges]

    def d_no_self(edges):
        # 排除与突变残基相连的边
        return [min(np.linalg.norm(wcoord[i] - mc), np.linalg.norm(wcoord[j] - mc))
                for (i, j) in edges if i != mi and j != mi]

    if not broken:
        continue
    da, ds = d_all(broken), d_no_self(broken)
    dk = d_all(kept)
    rows.append({
        'pdb_id': pid, 'mut_info': mut,
        'n_broken': len(broken), 'n_broken_nonself': len(ds),
        'mean_all': float(np.mean(da)) if da else np.nan,
        'mean_nonself': float(np.mean(ds)) if ds else np.nan,
        'mean_kept': float(np.mean(dk)) if dk else np.nan,
        'median_nonself': float(np.median(ds)) if ds else np.nan,
        'frac_nonself_lt5': float(np.mean(np.array(ds) < 5)) if ds else np.nan,
        'n_self_edges': len(broken) - len(ds),
    })

m = pd.DataFrame(rows)
m.to_csv(os.path.join(D, 'locality_nonself.csv'), index=False)

print(f'突变数: {len(m)}')
print(f'含自身边的突变: {(m.n_self_edges > 0).sum()} ({100*(m.n_self_edges>0).mean():.1f}%)')
print(f'自身断边总数: {int(m.n_self_edges.sum())} / {int(m.n_broken.sum())} 条断边')
print()
print('=== A) 原口径（全部断边, 逐突变均值）===')
o, lo_, hi_ = boot_ci(m.mean_all.to_numpy(), m.pdb_id.to_numpy())
print(f'  断裂 {o:.2f} Å  [{lo_:.2f}, {hi_:.2f}]')
o2, lo2, hi2 = boot_ci(m.mean_kept.to_numpy(), m.pdb_id.to_numpy())
print(f'  保持 {o2:.2f} Å  [{lo2:.2f}, {hi2:.2f}]')
print()
print('=== B) 排除自身边（应报告的"其他接触"）===')
sub = m.dropna(subset=['mean_nonself'])
o3, lo3, hi3 = boot_ci(sub.mean_nonself.to_numpy(), sub.pdb_id.to_numpy())
print(f'  断裂 {o3:.2f} Å  [{lo3:.2f}, {hi3:.2f}]  (n={len(sub)})')
o4, lo4, hi4 = boot_ci(sub.mean_kept.to_numpy(), sub.pdb_id.to_numpy())
print(f'  保持 {o4:.2f} Å  [{lo4:.2f}, {hi4:.2f}]')
# 逐边加权
w = np.average(m.dropna(subset=['mean_nonself']).mean_nonself,
               weights=m.dropna(subset=['mean_nonself']).n_broken_nonself)
print(f'  逐边加权断裂距离: {w:.2f} Å')
print(f'  <5 Å 比例: {sub.frac_nonself_lt5.mean()*100:.1f}%')
print(f'  断裂距离中位数: {sub.median_nonself.median():.2f} Å')
