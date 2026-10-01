# -*- coding: utf-8 -*-
"""
directional_and_type_analysis.py — 方向性分析 + 接触类型分解
==============================================================
回应审稿建议：
  (A) 方向性：net contact change ΔC = n_formed − n_broken、
      ratio n_broken/n_formed，与 ΔΔG 的关系是否不同于"总数"
  (B) 接触类型分解：broken / formed 按 hydro / elec / other 分层，
      各自与 ΔΔG 的关系

做法：对每个突变，从 WT 与 FoldX 突变体结构重建质心接触图（排氢、8 Å），
      得到 broken/formed 边集，再按边特征分型。
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from paths import DATA
from contact_graph_defs import extract_residues, repr_coords, AA_PROPERTIES

OUT = str(DATA)
TH = 8.0
BACKBONE = {'N', 'CA', 'C', 'O', 'OXT'}
AA20 = 'ACDEFGHIKLMNPQRSTVWY'
CACHE = os.path.join(OUT, '_dirtype_cache.npz')


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def cluster_ci(y, p, pid, B=1000, seed=0):
    rng = np.random.default_rng(seed)
    uniq, inv = np.unique(np.array(pid), return_inverse=True)
    pidx = [np.where(inv == i)[0] for i in range(len(uniq))]
    st = []
    for _ in range(B):
        idx = np.concatenate([pidx[i] for i in rng.integers(0, len(uniq), len(uniq))])
        if len(idx) >= 3:
            v = pearson(y[idx], p[idx])
            if not np.isnan(v):
                st.append(v)
    return (float(np.percentile(st, 2.5)), float(np.percentile(st, 97.5))) if st else (np.nan, np.nan)


def edges_of(pdb, chain, off, L):
    residues, seq = extract_residues(pdb, chain, off, L)
    if residues is None or len(residues) < 2:
        return None, None
    coords = repr_coords(residues, 'centroid')
    from scipy.spatial import cKDTree
    pairs = cKDTree(coords).query_pairs(r=TH, output_type='ndarray')
    out = {}
    for a, b in pairs:
        i, j = int(a), int(b)
        pi = AA_PROPERTIES.get(seq[i], [0, 0, 0])
        pj = AA_PROPERTIES.get(seq[j], [0, 0, 0])
        hy = 1.0 if (pi[0] > 1.5 and pj[0] > 1.5) else 0.0
        el = 1.0 if (pi[2] * pj[2] < 0) else 0.0
        out[(min(i, j), max(i, j))] = (hy, el)
    return out, seq


def main():
    df = pd.read_csv(os.path.join(OUT, 'benchmarks_s669_clean.csv'))
    align = pd.read_csv(os.path.join(OUT, 'alignment_map_s669.csv'))
    amap = {r['pdb_id']: r for _, r in align.iterrows()}
    if '_node_idx_ok' in df.columns:
        df = df[df['_node_idx_ok']].reset_index(drop=True)

    rows = []
    wt_cache = {}
    for _, r in df.iterrows():
        pid, mut = r['pdb_id'], r['mut_info']
        a = amap.get(pid)
        if a is None:
            continue
        L, chain, off = len(r['wt_seq']), a['chain'], int(a['offset'])
        wt_p = os.path.join(OUT, 'structures', f'pdb{pid.lower()}.ent')
        mt_p = os.path.join(OUT, 'mutant_structures_s669', pid, mut, 'work', 'out',
                            f'pdb{pid.lower()}_1.pdb')
        if not (os.path.exists(wt_p) and os.path.exists(mt_p)):
            continue
        if pid not in wt_cache:
            wt_cache[pid] = edges_of(wt_p, chain, off, L)
        W, sw = wt_cache[pid]
        M, sm = edges_of(mt_p, chain, off, L)
        if W is None or M is None or len(sw) != len(sm):
            continue
        broken = set(W) - set(M)
        formed = set(M) - set(W)

        def split(es, src):
            hy = sum(1 for e in es if src[e][0] > 0.5)
            el = sum(1 for e in es if src[e][1] > 0.5)
            return hy, el, len(es) - hy - el

        bh, be, bo = split(broken, W)
        fh, fe, fo = split(formed, M)
        nb, nf = len(broken), len(formed)
        rows.append({
            'pdb_id': pid, 'mut_info': mut, 'ddg': float(r['ddg']),
            'n_broken': nb, 'n_formed': nf, 'n_edit': nb + nf,
            'net_change': nf - nb,
            'ratio': (nb / nf) if nf > 0 else np.nan,
            'broken_hydro': bh, 'broken_elec': be, 'broken_other': bo,
            'formed_hydro': fh, 'formed_elec': fe, 'formed_other': fo,
        })

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT, 'directional_type_analysis.csv'), index=False)
    print(f'有效突变 {len(out)}')

    y = out['ddg'].values
    pid = out['pdb_id'].values
    print('\n=== (A) 方向性：与 ΔΔG 的相关 ===')
    print(f'{"特征":<22} {"r":>8} {"95% CI":>20} {"rho":>8}')
    from scipy.stats import spearmanr
    for col in ['n_broken', 'n_formed', 'n_edit', 'net_change', 'ratio']:
        v = out[col].values.astype(float)
        m = ~np.isnan(v)
        r = pearson(y[m], v[m])
        lo, hi = cluster_ci(y[m], v[m], pid[m])
        rho = spearmanr(y[m], v[m]).statistic
        sig = '显著' if (lo > 0 or hi < 0) else '不显著'
        print(f'  {col:<22} {r:>+8.4f} [{lo:>+6.3f},{hi:>+6.3f}] {rho:>+8.4f}  {sig}')

    print('\n=== (B) 接触类型分解：与 ΔΔG 的相关 ===')
    for col in ['broken_hydro', 'broken_elec', 'broken_other',
                'formed_hydro', 'formed_elec', 'formed_other']:
        v = out[col].values.astype(float)
        m = ~np.isnan(v)
        r = pearson(y[m], v[m])
        lo, hi = cluster_ci(y[m], v[m], pid[m])
        sig = '显著' if (lo > 0 or hi < 0) else '不显著'
        print(f'  {col:<22} {r:>+8.4f} [{lo:>+6.3f},{hi:>+6.3f}]  {sig}')

    print('\n=== 净变化 vs 总量：是否等价 ===')
    vb = out['n_broken'].values.astype(float)
    vn = out['net_change'].values.astype(float)
    print(f'  r(n_broken, ΔΔG) = {pearson(y, vb):+.4f}')
    print(f'  r(net_change, ΔΔG) = {pearson(y, vn):+.4f}')
    print(f'  r(n_broken, n_formed) = {pearson(vb, out["n_formed"].values):+.4f}')
    print(f'\n已保存 directional_type_analysis.csv')


if __name__ == '__main__':
    main()
