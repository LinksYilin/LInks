"""
mechanism_analysis_corrected.py — 修正后的机制分析（论文核心论点）
==================================================================
基于 data/edits_corrected.csv（一致氢排除 + 质控后 505 个突变对）。
计算：
  (1) 编辑特征与 ΔΔG 的关联（Pearson/Spearman + 蛋白簇 bootstrap CI）
  (2) 控制理化特征后的偏相关
  (3) 编辑特征与理化特征的冗余（相关矩阵 + VIF）
  (4) 按体积变化幅度分层
  (5) 断裂接触距突变位点的距离依赖
"""
import os
import sys

import numpy as np
import pandas as pd
from numpy.linalg import lstsq
from scipy.stats import spearmanr

from paths import DATA, FIGURES, PUB_FIGURES, ROOT, SRC, TOOLS  # 路径集中管理

SRC_PATH = str(SRC)
DATA_PATH = str(DATA)
ROOT_PATH = str(ROOT)
FIG_PATH = str(FIGURES)
PUB_PATH = str(PUB_FIGURES)
TOOLS_PATH = str(TOOLS)


DATA = DATA_PATH  # 来自 paths.py，可用 GED_ROOT 环境变量覆盖

KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,
      'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
VOL = {'A':88.6,'R':173.4,'N':114.1,'D':111.1,'C':108.5,'Q':143.8,'E':138.4,'G':60.1,'H':153.2,
       'I':166.7,'L':166.7,'K':168.6,'M':162.9,'F':189.9,'P':112.7,'S':89.0,'T':116.1,'W':227.8,'Y':193.6,'V':140.0}
CHG = {'D':-1,'E':-1,'K':1,'R':1}

def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])

def partial_corr(x, y, Z):
    x, y, Z = np.asarray(x, float), np.asarray(y, float), np.asarray(Z, float)
    Z1 = np.column_stack([np.ones(len(Z)), Z])
    rx = x - Z1 @ lstsq(Z1, x, rcond=None)[0]
    ry = y - Z1 @ lstsq(Z1, y, rcond=None)[0]
    return pearson(rx, ry)

def boot_ci(y, x, pid, metric, B=1000, seed=0):
    rng = np.random.default_rng(seed)
    prot = np.array(pid)
    uniq, inv = np.unique(prot, return_inverse=True)
    pidx = [np.where(inv == i)[0] for i in range(len(uniq))]
    st = []
    for _ in range(B):
        idx = np.concatenate([pidx[i] for i in rng.integers(0, len(uniq), len(uniq))])
        if len(idx) < 3:
            continue
        v = metric(y[idx], x[idx])
        if not np.isnan(v):
            st.append(v)
    if not st:
        return float('nan'), float('nan')
    return float(np.percentile(st, 2.5)), float(np.percentile(st, 97.5))

def vif(X):
    """方差膨胀因子：VIF_j = 1/(1-R2_j)，R2_j 为第 j 列对其余列回归的 R2。"""
    out = []
    for j in range(X.shape[1]):
        others = np.delete(X, j, axis=1)
        Z = np.column_stack([np.ones(len(others)), others])
        beta = lstsq(Z, X[:, j], rcond=None)[0]
        pred = Z @ beta
        ss_res = np.sum((X[:, j] - pred) ** 2)
        ss_tot = np.sum((X[:, j] - X[:, j].mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
        out.append(1.0 / (1.0 - r2) if r2 < 1 else np.inf)
    return out

def main():
    ed = pd.read_csv(os.path.join(DATA, 'edits_corrected.csv'))
    ed = ed[(ed['atom_def'] == 'centroid') & (ed['threshold'] == 8.0)]
    # edits_corrected.csv 已含 ddg，直接使用（无需再 merge）
    m = ed.dropna(subset=['ddg']).reset_index(drop=True)
    print(f'分析样本: {len(m)} 个突变对（质心图, 8 Å, 一致氢排除）')

    m['wt'] = m['mut_info'].str[0]
    m['mt'] = m['mut_info'].str[-1]
    m['d_vol'] = m.apply(lambda r: VOL.get(r['wt'], 140) - VOL.get(r['mt'], 140), axis=1)
    m['d_kd'] = m.apply(lambda r: KD.get(r['wt'], 0) - KD.get(r['mt'], 0), axis=1)
    m['vol_mt'] = m['mt'].map(VOL).fillna(140)
    m['kd_mt'] = m['mt'].map(KD).fillna(0)
    m['chg_mt'] = m['mt'].map(CHG).fillna(0)
    m['abs_dvol'] = m['d_vol'].abs()

    y = m['ddg'].values
    pid = m['pdb_id'].values
    phys_cols = ['d_vol', 'd_kd', 'vol_mt', 'kd_mt', 'chg_mt']
    P = m[phys_cols].values

    print('\n=== (1) 编辑特征 vs ΔΔG（Pearson + 蛋白簇 bootstrap 95% CI）===')
    print(f'{"特征":<12} {"Pearson":>9} {"95% CI":>19} {"Spearman":>9} {"偏相关|理化":>12}')
    for name in ['n_broken', 'n_formed', 'n_edit']:
        x = m[name].values
        r = pearson(x, y)
        lo, hi = boot_ci(y, x, pid, lambda a, b: pearson(b, a))
        rho = float(spearmanr(x, y)[0])
        pr = partial_corr(x, y, P)
        print(f'{name:<12} {r:>9.3f} [{lo:>7.3f},{hi:>7.3f}] {rho:>9.3f} {pr:>12.3f}')

    print('\n=== (2) 编辑特征与理化特征的最大 |相关|（冗余性）===')
    for name in ['n_broken', 'n_formed', 'n_edit']:
        cors = {c: pearson(m[name].values, m[c].values) for c in phys_cols}
        mx = max(cors.items(), key=lambda kv: abs(kv[1]))
        print(f'  {name}: 最大 |r| = {abs(mx[1]):.3f} (与 {mx[0]})')

    print('\n=== (3) VIF（编辑特征 + 理化特征联合）===')
    X = np.column_stack([m['n_broken'].values, m['n_formed'].values, P])
    names = ['n_broken', 'n_formed'] + phys_cols
    for nm, v in zip(names, vif(X)):
        print(f'  {nm:<10} VIF = {v:.2f}')

    print('\n=== (4) 按体积变化幅度分层（|Δvol| 三分位）===')
    q1, q2 = m['abs_dvol'].quantile([1/3, 2/3])
    for lab, sub in [('small', m[m['abs_dvol'] <= q1]),
                     ('medium', m[(m['abs_dvol'] > q1) & (m['abs_dvol'] <= q2)]),
                     ('large', m[m['abs_dvol'] > q2])]:
        if len(sub) < 10:
            continue
        r = pearson(sub['n_edit'].values, sub['ddg'].values)
        rho = float(spearmanr(sub['n_edit'], sub['ddg'])[0])
        print(f'  {lab:<7} n={len(sub):>3}  |Δvol| [{sub["abs_dvol"].min():>5.1f},{sub["abs_dvol"].max():>5.1f}]  '
              f'r={r:>6.3f}  rho={rho:>6.3f}')

    # (5) 距离依赖：需从 broken_edges 重算（此文件未存边列表，改用 3D 重算）
    print('\n=== (5) 断裂接触距突变位点距离 vs ΔΔG ===')

    sys.path.insert(0, os.path.dirname(__file__))
    from contact_graph_defs import extract_residues, pairwise_within, repr_coords

    align = pd.read_csv(os.path.join(DATA, 'alignment_map_s669.csv'))
    amap = {r['pdb_id']: r for _, r in align.iterrows()}
    ss = pd.read_csv(os.path.join(DATA, 'benchmarks_s669_clean.csv'))
    ssmap = {(r['pdb_id'], r['mut_info']): r for _, r in ss.iterrows()}

    dists, ddgs, pids, neds = [], [], [], []
    cache = {}

    def edges_from(pairs, th=8.0):
        return {k for k, d in pairs.items() if d < th}

    for _, row in m.iterrows():
        pid_, mut = row['pdb_id'], row['mut_info']
        a = amap.get(pid_)
        src = ssmap.get((pid_, mut))
        if a is None or src is None:
            continue
        base = f'pdb{pid_.lower()}'
        mt_path = os.path.join(DATA, 'mutant_structures_s669', pid_, mut, 'work', 'out', f'{base}_1.pdb')
        wt_path = os.path.join(DATA, 'structures', f'{base}.ent')
        if not (os.path.exists(mt_path) and os.path.exists(wt_path)):
            continue
        L = len(src['wt_seq'])
        ch, off = a['chain'], int(a['offset'])
        if pid_ not in cache:
            res, _ = extract_residues(wt_path, ch, off, L)
            cache[pid_] = (res, repr_coords(res, 'centroid'))
        res, coords = cache[pid_]
        if cache[pid_][0] is None:
            continue
        mt_res, _ = extract_residues(mt_path, ch, off, L)
        if mt_res is None or len(mt_res) != len(res):
            continue
        W = set()
        Wp = pairwise_within(coords, 10.0)
        W = edges_from(Wp)
        Mp = pairwise_within(repr_coords(mt_res, 'centroid'), 10.0)
        Mset = edges_from(Mp)
        broken = W - Mset
        mi = int(src['_node_idx'])  # 修正：用切片后索引
        if mi >= len(coords) or not broken:
            continue
        mc = coords[mi]
        ds = []
        for (i, j) in broken:
            if i < len(coords) and j < len(coords):
                ds.append(min(np.linalg.norm(coords[i] - mc), np.linalg.norm(coords[j] - mc)))
        if ds:
            dists.append(np.mean(ds))
            ddgs.append(row['ddg'])
            pids.append(pid_)
            neds.append(row['n_edit'])

    dists = np.array(dists); ddgs = np.array(ddgs); neds = np.array(neds)
    print(f'  有效样本 {len(dists)}')
    print(f'  平均断边距离 vs ΔΔG: r = {pearson(dists, ddgs):.3f}')
    print(f'  平均断边距离 vs 编辑总数: r = {pearson(dists, neds):.3f}')
    for lo_, hi_ in [(0, 5), (5, 10), (10, 15), (15, 20), (20, 100)]:
        msk = (dists >= lo_) & (dists < hi_)
        if msk.sum() >= 10:
            print(f'    [{lo_:>2},{hi_:>3}) Å: n={msk.sum():>3}, mean ΔΔG={ddgs[msk].mean():>5.2f}, '
                  f'r(dist,ddg)={pearson(dists[msk], ddgs[msk]):>6.3f}')

    # 保存
    m.to_csv(os.path.join(DATA, 'mechanism_corrected.csv'), index=False)
    print(f'\n已保存 {os.path.join(DATA, "mechanism_corrected.csv")}')

if __name__ == '__main__':
    main()
