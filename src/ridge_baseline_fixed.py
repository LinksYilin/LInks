# -*- coding: utf-8 -*-
"""
ridge_baseline_fixed.py — 在修复后的训练数据上重训理化岭回归基线
====================================================================
论文旧 Table 1 的 Ridge r=0.392 基于含氢 bug 的训练数据，必须重算。
本脚本用与其他模型**完全相同**的训练集划分与评估交集。
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from paths import DATA

OUT = str(DATA)
AA20 = 'ACDEFGHIKLMNPQRSTVWY'
KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4,
          H=-3.2, I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8,
          T=-0.7, W=-0.9, Y=-1.3, V=4.2)
VOL = dict(A=88.6, R=173.4, N=114.1, D=111.1, C=108.5, Q=143.8, E=138.4,
           G=60.1, H=153.2, I=166.7, L=166.7, K=168.6, M=162.9, F=189.9,
           P=112.7, S=89.0, T=116.1, W=227.8, Y=193.6, V=140.0)
CHG = dict(D=-1, E=-1, K=1, R=1)


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


def feats(wt, mt):
    return [KD.get(mt, 0) - KD.get(wt, 0), VOL.get(mt, 0) - VOL.get(wt, 0),
            KD.get(mt, 0), VOL.get(mt, 0), CHG.get(mt, 0) - CHG.get(wt, 0)]


def wt_from_graph(z, mi):
    return AA20[int(np.argmax(z['nodes'][mi, :20]))]


def main():
    from ladder_common import graph_dir
    from sklearn.linear_model import RidgeCV

    tr = pd.read_csv(os.path.join(OUT, 'training_merged_noleak_sc.csv'))
    ms, tm = tr[tr['source'] == 'megascale'], tr[tr['source'] == 'thermomutdb']
    n = min(len(ms), len(tm))
    tr = pd.concat([ms.sample(n=n, random_state=42), tm], ignore_index=True)

    Xtr, ytr = [], []
    for _, r in tr.iterrows():
        src, prot = r['source'], str(r['protein'])
        p = os.path.join(graph_dir(src, 'centroid'), f'{prot}.npz')
        if not os.path.exists(p):
            continue
        z = np.load(p, allow_pickle=True)
        mi = int(r['mut_idx'])
        if not (0 <= mi < z['nodes'].shape[0]):
            continue
        wt = wt_from_graph(z, mi)
        Xtr.append(feats(wt, str(r['mt_aa'])))
        ytr.append(float(r['ddg']))
    Xtr, ytr = np.array(Xtr), np.array(ytr)
    print(f'训练 {Xtr.shape}')

    m = RidgeCV(alphas=np.logspace(-3, 3, 30)).fit(Xtr, ytr)

    rows = []
    for bench in ['s669', 'ssym']:
        label = f'benchmarks_{bench}_clean.csv'
        df = pd.read_csv(os.path.join(OUT, label))
        if '_node_idx_ok' in df.columns:
            df = df[df['_node_idx_ok']]
        Xte, yte, pidv = [], [], []
        for _, r in df.iterrows():
            pid = r['pdb_id']
            p = os.path.join(graph_dir(bench, 'centroid'), f'{pid}.npz')
            if not os.path.exists(p):
                continue
            z = np.load(p, allow_pickle=True)
            mi = int(r['_node_idx'])
            if not (0 <= mi < z['nodes'].shape[0]):
                continue
            wt = wt_from_graph(z, mi)
            Xte.append(feats(wt, str(r['mut_info'])[-1]))
            yte.append(float(r['ddg']))
            pidv.append(pid)
        Xte, yte = np.array(Xte), np.array(yte)
        p = m.predict(Xte)
        r = pearson(yte, p)
        lo, hi = cluster_ci(yte, p, pidv)
        print(f'  {bench}: n={len(yte)}  r={r:.4f} [{lo:.3f},{hi:.3f}]  '
              f'alpha={m.alpha_:.3g}')
        rows.append({'benchmark': bench, 'model': 'ridge_physchem', 'n': len(yte),
                     'r': r, 'ci_low': lo, 'ci_high': hi,
                     'alpha': float(m.alpha_)})
        for i in range(len(yte)):
            rows.append({'benchmark': bench, 'model': 'ridge_physchem', 'n': len(yte),
                         'r': np.nan, 'ci_low': np.nan, 'ci_high': np.nan,
                         'alpha': np.nan})
    pd.DataFrame(rows[:len(rows)]).to_csv(os.path.join(OUT, 'ridge_fixed_results.csv'), index=False)
    print(f'\n已保存 ridge_fixed_results.csv')


if __name__ == '__main__':
    main()
