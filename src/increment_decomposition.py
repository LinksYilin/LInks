# -*- coding: utf-8 -*-
"""
increment_decomposition.py — 增量信息分解（WS2d，修正版）
============================================================
★ 修正：原版在评估集上拟合（样本内），导致 ESM-2 得到虚高的 r=0.91。
        本版严格用训练集拟合、在 S669/ssym 上评估。

信息源：
  P) 理化特征（5 个）
  S) 结构接触特征（突变位点度数、归一化度数、最近邻距离、平均邻距）
  E) ESM-2 序列嵌入（冻结，已缓存）

对 7 个组合分别拟合岭回归（在训练集上），在测试基准上报告 r 与
蛋白簇 bootstrap CI，并给出每个信息源的**边际增量**。
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from paths import DATA

OUT = str(DATA)
CACHE = os.path.join(OUT, 'esm_emb_cache')
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


def phys_feats(wt, mt):
    return [KD.get(mt, 0) - KD.get(wt, 0), VOL.get(mt, 0) - VOL.get(wt, 0),
            KD.get(mt, 0), VOL.get(mt, 0), CHG.get(mt, 0) - CHG.get(wt, 0)]


def struct_feats(z, mi):
    """从接触图取突变位点结构特征（与图定义无关的通用量）。"""
    ei = z['edge_index']
    n = z['nodes'].shape[0]
    mask = (ei[0] == mi)
    deg = int(mask.sum())
    if 'coords' in z.files and deg > 0:
        c = z['coords']
        d = np.linalg.norm(c[ei[1][mask]] - c[mi], axis=1)
        dmin, dmean, dmax = float(d.min()), float(d.mean()), float(d.max())
    else:
        dmin = dmean = dmax = np.nan
    return [deg, deg / max(n, 1), dmin, dmean, dmax]


def wt_from_graph(z, mi):
    return AA20[int(np.argmax(z['nodes'][mi, :20]))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', default='esm2_650m')
    ap.add_argument('--defs', default='centroid')
    args = ap.parse_args()

    # ---- 嵌入缓存 ----
    z = np.load(os.path.join(CACHE, f'emb_{args.model}.npz'), allow_pickle=True)
    E_all = z['emb']
    key2i = {k: i for i, k in enumerate(list(z['key']))}
    print(f'嵌入库 {E_all.shape}')

    # ---- 训练特征 ----
    from ladder_common import graph_dir

    tr = pd.read_csv(os.path.join(OUT, 'training_merged_noleak_sc.csv'))
    ms, tm = tr[tr['source'] == 'megascale'], tr[tr['source'] == 'thermomutdb']
    n = min(len(ms), len(tm))
    tr = pd.concat([ms.sample(n=n, random_state=42), tm], ignore_index=True)
    gcache = {}
    for src, gd in [('megascale', graph_dir('megascale', args.defs)),
                    ('thermomutdb', graph_dir('thermomutdb', args.defs))]:
        gcache[src] = gd

    Ptr, Str, Etr, ytr = [], [], [], []
    for _, r in tr.iterrows():
        src, prot = r['source'], str(r['protein'])
        p = os.path.join(gcache[src], f'{prot}.npz')
        k = f'{src}:{prot}:{int(r["mut_idx"])}:{r["mt_aa"]}'
        if not os.path.exists(p) or k not in key2i:
            continue
        zz = np.load(p, allow_pickle=True)
        mi = int(r['mut_idx'])
        if not (0 <= mi < zz['nodes'].shape[0]):
            continue
        wt = wt_from_graph(zz, mi)
        Ptr.append(phys_feats(wt, str(r['mt_aa'])))
        Str.append(struct_feats(zz, mi))
        Etr.append(E_all[key2i[k]])
        ytr.append(float(r['ddg']))
    Ptr, Str, Etr = np.array(Ptr), np.array(Str), np.array(Etr)
    ytr = np.array(ytr)
    print(f'训练: n={len(ytr)}, P={Ptr.shape}, S={Str.shape}, E={Etr.shape}')

    from sklearn.linear_model import RidgeCV
    from sklearn.preprocessing import StandardScaler

    combos = [('P（仅理化）', ['P']), ('S（仅结构）', ['S']), ('E（仅 ESM-2）', ['E']),
              ('P+S', ['P', 'S']), ('P+E', ['P', 'E']), ('S+E', ['S', 'E']),
              ('P+S+E', ['P', 'S', 'E'])]
    blocks = {'P': Ptr, 'S': Str, 'E': Etr}

    rows = []
    for bench in ['s669', 'ssym']:
        label = f'benchmarks_{bench}_clean.csv'
        df = pd.read_csv(os.path.join(OUT, label))
        if '_node_idx_ok' in df.columns:
            df = df[df['_node_idx_ok']].reset_index(drop=True)
        gd = graph_dir(bench, args.defs)

        Pte, Ste, Ete, yte, pidte = [], [], [], [], []
        for _, r in df.iterrows():
            pid = r['pdb_id']
            p = os.path.join(gd, f'{pid}.npz')
            k = f'{bench}:{pid}:{r["mut_info"]}'
            if not os.path.exists(p) or k not in key2i:
                continue
            zz = np.load(p, allow_pickle=True)
            mi = int(r['_node_idx'])
            if not (0 <= mi < zz['nodes'].shape[0]):
                continue
            wt = wt_from_graph(zz, mi)
            Pte.append(phys_feats(wt, str(r['mut_info'])[-1]))
            Ste.append(struct_feats(zz, mi))
            Ete.append(E_all[key2i[k]])
            yte.append(float(r['ddg']))
            pidte.append(pid)
        Pte, Ste, Ete = np.array(Pte), np.array(Ste), np.array(Ete)
        yte = np.array(yte)
        blocks_te = {'P': Pte, 'S': Ste, 'E': Ete}
        print(f'\n=== {bench}: n={len(yte)} ===')

        for name, ks in combos:
            Xtr = np.concatenate([blocks[k] for k in ks], axis=1)
            Xte = np.concatenate([blocks_te[k] for k in ks], axis=1)
            # 用训练集中位数插补缺失（结构特征在无接触时为 NaN）
            med = np.nanmedian(Xtr, axis=0)
            med = np.where(np.isnan(med), 0.0, med)
            Xtr = np.where(np.isnan(Xtr), med, Xtr)
            Xte = np.where(np.isnan(Xte), med, Xte)
            sc = StandardScaler().fit(Xtr)
            m = RidgeCV(alphas=np.logspace(-2, 3, 20)).fit(sc.transform(Xtr), ytr)
            p = m.predict(sc.transform(Xte))
            r = pearson(yte, p)
            lo, hi = cluster_ci(yte, p, pidte)
            rows.append({'benchmark': bench, 'features': name, 'n_features': Xtr.shape[1],
                         'n_test': len(yte), 'r': r, 'ci_low': lo, 'ci_high': hi})
            print(f'  {name:<16} 维度 {Xtr.shape[1]:>5}  r={r:.4f} [{lo:.3f},{hi:.3f}]')

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT, 'increment_decomposition.csv'), index=False)

    print('\n=== 边际增量（S669）===')
    s = out[out['benchmark'] == 's669'].set_index('features')['r']

    def g(k):
        return float(s.get(k, np.nan))
    print(f'  结构在理化之上  (P+S − P):     {g("P+S") - g("P（仅理化）"):+.4f}')
    print(f'  ESM 在理化之上  (P+E − P):     {g("P+E") - g("P（仅理化）"):+.4f}')
    print(f'  结构在 ESM 之上 (S+E − E):     {g("S+E") - g("E（仅 ESM-2）"):+.4f}  ← 关键')
    print(f'  结构在 P+E 之上 (P+S+E − P+E): {g("P+S+E") - g("P+E"):+.4f}  ← 关键')
    print(f'  全模型 (P+S+E):                {g("P+S+E"):.4f}')
    print('\n已保存 increment_decomposition.csv')


if __name__ == '__main__':
    main()
