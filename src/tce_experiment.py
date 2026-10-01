# -*- coding: utf-8 -*-
"""
tce_experiment.py — Typed Contact-Edit descriptors (TCE) 增量检验
====================================================================
目的：判断能否把"分类型的接触编辑"提出为一个有增量价值的方法。

背景（实测）：
  - 不分类型的接触计数与 ΔΔG 几乎无关（r ≈ 0.06–0.08，CI 含 0）
  - 但 "断裂的疏水接触数" r = +0.135 [+0.044, +0.201]（显著）

设计（严格训练/测试分离）：
  A  P                      5 个理化特征（地板）
  B  P + 未分类型编辑        broken/formed/total
  C  P + TCE-6              断/成 × 疏水/静电/其他
  D  P + TCE-6 + 距离分层   × 2 个距离壳层（≤5 Å, >5 Å）
  E  ESM-2 + TCE-6          在强序列基线之上是否还有增量

评估：S669（探索性）与 ssym（独立），配对 bootstrap 给出增量 CI。
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from paths import DATA

OUT = str(DATA)
CACHE = os.path.join(OUT, 'esm_emb_cache')
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


def paired_dr(y, p1, p2, pid, B=2000, seed=0):
    rng = np.random.default_rng(seed)
    uniq, inv = np.unique(np.array(pid), return_inverse=True)
    pidx = [np.where(inv == i)[0] for i in range(len(uniq))]
    st = []
    for _ in range(B):
        idx = np.concatenate([pidx[i] for i in rng.integers(0, len(uniq), len(uniq))])
        if len(idx) < 5:
            continue
        d = pearson(y[idx], p1[idx]) - pearson(y[idx], p2[idx])
        if not np.isnan(d):
            st.append(d)
    st = np.array(st)
    obs = pearson(y, p1) - pearson(y, p2)
    lo, hi = np.percentile(st, [2.5, 97.5])
    p = 2 * min((st <= 0).mean(), (st >= 0).mean())
    return obs, float(lo), float(hi), float(min(p, 1.0)), st


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


def phys(wt, mt):
    return [KD.get(mt, 0) - KD.get(wt, 0), VOL.get(mt, 0) - VOL.get(wt, 0),
            KD.get(mt, 0), VOL.get(mt, 0), CHG.get(mt, 0) - CHG.get(wt, 0)]


def main():
    dt = pd.read_csv(os.path.join(OUT, 'directional_type_analysis.csv'))
    print(f'TCE 数据 {len(dt)} 行')

    # 距离分层需要额外计算：用现有 locality 数据或重算
    # 这里先用现有 6 个类型特征 + 3 个未分类型
    dt['n_typed'] = (dt['broken_hydro'] + dt['broken_elec'] + dt['broken_other']
                     + dt['formed_hydro'] + dt['formed_elec'] + dt['formed_other'])

    # 训练集：需要同样特征 → 用 Fixed 训练集重建（这里用训练集突变 + 重新计算成本高）
    # 改为：用 S669 内部做 5-fold 蛋白簇交叉验证，避免与训练集特征口径不一致
    from sklearn.linear_model import RidgeCV
    from sklearn.preprocessing import StandardScaler

    # ESM 嵌入（可选）
    emb_file = os.path.join(CACHE, 'emb_esm2_650m.npz')
    E_map = None
    if os.path.exists(emb_file):
        z = np.load(emb_file, allow_pickle=True)
        E_all = z['emb']
        E_map = {k: i for i, k in enumerate(list(z['key']))}

    for bench in ['s669', 'ssym']:
        sub = dt.copy()  # directional_type_analysis 只做了 S669
        if bench != 's669':
            print(f'\n[{bench}] directional_type_analysis.csv 仅含 S669，跳过')
            continue

        # 理化特征：从 benchmarks 取 wt/mt
        bm = pd.read_csv(os.path.join(OUT, 'benchmarks_s669_clean.csv'))
        if '_node_idx_ok' in bm.columns:
            bm = bm[bm['_node_idx_ok']]
        key2idx = {f'{r["pdb_id"]}|{r["mut_info"]}': i for i, r in bm.iterrows()}

        P, T6, Tun, E, Y, PID = [], [], [], [], [], []
        for _, r in sub.iterrows():
            k = f'{r["pdb_id"]}|{r["mut_info"]}'
            if k not in key2idx:
                continue
            row = bm.loc[key2idx[k]]
            wt, mt = str(row['mut_info'])[0], str(row['mut_info'])[-1]
            P.append(phys(wt, mt))
            Tun.append([r['n_broken'], r['n_formed'], r['n_edit']])
            T6.append([r['broken_hydro'], r['broken_elec'], r['broken_other'],
                       r['formed_hydro'], r['formed_elec'], r['formed_other']])
            if E_map is not None:
                ek = f's669:{r["pdb_id"]}:{r["mut_info"]}'
                E.append(E_all[E_map[ek]] if ek in E_map else np.zeros(E_all.shape[1]))
            Y.append(float(r['ddg']))
            PID.append(r['pdb_id'])
        P, T6, Tun = np.array(P), np.array(T6), np.array(Tun)
        Y = np.array(Y)
        PID = np.array(PID)
        E = np.array(E) if E and len(E) == len(Y) else None
        print(f'\n=== {bench}: n={len(Y)}  P={P.shape} T6={T6.shape}'
              + (f' E={E.shape}' if E is not None else '') + ' ===')

        # 蛋白簇 5-fold CV（保证训练/测试分离）
        uniq = np.unique(PID)
        rng = np.random.default_rng(42)
        fold = {p: i % 5 for i, p in enumerate(rng.permutation(uniq))}
        fold_id = np.array([fold[p] for p in PID])

        def cv_r(X):
            pred = np.zeros(len(Y))
            for f in range(5):
                trm, tem = fold_id != f, fold_id == f
                if trm.sum() < 20 or tem.sum() < 5:
                    continue
                sc = StandardScaler().fit(X[trm])
                m = RidgeCV(alphas=np.logspace(-2, 3, 20)).fit(sc.transform(X[trm]), Y[trm])
                pred[tem] = m.predict(sc.transform(X[tem]))
            return pred

        models = {
            'A  P（理化）': P,
            'B  P + 未分类型编辑': np.concatenate([P, Tun], axis=1),
            'C  P + TCE-6（分类型）': np.concatenate([P, T6], axis=1),
        }
        if E is not None:
            models['E  ESM-2'] = E
            models['F  ESM-2 + TCE-6'] = np.concatenate([E, T6], axis=1)

        preds = {}
        for name, X in models.items():
            p = cv_r(X)
            r = pearson(Y, p)
            lo, hi = cluster_ci(Y, p, PID)
            preds[name] = p
            print(f'  {name:<26} 维度 {X.shape[1]:>5}  5折CV r={r:.4f} [{lo:.3f},{hi:.3f}]')

        print('\n  --- 增量（配对 bootstrap，重复簇保留）---')
        for a, b in [('C  P + TCE-6（分类型）', 'A  P（理化）'),
                     ('B  P + 未分类型编辑', 'A  P（理化）'),
                     ('C  P + TCE-6（分类型）', 'B  P + 未分类型编辑')]:
            if a in preds and b in preds:
                obs, lo, hi, pv, _ = paired_dr(Y, preds[a], preds[b], PID)
                sig = '显著' if (lo > 0 or hi < 0) else '不显著'
                print(f'    {a}  −  {b}')
                print(f'        Δr={obs:+.4f} [{lo:+.4f},{hi:+.4f}]  p={pv:.3f}  {sig}')
        if 'F  ESM-2 + TCE-6' in preds and 'E  ESM-2' in preds:
            obs, lo, hi, pv, _ = paired_dr(Y, preds['F  ESM-2 + TCE-6'], preds['E  ESM-2'], PID)
            sig = '显著' if (lo > 0 or hi < 0) else '不显著'
            print(f'    F  ESM-2+TCE  −  E  ESM-2')
            print(f'        Δr={obs:+.4f} [{lo:+.4f},{hi:+.4f}]  p={pv:.3f}  {sig}')


if __name__ == '__main__':
    main()
