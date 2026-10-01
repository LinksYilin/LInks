# -*- coding: utf-8 -*-
"""
analyze_ladder.py — 表示效应 vs 模型能力（WS1 核心产出）
==========================================================
回答核心问题：
  "换图定义带来的 Δr，是否随模型能力增强而变大？"
  若整条曲线贴近 0 → 有力证明是"表示无关"，而非"模型太弱"。

产出：
  1. 每个 骨架 × 定义 的 r（种子均值 ± 范围）
  2. 每个 骨架 内部的表示效应：两两定义的配对 Δr + 蛋白簇 bootstrap CI
  3. 表示效应 vs 参数量 的汇总表（供作图）
  4. 交互项检验：表示效应是否依赖骨架（置换检验）
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from paths import DATA

OUT = str(DATA)
DEFS = ['ca', 'cb', 'centroid', 'allatom']
ORDER = ['gnn_global', 'gnn_local', 'gnn_edge', 'deep_gine', 'egnn']


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def clusters(pid):
    prot = np.array(pid)
    uniq, inv = np.unique(prot, return_inverse=True)
    return [np.where(inv == i)[0] for i in range(len(uniq))], len(uniq)


def paired_dr(y, p1, p2, pidx, n_uniq, B=2000, seed=0):
    rng = np.random.default_rng(seed)
    st = []
    for _ in range(B):
        idx = np.concatenate([pidx[i] for i in rng.integers(0, n_uniq, n_uniq)])
        if len(idx) < 5:
            continue
        d = pearson(y[idx], p1[idx]) - pearson(y[idx], p2[idx])
        if not np.isnan(d):
            st.append(d)
    st = np.array(st)
    obs = pearson(y, p1) - pearson(y, p2)
    lo, hi = np.percentile(st, [2.5, 97.5])
    p = 2 * min((st <= 0).mean(), (st >= 0).mean())
    # 等价性检验：可排除的效应上界
    return obs, float(lo), float(hi), float(min(p, 1.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tag', default='ladder')
    ap.add_argument('--bench', default='s669')
    args = ap.parse_args()

    res = pd.read_csv(os.path.join(OUT, f'{args.tag}_results.csv'))
    pred = pd.read_csv(os.path.join(OUT, f'{args.tag}_predictions.csv'))
    res = res[res['benchmark'] == args.bench]
    pred = pred[pred['benchmark'] == args.bench]
    models = [m for m in ORDER if m in set(res['model'])]
    print(f'骨架 {len(models)} 个, 定义 {sorted(set(res["atom_def"]))}, 基准 {args.bench}')

    # ---- 1) r 表 ----
    piv = res.pivot_table(index='model', columns='atom_def', values='r', aggfunc='mean')
    piv = piv.reindex(index=models, columns=[d for d in DEFS if d in piv.columns])
    npar = res.groupby('model')['n_params'].first()
    print('\n=== 各骨架 × 定义 的 r（种子均值）===')
    out = piv.copy()
    out.insert(0, 'n_params', npar.reindex(models))
    print(out.round(4).to_string())

    # ---- 2) 每个骨架内部的表示效应 ----
    rows = []
    y = None
    pidx = n_uniq = None
    for model in models:
        sub = pred[pred['model'] == model]
        piv_p = sub.pivot_table(index=['protein_id', 'mutation_id'],
                                columns=['atom_def', 'seed'], values='y_pred', aggfunc='mean')
        # 每个定义取种子均值
        defs_avail = sorted({c[0] for c in piv_p.columns})
        def arr(d):
            cols = [c for c in piv_p.columns if c[0] == d]
            return piv_p[cols].mean(axis=1).values
        yy = sub.groupby(['protein_id', 'mutation_id'])['y_true'].first().loc[piv_p.index].values
        pidv = piv_p.index.get_level_values(0).values
        pidx, n_uniq = clusters(pidv)
        y = yy
        if 'centroid' not in defs_avail:
            continue
        pc = arr('centroid')
        for d in defs_avail:
            if d == 'centroid':
                continue
            pa = arr(d)
            obs, lo, hi, pv = paired_dr(yy, pa, pc, pidx, n_uniq)
            rows.append({'model': model, 'n_params': int(npar[model]),
                         'comparison': f'{d} vs centroid',
                         'r_alt': pearson(yy, pa), 'r_centroid': pearson(yy, pc),
                         'delta_r': obs, 'ci_low': lo, 'ci_high': hi, 'p': pv,
                         'excludes_zero': (lo > 0 or hi < 0),
                         'max_excludable_effect': max(abs(lo), abs(hi))})
    eff = pd.DataFrame(rows)
    print('\n=== 表示效应（各定义 vs 侧链质心，骨架内配对）===')
    if len(eff):
        print(eff[['model', 'n_params', 'comparison', 'delta_r', 'ci_low', 'ci_high',
                   'p', 'excludes_zero']].round(4).to_string(index=False))

    # ---- 3) 汇总：表示效应幅度 vs 能力 ----
    if len(eff):
        summ = eff.groupby(['model', 'n_params']).agg(
            mean_abs_delta_r=('delta_r', lambda s: float(np.mean(np.abs(s)))),
            max_abs_delta_r=('delta_r', lambda s: float(np.max(np.abs(s)))),
            n_comparisons=('delta_r', 'size'),
            n_significant=('excludes_zero', 'sum'),
        ).reset_index().sort_values('n_params')
        print('\n=== 表示效应幅度 vs 模型能力 ===')
        print(summ.round(4).to_string(index=False))
        summ.to_csv(os.path.join(OUT, f'{args.tag}_representation_effect.csv'), index=False)
    eff.to_csv(os.path.join(OUT, f'{args.tag}_paired_effects.csv'), index=False)
    out.to_csv(os.path.join(OUT, f'{args.tag}_summary.csv'))

    # ---- 4) 结论提示 ----
    print('\n=== 判读 ===')
    if len(eff):
        any_sig = eff['excludes_zero'].any()
        mx = eff['delta_r'].abs().max()
        print(f'  表示效应绝对值最大 {mx:.3f}')
        print(f'  是否有任何（骨架, 定义对）达到显著: {"是" if any_sig else "否"}')
        if not any_sig:
            bound = eff['max_excludable_effect'].max()
            print(f'  → 在所有测试的骨架下，均可排除 |Δr| > {bound:.3f} 的表示效应')
            print('  → 支持"表示无关"，而非"模型太弱"')
        else:
            sig = eff[eff['excludes_zero']]
            print(f'  → 有 {len(sig)} 个比较达到显著，需如实报告：')
            print(sig[['model', 'comparison', 'delta_r', 'ci_low', 'ci_high']].round(3).to_string(index=False))


if __name__ == '__main__':
    main()
