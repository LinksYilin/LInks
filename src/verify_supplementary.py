# -*- coding: utf-8 -*-
"""verify_supplementary.py — 核对补充材料中的数字与结果文件"""
import os
import re

import numpy as np
import pandas as pd

D = r'D:\GED_mutation\data'
SUP = r'D:\GED_mutation\补充材料_Supplementary.md'


def main():
    text = open(SUP, encoding='utf-8').read()
    results = []

    def chk(label, expected, actual, tol=0.0015):
        ok = abs(float(expected) - float(actual)) <= tol
        results.append(ok)
        print(f'  {"✅" if ok else "❌"} {label:<44} 补充材料 {expected:<10} 实测 {actual}')

    # S5：13 项比较
    holm = pd.read_csv(os.path.join(D, 'ladder_paired_effects_audited_holm.csv'))
    chk('S5 比较总数 13', 13, len(holm))
    chk('S5 Holm 后显著数', 0, int((holm.holm_p < 0.05).sum()))
    chk('S5 未校正显著数', 1,
        int((((holm.fixed_lo > 0) | (holm.fixed_hi < 0))).sum()))
    chk('S5 最小 Holm P', 0.585, round(float(holm.holm_p.min()), 3), tol=0.002)
    # 抽查三行
    row = holm[(holm.model == 'gnn_global') & (holm.comparison == 'allatom vs centroid')].iloc[0]
    chk('S5 gnn_global allatom Δr', -0.101, round(float(row.delta_r), 3), tol=0.002)
    chk('S5 gnn_global allatom 未校正 P', 0.045, round(float(row.seed_aware_p), 3), tol=0.002)
    row = holm[(holm.model == 'gnn_local') & (holm.comparison == 'ca vs centroid')].iloc[0]
    chk('S5 gnn_local ca Δr', -0.064, round(float(row.delta_r), 3), tol=0.002)

    # S5：趋势
    tr = pd.read_csv(os.path.join(D, 'ladder_capacity_trend_audited.csv'))
    chk('S5 趋势 rho', -0.800, round(float(tr.rho.iloc[0]), 3), tol=0.002)
    chk('S5 趋势精确 P', 0.333, round(float(tr.exact_p.iloc[0]), 3), tol=0.002)
    vals = sorted(tr.mean_abs_dr.round(3).tolist())
    exp = [0.021, 0.043, 0.048, 0.048]
    ok = all(abs(a - b) <= 0.002 for a, b in zip(vals, exp))
    results.append(ok)
    print(f'  {"✅" if ok else "❌"} {"S5 四档平均|Δr|":<44} 补充材料 {exp} 实测 {vals}')

    # S6：阈值
    th = pd.read_csv(os.path.join(D, 'threshold_sensitivity_prediction.csv'))
    for bm, exp in [('S669', [0.410, 0.386, 0.382, 0.353, 0.353]),
                    ('ssym', [0.381, 0.371, 0.410, 0.359, 0.348])]:
        s = th[th.benchmark == bm].sort_values('threshold')['r'].tolist()
        for i, (v, e) in enumerate(zip(s, exp)):
            chk(f'S6 {bm} {[6,7,8,9,10][i]}Å', f'{e:.3f}', f'{v:.3f}')
    # 边界数
    e = pd.read_csv(os.path.join(D, 'edits_corrected.csv'))
    for t, exp in [(6.0, 140), (7.0, 237), (8.0, 338), (9.0, 465), (10.0, 615)]:
        sub = e[(e.threshold == t) & (e.atom_def == 'centroid')]
        chk(f'S6 {t:g}Å 平均接触数', exp, round(float(sub.wt_edges.mean()), 0), tol=1.5)

    # S7：EGNN
    egnn = pd.read_csv(os.path.join(D, 'ladder_egnn_legacy_s2024_results.csv'))
    lad = pd.read_csv(os.path.join(D, 'ladder_results.csv'))
    seeds = lad[(lad.model == 'egnn') & (lad.atom_def == 'centroid')]['r'].tolist()
    seeds += egnn[egnn.atom_def == 'centroid']['r'].tolist()
    chk('S7 三种子均值', 0.131, round(float(np.mean(seeds)), 3), tol=0.002)
    chk('S7 样本 SD', 0.241, round(float(np.std(seeds, ddof=1)), 3), tol=0.002)
    chk('S7 Cα 三种子最小值', -0.053,
        round(float(min(lad[(lad.model == 'egnn') & (lad.atom_def == 'ca')]['r'].tolist() +
                         egnn[egnn.atom_def == 'ca']['r'].tolist())), 3), tol=0.002)

    # S2：样本量
    for val in ['543', '512', '511', '508', '505', '342', '792', '905']:
        pass
    chk('S2 训练蛋白数 420', 420,
        len(set(pd.read_csv(os.path.join(D, 'training_merged_noleak_sc.csv')).protein)))

    print()
    print('=' * 74)
    n_ok = sum(results)
    print(f'补充材料核验：{n_ok}/{len(results)} 通过')
    print('=' * 74)
    return 0 if n_ok == len(results) else 1


if __name__ == '__main__':
    raise SystemExit(main())
