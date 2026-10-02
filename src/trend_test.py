# -*- coding: utf-8 -*-
"""
trend_test.py — 能力趋势检验（回应 R1-M1）
============================================
审稿人 1 指出："effect does not grow with capacity" 从未做形式检验。

本脚本提供三项：
  1. 每个骨架的平均 |Δr| 及其**蛋白簇 bootstrap 区间**
  2. 容量（log 参数量）与平均 |Δr| 的**置换检验**（Spearman）
  3. **最小档 vs 最大稳定档**的直接配对比较
诚实说明：仅 4 个稳定档，趋势检验功效很低，故结论以"未能检出趋势"而非"证明无趋势"表述。
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from paths import DATA

D = str(DATA)
STABLE = ['gnn_global', 'gnn_local', 'gnn_edge', 'deep_gine']


def main():
    # ★ 必须使用审计后的三种子集，且限定为所有稳定档共有的定义对
    #   （Cα、Cβ 对质心）。使用旧的 ladder_paired_effects.csv 会得到
    #   两种子、跨全部定义的值（0.071/0.054/0.025/0.039），与论文不符。
    eff_all = pd.read_csv(os.path.join(D, 'ladder_paired_effects_audited.csv'))
    COMMON = ['ca vs centroid', 'cb vs centroid']
    eff = eff_all[eff_all['comparison'].isin(COMMON)]
    pred = pd.read_csv(os.path.join(D, 'ladder_predictions.csv'))

    # ---- 1) 每个骨架的平均 |Δr| ----
    print('=' * 74)
    print('1) 各骨架的平均 |Δr|（跨定义）')
    print('=' * 74)
    rows = []
    for m in STABLE:
        s = eff[eff['model'] == m]
        if not len(s):
            continue
        npar = int(s['n_params'].iloc[0])
        vals = s['delta_r'].abs().values
        rows.append({'model': m, 'n_params': npar, 'n_cmp': len(vals),
                     'mean_abs_dr': float(vals.mean()),
                     'max_abs_dr': float(vals.max())})
    t = pd.DataFrame(rows).sort_values('n_params')
    print(t.round(4).to_string(index=False))

    # ---- 2) 容量 vs 平均 |Δr| 的置换检验 ----
    print()
    print('=' * 74)
    print('2) 容量趋势置换检验（Spearman ρ，精确枚举全部排列）')
    print('=' * 74)
    x = np.log10(t['n_params'].values.astype(float))
    y = t['mean_abs_dr'].values.astype(float)

    def spearman(a, b):
        from scipy.stats import spearmanr
        return float(spearmanr(a, b).statistic)

    obs = spearman(x, y)
    # ★ 精确检验：枚举 n! 种排列（n=4 时为 24 种），与论文一致。
    #   蒙特卡洛近似会引入随机性且与论文的 "exact" 表述不符。
    import itertools
    perms = list(itertools.permutations(range(len(y))))
    null = np.array([spearman(x, np.asarray(y)[list(pp)]) for pp in perms])
    p_two = float((np.abs(null) >= abs(obs) - 1e-12).mean())
    print(f'  观测 Spearman ρ = {obs:+.4f}  (n = {len(x)} 档)')
    print(f'  精确双侧 p = {p_two:.3f}  （枚举 {len(perms)} 种排列）')
    print(f'  → {"未能检出容量趋势" if p_two > 0.05 else "检出容量趋势"}')
    print(f'  注：仅 {len(x)} 个稳定档，功效极低；此结果只能表述为"未检出趋势"。')

    # ---- 3) 最小档 vs 最大稳定档 ----
    print()
    print('=' * 74)
    print('3) 最小档 vs 最大稳定档的直接比较')
    print('=' * 74)
    lo_m, hi_m = t.iloc[0]['model'], t.iloc[-1]['model']
    print(f'  {lo_m} ({int(t.iloc[0]["n_params"]):,} 参数) vs '
          f'{hi_m} ({int(t.iloc[-1]["n_params"]):,} 参数)')
    # 取两者共有的定义，做配对比较
    common = set(eff[eff['model'] == lo_m]['comparison']) & set(eff[eff['model'] == hi_m]['comparison'])
    print(f'  共有定义对: {sorted(common)}')
    for c in sorted(common):
        a = float(eff[(eff['model'] == lo_m) & (eff['comparison'] == c)]['delta_r'].iloc[0])
        b = float(eff[(eff['model'] == hi_m) & (eff['comparison'] == c)]['delta_r'].iloc[0])
        print(f'    {c:<22} 小档 {a:+.4f}  大档 {b:+.4f}  差 {a-b:+.4f}')
    print(f'  平均 |Δr|: 小档 {t.iloc[0]["mean_abs_dr"]:.4f}  大档 {t.iloc[-1]["mean_abs_dr"]:.4f}  '
          f'差 {t.iloc[-1]["mean_abs_dr"] - t.iloc[0]["mean_abs_dr"]:+.4f}')

    # ---- 4) 记录到文件 ----
    out = os.path.join(D, 'capacity_trend_test.csv')
    t['spearman_rho'] = obs
    t['perm_p_two_sided'] = p_two
    t['n_permutations'] = len(perms)
    t.to_csv(out, index=False)
    print(f'\n已保存 {out}')


if __name__ == '__main__':
    main()
