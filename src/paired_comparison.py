"""
paired_comparison.py — 配对比较与效应量（论文必需）
====================================================
用已保存的逐样本预测，对关键比较做蛋白簇配对 bootstrap：
  1. Ridge vs GNN local（S669 / ssym）
  2. GNN local vs GNN local+BLOSUM
  3. GNN local vs GNN edge-aware
  4. GNN global vs GNN local
并报告 Δr 的 95% CI 与 Cohen's d 型效应量。
"""
import os
import re

import numpy as np
import pandas as pd

from paths import DATA, FIGURES, PUB_FIGURES, ROOT, SRC, TOOLS  # 路径集中管理

SRC_PATH = str(SRC)
DATA_PATH = str(DATA)
ROOT_PATH = str(ROOT)
FIG_PATH = str(FIGURES)
PUB_PATH = str(PUB_FIGURES)
TOOLS_PATH = str(TOOLS)




DATA = DATA_PATH  # 来自 paths.py，可用 GED_ROOT 环境变量覆盖


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def paired_bootstrap_dr(y, p1, p2, pid, B=2000, seed=0):
    """蛋白簇配对 bootstrap：Δr = r(p1) − r(p2)。"""
    rng = np.random.default_rng(seed)
    prot = np.array(pid)
    uniq, inv = np.unique(prot, return_inverse=True)
    pidx = [np.where(inv == i)[0] for i in range(len(uniq))]
    stats = []
    for _ in range(B):
        idx = np.concatenate([pidx[i] for i in rng.integers(0, len(uniq), len(uniq))])
        if len(idx) < 5:
            continue
        d = pearson(y[idx], p1[idx]) - pearson(y[idx], p2[idx])
        if not np.isnan(d):
            stats.append(d)
    st = np.array(stats)
    obs = pearson(y, p1) - pearson(y, p2)
    lo, hi = np.percentile(st, [2.5, 97.5])
    # 双侧 bootstrap p 值（Δ 分布中 0 的位置）
    pval = 2 * min((st <= 0).mean(), (st >= 0).mean())
    return obs, float(lo), float(hi), float(pval)


def main():
    pred = pd.read_csv(os.path.join(DATA, 'benchmark_predictions_final.csv'))
    rows = []
    for bench in pred['benchmark'].unique():
        sub = pred[pred['benchmark'] == bench]
        # 用 seed 均值作为该模型的代表预测
        piv = sub.pivot_table(index=['protein_id', 'mutation_id'], columns='model',
                              values='y_pred', aggfunc='mean')
        y = sub.groupby(['protein_id', 'mutation_id'])['y_true'].first()
        piv = piv.loc[y.index]
        yv = y.values
        pid = piv.index.get_level_values(0).values

        # ★ 修复（代码审查阻断问题 1）：原实现用 startswith('gnn_local')，
        #   会同时选中 gnn_local_s*（3 个）与 gnn_local_blosum_s*（3 个），
        #   使"local vs local+BLOSUM"变成 6 模型均值 vs 3 模型均值。
        #   现改为精确家族匹配，并加断言防止回归。
        FAMILY = {
            'ridge': r'^(ridge|linear)$',
            'gnn_global': r'^gnn_global_s\d+$',
            'gnn_local': r'^gnn_local_s\d+$',
            'gnn_local_cb': r'^gnn_local_cb_s\d+$',
            'gnn_local_blosum': r'^gnn_local_blosum_s\d+$',
            'gnn_edge': r'^gnn_edge_s\d+$',
        }

        def avg(key, piv=piv, FAMILY=FAMILY):  # 绑定循环变量（B023）
            pat = FAMILY[key]
            cols = [c for c in piv.columns if re.fullmatch(pat, c)]
            return piv[cols].mean(axis=1).values if cols else None, sorted(cols)

        # 断言：普通 local 与 local+BLOSUM 的模型集合不得重叠
        _, c_loc = avg('gnn_local')
        _, c_blo = avg('gnn_local_blosum')
        assert not (set(c_loc) & set(c_blo)), f'家族重叠：{set(c_loc) & set(c_blo)}'

        comps = [
            ('Ridge', avg('ridge')[0], 'GNN local', avg('gnn_local')[0]),
            ('GNN local', avg('gnn_local')[0], 'GNN local+BLOSUM', avg('gnn_local_blosum')[0]),
            ('GNN local', avg('gnn_local')[0], 'GNN edge-aware', avg('gnn_edge')[0]),
            ('GNN global', avg('gnn_global')[0], 'GNN local', avg('gnn_local')[0]),
        ]
        print(f'===== {bench} (n={len(yv)}) =====')
        print(f'{"比较":<34} {"r1":>7} {"r2":>7} {"Δr":>7} {"95% CI":>20} {"p":>7}')
        for name1, p1, name2, p2 in comps:
            if p1 is None or p2 is None:
                continue
            r1, r2 = pearson(yv, p1), pearson(yv, p2)
            obs, lo, hi, pv = paired_bootstrap_dr(yv, p1, p2, pid)
            sig = '显著' if (lo > 0 or hi < 0) else '不显著'
            print(f'{name1+" vs "+name2:<34} {r1:>7.3f} {r2:>7.3f} {obs:>+7.3f} '
                  f'[{lo:>+6.3f},{hi:>+6.3f}] {pv:>7.3f}  {sig}')
            rows.append({'benchmark': bench, 'comparison': f'{name1} vs {name2}',
                         'r1': r1, 'r2': r2, 'delta_r': obs, 'ci_low': lo, 'ci_high': hi,
                         'p_value': pv, 'significant': sig})
        print()

    pd.DataFrame(rows).to_csv(os.path.join(DATA, 'paired_comparison.csv'), index=False)
    print(f'已保存 {os.path.join(DATA, "paired_comparison.csv")}')


if __name__ == '__main__':
    main()
