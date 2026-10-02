# -*- coding: utf-8 -*-
"""within_protein_metric.py — 对四种图定义计算蛋白内 Fisher-z 相关

审稿人 R1-M4（Blocking）：论文报告过一次蛋白内 Fisher-z 指标，结果显示
架构排序**反转**（global pooling 蛋白内最高 0.420-0.442，但突变级最低），
却从未把该指标用于四种接触图定义之间的比较。因此无法判断表示效应零结果
是否只是**汇总方式**的产物。

此处对每个编码器 × 定义计算：蛋白内 Pearson（每个蛋白 ≥5 个突变），
用 Fisher z 无权重平均汇总，并按蛋白簇成对 bootstrap 比较定义差异。
"""
import os

import numpy as np
import pandas as pd

R = r'D:\GED_mutation\release\results'
OUT = r'D:\GED_mutation\data'
P = pd.read_csv(os.path.join(R, 'ladder_predictions.csv'))

# 只看 S669（与论文主分析一致）
P = P[P.benchmark == 's669']

DEFS = ['ca', 'cb', 'centroid', 'allatom']
MODELS = ['gnn_global', 'gnn_local', 'gnn_edge', 'deep_gine']
REF = 'centroid'
B = 2000
rng = np.random.default_rng(42)


def within_protein_fisher(df):
    """返回 (汇总 z, 使用的蛋白数)。每蛋白需 >=5 个突变。"""
    zs = []
    for _, g in df.groupby('protein_id'):
        if len(g) < 5:
            continue
        if g.y_true.std() == 0 or g.y_pred.std() == 0:
            continue
        r = float(np.corrcoef(g.y_true, g.y_pred)[0, 1])
        r = max(min(r, 0.9999), -0.9999)
        zs.append(np.arctanh(r))
    if not zs:
        return np.nan, 0
    return float(np.mean(zs)), len(zs)


print('=' * 78)
print('蛋白内 Fisher-z 相关（每蛋白 >=5 个突变，按种子平均预测后计算）')
print('=' * 78)
print(f'{"encoder":<12}{"def":<10}{"n_prot":>7}{"within_z":>11}{"as_r":>9}')
print('-' * 50)

rows = []
for m in MODELS:
    for d in DEFS:
        sub = P[(P.model == m) & (P.atom_def == d)]
        if sub.empty:
            continue
        # 先按蛋白×突变对多种子求平均预测
        avg = (sub.groupby(['protein_id', 'mutation_id'])
               .agg(y_true=('y_true', 'first'), y_pred=('y_pred', 'mean'))
               .reset_index())
        z, n = within_protein_fisher(avg)
        rows.append({'model': m, 'atom_def': d, 'n_prot': n,
                     'within_z': z, 'within_r': float(np.tanh(z)) if not np.isnan(z) else np.nan})
        print(f'{m:<12}{d:<10}{n:>7}{z:>11.4f}{np.tanh(z):>9.4f}')

res = pd.DataFrame(rows)

print()
print('=' * 78)
print(f'与参考定义（{REF}）的配对差异（蛋白簇 bootstrap，B = {B}）')
print('=' * 78)

# 预计算每个 (model, def, protein) 的 z
zmap = {}
for m in MODELS:
    for d in DEFS:
        sub = P[(P.model == m) & (P.atom_def == d)]
        if sub.empty:
            continue
        avg = (sub.groupby(['protein_id', 'mutation_id'])
               .agg(y_true=('y_true', 'first'), y_pred=('y_pred', 'mean'))
               .reset_index())
        per = {}
        for pid, g in avg.groupby('protein_id'):
            if len(g) < 5 or g.y_true.std() == 0 or g.y_pred.std() == 0:
                continue
            r = float(np.corrcoef(g.y_true, g.y_pred)[0, 1])
            per[pid] = np.arctanh(max(min(r, 0.9999), -0.9999))
        zmap[(m, d)] = per

pairs = []
for m in MODELS:
    common_prots = set(zmap.get((m, REF), {}))
    for d in DEFS:
        if d == REF:
            continue
        common = sorted(common_prots & set(zmap.get((m, d), {})))
        if len(common) < 5:
            continue
        diff_obs = float(np.mean([zmap[(m, d)][p] for p in common])
                         - np.mean([zmap[(m, REF)][p] for p in common]))
        bs = np.empty(B)
        for b in range(B):
            pick = rng.choice(common, size=len(common), replace=True)
            bs[b] = float(np.mean([zmap[(m, d)][p] for p in pick])
                          - np.mean([zmap[(m, REF)][p] for p in pick]))
        lo, hi = float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))
        ptail = min(float(np.mean(bs <= 0)), float(np.mean(bs >= 0)))
        pairs.append({'model': m, 'comparison': f'{d} vs {REF}', 'n_prot': len(common),
                      'delta_z': diff_obs, 'ci_lo': lo, 'ci_hi': hi,
                      'p': min(2 * ptail, 1.0),
                      'excludes_zero': bool(lo > 0 or hi < 0)})

pr = pd.DataFrame(pairs)
if not pr.empty:
    print(pr.to_string(index=False, float_format=lambda v: f'{v:.4f}'))
    print()
    print(f'区间排除零: {int(pr["excludes_zero"].sum())} / {len(pr)}')
    print(f'最小未校正 P: {pr["p"].min():.4f}')
    # Holm
    pv = pr['p'].to_numpy()
    m_ = len(pv)
    order = np.argsort(pv)
    holm = np.empty_like(pv)
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, (m_ - rank) * pv[i])
        holm[i] = min(run, 1.0)
    pr['p_holm'] = holm
    print(f'Holm 校正后 P < 0.05: {int((pr["p_holm"] < 0.05).sum())} 项')
else:
    print('（数据不足，无法配对）')

os.makedirs(OUT, exist_ok=True)
res.to_csv(os.path.join(OUT, 'within_protein_by_definition.csv'), index=False)
if not pr.empty:
    pr.to_csv(os.path.join(OUT, 'within_protein_paired.csv'), index=False)
print(f'\n已保存至 {OUT}')
