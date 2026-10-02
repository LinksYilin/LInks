# -*- coding: utf-8 -*-
"""cluster_level_typed_counts.py — 用蛋白簇 bootstrap 重算六类接触计数的 P 值

审稿人 R1-M2（Blocking）：论文用突变级未校正 P = 0.0024 判定破碎疏水接触
"通过多重比较校正"，但全部分析的独立单位是**蛋白簇**（S669 仅 88 个蛋白），
且论文自己在别处都用蛋白簇 bootstrap。

此处按蛋白簇重采样（整簇抽样，保留簇内结构）计算：
  - 蛋白簇 bootstrap 95% 百分位区间
  - 双侧 bootstrap 符号概率（与论文 §3.6 声明的约定一致）
  - Holm 与 Bonferroni 校正
"""
import os

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths
R = str(paths.DATA)
OUT = str(paths.DATA)
D = pd.read_csv(os.path.join(R, 'directional_type_analysis.csv'))

TYPED = ['broken_hydro', 'broken_elec', 'broken_other',
         'formed_hydro', 'formed_elec', 'formed_other']

y = D['ddg'].to_numpy(dtype=float)
pid = D['pdb_id'].to_numpy()
proteins = np.unique(pid)
idx_by_prot = {p: np.where(pid == p)[0] for p in proteins}
B = 5000
rng = np.random.default_rng(42)


def pearson(a, b):
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


print(f'数据: {len(D)} 对突变, {len(proteins)} 个蛋白簇, B = {B}')
print()

rows = []
for c in TYPED:
    x = D[c].to_numpy(dtype=float)
    obs = pearson(x, y)
    bs = np.empty(B)
    for b in range(B):
        pick = rng.choice(proteins, size=len(proteins), replace=True)
        ii = np.concatenate([idx_by_prot[p] for p in pick])
        bs[b] = pearson(x[ii], y[ii])
    bs = bs[~np.isnan(bs)]
    lo, hi = float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))
    # 双侧 bootstrap 符号概率
    ptail = min(float(np.mean(bs <= 0)), float(np.mean(bs >= 0)))
    p_two = min(2 * ptail, 1.0)
    rows.append({'count': c, 'r': obs, 'ci_lo': lo, 'ci_hi': hi,
                 'p_cluster': p_two, 'excludes_zero': bool(lo > 0 or hi < 0)})

res = pd.DataFrame(rows)
p = res['p_cluster'].to_numpy()
m = len(p)
order = np.argsort(p)
holm = np.empty_like(p)
running = 0.0
for rank, i in enumerate(order):
    running = max(running, (m - rank) * p[i])
    holm[i] = min(running, 1.0)
res['p_holm'] = holm
res['p_bonferroni'] = np.minimum(p * m, 1.0)

print(res.to_string(index=False, float_format=lambda v: f'{v:.4f}'))
print()
print(f'区间排除零: {int(res["excludes_zero"].sum())} / {m}')
print(f'Holm 校正后 P < 0.05: {int((res["p_holm"] < 0.05).sum())} 项')
print()

bh = res[res['count'] == 'broken_hydro'].iloc[0]
print('关键对比（破碎疏水接触）:')
print(f'  论文当前引用：突变级未校正 P = 0.0024，据此称"通过多重比较校正"')
print(f'  蛋白簇区间  : [{bh.ci_lo:+.4f}, {bh.ci_hi:+.4f}]  '
      f'({"排除零" if bh.excludes_zero else "含零"})')
print(f'  蛋白簇 P    : {bh.p_cluster:.4f}')
print(f'  Holm 校正 P : {bh.p_holm:.4f}     Bonferroni P: {bh.p_bonferroni:.4f}')
print()

os.makedirs(OUT, exist_ok=True)
res.to_csv(os.path.join(OUT, 'cluster_level_typed_counts.csv'), index=False)
print(f'已保存 {os.path.join(OUT, "cluster_level_typed_counts.csv")}')
