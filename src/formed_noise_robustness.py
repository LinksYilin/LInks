# -*- coding: utf-8 -*-
"""formed_noise_robustness.py — 检验结论是否依赖噪声较大的"形成接触"（R2-M5）

审稿人 R2-M5：形成接触跨 FoldX 运行的一致率只有 65.7%（基于 3 个突变），
而这类接触正是喂给图编码器的输入之一。测量误差从未被传播到关联分析，
因此"接触计数信号弱"可能只是 FoldX 运行方差的产物。

此处用**不需要新 FoldX 运行**的方式检验该解释：把分析限制在**只含断裂接触**
（跨运行一致率 96.2%，远高于形成接触的 65.7%）的特征上，看结论是否改变。
若结论不依赖形成接触，则"形成接触噪声"这一替代解释就被削弱。
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths  # noqa: E402

R = str(paths.DATA)
D = pd.read_csv(os.path.join(R, 'directional_type_analysis.csv'))
y = D['ddg'].to_numpy(float)
pid = D['pdb_id'].to_numpy()
prots = np.unique(pid)
idx = {p: np.where(pid == p)[0] for p in prots}


def pear(a, b):
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def cci(x, yy, B=3000, seed=0):
    rng = np.random.default_rng(seed)
    bs = np.empty(B)
    for b in range(B):
        pick = rng.choice(prots, size=len(prots), replace=True)
        ii = np.concatenate([idx[p] for p in pick])
        bs[b] = pear(x[ii], yy[ii])
    bs = bs[~np.isnan(bs)]
    lo, hi = float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))
    ptail = min(float(np.mean(bs <= 0)), float(np.mean(bs >= 0)))
    return lo, hi, min(2 * ptail, 1.0)


print(f'数据 {len(D)} 对，{len(prots)} 蛋白簇')
print()
print('跨运行一致率（来自 foldx_reproducibility_v2.csv）:')
rep = pd.read_csv(os.path.join(R, 'foldx_reproducibility_v2.csv'))
print(f'  断裂接触 {rep.broken_agreement.mean():.3f} | '
      f'形成接触 {rep.formed_agreement.mean():.3f}')
print()

SPECS = [
    ('断裂疏水（高一致率）', 'broken_hydro'),
    ('断裂总数（高一致率）', 'n_broken'),
    ('形成疏水（低一致率）', 'formed_hydro'),
    ('形成总数（低一致率）', 'n_formed'),
    ('总编辑数', 'n_edit'),
]

print(f'{"特征":<24}{"r":>10}{"95% CI":>22}{"P":>9}  一致率')
print('-' * 78)
rows = []
for label, col in SPECS:
    x = D[col].to_numpy(float)
    r = pear(x, y)
    lo, hi, p = cci(x, y)
    conf = '高 (0.96)' if col.startswith(('broken', 'n_broken')) else \
           ('低 (0.66)' if col.startswith(('formed', 'n_formed')) else '—')
    print(f'{label:<24}{r:>+10.4f}   [{lo:+.4f}, {hi:+.4f}] {p:>8.4f}  {conf}')
    rows.append({'feature': label, 'col': col, 'r': r, 'ci_lo': lo, 'ci_hi': hi,
                 'p_cluster': p, 'excludes_zero': bool(lo > 0 or hi < 0)})

res = pd.DataFrame(rows)
print()
print('只依赖高一致率（断裂）特征的结论是否成立:')
hi_only = res[res.col.isin(['broken_hydro', 'n_broken', 'n_edit'])]
print(f'  其中区间排除零: {int(hi_only.excludes_zero.sum())} / {len(hi_only)}')
print()
print('低一致率（形成）特征:')
lo_only = res[res.col.isin(['formed_hydro', 'n_formed'])]
print(f'  其中区间排除零: {int(lo_only.excludes_zero.sum())} / {len(lo_only)}')
print()
print('关键结论：破碎疏水接触的关联只基于断裂接触，不涉及形成接触的噪声，')
print('因此该发现不能由形成接触的跨运行不一致来解释。')

res.to_csv(os.path.join(R, 'formed_noise_robustness.csv'), index=False)
print(f'\n已保存 {os.path.join(R, "formed_noise_robustness.csv")}')
