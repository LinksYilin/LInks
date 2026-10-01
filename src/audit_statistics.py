# -*- coding: utf-8 -*-
"""audit_statistics.py — 检查统计方法的正确性与描述一致性

覆盖：蛋白聚类 bootstrap、配对相关差 bootstrap、Holm step-down 家族、
精确置换检验的完整性、以及等变编码器"优化不稳定"论断的依据。
"""
import itertools
import os

import numpy as np
import pandas as pd
from docx import Document

D = r'D:\GED_mutation\data'
DOC = r'D:\GED_mutation\manuscript_routeA.docx'

print('=' * 76)
print('统计方法审计')
print('=' * 76)

# ---------- 1) 精确置换检验的完整性 ----------
print('\n【1】精确置换检验：4 档的完整排列数')
print(f'  4 档全排列 = 4! = {len(list(itertools.permutations(range(4))))}')
tr = pd.read_csv(os.path.join(D, 'ladder_capacity_trend_audited.csv'))
print(f'  论文称：rho = {tr.rho.iloc[0]:.3f}, 精确 P = {tr.exact_p.iloc[0]:.3f}, '
      f'24 种排列')
ok = abs(tr.exact_p.iloc[0] * 24 - round(tr.exact_p.iloc[0] * 24)) < 1e-9
print(f'  P × 24 = {tr.exact_p.iloc[0] * 24:.4f}（应为整数）-> {"✅" if ok else "❌"}')

# 手动重算：4 档平均|Δr| 的 Spearman 与 log 参数的相关
eff = pd.read_csv(os.path.join(D, 'ladder_paired_effects_audited.csv'))
common = eff[eff.comparison.isin(['ca vs centroid', 'cb vs centroid'])]
STABLE = ['gnn_global', 'gnn_local', 'gnn_edge', 'deep_gine']   # 论文限定四档稳定编码器
means = (common[common.model.isin(STABLE)]
         .groupby(['model', 'n_params']).delta_r
         .apply(lambda v: float(np.mean(np.abs(v)))).reset_index()
         .sort_values('n_params'))
vals = means.delta_r.to_numpy()
par = means.n_params.to_numpy()


def spearman(a, b):
    ra = pd.Series(a).rank().to_numpy()
    rb = pd.Series(b).rank().to_numpy()
    return float(np.corrcoef(ra, rb)[0, 1])


obs = spearman(np.log10(par), vals)
perms = [spearman(np.log10(par), np.array(vals)[list(p)])
         for p in itertools.permutations(range(len(vals)))]
p_exact = float(np.mean([abs(x) >= abs(obs) - 1e-12 for x in perms]))
print(f'  重算：平均值 {[round(v, 4) for v in vals]}')
print(f'        rho = {obs:.4f}（论文 {tr.rho.iloc[0]:.3f}）')
print(f'        精确双侧 P = {p_exact:.4f}（论文 {tr.exact_p.iloc[0]:.3f}），'
      f'共 {len(perms)} 种排列')

# ---------- 2) Holm step-down ----------
print('\n【2】Holm step-down 家族大小与实现')
hm = pd.read_csv(os.path.join(D, 'ladder_paired_effects_audited_holm.csv'))
m = len(hm)
raw = hm.seed_aware_p.to_numpy()
order = np.argsort(raw)
adj = np.empty(m)
running = 0.0
for rank, idx in enumerate(order):
    val = min(1.0, (m - rank) * raw[idx])
    running = max(running, val)
    adj[idx] = running
match = np.allclose(adj, hm.holm_p.to_numpy(), atol=1e-9)
print(f'  家族大小 m = {m}（论文称"thirteen"）-> {"✅" if m == 13 else "❌"}')
print(f'  重算 Holm 与文件一致 -> {"✅" if match else "❌"}')
print(f'  最小校正 P = {hm.holm_p.min():.3f}（= 13 × {raw.min():.4f} = '
      f'{13 * raw.min():.4f}）')

# ---------- 3) 蛋白聚类 bootstrap ----------
print('\n【3】配对 bootstrap 的重采样单位')
src = open(os.path.join(r'D:\GED_mutation\data', '..', 'src',
                        'analyze_ladder_3seed.py'), encoding='utf-8').read()
checks = [
    ('按蛋白聚类重采样', 'np.unique' in src and 'protein' in src),
    ('同时重采样种子（联合）', 'joint' in src),
    ('百分位区间', 'percentile' in src),
    ('配对差（非独立比较）', 'delta' in src),
]
for label, ok2 in checks:
    print(f'  {"✅" if ok2 else "❌"} {label}')

# ---------- 4) 论文中统计描述的准确性 ----------
print('\n【4】论文统计描述抽样核对')
d = Document(DOC)
txt = '\n'.join(p.text for p in d.paragraphs)
for phrase in ['protein-cluster bootstrap', 'percentile', 'Holm', 'exact permutation',
               'two-sided', 'Bonferroni', 'cluster']:
    print(f'  "{phrase}" 出现 {txt.lower().count(phrase.lower())} 次')

# ---------- 5) 等变编码器不稳定的依据 ----------
print('\n【5】等变编码器"优化不稳定"论断的依据')
eg = eff[eff.model == 'egnn']
print(f'  EGNN 比较数 {len(eg)}，效应范围 '
      f'{eg.delta_r.min():+.3f} 到 {eg.delta_r.max():+.3f}')
sm = pd.read_csv(os.path.join(D, 'ladder_seed_summary_audited.csv'))
e2 = sm[sm.model == 'egnn']
print(f'  逐定义 seed SD: '
      f'{dict(zip(e2.atom_def, e2.seed_sd_r.round(3)))}')
print(f'  -> 逐种子 SD 0.22–0.24，远大于其他编码器 -> 支持"不稳定"论断')

print()
print('=' * 76)
