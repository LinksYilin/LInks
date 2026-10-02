# -*- coding: utf-8 -*-
"""foldx_bottleneck_stratified.py — 检验结构特征的零关联是否只是 FoldX 误差的产物

审稿人 R2-M1（Blocking）：所有结构量都来自 FoldX，而 FoldX 与实验标签仅
r = 0.26。因此"结构增量不显著"既可能是表示层无信息，也可能是结构流无信息。

R2 自己给出的替代方案：按 FoldX 误差分层，看特征在预测良好的层内是否仍为零关联。

此处执行三项检验：
 (1) 接触特征与 **实验** ΔΔG 的相关 vs 与 **FoldX 预测** ΔΔG 的相关
     —— 若前者显著弱于后者，说明特征部分在追踪 FoldX 的误差
 (2) 按 FoldX 绝对误差分层（预测良好的下/上半），比较关联强度
 (3) 在 FoldX 预测良好的子集内，检查破碎疏水接触是否恢复关联
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths  # noqa: E402

R = str(paths.DATA)
E = pd.read_csv(os.path.join(R, 'directional_type_analysis.csv'))
F = pd.read_csv(os.path.join(R, 'foldx_parsed_energy.csv'))

# 键对齐：mut_info vs mutation_id
E = E.rename(columns={'mut_info': 'mutation_id'})
M = E.merge(F[['pdb_id', 'mutation_id', 'foldx_ddg', 'experimental_ddg']],
            on=['pdb_id', 'mutation_id'], how='inner')
print(f'合并后可分析: {len(M)} 对（接触表 {len(E)}，FoldX 表 {len(F)}）')
print()

FEATS = ['broken_hydro', 'n_broken', 'n_formed', 'n_edit',
         'broken_elec', 'broken_other', 'formed_hydro']


def pear(a, b):
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def cluster_ci(x, y, pid, B=2000, seed=0):
    rng = np.random.default_rng(seed)
    prots = np.unique(pid)
    idx = {p: np.where(pid == p)[0] for p in prots}
    bs = np.empty(B)
    for b in range(B):
        pick = rng.choice(prots, size=len(prots), replace=True)
        ii = np.concatenate([idx[p] for p in pick])
        bs[b] = pear(x[ii], y[ii])
    bs = bs[~np.isnan(bs)]
    return float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


print('=' * 88)
print('(1) 特征与实验标签 vs 与 FoldX 预测的相关')
print('=' * 88)
print(f'{"feature":<16}{"r vs experimental":>20}{"r vs FoldX":>14}{"差":>10}')
print('-' * 62)
rows = []
for f in FEATS:
    r_exp = pear(M[f].to_numpy(float), M['experimental_ddg'].to_numpy(float))
    r_fx = pear(M[f].to_numpy(float), M['foldx_ddg'].to_numpy(float))
    rows.append({'feature': f, 'r_experimental': r_exp, 'r_foldx': r_fx,
                 'diff': r_fx - r_exp})
    print(f'{f:<16}{r_exp:>20.4f}{r_fx:>14.4f}{r_fx - r_exp:>+10.4f}')

print()
print('FoldX 预测自身与实验标签的相关: '
      f"{pear(M['foldx_ddg'].to_numpy(float), M['experimental_ddg'].to_numpy(float)):.4f}")
print()

print('=' * 88)
print('(2) 按 FoldX 绝对误差分层（中位数切分）')
print('=' * 88)
M['foldx_err'] = (M['foldx_ddg'] - M['experimental_ddg']).abs()
med = M['foldx_err'].median()
well = M[M['foldx_err'] <= med]
poor = M[M['foldx_err'] > med]
print(f'中位绝对误差 = {med:.3f} kcal/mol')
print(f'预测良好的下 {len(well)} 对（{well.pdb_id.nunique()} 蛋白）'
      f'｜预测较差的上 {len(poor)} 对（{poor.pdb_id.nunique()} 蛋白）')
print()
print(f'{"feature":<16}{"良好层 r":>12}{"较差层 r":>12}')
print('-' * 42)
for f in FEATS:
    a = pear(well[f].to_numpy(float), well['experimental_ddg'].to_numpy(float))
    b = pear(poor[f].to_numpy(float), poor['experimental_ddg'].to_numpy(float))
    print(f'{f:<16}{a:>12.4f}{b:>12.4f}')

print()
print('=' * 88)
print('(3) 破碎疏水接触：全样本 vs FoldX 预测良好子集（蛋白簇 CI）')
print('=' * 88)
for label, sub in [('全样本', M), ('FoldX 预测良好（下半）', well),
                   ('FoldX 预测较差（上半）', poor)]:
    x = sub['broken_hydro'].to_numpy(float)
    y = sub['experimental_ddg'].to_numpy(float)
    r = pear(x, y)
    lo, hi = cluster_ci(x, y, sub['pdb_id'].to_numpy())
    excl = '排除零' if (lo > 0 or hi < 0) else '含零'
    print(f'{label:<24} n={len(sub):>4}  r={r:+.4f}  '
          f'CI [{lo:+.4f}, {hi:+.4f}]  {excl}')

print()
out = pd.DataFrame(rows)
out.to_csv(os.path.join(R, 'foldx_bottleneck_stratified.csv'), index=False)
print(f'已保存 {os.path.join(R, "foldx_bottleneck_stratified.csv")}')
