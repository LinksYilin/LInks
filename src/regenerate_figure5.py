# -*- coding: utf-8 -*-
"""
regenerate_figure5.py — 重新生成 Figure 5（旧 panel a 用了失效的 BLOSUM 数据）
================================================================================
新面板：
  a) 质心接触编辑数 vs 实验 ΔΔG（n = 505，r = 0.081）
  b) 分类型接触编辑与 ΔΔG 的相关（新分析，含蛋白簇 bootstrap CI）
  c) FoldX 预测 vs 实验 ΔΔG（n = 538，r = 0.26）
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from paths import DATA, PUB_FIGURES

D = str(DATA)
OUT = str(PUB_FIGURES)
plt.rcParams.update({'font.size': 7, 'axes.linewidth': 0.6, 'font.family': 'DejaVu Sans'})


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else float('nan')


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


def main():
    edits = pd.read_csv(os.path.join(D, 'edits_corrected.csv'))
    e8 = edits[(edits['threshold'] == 8.0) & (edits['atom_def'] == 'centroid')].copy()
    dt = pd.read_csv(os.path.join(D, 'directional_type_analysis.csv'))
    fx = pd.read_csv(os.path.join(D, 'foldx_parsed_energy.csv'))

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6), constrained_layout=True,
                             gridspec_kw={'width_ratios': [1, 1.05, 1]})

    # ---- a) 编辑数 vs ΔΔG ----
    ax = axes[0]
    x, y = e8['n_edit'].values.astype(float), e8['ddg'].values.astype(float)
    ax.scatter(x, y, s=5, alpha=0.45, color='#00876C', edgecolor='none')
    b, a = np.polyfit(x, y, 1)
    xs = np.linspace(0, np.percentile(x, 99), 50)
    ax.plot(xs, a + b * xs, color='#005544', lw=1.0)
    r = pearson(y, x)
    lo, hi = cluster_ci(y, x, e8['pdb_id'].values)
    ax.text(0.97, 0.95, f'r = {r:.3f}\n95% CI {lo:.3f} to {hi:.3f}\nn = {len(y)}',
            transform=ax.transAxes, ha='right', va='top', fontsize=5.6)
    ax.set_xlabel('Total side-chain-centroid contact edits')
    ax.set_ylabel('Experimental \u0394\u0394G (kcal mol\u207b\u00b9)')
    ax.set_title('a  Unaggregated edit counts', loc='left', fontsize=7)

    # ---- b) 分类型接触编辑相关 ----
    ax = axes[1]
    cols = [('broken_hydro', 'broken\nhydro-\nphobic'), ('broken_elec', 'broken\nelectro-\nstatic'),
            ('broken_other', 'broken\nother'), ('formed_hydro', 'formed\nhydro-\nphobic'),
            ('formed_elec', 'formed\nelectro-\nstatic'), ('formed_other', 'formed\nother')]
    yy = dt['ddg'].values
    pid = dt['pdb_id'].values
    rs, los, his, labs = [], [], [], []
    for c, lab in cols:
        v = dt[c].values.astype(float)
        rs.append(pearson(yy, v))
        l, h = cluster_ci(yy, v, pid)
        los.append(l); his.append(h); labs.append(lab)
    xs_ = np.arange(len(cols))
    sig = [i for i, (l, h) in enumerate(zip(los, his)) if l > 0 or h < 0]
    ax.errorbar(xs_, rs, yerr=[np.array(rs) - np.array(los), np.array(his) - np.array(rs)],
                fmt='none', ecolor='#999999', elinewidth=0.8, capsize=2, zorder=2)
    ax.scatter(xs_, rs, s=22, zorder=3,
               color=['#D55E00' if i in sig else '#0072B2' for i in range(len(cols))])
    ax.axhline(0, color='#CCCCCC', lw=0.6)
    ax.set_xticks(xs_, labs, fontsize=4.9)
    ax.set_ylim(-0.15, 0.26)
    ax.set_ylabel('Pearson $r$ with \u0394\u0394G')
    ax.set_title('b  Contact edits by type', loc='left', fontsize=7)
    if sig:
        ax.annotate('only broken\nhydrophobic\ncontacts exclude 0', (sig[0], rs[sig[0]]),
                    fontsize=4.9, color='#D55E00', ha='left', va='bottom',
                    xytext=(6, 4), textcoords='offset points')

    # ---- c) FoldX vs 实验 ----
    ax = axes[2]
    ax.scatter(fx['experimental_ddg'], fx['foldx_ddg'], s=5, alpha=0.45,
               color='#7B5EA7', edgecolor='none')
    b2, a2 = np.polyfit(fx['experimental_ddg'], fx['foldx_ddg'], 1)
    lim = [-10, 11]
    ax.plot(lim, lim, ls='--', color='#AAAAAA', lw=0.8, label='identity')
    ax.plot(np.array(lim), a2 + b2 * np.array(lim), color='#4B2E83', lw=1.0,
            label='linear fit')
    rf = pearson(fx['experimental_ddg'].values, fx['foldx_ddg'].values)
    mae = float(np.mean(np.abs(fx['experimental_ddg'] - fx['foldx_ddg'])))
    ax.text(0.03, 0.97, f'r = {rf:.2f}\nMAE = {mae:.2f} kcal mol\u207b\u00b9\nn = {len(fx)}',
            transform=ax.transAxes, ha='left', va='top', fontsize=5.6)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel('Experimental \u0394\u0394G')
    ax.set_ylabel('FoldX \u0394\u0394G (kcal mol\u207b\u00b9)')
    ax.set_title('c  FoldX energy signal', loc='left', fontsize=7)
    ax.legend(fontsize=5.2, frameon=False, loc='lower right')

    for ext in ['png', 'pdf']:
        fig.savefig(os.path.join(OUT, f'figure5_contact_diagnostics.{ext}'), dpi=400)
    plt.close(fig)
    print('  figure5_contact_diagnostics.png/pdf')
    print(f'  a: r = {r:.3f} [{lo:.3f}, {hi:.3f}] n = {len(y)}')
    print(f'  b: 显著项 = {[cols[i][0] for i in sig]}')
    print(f'  c: r = {rf:.3f}, MAE = {mae:.2f}, n = {len(fx)}')


if __name__ == '__main__':
    main()
