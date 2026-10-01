# -*- coding: utf-8 -*-
"""
routeA_figures.py — 路线 A 新增图：Figure 3（能力阶梯）+ Figure 4（信息预算）
==============================================================================
Figure 3: 表示效应 vs 模型能力（回答"是不是模型太弱"）
Figure 4: 信息预算（P / S / E / 融合 的 r，含增量瀑布）

数据来源：
  data/ladder_results.csv            骨架 × 定义 × 种子
  data/ladder_paired_effects.csv     骨架内表示效应配对检验
  data/increment_decomposition.csv   特征组合的样本外 r
  data/esm2_fusion_results.csv       序列基线 + 融合
  data/ridge_fixed_results.csv       理化岭回归
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

OUT = str(PUB_FIGURES)
os.makedirs(OUT, exist_ok=True)

COL = {'ca': '#0072B2', 'cb': '#D55E00', 'centroid': '#009E73', 'allatom': '#CC79A7'}
LAB = {'ca': 'Cα', 'cb': 'Cβ', 'centroid': 'side-chain\ncentroid', 'allatom': 'all-atom'}
ORDER = ['gnn_global', 'gnn_local', 'gnn_edge', 'deep_gine', 'egnn']
NICE = {'gnn_global': 'GCN\nglobal', 'gnn_local': 'GCN\nlocal', 'gnn_edge': 'GINE\nedge',
        'deep_gine': 'deep GINE\nattention', 'egnn': 'EGNN'}

plt.rcParams.update({'font.size': 7, 'axes.linewidth': 0.6, 'font.family': 'DejaVu Sans'})


def fig3():
    res = pd.read_csv(os.path.join(str(DATA), 'ladder_seed_summary_audited.csv'))
    eff = pd.read_csv(os.path.join(str(DATA), 'ladder_paired_effects_audited.csv'))
    # Chart seed-mean r (not a mixture of seed and graph-definition ranges).
    assert res.model.nunique() == 5 and eff.shape[0] == 13

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9), constrained_layout=True,
                             gridspec_kw={'width_ratios': [1.15, 1]})

    # ---- 面板 a：各骨架 × 定义 的 r ----
    ax = axes[0]
    defs = [d for d in ['ca', 'cb', 'centroid', 'allatom'] if d in set(res['atom_def'])]
    piv = res.pivot_table(index='model', columns='atom_def', values='seed_mean_r', aggfunc='mean')
    npar = res.groupby('model')['n_params'].first()
    models = [m for m in ORDER if m in piv.index]
    xs = np.arange(len(models))
    for k, d in enumerate(defs):
        if d not in piv.columns:
            continue
        off = (k - (len(defs) - 1) / 2) * 0.16
        ax.scatter(xs + off, piv[d].reindex(models), s=26, marker='o',
                   color=COL[d], label=LAB[d].replace('\n', ' '), zorder=3,
                   edgecolor='white', linewidth=0.4)
    ax.set_xticks(xs, [f'{NICE[m]}\n({npar[m]/1000:.1f}K)' for m in models], fontsize=5.6)
    ax.set_ylabel('S669 Pearson $r$')
    ax.set_title('a  Model performance by graph definition', loc='left', fontsize=7)
    ax.grid(axis='y', color='#EEEEEE', lw=0.5, zorder=0)
    ax.legend(fontsize=5.2, frameon=False, ncol=4, loc='upper center',
              bbox_to_anchor=(0.5, -0.16), columnspacing=0.8, handletextpad=0.3)
    ax.set_ylim(0, 0.46)

    # ---- 面板 b：表示效应 vs 参数量（仅用所有稳定档共有的定义对）----
    ax = axes[1]
    if len(eff):
        COMMON = ['ca vs centroid', 'cb vs centroid']
        common = eff[eff['comparison'].isin(COMMON)]
        s = common.groupby(['model', 'n_params'])['delta_r'].apply(
            lambda v: float(np.mean(np.abs(v)))).reset_index()
        s = s.sort_values('n_params')
        stable = s[s['model'] != 'egnn']
        unstable = s[s['model'] == 'egnn']
        ax.plot(stable['n_params'], stable['delta_r'] * 1000, 'o-', color='#333333',
                lw=1.0, ms=4.5, zorder=3, label='stable encoders')
        if len(unstable):
            ax.plot(unstable['n_params'], unstable['delta_r'] * 1000, 'o', mfc='white',
                    mec='#CC3311', mew=1.2, ms=6, zorder=4, label='EGNN (unstable)')
            ax.plot([stable['n_params'].iloc[-1], unstable['n_params'].iloc[0]],
                    [stable['delta_r'].iloc[-1] * 1000, unstable['delta_r'].iloc[0] * 1000],
                    ls=':', color='#CC3311', lw=0.9, zorder=2)
            ax.annotate('EGNN: seed-to-seed\ncollapse, not a\nrepresentation effect',
                        (float(unstable['n_params'].iloc[0]), float(unstable['delta_r'].iloc[0]) * 1000),
                        fontsize=5.0, color='#CC3311', ha='right', va='top',
                        xytext=(-6, -6), textcoords='offset points')
        # 标注未校正时若落在共有定义集内的显著比较
        sig = common[(common['fixed_lo'] > 0) | (common['fixed_hi'] < 0)]
        for _, r in sig.iterrows():
            xs_ = float(s[s['model'] == r['model']]['n_params'].iloc[0])
            ys_ = float(s[s['model'] == r['model']]['delta_r'].iloc[0]) * 1000
            ax.annotate('1/13 uncorrected;\n0/13 after Holm', (xs_, ys_),
                        fontsize=4.8, color='#0072B2', ha='left', va='center',
                        xytext=(14, -14), textcoords='offset points',
                        arrowprops=dict(arrowstyle='-', lw=0.6, color='#0072B2'))
        seen = []
        for _, r in stable.iterrows():
            xy = (float(r['n_params']), float(r['delta_r']) * 1000)
            dy = 6
            for (px, py) in seen:
                if abs(np.log10(xy[0]) - np.log10(px)) < 0.06 and abs(xy[1] - py) < 12:
                    dy = -12
            seen.append(xy)
            ax.annotate(NICE.get(r['model'], r['model']).replace('\n', ' '),
                        xy, fontsize=5.0, xytext=(4, dy), textcoords='offset points')
        ax.axhline(0, color='#999999', lw=0.6, ls='--')
        ax.set_xscale('log')
        ax.set_xlabel('encoder parameters')
        ax.set_ylabel(r'mean $|\Delta r|$ (C$\alpha$, C$\beta$ vs centroid, $\times 10^{-3}$)')
        ax.set_title('b  Representation effect vs capacity', loc='left', fontsize=7)
        ax.grid(color='#EEEEEE', lw=0.5, zorder=0)
        ax.legend(fontsize=5.0, frameon=False, loc='lower right')

    for ext in ['png', 'pdf']:
        fig.savefig(os.path.join(OUT, f'figure3_capacity_ladder.{ext}'), dpi=400)
    plt.close(fig)
    print('  figure3_capacity_ladder.png/pdf')


def fig4():
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9), constrained_layout=True)

    # ---- 面板 a：信息源 r ----
    ax = axes[0]
    rows = []
    fd = os.path.join(str(DATA), 'increment_decomposition.csv')
    if os.path.exists(fd):
        d = pd.read_csv(fd)
        for _, r in d[d['benchmark'] == 's669'].iterrows():
            rows.append((r['features'], r['r'], r['ci_low'], r['ci_high']))
    if rows:
        # 特征名映射为英文（图中无法渲染中文）
        EN = {'P（仅理化）': 'physicochemical (P)', 'S（仅结构）': 'structural (S)',
              'E（仅 ESM-2）': 'ESM-2 (E)', 'P+S': 'P + S', 'P+E': 'P + E',
              'S+E': 'S + E', 'P+S+E': 'P + S + E'}
        names = [EN.get(r[0], r[0]) for r in rows]
        rs = [r[1] for r in rows]
        lo = [r[2] for r in rows]
        hi = [r[3] for r in rows]
        ys = np.arange(len(rows))[::-1]
        ax.errorbar(rs, ys, xerr=[np.array(rs) - np.array(lo), np.array(hi) - np.array(rs)],
                    fmt='o', ms=4, color='#0072B2', ecolor='#888888', elinewidth=0.8,
                    capsize=2, zorder=3)
        ax.set_yticks(ys, names, fontsize=5.6)
        ax.axvline(0, color='#CCCCCC', lw=0.6)
        ax.set_xlabel('S669 Pearson $r$ (out-of-sample)')
        ax.set_title('a  Information sources', loc='left', fontsize=7)
        ax.grid(axis='x', color='#EEEEEE', lw=0.5, zorder=0)

    # ---- 面板 b：序列基线 + 融合 ----
    ax = axes[1]
    fus = os.path.join(str(DATA), 'esm2_fusion_results.csv')
    if os.path.exists(fus):
        f = pd.read_csv(fus)
        f = f[~f['model'].astype(str).str.startswith('increment')]
        benches = ['s669', 'ssym']
        models = ['ESM-only', 'Graph-only', 'ESM+Graph']
        xs = np.arange(len(benches))
        mk = {'ESM-only': 's', 'Graph-only': '^', 'ESM+Graph': 'D'}
        cl = {'ESM-only': '#0072B2', 'Graph-only': '#D55E00', 'ESM+Graph': '#009E73'}
        for m in models:
            v = [float(f[(f['benchmark'] == b) & (f['model'] == m)]['r'].iloc[0])
                 if len(f[(f['benchmark'] == b) & (f['model'] == m)]) else np.nan
                 for b in benches]
            ax.plot(xs, v, marker=mk[m], color=cl[m], lw=1.0, ms=5, label=m, zorder=3)
        ax.set_xticks(xs, ['S669\n(n=508)', 'ssym\n(n=342)'], fontsize=6)
        ax.set_ylabel('Pearson $r$')
        ax.set_title('b  Structure on top of a sequence baseline', loc='left', fontsize=7)
        ax.grid(axis='y', color='#EEEEEE', lw=0.5, zorder=0)
        ax.legend(fontsize=5.6, frameon=False, loc='upper left')
        # 增量标注（放在可见范围内，靠近对应基准的右侧）
        inc = pd.read_csv(fus)
        inc = inc[inc['model'].astype(str).str.startswith('increment')]
        ylo, yhi = ax.get_ylim()
        for k, b in enumerate(benches):
            r0 = inc[inc['benchmark'] == b]
            if len(r0):
                dv = float(r0['r'].iloc[0])
                lo_v = float(r0['ci_low'].iloc[0]); hi_v = float(r0['ci_high'].iloc[0])
                ax.text(k, ylo + 0.03 * (yhi - ylo),
                        f'$\\Delta r$ = {dv:+.3f}\n[{lo_v:+.3f}, {hi_v:+.3f}]',
                        fontsize=5.0, ha='center', va='bottom', color='#555555')

    for ext in ['png', 'pdf']:
        fig.savefig(os.path.join(OUT, f'figure4_information_budget.{ext}'), dpi=400)
    plt.close(fig)
    print('  figure4_information_budget.png/pdf')


if __name__ == '__main__':
    print('生成路线 A 新增图...')
    fig3()
    fig4()
    print(f'输出目录: {OUT}')
