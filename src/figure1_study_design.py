# -*- coding: utf-8 -*-
"""
figure1_study_design.py — Figure 1：研究设计与评估协议流程图（v2，版面修正）
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = r'D:\GED_mutation\figures\publication'
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.size': 6.2, 'font.family': 'DejaVu Sans'})

STY = dict(boxstyle='round,pad=0.40', linewidth=0.7, edgecolor='#333333')
CLR = {'data': '#DCE9F5', 'proc': '#F5E6DC', 'bench': '#DFF0E0', 'enc': '#EFE3F0'}


def box(ax, x, y, w, h, text, fc, fs=6.0, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, facecolor=fc, **STY))
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=fs,
            fontweight='bold' if bold else 'normal', linespacing=1.4)


def arrow(ax, p1, p2, color='#555555', ls='-'):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle='-|>', mutation_scale=6.5,
                                 linewidth=0.7, linestyle=ls, color=color,
                                 shrinkA=1, shrinkB=1))


def main():
    fig, ax = plt.subplots(figsize=(7.2, 3.9), constrained_layout=True)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 60)
    ax.axis('off')

    # ================= 第 1 行：训练数据 =================
    box(ax, 0.5, 48, 20, 8.5, 'MegaScale\ntrain_s669 split', CLR['data'])
    box(ax, 0.5, 37.5, 20, 8.5, 'ThermoMutDB', CLR['data'])
    box(ax, 25, 42, 20, 10, 'BLAST leakage\nfilter vs S669 / ssym', CLR['proc'])
    arrow(ax, (20.5, 52.2), (25, 48.5))
    arrow(ax, (20.5, 41.7), (25, 45.5))
    box(ax, 49, 42, 21, 10, 'Training set\n7,905 examples\n(1:1 sampling)', CLR['data'], bold=True)
    arrow(ax, (45, 47), (49, 47))

    # ================= 第 2 行：基准 =================
    box(ax, 25, 25, 20, 9, 'S669  (exploratory)\n511 mutations, 88 proteins', CLR['bench'])
    box(ax, 49, 25, 21, 9, 'ssym  (development-\nindependent)\n342 mutations, 15 proteins',
        CLR['bench'])
    # 训练集 → 两个基准（泄漏过滤方向）
    arrow(ax, (74.5, 42), (79, 34), ls=':')
    arrow(ax, (79, 34), (35, 34), ls=':')
    arrow(ax, (35, 34), (35, 34.2), ls=':')

    # ================= 第 3 行：图定义 =================
    box(ax, 0.5, 11, 98, 10,
        'Four residue contact-graph definitions  —  wild-type and FoldX-modelled mutant '
        'structures processed identically\n'
        r'C$\alpha$          C$\beta$          side-chain centroid (H excluded)          all-atom',
        CLR['proc'], fs=6.0)
    arrow(ax, (35, 25), (35, 21))
    arrow(ax, (59.5, 25), (59.5, 21))

    # ================= 第 4 行：模型 =================
    box(ax, 0.5, 0.5, 47, 8,
        'Five encoder capacities, 5.8 k \u2192 509 k parameters\n'
        'mean-pooled GCN | mutation-site GCN | edge-aware GINE |\n'
        'attention-pooled deep GINE | E(3)-equivariant network',
        CLR['enc'], fs=5.5)
    box(ax, 51.5, 0.5, 47, 8,
        'Reference models\n'
        'physicochemical ridge baseline | ESM-2 650M with learned head |\n'
        'sequence\u2013structure fusion',
        CLR['enc'], fs=5.5)
    arrow(ax, (24, 11), (24, 8.5))
    arrow(ax, (75, 11), (75, 8.5))

    # 底部统计说明
    ax.text(50, 58.6, 'All paired comparisons use 95% protein-cluster bootstrap intervals; '
                      'representation effects are paired within each encoder.',
            ha='center', va='center', fontsize=5.5, style='italic', color='#555555')

    for ext in ['png', 'pdf']:
        fig.savefig(os.path.join(OUT, f'figure1_study_design.{ext}'), dpi=400)
    plt.close(fig)
    print('  figure1_study_design.png/pdf')


if __name__ == '__main__':
    main()
