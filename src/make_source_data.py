# -*- coding: utf-8 -*-
"""make_source_data.py — 生成与论文图号对应的源数据文件（回应 E-3）

现状问题：
  figure1_model_metrics_source_data.csv 支持的是未使用的模型对比图，不是 Figure 1
  figure3_contact_signal_/figure3_foldx_ 装的其实是 Figure 5 的面板数据
  Figure 1（研究设计）与 Figure 3（能力阶梯）没有任何源数据文件
"""
import os
import shutil

import numpy as np
import pandas as pd

D = r'D:\GED_mutation\data'
FIG = r'D:\GED_mutation\figures\publication'

# ---------- Figure 3a: 每个编码器 × 定义 的三种子均值 ----------
sm = pd.read_csv(os.path.join(D, 'ladder_seed_summary_audited.csv'))
sm[['model', 'atom_def', 'n_params', 'n_seeds', 'seed_mean_r', 'seed_sd_r']].to_csv(
    os.path.join(FIG, 'figure3a_capacity_ladder_source_data.csv'), index=False)

# ---------- Figure 3b: 共有定义对的平均 |Δr| ----------
eff = pd.read_csv(os.path.join(D, 'ladder_paired_effects_audited.csv'))
common = eff[eff.comparison.isin(['ca vs centroid', 'cb vs centroid'])]
(common.groupby(['model', 'n_params'])
      .agg(mean_abs_delta_r=('delta_r', lambda v: float(np.mean(np.abs(v)))),
           n_comparisons=('delta_r', 'size'))
      .reset_index()
      .sort_values('n_params')
      .to_csv(os.path.join(FIG, 'figure3b_capacity_trend_source_data.csv'), index=False))

# ---------- Figure 5b: 分类接触计数的相关 ----------
t = pd.read_csv(os.path.join(D, 'directional_type_analysis.csv'))
t.to_csv(os.path.join(FIG, 'figure5b_contact_type_source_data.csv'), index=False)

# ---------- 保留既有文件名但重命名为 Figure 5 对应 ----------
mv = [
    ('figure3_contact_signal_source_data.csv', 'figure5a_contact_edits_source_data.csv'),
    ('figure3_foldx_source_data.csv', 'figure5c_foldx_energy_source_data.csv'),
]
for old, new in mv:
    src, dst = os.path.join(FIG, old), os.path.join(FIG, new)
    if os.path.exists(src) and not os.path.exists(dst):
        shutil.copy2(src, dst)
        os.remove(src)
        print(f'  重命名 {old} -> {new}')

# figure1_model_metrics 支持的是未使用的图，改名以免误导
old1 = os.path.join(FIG, 'figure1_model_metrics_source_data.csv')
new1 = os.path.join(FIG, 'unused_model_comparison_source_data.csv')
if os.path.exists(old1) and not os.path.exists(new1):
    shutil.copy2(old1, new1)
    os.remove(old1)
    print('  重命名 figure1_model_metrics -> unused_model_comparison（该图未在论文中使用）')

# ---------- Figure 1: 研究设计是示意图，无数据；建立说明文件 ----------
with open(os.path.join(FIG, 'figure1_study_design_README.txt'), 'w', encoding='utf-8') as f:
    f.write(
        'Figure 1 is a study-design schematic and has no data source.\n'
        'The quantities annotated on it are: 7,905 training examples from 420 proteins;\n'
        '511 S669 mutations on 88 proteins; 342 ssym mutations on 15 proteins;\n'
        'and an encoder capacity range of 5.8 k to 509 k parameters.\n')
print('  建立 figure1_study_design_README.txt（示意图无数据源）')

# ---------- 映射表 ----------
with open(os.path.join(FIG, 'SOURCE_DATA_MAP.md'), 'w', encoding='utf-8') as f:
    f.write('# Figure to source-data mapping\n\n')
    f.write('| Figure | Panel | Source data file |\n|---|---|---|\n')
    f.write('| Figure 1 | — (schematic) | `figure1_study_design_README.txt` (no data; '
            'annotated values listed there) |\n')
    f.write('| Figure 2 | a, b | `figure2_contact_edits_source_data.csv` |\n')
    f.write('| Figure 2 | c | `figure2_repeatability_source_data.csv` |\n')
    f.write('| Figure 3 | a | `figure3a_capacity_ladder_source_data.csv` |\n')
    f.write('| Figure 3 | b | `figure3b_capacity_trend_source_data.csv` |\n')
    f.write('| Figure 4 | a | `increment_decomposition.csv` (`data/`) |\n')
    f.write('| Figure 4 | b | `esm2_fusion_results.csv` and `esm2_fusion_predictions.csv` '
            '(`data/`) |\n')
    f.write('| Figure 5 | a | `figure5a_contact_edits_source_data.csv` |\n')
    f.write('| Figure 5 | b | `figure5b_contact_type_source_data.csv` |\n')
    f.write('| Figure 5 | c | `figure5c_foldx_energy_source_data.csv` |\n')
    f.write('\n`unused_model_comparison_source_data.csv` supports a panel that is not in '
            'the manuscript.\n')

print()
print('现有源数据文件:')
for f_ in sorted(os.listdir(FIG)):
    if 'source_data' in f_ or 'SOURCE_DATA' in f_ or f_.endswith('README.txt'):
        print(f'  {f_}')
