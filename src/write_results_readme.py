# -*- coding: utf-8 -*-
"""write_results_readme.py — 为 release/results/ 生成文件说明，并修正顶层 README"""
import os
from datetime import datetime

import pandas as pd

ROOT = r'D:\GED_mutation'
REL = os.path.join(ROOT, 'release')
RES = os.path.join(REL, 'results')

# 按用途分组的说明
GROUPS = {
    '基准与样本': [
        ('benchmarks_s669_clean.csv', 'S669 基准（543 条），含 _node_idx 与通过标记'),
        ('benchmarks_ssym_clean.csv', 'ssym 基准（342 条）'),
        ('benchmarks_ssym_excluded.csv', '记录 ssym 无任何排除（单行说明）'),
        ('data_flow_skip_log.csv', '逐条排除记录，含是否进入最终集合的布尔列'),
        ('training_merged_noleak_sc.csv', '去泄漏后的训练集（MegaScale + ThermoMutDB）'),
    ],
    '接触变化（论文核心）': [
        ('edits_corrected.csv',
         '★ 505 对 × 4 定义 × 5 阈值 = 10,100 行；论文所有接触变化率的唯一来源'),
        ('edits_corrected_summary.csv', '上表的聚合摘要'),
        ('true_edits_s669_sc.csv',
         '⚠ 命名有歧义：本文件是 side-chain centroid 定义（mean broken 2.11）'),
        ('true_edits_s669.csv',
         '⚠ 命名有歧义：本文件其实是 Cβ 定义（mean broken 0.097）。见 prepare_figure_inputs.py'),
        ('directional_type_analysis.csv', '逐突变的分类接触计数（505 行），论文 §4.5 的来源'),
        ('locality_corrected.csv', '距离统计，含自身接触的旧口径（4.03 Å）'),
        ('locality_nonself.csv',
         '排除自身接触后的口径（7.93 Å），论文 §3.3 使用的版本'),
    ],
    '能力阶梯': [
        ('ladder_results.csv', '逐种子结果，seeds 42/123（EGNN 为 508,934 配置）'),
        ('ladder_s2024_results.csv',
         '⚠ seed 2024 的 508,938 参数对照 run（坐标缩放变体），不进入审计汇总'),
        ('ladder_egnn_legacy_s2024_results.csv',
         '★ architecture-matched 的 seed 2024 EGNN（508,934），审计汇总使用这个'),
        ('ladder_seed_summary_audited.csv', '★ 逐编码器 × 定义的三种子均值与 SD'),
        ('ladder_paired_effects_audited.csv', '★ 13 项配对比较，含联合种子+蛋白 bootstrap'),
        ('ladder_paired_effects_audited_holm.csv', '★ 上表加 Holm 校正列'),
        ('ladder_capacity_trend_audited.csv', '共有定义对的趋势检验（ρ=−0.80, exact P=0.333）'),
        ('ladder_ssym_results.csv', 'ssym 上的阶梯（两个种子）'),
        ('ladder_predictions.csv', '逐样本预测，seeds 42/123'),
        ('ladder_s2024_predictions.csv', '逐样本预测，seed 2024（含 508,938 EGNN）'),
        ('ladder_egnn_legacy_s2024_predictions.csv', '逐样本预测，seed 2024（508,934 EGNN）'),
        ('ladder_paired_effects.csv',
         '⚠ 已被取代的两种子配对效应，不要用于论文数字'),
    ],
    '序列基线与融合': [
        ('esm2_fusion_results.csv', '★ 结构增量：ESM-only / Graph-only / ESM+Graph'),
        ('esm2_fusion_predictions.csv', '逐样本预测，用于复算增量区间'),
        ('esm2_zeroshot_s669_esm2_650m.csv', 'ESM-2 零样本打分（512 条）'),
        ('esm2_zeroshot_ssym_esm2_650m.csv', 'ssym 上的零样本打分'),
        ('ridge_fixed_results.csv', '五种理化特征的岭回归基线'),
        ('increment_decomposition.csv', '特征块分解（P/S/E 及其组合）'),
    ],
    '支撑分析': [
        ('foldx_reproducibility_v2.csv', '三次 FoldX 运行的一致率（96.2% 断裂 / 65.7% 形成）'),
        ('edits_scwrl4_local8A.csv', 'SCWRL4 对照（507 对）'),
        ('scwrl4_cb_displacement.csv', 'SCWRL4 的 Cβ 位移（11 例，Cα 恒为 0）'),
        ('threshold_sensitivity_prediction.csv', '阈值敏感性的预测值'),
        ('hydrogen_bias_audit.csv', '氢修复前后对比（538 对）'),
        ('equivalence_analysis.csv', '可排除效应上界汇总'),
        ('symmetry_test.csv', '方向对称性检验'),
        ('paired_comparison.csv', '模型间配对比较'),
    ],
}


def main():
    files = sorted(f for f in os.listdir(RES) if f.endswith('.csv'))
    described = {n for v in GROUPS.values() for n, _ in v}

    lines = [
        '# Released result tables',
        '',
        f'{len(files)} CSV files. Generated {datetime.now():%Y-%m-%d}.',
        '',
        'These are the derived tables behind every number and figure in the manuscript.',
        'The repository does **not** ship the raw inputs (PDB structures, FoldX mutant',
        'structures, contact-graph `.npz` files, ESM embedding caches) because they total',
        'roughly 21 GB; the scripts that rebuild them from public sources are in `src/`.',
        '',
        '## How paths resolve',
        '',
        '`src/paths.py` sets `DATA` to `data/` if that directory exists and otherwise to',
        '`results/`, so scripts run unchanged in both the full working copy and this release.',
        'Set `GED_DATA` to override.',
        '',
        '## File guide',
        '',
        'Files marked ★ are the ones the manuscript numbers come from. Files marked ⚠ are',
        'superseded or have a name that does not match their content — do not use them for',
        'reported values.',
        '',
    ]
    for group, items in GROUPS.items():
        lines += [f'### {group}', '', '| File | Contents |', '|---|---|']
        for name, desc in items:
            mark = 'yes' if name in files else 'MISSING'
            lines.append(f'| `{name}` | {desc} |' if mark == 'yes'
                         else f'| `{name}` | {desc} — **not present** |')
        lines.append('')

    other = [f for f in files if f not in described]
    if other:
        lines += ['### Other tables', '', '| File |', '|---|']
        lines += [f'| `{f}` |' for f in other]
        lines.append('')

    out = os.path.join(RES, 'README.md')
    # 同时写到项目根，便于 build_release 复制到发布包（release/results 会被重建清空）
    open(os.path.join(ROOT, 'results_README.md'), 'w', encoding='utf-8').write('\n'.join(lines))
    open(out, 'w', encoding='utf-8').write('\n'.join(lines))
    print(f'已写 {out}')
    print(f'  列出 {len(files)} 个 CSV，其中 {len([f for f in files if f in described])} 个有说明')
    print(f'  未说明 {len(other)} 个')
    missing = sorted(described - set(files))
    if missing:
        print(f'  ⚠ 说明里提到但文件不存在: {missing}')


if __name__ == '__main__':
    main()
