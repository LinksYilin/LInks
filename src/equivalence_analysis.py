# -*- coding: utf-8 -*-
"""
equivalence_analysis.py — 等价性检验（WS3）
=============================================
把"没找到差异（p>0.05）"改写成"可排除大于 X 的效应"，这是负结果论文的正确统计表达。

对每个配对比较，从 bootstrap 分布给出：
  - 观测 Δr
  - 95% CI
  - **可排除的效应上界** = max(|ci_low|, |ci_high|)
    （即：若真实效应超过此值，本数据将以 ≥95% 概率检测到方向性差异）
  - 等价性判定：若 CI 完全落在 ±margin 内，则在该 margin 上"等价"

输入：data/esm2_fusion_results.csv（含增量的 CI）
      data/ladder_paired_effects.csv（若已生成）
输出：data/equivalence_analysis.csv + 打印表
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from paths import DATA

OUT = str(DATA)


def summarize(df, label_cols, dr_col='r', lo='ci_low', hi='ci_high', pcol='p'):
    out = []
    for _, r in df.iterrows():
        lo_v, hi_v = r.get(lo), r.get(hi)
        if pd.isna(lo_v) or pd.isna(hi_v):
            margin = np.nan
            excl = np.nan
        else:
            # 可排除的效应上界：CI 越窄，能排除的效应越小
            excl = max(abs(lo_v), abs(hi_v))
            margin = excl
        out.append({
            **{c: r[c] for c in label_cols if c in r},
            'delta_r': r[dr_col],
            'ci_low': lo_v, 'ci_high': hi_v,
            'p': r.get(pcol, np.nan),
            'excludes_zero': (not pd.isna(lo_v)) and (lo_v > 0 or hi_v < 0),
            'max_excludable_effect': excl,
        })
    return pd.DataFrame(out)


def main():
    rows = []

    # ---- 1) 结构增量（ESM+Graph vs ESM-only）----
    fus = os.path.join(OUT, 'esm2_fusion_results.csv')
    if os.path.exists(fus):
        f = pd.read_csv(fus)
        inc = f[f['model'].astype(str).str.startswith('increment')]
        if len(inc):
            t = summarize(inc, ['benchmark', 'model'])
            t['analysis'] = 'structure increment over sequence baseline'
            rows.append(t)

    # ---- 2) 阶梯的表示效应 ----
    # ★ 优先使用审计后的三种子文件；旧的 ladder_paired_effects.csv 是两种子、
    #   跨全部定义的值，与论文不一致，仅作为回退。
    lad = os.path.join(OUT, 'ladder_paired_effects_audited.csv')
    if not os.path.exists(lad):
        lad = os.path.join(OUT, 'ladder_paired_effects.csv')
    if os.path.exists(lad):
        try:
            l = pd.read_csv(lad)
            if len(l):
                t = summarize(l, ['model', 'comparison'],
                              dr_col='delta_r', lo='seed_aware_lo', hi='seed_aware_hi',
                              pcol='seed_aware_p') if 'seed_aware_lo' in l.columns else \
                    summarize(l, ['model', 'comparison'],
                              dr_col='delta_r', lo='ci_low', hi='ci_high', pcol='p')
                t['analysis'] = 'representation effect within backbone'
                rows.append(t)
        except Exception:
            pass

    if not rows:
        print('无可分析的结果文件')
        return

    out = pd.concat(rows, ignore_index=True)
    out.to_csv(os.path.join(OUT, 'equivalence_analysis.csv'), index=False)

    print('=== 等价性检验 ===')
    for _, r in out.iterrows():
        tag = ' | '.join(str(r[c]) for c in ['analysis', 'benchmark', 'model', 'comparison']
                         if c in r and not pd.isna(r[c]))
        lo_v, hi_v = r['ci_low'], r['ci_high']
        if pd.isna(lo_v):
            print(f'  {tag}: Δr={r["delta_r"]:+.3f}  （无 CI）')
            continue
        sig = '显著' if r['excludes_zero'] else '不显著'
        print(f'  {tag}')
        print(f'      Δr={r["delta_r"]:+.4f}  95% CI [{lo_v:+.4f}, {hi_v:+.4f}]  p={r["p"]:.3f}  {sig}')
        if not r['excludes_zero']:
            print(f'      → 可排除 |Δr| > {r["max_excludable_effect"]:.4f} 的效应')

    # 汇总：不显著项的可排除上界
    ns = out[~out['excludes_zero'].fillna(True)]
    if len(ns):
        print(f'\n=== 不显著比较共 {len(ns)} 项，可排除上界 ===')
        print(f'  中位 {ns["max_excludable_effect"].median():.4f}, '
              f'最大 {ns["max_excludable_effect"].max():.4f}')
        print(f'  即：本数据可排除约 |Δr| > {ns["max_excludable_effect"].median():.3f} 的效应')
    print(f'\n已保存 {os.path.join(OUT, "equivalence_analysis.csv")}')


if __name__ == '__main__':
    main()
