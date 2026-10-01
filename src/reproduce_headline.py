# -*- coding: utf-8 -*-
"""reproduce_headline.py — 仅用发布包 results/ 的数据复算论文头条数字

验证论文的可用性声明：代码与结果表足以复现报告的数字。
只读取 release/results/*.csv，不使用论文或项目中间产物。
"""
import os

import numpy as np
import pandas as pd

REL = r'D:\GED_mutation\release\results'
RESULTS = []


def rd(name):
    p = os.path.join(REL, name)
    return pd.read_csv(p) if os.path.exists(p) else None


def chk(label, paper, got, tol=0.0015):
    try:
        ok = abs(float(paper) - float(got)) <= tol
    except (TypeError, ValueError):
        ok = str(paper) == str(got)
    RESULTS.append(ok)
    print(f'  {"OK " if ok else "FAIL"} {label:<42} paper {paper:<11} recomputed {got}')
    return ok


def main():
    print('=' * 76)
    print('Recomputing the manuscript headline numbers from released data only')
    print('=' * 76)

    # 1) contact-change rates
    ed = rd('edits_corrected.csv')
    if ed is not None:
        print('\n[contact-change rates, 505 pairs at 8 A]')
        s = ed[ed.threshold == 8.0]
        for ad, exp in [('ca', 0.0), ('cb', 9.9), ('centroid', 90.7), ('allatom', 81.6)]:
            sub = s[s.atom_def == ad]
            chk(f'{ad} rate %', f'{exp:.1f}', f'{100 * (sub.n_edit > 0).mean():.1f}', tol=0.06)
        chk('centroid pair count', 505, len(s[s.atom_def == 'centroid']))

    # 2) FoldX reproducibility
    ro = rd('foldx_reproducibility_v2.csv')
    if ro is not None:
        print('\n[FoldX reproducibility, three runs]')
        chk('broken agreement %', '96.2', f'{ro.broken_agreement.mean() * 100:.1f}', tol=0.1)
        chk('formed agreement %', '65.7', f'{ro.formed_agreement.mean() * 100:.1f}', tol=0.1)

    # 3) capacity ladder and multiplicity
    hm = rd('ladder_paired_effects_audited_holm.csv')
    if hm is not None:
        print('\n[capacity ladder and multiplicity]')
        chk('comparison count', 13, len(hm))
        chk('uncorrected significant', 1,
            int((((hm.fixed_lo > 0) | (hm.fixed_hi < 0))).sum()))
        chk('significant after Holm', 0, int((hm.holm_p < 0.05).sum()))
        chk('smallest Holm P', 0.585, f'{hm.holm_p.min():.3f}', tol=0.002)

    # 4) capacity trend
    tr = rd('ladder_capacity_trend_audited.csv')
    if tr is not None:
        print('\n[capacity trend, common definition pairs]')
        chk('Spearman rho', -0.800, f'{tr.rho.iloc[0]:.3f}', tol=0.002)
        chk('exact two-sided P', 0.333, f'{tr.exact_p.iloc[0]:.3f}', tol=0.002)

    # 5) structure increment
    fu = rd('esm2_fusion_results.csv')
    if fu is not None:
        print('\n[structure increment over the sequence baseline]')
        for bm, er, elo, ehi in [('s669', -0.013, -0.036, 0.009),
                                 ('ssym', 0.025, -0.007, 0.070)]:
            row = fu[(fu.benchmark == bm) & (fu.model.str.contains('increment'))]
            if len(row):
                r0 = row.iloc[0]
                chk(f'{bm} delta r', f'{er:.3f}', f'{r0.r:.3f}', tol=0.002)
                chk(f'{bm} 95% low', f'{elo:.3f}', f'{r0.ci_low:.3f}', tol=0.002)
                chk(f'{bm} 95% high', f'{ehi:.3f}', f'{r0.ci_high:.3f}', tol=0.002)

    # 6) typed contact correlations, recomputed from per-mutation rows
    tp = rd('directional_type_analysis.csv')
    if tp is not None:
        print('\n[typed contact correlations, recomputed from per-mutation rows]')
        for col, exp in [('n_broken', 0.061), ('n_formed', 0.083), ('n_edit', 0.081),
                         ('net_change', -0.013), ('broken_hydro', 0.135),
                         ('broken_elec', 0.047), ('broken_other', 0.027),
                         ('formed_hydro', 0.096), ('formed_elec', -0.046),
                         ('formed_other', 0.063)]:
            if col in tp.columns:
                r = np.corrcoef(tp.ddg, tp[col])[0, 1]
                chk(f'{col} r', f'{exp:.3f}', f'{r:.3f}', tol=0.002)

    print()
    print('=' * 76)
    print(f'Recomputed: {sum(RESULTS)}/{len(RESULTS)} pass')
    print('=' * 76)
    return 0 if all(RESULTS) else 1


if __name__ == '__main__':
    raise SystemExit(main())
