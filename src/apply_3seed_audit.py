# -*- coding: utf-8 -*-
"""
apply_3seed_audit.py — 把审计后的三种子结果写入论文
======================================================
旧稿的容量阶梯基于两种子，且未做多重比较校正。
本脚本按 ladder_seed_summary_audited.csv 与 ladder_paired_effects_audited.csv
更新：摘要 / Results 4.2 / Table 1 / Figure 3 图注 / Discussion / Conclusions。
"""
import os

import pandas as pd
from docx import Document

DOC = r'D:\GED_mutation\manuscript_routeA.docx'
D = r'D:\GED_mutation\data'
DEL = '\u0394'
MINUS = '\u2212'
RHO = '\u03c1'
ALPHA = '\u03b1'
BETA = '\u03b2'

summ = pd.read_csv(os.path.join(D, 'ladder_seed_summary_audited.csv'))
eff = pd.read_csv(os.path.join(D, 'ladder_paired_effects_audited.csv'))
holm = pd.read_csv(os.path.join(D, 'ladder_paired_effects_audited_holm.csv'))
trend = pd.read_csv(os.path.join(D, 'ladder_capacity_trend_audited.csv'))

STABLE = ['gnn_global', 'gnn_local', 'gnn_edge', 'deep_gine']


def rng(model):
    s = summ[summ.model == model]
    return float(s.seed_mean_r.min()), float(s.seed_mean_r.max())


def main():
    n_sig_raw = int(((eff.fixed_lo > 0) | (eff.fixed_hi < 0)).sum())
    n_sig_holm = int((holm.holm_p < 0.05).sum())
    g, G = rng('gnn_global')
    l, L = rng('gnn_local')
    e, E = rng('gnn_edge')
    dg, DG = rng('deep_gine')

    repl = []

    # ---------------- 摘要 ----------------
    repl.append((
        'Prediction varied little across definitions: on a five-encoder capacity ladder spanning '
        '5.8 k to 509 k parameters, twelve of thirteen paired comparisons were non-significant, and '
        'the single significant one occurred at the smallest encoder, where the all-atom definition '
        'performed worse (' + DEL + 'r = ' + MINUS + '0.109, 95% CI ' + MINUS + '0.187 to '
        + MINUS + '0.017).',
        'Prediction varied little across definitions: on a five-encoder capacity ladder spanning '
        '5.8 k to 509 k parameters, one of thirteen paired comparisons was significant before '
        'correction and none survived Holm correction across the thirteen; the uncorrected case was '
        'the all-atom definition at the smallest encoder (' + DEL + 'r = ' + MINUS + '0.101, 95% CI '
        + MINUS + '0.186 to ' + MINUS + '0.002).'))

    # ---------------- Results 4.2 ----------------
    old_open = ('A single encoder family cannot distinguish \u201cthe representation does not '
                'matter\u201d from \u201cthis encoder cannot exploit it\u201d.')
    new_open = (
        'A single encoder family cannot distinguish \u201cthe representation does not matter\u201d '
        'from \u201cthis encoder cannot exploit it\u201d. We repeated the four-definition '
        'comparison at five encoder capacities, using three seeds for the four stable encoders and '
        'an architecture-matched three-seed set for the equivariant network.')
    repl.append((old_open, new_open))

    # 替换 4.2 的均值句（旧：0.071/0.054/0.025/0.039 + EGNN 0.061）
    repl.append((
        'Across the four stable encoders (5.8 k to 460 k parameters) the eleven paired definition '
        'comparisons yielded one significant result, and the mean absolute effect was 0.071, 0.054, '
        '0.025 and 0.039 respectively (Figure 3b).',
        'Restricting attention to the definition pairs available at every rung (C' + ALPHA + ' and '
        'C' + BETA + ' against the side-chain centroid), the mean absolute effect was '
        f'{(holm[holm.model=="gnn_global"].delta_r.abs().mean()):.3f}, '
        f'{(holm[holm.model=="gnn_local"].delta_r.abs().mean()):.3f}, '
        f'{(holm[holm.model=="gnn_edge"].delta_r.abs().mean()):.3f} and '
        f'{(holm[holm.model=="deep_gine"].delta_r.abs().mean()):.3f} across the four stable rungs, '
        'and an exact permutation test on the Spearman correlation with log parameter count gave '
        + RHO + ' = ' + MINUS + '0.80 with a two-sided exact P = 0.333 over all 24 permutations '
        '(Figure 3b).'))

    # 替换单显著句
    repl.append((
        'The single significant comparison occurred at the smallest encoder, where the all-atom '
        'definition performed worse than the side-chain centroid (' + DEL + 'r = ' + MINUS +
        '0.109, 95% CI ' + MINUS + '0.187 to ' + MINUS + '0.017; P = 0.019).',
        'One comparison was significant before correction: at the smallest encoder the all-atom '
        'definition performed worse than the side-chain centroid (' + DEL + 'r = ' + MINUS +
        '0.101, joint seed and protein bootstrap 95% CI ' + MINUS + '0.186 to ' + MINUS + '0.002; '
        'uncorrected P = 0.045). No comparison survived Holm correction across the thirteen '
        '(smallest corrected P = 0.585).'))

    # 替换 EGNN 不稳定句（三种子同架构）
    repl.append((
        'The 509 k-parameter equivariant network gave a mean absolute effect of 0.061 with no '
        'significant comparison, but it collapsed on a subset of seeds in every definition tested '
        '(for the side-chain-centroid graphs, r = 0.393 and ' + MINUS + '0.081 for two seeds), so '
        'its comparisons are reported as an optimisation instability rather than as representation '
        'effects.',
        'The 509 k-parameter equivariant network gave no significant comparison, but it collapsed '
        'on a subset of seeds in every definition tested. Its three architecture-matched seeds on '
        'the side-chain-centroid graphs gave r = 0.393, ' + MINUS + '0.081 and 0.079 (mean 0.131, '
        'SD 0.197), and on C' + BETA + ' r = ' + MINUS + '0.072, 0.358 and ' + MINUS + '0.014, so '
        'its comparisons are reported as an optimisation instability rather than as representation '
        'effects.'))

    # 替换单调性句（数值更新）
    repl.append((
        'The 50 k-parameter GINE encoder reached r = 0.402 (two-seed mean) on the C' + BETA +
        ' definition, the highest value among the graph encoders and above the physicochemical '
        'ridge baseline (r = 0.390), whereas the 460 k-parameter GINE reached 0.322 on the same '
        'definition.',
        f'The 50 k-parameter GINE encoder reached r = {E:.3f} (three-seed ensemble) on the C{BETA} '
        f'definition, the highest value among the graph encoders and close to the physicochemical '
        f'ridge baseline (r = 0.390), whereas the 460 k-parameter GINE reached {DG:.3f} on the same '
        'definition.'))

    # ---------------- Discussion 5.4 ----------------
    repl.append((
        'Across the four stable encoders, ten of the eleven paired definition comparisons were '
        'non-significant, and the mean absolute effect did not grow with capacity (0.071, 0.054, '
        '0.025 and 0.039).',
        'Across the four stable encoders, none of the thirteen paired definition comparisons '
        'survived Holm correction, and the mean absolute effect over the definition pairs common to '
        'every rung did not grow with capacity.'))

    repl.append((
        'we repeated the four-definition comparison at five encoder capacities spanning 5.8 k to '
        '509 k parameters. Across the four stable encoders, one of the thirteen paired comparisons '
        'reached significance: all-atom versus side-chain centroid at the smallest (5.8 k) encoder, '
        + DEL + 'r = ' + MINUS + '0.109 (95% CI ' + MINUS + '0.187 to ' + MINUS + '0.017; '
        'P = 0.019), where the all-atom definition performed worse. The remaining twelve '
        'comparisons were non-significant.',
        'we repeated the four-definition comparison at five encoder capacities spanning 5.8 k to '
        '509 k parameters. One of the thirteen paired comparisons was significant before '
        'correction, all-atom versus side-chain centroid at the smallest (5.8 k) encoder, where '
        'the all-atom definition performed worse (' + DEL + 'r = ' + MINUS + '0.101; uncorrected '
        'P = 0.045), and none survived Holm correction across the thirteen. The remaining twelve '
        'comparisons were non-significant before correction.'))

    # ---------------- Conclusions ----------------
    repl.append((
        'twelve of the thirteen paired definition comparisons were non-significant, and adding a '
        'contact-graph branch to an ESM-2 650M sequence baseline produced no increment larger than '
        '|' + DEL + 'r| \u2248 0.04 on S669 and \u2248 0.07 on ssym.',
        'no paired definition comparison survived Holm correction across the thirteen, and adding '
        'a contact-graph branch to an ESM-2 650M sequence baseline produced no increment larger '
        'than about +0.01 on S669 and +0.07 on ssym.'))

    doc = Document(DOC)
    applied = 0
    for p in doc.paragraphs:
        if not p.runs:
            continue
        full = ''.join(r.text for r in p.runs)
        new = full
        for a, b in repl:
            if a in new:
                new = new.replace(a, b)
                applied += 1
        if new != full:
            p.runs[0].text = new
            for r in p.runs[1:]:
                r.text = ''
    doc.save(DOC)
    print(f'应用 {applied} 处文本修订')
    print(f'  未校正显著 {n_sig_raw}/13；Holm 后显著 {n_sig_holm}/13')


if __name__ == '__main__':
    main()
