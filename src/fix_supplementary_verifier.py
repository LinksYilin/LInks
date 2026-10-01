# -*- coding: utf-8 -*-
"""fix_supplementary_verifier.py — 同步补充材料到修正后的口径"""
import os
import re

SUP = r'D:\GED_mutation\补充材料_Supplementary.md'
ROOT_SUP = os.path.join(r'D:\GED_mutation', '补充材料_Supplementary.md')

RULES = [
    # S2：嵌套声明
    ('The evaluation sets are nested subsets of the 543-mutation S669 benchmark.',
     'The S669 evaluation sets were built by successive filters from the 543-mutation benchmark. '
     'The filters are not strictly nested: the 505-pair contact-change set was constructed from '
     'the 538 modelled mutant structures independently of the 511- and 508-mutation predictive '
     'sets, so it contains four pairs that the predictive sets lack and omits seven that they '
     'contain.'),
    ('| Contact graph available | 511 | 1 | graph could not be built (1G3P) |',
     '| Contact graph available | 511 | 1 | graph could not be built (1G3P) |'),
    ('| Matching FoldX mutant structure | 505 | 3 | modelled mutant residue count differed |',
     '| FoldX output with matching residue count | 505 | 33 removed, 2 re-admitted | incomplete or '
     'mismatched FoldX output, measured against the 538 modelled pairs |'),
    ('Consequently: **505** mutations for contact-change statistics, **511** for predictive '
     'comparisons without sequence features, **508** for sequence and fusion comparisons.',
     'Consequently: **505** mutations for contact-change statistics, **511** for predictive '
     'comparisons without sequence features, **508** for sequence and fusion comparisons. The '
     '505-pair set was built from the 538 modelled pairs (538 \u2212 33 + 2 = 507 candidate pairs, '
     'of which 505 carry a complete graph).'),
    ('Training data came from MegaScale plus ThermoMutDB after BLAST leakage filtering, sampled '
     '1:1 to 7,905 examples drawn from **420 proteins** (232 MegaScale, 188 ThermoMutDB).',
     'Training data came from MegaScale plus ThermoMutDB after BLAST leakage filtering, sampled '
     '1:1 to 7,905 examples drawn from **420 distinct proteins** (232 sampled from MegaScale and '
     '191 from ThermoMutDB, with three proteins present in both sources).'),

    # S7：EGNN 三种子说明已正确，补一句关于 ensemble 标签
    ('Rescaling the coordinate update did not stabilise the encoder.',
     'All three seeds share the same 508,934-parameter architecture. Rescaling the coordinate '
     'update did not stabilise the encoder.'),

    # S9：ensemble 声明修正
    ('Every reported correlation is computed on predictions averaged over the stated seeds, and '
     'every paired interval is computed on that seed-averaged prediction.',
     'Two aggregations are used and are distinguished throughout. The correlation reported for '
     'each encoder and graph definition is the unweighted mean of the per-seed correlations '
     '(`seed_mean_r`). The paired definition effects reported in S5 are computed from the joint '
     'seed and protein bootstrap, in which the seed-averaged prediction enters the resampling. '
     'Averaging per-seed correlations and correlating seed-averaged predictions give different '
     'values; for the four stable encoders the difference is at most 0.033, and for the '
     'equivariant network it reaches 0.216 because its per-seed correlations are bimodal.'),

    # S6：阈值表的种子说明
    ('(GCN with mutation-site pooling, seeds 42 and 123, predictions averaged over the two seeds '
     'before correlating)',
     '(GCN with mutation-site pooling, seeds 42 and 123, predictions averaged over the two seeds '
     'before correlating; the ladder in S5 uses three seeds for the same encoder)'),
]

# S1.1：before/after 表中的过时范围
RULES.append(
    ('| Edge-aware GINE, range across definitions (S669) | 0.289\u20130.354 (three seeds) | '
     '0.340\u20130.373 (three seeds) |',
     '| Edge-aware GINE, range across definitions (S669) | 0.289\u20130.354 (three seeds, '
     'H-included training graphs) | 0.340\u20130.373 (means of three per-seed correlations) |'))


def main():
    for path in {SUP, ROOT_SUP}:
        if not os.path.exists(path):
            continue
        t = open(path, encoding='utf-8').read()
        orig = t
        n = 0
        for a, b in RULES:
            if a in t:
                t = t.replace(a, b)
                n += 1
        if t != orig:
            open(path, 'w', encoding='utf-8').write(t)
            print(f'{path}: 修改 {n} 处')

    t = open(SUP, encoding='utf-8').read()
    print()
    for pat, name in [('nested subsets', '残留 nested subsets'),
                      ('188 ThermoMutDB', '残留 188'),
                      ('191 from ThermoMutDB', '新 191'),
                      ('not strictly nested', '新嵌套说明'),
                      ('unweighted mean of the per-seed correlations', '新 ensemble 说明')]:
        print(f'  {name}: {t.count(pat)}')


if __name__ == '__main__':
    main()
