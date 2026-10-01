# -*- coding: utf-8 -*-
"""fix_verifier_5.py — 修正 §4.1 的嵌套声明与两种子数值，以及"ensemble"误标"""
from docx import Document

DOC = r'D:\GED_mutation\manuscript_routeA.docx'
DEL = '\u0394'
A = '\u03b1'
B = '\u03b2'
MIN = '\u2212'
ALPHA = A

RULES = [
    # 9) 嵌套声明 → 如实说明各阶段计数与集合关系
    ('The evaluation sets are nested subsets of the 543-mutation S669 benchmark: 512 mutations '
     'passed the residue-index verification, 511 had a contact graph for the required definition, '
     '508 additionally had a cached ESM-2 embedding, and 505 additionally had a FoldX mutant '
     'structure whose residue count matched the wild-type. Contact-change statistics therefore '
     'use the 505 quality-controlled pairs, predictive comparisons that do not involve sequence '
     'features use all 511 mutations, and the sequence and fusion comparisons use 508.',
     'The S669 evaluation sets were built from the 543-mutation benchmark by successive filters, '
     'with the counts measured directly at each stage: 512 mutations passed the residue-index '
     'verification (31 removed, all in 3DV0), 511 had a contact graph for the required definition '
     '(1 removed, 1G3P), and 508 additionally had a cached ESM-2 embedding (3 removed, in 1O6X '
     'and 2HBB). The 505-pair contact-change set was built from the 538 modelled mutant structures '
     'by requiring a complete FoldX output whose modelled mutant matched the wild-type residue '
     'count; this filter removes 33 pairs and re-admits two, so the contact-change set is not a '
     'subset of the 511- or 508-mutation predictive sets. Contact-change statistics therefore use '
     'the 505 quality-controlled pairs, predictive comparisons that do not involve sequence '
     'features use the 511 mutations, and the sequence and fusion comparisons use 508.'),

    # 11) §4.1 两种子 → 三种子（mutation-site GCN）
    ('point estimates of r = 0.380 for the side-chain centroid, 0.344 for C\u03b2, 0.323 for '
     'all-atom minimum distance, and 0.312 for C' + ALPHA + '. No definition differed '
     'significantly from the centroid in a protein-cluster paired bootstrap (C' + ALPHA + ' ' + DEL +
     'r = ' + MIN + '0.068, 95% CI ' + MIN + '0.129 to +0.011, P = 0.096; C' + B + ' ' + DEL +
     'r = ' + MIN + '0.037, ' + MIN + '0.088 to +0.018, P = 0.17; all-atom ' + DEL + 'r = ' +
     MIN + '0.058, ' + MIN + '0.115 to +0.016, P = 0.10).',
     'point estimates of r = 0.378 for the side-chain centroid, 0.347 for C' + B + ', 0.331 for '
     'all-atom minimum distance, and 0.314 for C' + ALPHA + ' (means of three per-seed '
     'correlations). No definition differed significantly from the centroid in a joint seed and '
     'protein-cluster paired bootstrap (C' + ALPHA + ' ' + DEL + 'r = ' + MIN + '0.064, 95% CI ' +
     MIN + '0.128 to +0.016, P = 0.11; C' + B + ' ' + DEL + 'r = ' + MIN + '0.031, ' + MIN +
     '0.085 to +0.029, P = 0.28; all-atom ' + DEL + 'r = ' + MIN + '0.049, ' + MIN + '0.107 to '
     '+0.031, P = 0.18).'),

    # 15) "three-seed ensemble" → 正确标签（Table 1 脚注）
    ('the capacity ladder reports the three-seed ensemble value for each encoders range across '
     'graph definitions on S669',
     'the capacity ladder reports, for each encoder, the mean of the three per-seed correlations '
     'across graph definitions on S669'),
    ('Encoder results on S669 are reported as the range of three-seed ensemble correlations across',
     'Encoder results on S669 are reported as the range of three-seed mean-of-per-seed '
     'correlations across'),

    # 15) 正文 [123]
    ('(all three-seed ensembles), and the 509 k-parameter equivariant network',
     '(all means of three per-seed correlations), and the 509 k-parameter equivariant network'),
    ('The highest value among graph encoders was r = 0.373 (three-seed ensemble) on the C' + B +
     ' definition',
     'The highest value among graph encoders was r = 0.373 (mean of three per-seed correlations) '
     'on the C' + B + ' definition'),

    # 15) Figure 3 图注
    ('Points are three-seed ensembles; the equivariant network uses three architecture-matched '
     'seeds',
     'Points are means of three per-seed correlations; the equivariant network uses three '
     'architecture-matched seeds'),
]


def main():
    d = Document(DOC)
    n = 0
    for p in d.paragraphs:
        if not p.runs:
            continue
        full = ''.join(r.text for r in p.runs)
        new = full
        for a, b in RULES:
            if a in new:
                new = new.replace(a, b)
                n += 1
        if new != full:
            p.runs[0].text = new
            for r in p.runs[1:]:
                r.text = ''
    d.save(DOC)
    print(f'修改 {n} 处')
    d2 = Document(DOC)
    txt = '\n'.join(p.text for p in d2.paragraphs)
    print()
    for pat, name in [('nested subsets', '残留 nested subsets'),
                      ('r = 0.380 for the side-chain', '残留 0.380'),
                      ('three-seed ensemble', '残留 ensemble 误标'),
                      ('not a subset of the 511', '新嵌套说明'),
                      ('mean of three per-seed', '新标签')]:
        print(f'  {name}: {txt.count(pat)}')


if __name__ == '__main__':
    main()
