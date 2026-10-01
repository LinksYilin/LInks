# -*- coding: utf-8 -*-
"""
fix_contrib_and_framing.py — 回应 R2-m2（贡献列表改为发现）与 R2-M3（信息预算措辞）
"""
from docx import Document

DOC = r'D:\GED_mutation\manuscript_routeA.docx'
DEL = '\u0394'
ALPHA = '\u03b1'
BETA = '\u03b2'
MINUS = '\u2212'

NEW = [
    # R2-m2：贡献列表改成"发现"
    ('This study makes four contributions to the extant literature.',
     'This study reports four findings.'),

    ('1.    A representative comparison. We show that C' + BETA + ' graphs remain unchanged '
     'under fixed-backbone side-chain repacking and compare the side-chain centroid graphs '
     'against the null representation, quantifying the response of each to side-chain repacking.',
     '1.    What a definition reports is set by the definition, not by the mutation. Processing '
     'wild-type and modelled mutant structures identically, C' + ALPHA + ' graphs recorded no '
     'contact change at the reference cutoff whereas C' + BETA + ', side-chain-centroid and '
     'all-atom graphs changed in 9.9%, 90.7% and 81.6% of 505 quality-controlled pairs. A '
     'representation can therefore be vacuous by construction, and the choice is a substantive '
     'modelling decision rather than an implementation detail.'),

    ('2.    A FoldX contact-change description. We compare wild-type and modelled mutant centroid '
     'graphs, quantify the reproducibility of broken and formed contacts across FoldX runs, and '
     'treat the result as an observed structural description rather than a causal attribution.',
     '2.    The contact-change description is reproducible but engine-dependent, and local. '
     'Repeated FoldX runs agreed on 96.2% of broken-contact sets but only 65.7% of formed-contact '
     'sets across three case-study mutations, and broken contacts sat on average 4.03 ' + '\u00c5' +
     ' from the mutated residue versus 15.70 ' + '\u00c5' + ' for contacts that persisted. The '
     'C' + ALPHA + ' result reproduced under SCWRL4 while the C' + BETA + ' contact-change rate did '
     'not (9.9% versus 97.8%), so absolute contact counts are a property of the modelling engine.'),

    ('3.    A controlled evaluation was performed. Under a fixed protocol and leakage control, we '
     'tested whether each representation and the structural descriptors derived from it improved '
     'the ' + DEL + DEL + 'G prediction beyond a physicochemical baseline on the exploratory S669 '
     'benchmark.',
     '3.    What a definition predicts does not follow from what it reports. Across a five-encoder '
     'capacity ladder from 5.8 k to 509 k parameters, one of thirteen paired definition comparisons '
     'was significant before correction and none survived Holm correction, and the mean absolute '
     'effect did not grow with capacity (exact permutation test over four rungs, ' + '\u03c1' + ' = '
     + MINUS + '0.80, P = 0.333).'),

    ('4.    A mechanistic analysis. We examined why representational sensitivity does not translate '
     'into predictive gain by quantifying the relationship between contact-based descriptors and '
     'physicochemical features, the strength of the association between individual contact changes '
     'and ' + DEL + DEL + 'G, and how that association depends on the spatial relationship between '
     'the contact and the mutated residue.',
     '4.    Structure adds no reliable increment once sequence is available, and the one structural '
     'signal is redundant. Adding a side-chain-centroid graph to an ESM-2 650M baseline produced no '
     'increment above +0.010 on S669 or +0.070 on ssym, and the only contact-derived quantity that '
     'survived multiplicity control, the count of broken hydrophobic contacts (r = 0.135), added no '
     'increment above about +0.04 over five physicochemical features.'),

    # R2-M3：信息预算措辞降级为"特征块核算"
    ('4.4 An information budget: where the predictive signal is',
     '4.4 Feature-block accounting: where the predictive signal is'),

    ('Figure 4 | An information budget for ' + DEL + DEL + 'G prediction.',
     'Figure 4 | Feature-block accounting for ' + DEL + DEL + 'G prediction.'),
]

# 标题也含 "information-budget"，检查后单独处理
TITLE_OLD = ('Contact-graph definitions shape what a model describes but not what it predicts: '
             'a controlled comparison across five encoder capacities')


def main():
    doc = Document(DOC)
    n = 0
    for p in doc.paragraphs:
        if not p.runs:
            continue
        full = ''.join(r.text for r in p.runs)
        new = full
        for a, b in NEW:
            if a in new:
                new = new.replace(a, b)
                n += 1
        if new != full:
            p.runs[0].text = new
            for r in p.runs[1:]:
                r.text = ''
    doc.save(DOC)
    print(f'应用 {n} 处（贡献列表改发现 + 信息预算措辞）')

    # 验证
    d2 = Document(DOC)
    txt = '\n'.join(p.text for p in d2.paragraphs)
    for pat, name in [('information budget', '旧"information budget"'),
                      ('An information budget', '旧标题'),
                      ('reports four findings', '新贡献列表'),
                      ('Feature-block accounting', '新小节标题')]:
        print(f'  {name}: {txt.count(pat)}')


if __name__ == '__main__':
    main()
