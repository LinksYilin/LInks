# -*- coding: utf-8 -*-
"""fix_verifier_4.py — 修正 §4.9、§4.1、§4.7、§4.3、Figure 5 标题"""
from docx import Document

DOC = r'D:\GED_mutation\manuscript_routeA.docx'
DEL = '\u0394'
A = '\u03b1'
B = '\u03b2'
MIN = '\u2212'

RULES = [
    # 6) §4.9 案例研究距离（实测 19.41 / 17.76 / 6.78）
    ('The broken contacts were 6.5 \u00c5 from the mutated residue in each case.',
     'The mean distance from the mutated residue to a broken contact was 19.41 \u00c5 for R244E, '
     '17.76 \u00c5 for L32P and 6.78 \u00c5 for W36A, so the three cases span a wide range of '
     'contact geometry rather than sharing one local pattern.'),

    # 19) §4.7 SCWRL4 / FoldX Cβ 位移：改为实测值并限定范围
    ('FoldX largely retained C\u03b2 coordinates (mean displacement 1 \u00d7 10\u207b\u2074 \u00c5; '
     'maximum 0.012 \u00c5), whereas SCWRL4 rebuilt them (mean 0.096 \u00c5; maximum 0.300 \u00c5).',
     'FoldX BuildModel retains the wild-type side-chain template for unmutated residues, whereas '
     'SCWRL4 rebuilt them. Over the 11 mutations for which matched wild-type input and SCWRL4 '
     'output were compared, the SCWRL4 C' + B + ' displacement had a per-case mean of 0.045 to '
     '0.136 \u00c5 and a maximum of 0.425 \u00c5, while the C' + A + ' displacement was 0.0000 '
     '\u00c5 in every case.'),

    # 20) §4.3 ESM-2 zero-shot 的样本量
    ('Its weakness relative to the graph encoders on the same mutations indicates',
     'Its weakness relative to the graph encoders, evaluated on the 512 S669 mutations for which '
     'an ESM-2 embedding was available, indicates'),
    ('it scores below the graph encoders on the same mutations.',
     'it scores below the graph encoders over the mutations all three can score.'),

    # 14) Figure 5 标题去掉不存在的 cross-engine 面板
    ('Figure 5 | Contact-edit diagnostics and cross-engine consistency.',
     'Figure 5 | Contact-edit diagnostics for the side-chain-centroid definition.'),
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
    for pat, name in [('6.5 \u00c5 from the mutated', '残留 6.5 Å'),
                      ('maximum 0.300', '残留 0.300'),
                      ('mean 0.096', '残留 0.096'),
                      ('cross-engine consistency', '残留 cross-engine'),
                      ('19.41', '新 19.41'), ('0.425', '新 0.425')]:
        print(f'  {name}: {txt.count(pat)}')


if __name__ == '__main__':
    main()
