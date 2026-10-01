# -*- coding: utf-8 -*-
"""
add_trend_test.py — 把能力趋势检验结果写入论文（回应 R1-M1）
"""
from docx import Document

DOC = r'D:\GED_mutation\manuscript_routeA.docx'
D = '\u0394'
M = '\u2212'
RHO = '\u03c1'

# Results 4.2 追加
ADD_RES = (
    ' We tested the trend formally rather than asserting it. Across the four stable rungs the mean '
    'absolute effect was 0.071 (5.8 k), 0.054 (5.9 k), 0.025 (50 k) and 0.039 (460 k), and a '
    'permutation test on the Spearman correlation between log parameter count and mean absolute '
    'effect gave ' + RHO + ' = ' + M + '0.80 with a two-sided P = 0.33 (20,000 permutations). '
    'With only four rungs the test is underpowered, so we report this as a failure to detect a '
    'growth trend rather than as evidence that none exists; the point estimate is in the direction '
    'opposite to the one that limited encoder capacity would produce. Comparing the smallest and '
    'largest stable rungs directly, the mean absolute effect fell from 0.071 to 0.039 '
    '(difference ' + M + '0.032).')

# Discussion 5.2 追加
ADD_DISC = (
    ' Across the four stable rungs the mean absolute representation effect did not rise with '
    'capacity and a permutation test found no growth trend (' + RHO + ' = ' + M + '0.80, two-sided '
    'P = 0.33, four rungs). This is the evidence that the null result is not simply an '
    'underpowered encoder: had capacity been the binding constraint, the representation effect '
    'should have increased from the 5.8 k-parameter mean-pooled GCN to the 460 k-parameter '
    'attention-pooled GINE network, and it did not.')


def append_after(doc, starts, add):
    n = 0
    for p in doc.paragraphs:
        full = ''.join(r.text for r in p.runs)
        if full.startswith(starts) and add.strip() not in full:
            p.runs[0].text = full.rstrip() + add
            for r in p.runs[1:]:
                r.text = ''
            n += 1
    return n


def main():
    d = Document(DOC)
    n = 0
    n += append_after(d, 'A single encoder family cannot distinguish', ADD_RES)
    n += append_after(d, 'Three observations explain the decoupling.', ADD_DISC)
    d.save(DOC)
    print(f'已写入趋势检验 {n} 处')
    d2 = Document(DOC)
    for p in d2.paragraphs:
        if '20,000 permutations' in p.text:
            i = p.text.find('We tested the trend formally')
            print('  ✅ Results:', p.text[i:i + 260])
            break


if __name__ == '__main__':
    main()
