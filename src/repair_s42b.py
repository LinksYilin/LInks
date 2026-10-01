# -*- coding: utf-8 -*-
"""repair_s42b.py — 用 insert_paragraph_before 重建 §4.2"""
from docx import Document

DOC = r'D:\GED_mutation\manuscript_routeA.docx'
A = '\u03b1'
B = '\u03b2'
RHO = '\u03c1'
DEL = '\u0394'
MIN = '\u2212'

BODY = [
    ('A single encoder family cannot distinguish \u201cthe representation does not matter\u201d '
     'from \u201cthis encoder cannot exploit it\u201d. We therefore trained five encoders of '
     'increasing capacity on identical data and splits, and repeated the four-definition '
     'comparison at each rung: a mean-pooled GCN (5.8 k parameters), a mutation-site-pooled GCN '
     '(5.9 k), a GINE model with explicit edge features (50 k), a six-layer GINE with attention '
     'pooling and a mutation-site local readout (460 k), and an E(3)-equivariant network with '
     'coordinate updates (509 k). Every model was trained for 20 epochs with the same optimisation '
     'settings, using three seeds for the four stable encoders and an architecture-matched '
     'three-seed set for the equivariant network. Representation effects were estimated as '
     'within-encoder paired differences against the side-chain-centroid definition, using a joint '
     'seed and protein-cluster bootstrap (B = 2,000).'),
    ('Restricting attention to the definition pairs available at every rung (C' + A + ' and C' + B +
     ' against the side-chain centroid), the mean absolute effect was 0.048, 0.048, 0.021 and '
     '0.043 across the four stable rungs, and an exact permutation test on the Spearman correlation '
     'with log parameter count gave ' + RHO + ' = ' + MIN + '0.80 with a two-sided exact P = 0.333 '
     'over all 24 permutations (Figure 3b). One comparison was significant before correction: at '
     'the smallest encoder the all-atom definition performed worse than the side-chain centroid ('
     + DEL + 'r = ' + MIN + '0.101, joint seed and protein bootstrap 95% CI ' + MIN + '0.186 to '
     + MIN + '0.002; uncorrected P = 0.045). No comparison survived Holm correction across the '
     'thirteen (smallest corrected P = 0.585). The effect therefore does not increase with '
     'capacity, which is the pattern expected if encoder strength were the bottleneck. With only '
     'four rungs the trend test is underpowered, so we report this as a failure to detect a growth '
     'trend rather than as evidence that none exists; the point estimate lies in the direction '
     'opposite to the one that limited encoder capacity would produce.'),
    ('The 509 k-parameter equivariant network gave a mean absolute effect of 0.059 with no '
     'significant comparison, but it was unstable across seeds in every definition tested. On the '
     'side-chain-centroid graphs, three seeds gave r = 0.393, ' + MIN + '0.081 and 0.079 (mean '
     '0.131, sample SD 0.241); a variant whose coordinate-update step was rescaled from zero at '
     'initialisation did not reduce this spread (r = ' + MIN + '0.182, 0.346 and 0.076; mean 0.080, '
     'sample SD 0.264). Its comparisons are therefore reported as an optimisation instability '
     'rather than as representation effects.'),
    ('Model performance itself was not monotone in capacity. The 50 k-parameter GINE encoder '
     'reached r = 0.373, the highest value among the graph encoders on the C' + B + ' definition '
     'and close to the physicochemical ridge baseline (r = 0.390), whereas the 460 k-parameter '
     'GINE reached 0.331 on the same definition. These values are means of the three per-seed '
     'correlations. The same ladder on the development-independent ssym benchmark gave 0.147, '
     '0.409, 0.452 and 0.478 for the four stable encoders, so the 460 k-parameter attention-pooled '
     'GINE network exceeded both the physicochemical baseline (0.305) and the narrower sequence '
     'head (0.458), and fell below only the wider sequence head (0.516). The graph encoders '
     'therefore carry real predictive signal on the independent benchmark.'),
]


def main():
    d = Document(DOC)
    paras = list(d.paragraphs)
    i43 = next(i for i, p in enumerate(paras) if p.text.strip().startswith('4.3 '))
    anchor = paras[i43]
    print(f'§4.3 标题在段 {i43}: {anchor.text.strip()[:50]}')

    for text in BODY:
        anchor.insert_paragraph_before(text)

    d.save(DOC)
    print(f'已插入 {len(BODY)} 段')

    d2 = Document(DOC)
    ps = [p.text.strip() for p in d2.paragraphs]
    i42 = next(i for i, t in enumerate(ps) if t.startswith('4.2 '))
    i43b = next(i for i, t in enumerate(ps) if t.startswith('4.3 '))
    print(f'\n§4.2 现有 {i43b - i42 - 1} 段:')
    for i in range(i42, i43b):
        if ps[i]:
            print(f'  [{i}] {ps[i][:120]}')


if __name__ == '__main__':
    main()
