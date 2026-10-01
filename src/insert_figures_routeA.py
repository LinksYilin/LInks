# -*- coding: utf-8 -*-
"""
insert_figures_routeA.py — 按路线 A 重排图：插入新 Figure 3 / Figure 4，旧 Figure 3 → Figure 5
==============================================================================================
现有：Figure 1（旧数据）/ Figure 2（有效）/ Figure 3（panel a 失效，b/c 有效）
目标：Figure 1 / Figure 2 / **Figure 3 能力阶梯（新）** / **Figure 4 信息预算（新）** / Figure 5（旧 3）
"""
import os
from docx import Document
from docx.shared import Inches

DOC = r'D:\GED_mutation\manuscript_routeA.docx'
FIGDIR = r'D:\GED_mutation\figures\publication'

CAP3 = (
    'Figure 3 | A capacity ladder: the representation effect does not grow with encoder size. '
    'a, Mutation-level Pearson r on S669 for five encoders spanning 5.8 k to 509 k parameters, '
    'each trained on the audited training set under four graph definitions (C\u03b1, C\u03b2, '
    'side-chain centroid, all-atom). b, Mean absolute representation effect across definitions '
    '(within-encoder paired difference against the side-chain-centroid definition) as a function '
    'of encoder capacity. Across the four stable encoders none of the eight paired comparisons '
    'reached significance. The 509 k-parameter equivariant network (open marker) produced a larger '
    'apparent effect but collapsed on a subset of seeds in every definition tested, so its value '
    'reflects optimisation instability rather than a representation effect. Encoder results are '
    'reported as the range across definitions in Table 1; the corresponding ssym graphs were not '
    'constructed for these encoders.')

CAP4 = (
    'Figure 4 | An information budget for \u0394\u0394G prediction. a, Out-of-sample S669 Pearson '
    'r for physicochemical (P), graph-derived structural (S) and frozen ESM-2 (E) feature blocks '
    'and their combinations, fitted on the training set and evaluated on the benchmark. Error bars '
    'are 95% protein-cluster bootstrap intervals. b, Sequence baseline and structure increment. '
    'Three models were trained under an identical protocol on S669 and ssym: a frozen ESM-2 650M '
    'embedding with a learned head, the same embedding fused with a six-layer GINE graph branch, '
    'and the graph branch alone. Annotations give the paired increment of the fusion over the '
    'sequence-only model with its 95% protein-cluster bootstrap interval (B = 2,000); both '
    'intervals include zero. The graph branch alone performs close to the sequence baseline, so '
    'the null increment is not attributable to an under-powered graph branch.')

CAP5_OLD = (
    'Figure 5 | Contact-edit diagnostics and cross-engine consistency. a, Paired mutation-level '
    'Pearson correlations for local GNNs with and without additional sequence-derived features. '
    'b, Experimental \u0394\u0394G versus total side-chain-centroid contact edits for 505 '
    'quality-controlled mutations; the ordinary least-squares line is descriptive (r = 0.081, 95% '
    'CI \u22120.020 to +0.214). c, FoldX-predicted versus experimental \u0394\u0394G for 538 '
    'available pairs; the solid line is an ordinary least-squares fit and the dashed line is '
    'identity (r = 0.26). Panels b and c use separate available-case datasets.')


def insert_picture_before(par, image, width_in=6.2):
    doc = par._parent
    newp = doc.add_paragraph()
    newp.alignment = 1
    newp.add_run().add_picture(image, width=Inches(width_in))
    par._p.addprevious(newp._p)
    return newp


def insert_paragraph_before(par, text):
    doc = par._parent
    newp = doc.add_paragraph()
    newp.add_run(text)
    par._p.addprevious(newp._p)
    return newp


def main():
    doc = Document(DOC)
    ps = [p.text.strip() for p in doc.paragraphs]

    # 定位旧 Figure 3 caption
    i3 = next(i for i, t in enumerate(ps) if t.startswith('Figure 3 |'))

    # 在旧 Figure 3 之前插入新 Figure 3 与 Figure 4（图 + 图注）
    f3 = os.path.join(FIGDIR, 'figure3_capacity_ladder.png')
    f4 = os.path.join(FIGDIR, 'figure4_information_budget.png')
    anchor = doc.paragraphs[i3]      # ★ 只取一次：插入后索引会偏移
    for img, cap in [(f3, CAP3), (f4, CAP4)]:
        if not os.path.exists(img):
            print(f'  ⚠ 缺图 {img}')
            continue
        insert_picture_before(anchor, img)
        insert_paragraph_before(anchor, cap)
    print('  [1] 新 Figure 3 / Figure 4 已插入')
    doc.save(DOC)          # ★ 必须先保存，否则插入丢失

    # 旧 Figure 3 → Figure 5
    doc2 = Document(DOC)
    n = 0
    for p in doc2.paragraphs:
        if p.text.strip().startswith('Figure 3 | Fixed edit descriptors'):
            if p.runs:
                p.runs[0].text = CAP5_OLD
                for r in p.runs[1:]:
                    r.text = ''
            n += 1
    doc2.save(DOC)
    print(f'  [2] 旧 Figure 3 → Figure 5（{n} 处）')


if __name__ == '__main__':
    main()
