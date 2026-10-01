# -*- coding: utf-8 -*-
"""
rebuild_table1.py — 用修复后数据重建 Table 1
==============================================
旧表的所有 GNN 数字基于含氢 bug 的训练数据，已失效。
新表结构：编码器/基线 × (S669 r, ssym r)，并显式标注未跑的单元格。
"""
import os
import sys

import pandas as pd
from docx import Document
from docx.shared import Pt
from docx.oxml.ns import qn

sys.path.insert(0, os.path.dirname(__file__))
from paths import DATA

DOC = r'D:\GED_mutation\manuscript_routeA.docx'
OUT = str(DATA)

ROWS = [
    ['Method', 'S669 Pearson r', 'ssym Pearson r'],

    ['Ridge, 5 physicochemical features',
     '0.390 [0.264, 0.494]', '0.305 [0.141, 0.470]'],
    ['ESM-2 650M, zero-shot (no \u0394\u0394G supervision)',
     '0.317 [0.183, 0.446]', '0.266 [\u22120.008, 0.403]'],
    ['ESM-2 650M, frozen embeddings + learned head',
     '0.348\u20130.392 \u00a7', '0.458\u20130.516 \u00a7'],
    ['GCN, global pooling (5.8 k)',
     '0.052\u20130.150 \u2020', '0.147'],
    ['GCN, mutation-site pooling (5.9 k)',
     '0.314\u20130.378 \u2020', '0.409'],
    ['GINE, edge-aware (50 k)',
     '0.340\u20130.373 \u2020', '0.452'],
    ['Deep GINE, attention pooling (460 k)',
     '0.314\u20130.362 \u2020', '0.478'],
    ['EGNN, E(3)-equivariant (509 k) \u2021',
     '0.091\u20130.131 \u2020', '0.235'],
    ['ESM-2 650M + contact-graph fusion',
     '0.379', '0.540'],
    ['Structure increment (fusion \u2212 ESM-2 only)',
     '\u22120.013 [\u22120.036, +0.010]', '+0.025 [\u22120.007, +0.070]'],
]

FOOTNOTE = (
    'Table 1 | Predictive performance on the common S669 intersection and on ssym. All methods '
    'were evaluated on mutations shared with the training set after leakage filtering. Values for '
    'the ESM-2 baseline and the fusion are three-seed means; the capacity ladder reports the '
    'three-seed ensemble value for each encoders range across graph definitions on S669 '
    '(Figure 3) and the side-chain-centroid value on ssym. '
    '\u2020 Encoder results on S669 are reported as the range of three-seed ensemble correlations '
    'across the C\u03b1, C\u03b2, side-chain-centroid and all-atom definitions; the ssym column '
    'gives the side-chain-centroid definition, which is the reference definition throughout. '
    '\u2021 The equivariant encoder was run with the same architecture as the rest of the ladder '
    '(508,934 parameters) for all three seeds, and collapsed on a subset of seeds in every '
    'definition tested: on side-chain-centroid graphs its S669 seeds gave r = 0.393, \u22120.081 '
    'and 0.079 (mean 0.131, SD 0.197), and on C\u03b2 r = \u22120.072, 0.358 and \u22120.014. Its '
    'values reflect optimisation instability rather than a representation effect, and its ssym '
    'row uses two seeds. '
    'On the development-independent benchmark the 460 k-parameter attention-pooled GINE network '
    'reached the highest value of any graph encoder (0.478), above the physicochemical baseline '
    '(0.305) and the narrower sequence head (0.458), and below only the wider sequence head (0.516). '
    'The structure increment is the paired difference between the fusion and the sequence-only '
    'model, with a 95% protein-cluster bootstrap interval (B = 2,000). These intervals are '
    'asymmetric and are reported directionally: on S669 the increment is excluded below '
    '\u22120.036 and above +0.010; on ssym, below \u22120.007 and above +0.070. A single symmetric '
    'bound cannot be quoted for either benchmark. Of the thirteen paired graph-definition '
    'comparisons, one was significant before correction and none survived Holm correction. '
    '\u00a7 The sequence baseline was '
    'evaluated in two head widths over the same frozen embeddings. The narrower head (0.348\u20130.388 '
    'on S669, 0.458\u20130.493 on ssym) is a three-seed mean. The wider head (0.392 on S669, 0.516 '
    'on ssym) is the branch actually used inside the fusion comparison and is also a three-seed '
    'mean, computed by averaging the three per-seed predictions before correlating; it is therefore '
    'the appropriate reference for the structure increment.'
)


def main():
    doc = Document(DOC)
    old = doc.tables[0]

    # 用新表整体替换旧表（旧表可能有合并单元格，逐格操作不可靠）
    new = doc.add_table(rows=len(ROWS), cols=3)
    try:
        new.style = old.style
    except Exception:
        pass
    old._tbl.addprevious(new._tbl)
    old._tbl.getparent().remove(old._tbl)

    for ri, row in enumerate(ROWS):
        for ci, val in enumerate(row):
            cell = new.cell(ri, ci)
            cell.text = ''
            p = cell.paragraphs[0]
            run = p.add_run(val)
            run.font.size = Pt(7.5)
            if ri == 0:
                run.bold = True

    # 更新脚注
    n = 0
    for p in doc.paragraphs:
        if p.text.strip().startswith('Table 1 |'):
            if p.runs:
                p.runs[0].text = FOOTNOTE
                for r in p.runs[1:]:
                    r.text = ''
            n += 1
    doc.save(DOC)
    print(f'Table 1 重建完成：{len(ROWS)-1} 数据行 × 3 列；脚注更新 {n} 处')
    for r in ROWS:
        print('  | ' + ' | '.join(r))


if __name__ == '__main__':
    main()
