# -*- coding: utf-8 -*-
"""fix_verifier_1.py — 修正独立验证者报告的第 1-5 项（§4.2 核心数字）

1. 平均|Δr| 0.066/0.048/0.019/0.043 → 共有对 0.048/0.048/0.021/0.043
2. 删除过时的重复趋势段（0.071/0.054/0.025/0.039, 20000 次置换）
3. EGNN 平均效应 0.061 → 0.059
4. "两种子" → "三种子"
5. deep_gine "Cβ 0.362" → 0.331
"""
import os
import re

import pandas as pd
from docx import Document

DOC = r'D:\GED_mutation\manuscript_routeA.docx'
D = r'D:\GED_mutation\data'
DEL = '\u0394'
A = '\u03b1'
B = '\u03b2'
RHO = '\u03c1'

eff = pd.read_csv(os.path.join(D, 'ladder_paired_effects_audited.csv'))
COMMON = ['ca vs centroid', 'cb vs centroid']
means = (eff[eff.comparison.isin(COMMON)].groupby('model').delta_r
         .apply(lambda v: v.abs().mean()).to_dict())
print('共有对平均|Δr|:', {k: round(v, 4) for k, v in means.items()})

g, l, e, dg = (means['gnn_global'], means['gnn_local'],
               means['gnn_edge'], means['deep_gine'])
eg = float(eff[eff.model == 'egnn'].delta_r.abs().mean())

RULES = [
    # 1) 平均 |Δr| 数值
    ('the mean absolute effect was 0.066, 0.048, 0.019 and 0.043 across the four '
     'stable rungs,',
     f'the mean absolute effect was {g:.3f}, {l:.3f}, {e:.3f} and {dg:.3f} across the four '
     'stable rungs,'),
    ('the mean absolute effect was 0.066, 0.048, 0.019 and 0.043',
     f'the mean absolute effect was {g:.3f}, {l:.3f}, {e:.3f} and {dg:.3f}'),
    # 3) EGNN 平均效应
    ('The 509 k-parameter equivariant network gave a mean absolute effect of 0.061',
     f'The 509 k-parameter equivariant network gave a mean absolute effect of {eg:.3f}'),
    # 5) deep_gine Cβ
    ('whereas the 460 k-parameter GINE reached 0.362 on the same definition',
     'whereas the 460 k-parameter GINE reached 0.331 on the same definition'),
]

# 2) 要删除的整段（过时重复趋势）
STALE_MARK = '0.071 (5.8 k), 0.054 (5.9 k), 0.025 (50 k) and 0.039 (460 k)'


def main():
    doc = Document(DOC)
    n_num, n_del = 0, 0
    for p in doc.paragraphs:
        if not p.runs:
            continue
        full = ''.join(r.text for r in p.runs)
        # 2) 整段删除
        if STALE_MARK in full and 'permutation' in full.lower():
            for r in list(p.runs):
                r._r.getparent().remove(r._r)
            n_del += 1
            print('  删除过时趋势段:', full[:80] + '...')
            continue
        new = full
        for a, b in RULES:
            if a in new:
                new = new.replace(a, b)
                n_num += 1
        if new != full:
            p.runs[0].text = new
            for r in p.runs[1:]:
                r.text = ''
    doc.save(DOC)
    print(f'\n数值替换 {n_num} 处; 删除段落 {n_del} 处')

    # 4) 两种子 → 三种子
    d2 = Document(DOC)
    pat = re.compile(r'(capacity ladder used two seeds|same optimisation settings and two seeds'
                     r'|two seeds per rung|two-seed bounds)')
    hits = [(i, p.text) for i, p in enumerate(d2.paragraphs) if pat.search(p.text)]
    print(f'\n仍称"两种子"的段落: {len(hits)}')
    for i, t in hits:
        for m in re.finditer(pat, t):
            s = max(0, m.start() - 90)
            print(f'  [{i}] …{t[s:m.end()+90]}…')
            print()


if __name__ == '__main__':
    main()
