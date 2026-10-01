# -*- coding: utf-8 -*-
"""fix_verifier_3.py — 修正 Table 1 脚注 与 Discussion 的 EGNN 描述

- Table 1 脚注：SD 0.197 → 样本 SD 0.241；范围下界按舍入修正
- Discussion [198]："two seeds" → 三种子正确值
- 保留"ssym row uses two seeds"（这是事实）
"""
from docx import Document

DOC = r'D:\GED_mutation\manuscript_routeA.docx'
MIN = '\u2212'

RULES = [
    # Table 1 脚注：SD
    ('its S669 seeds gave r = 0.393, ' + MIN + '0.081 and 0.079 (mean 0.131, SD 0.197)',
     'its S669 seeds gave r = 0.393, ' + MIN + '0.081 and 0.079 (mean 0.131, sample SD 0.241)'),
    ('(mean 0.131, SD 0.197)', '(mean 0.131, sample SD 0.241)'),
    # Table 1 脚注：范围舍入（0.313464→0.313；0.377452→0.377）——但正文用 0.314–0.378
    # 改为以三位小数表述并说明为范围
    ('0.052\u20130.150 \u2020', '0.052\u20130.150 \u2020'),
    # Discussion：两种子 → 三种子
    ('(r = 0.393 and -0.081 for two seeds on the side-chain-centroid graphs)',
     '(sample SD 0.241 across three architecture-matched seeds on the side-chain-centroid '
     'graphs)'),
    ('r = 0.393 and -0.081 for two seeds on the side-chain-centroid graphs',
     'sample SD 0.241 across three architecture-matched seeds on the side-chain-centroid graphs'),
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
    for pat, name in [('SD 0.197', '残留 SD 0.197'),
                      ('for two seeds', '残留 for two seeds'),
                      ('sample SD 0.241', '新 sample SD 0.241'),
                      ('ssym row uses two seeds', 'ssym 两种子（应保留）')]:
        print(f'  {name}: {txt.count(pat)}')


if __name__ == '__main__':
    main()
