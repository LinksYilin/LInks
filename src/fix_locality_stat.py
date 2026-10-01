# -*- coding: utf-8 -*-
"""fix_locality_stat.py — 用排除自身接触后的正确局域性统计替换

原统计 mean_d_broken 取 min(端点到突变中心距离)，对与突变残基相连的边恒为 0，
358 个突变中 253 个（70.7%）含自身边（共 382/1068 条断边），把均值拉到 4.03 Å。

排除自身边后（n=234）：
  断裂 7.93 Å  [6.63, 10.11]；保持 15.93 Å [13.93, 18.86]
  逐边加权 11.93 Å；中位数 5.88 Å；<5 Å 占 30.1%
定性结论不变（断裂接触仍显著更靠近突变位点）。
"""
from docx import Document

DOC = r'D:\GED_mutation\manuscript_routeA.docx'

RULES = [
    # §3.3 方法描述：如实说明统计口径
    ('Across 358 mutation pairs with at least one broken contact, broken contacts lay at a mean '
     'distance of 4.03 \u00c5 from the mutated residue (95% protein-cluster bootstrap CI '
     '3.11\u20135.42 \u00c5), whereas unchanged contacts lay at 15.70 \u00c5 (13.95\u201318.14 '
     '\u00c5); the median distance for a broken contact was 2.70 \u00c5 and 65.8% lay within 5 '
     '\u00c5.',
     'For each mutation we measured the distance from the mutated residue to each contact it '
     'lost, taking the nearer endpoint of the contact. A contact incident on the mutated residue '
     'itself is therefore at distance zero; such self-contacts account for 382 of the 1,068 '
     'broken contacts and appear for 253 of the 358 mutation pairs. Excluding them, and averaging '
     'over the 234 mutation pairs that retain at least one non-self broken contact, broken '
     'contacts lay at a mean distance of 7.93 \u00c5 from the mutated residue (95% '
     'protein-cluster bootstrap CI 6.63\u201310.11 \u00c5), whereas unchanged contacts lay at '
     '15.93 \u00c5 (13.93\u201318.86 \u00c5); the median distance for a broken contact was 5.88 '
     '\u00c5 and 30.1% lay within 5 \u00c5. Counting each broken contact once rather than each '
     'mutation once gives 11.93 \u00c5.'),

    # 贡献列表
    ('broken contacts sat on average 4.03 \u00c5 from the mutated residue versus 15.70 \u00c5 for '
     'contacts that persisted',
     'broken contacts sat on average 7.9 \u00c5 from the mutated residue versus 15.9 \u00c5 for '
     'contacts that persisted, after excluding contacts incident on the mutated residue itself'),

    # Discussion 5.1
    ('Broken contacts were strongly localised (4.03 \u00c5 versus 15.70 \u00c5 for unchanged '
     'contacts)',
     'Broken contacts were strongly localised (7.9 \u00c5 versus 15.9 \u00c5 for unchanged '
     'contacts)'),
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
    for pat, name in [('4.03', '残留 4.03'), ('15.70', '残留 15.70'),
                      ('7.93', '新 7.93'), ('15.93', '新 15.93'),
                      ('382 of the 1,068', '自身边说明')]:
        print(f'  {name}: {txt.count(pat)}')


if __name__ == '__main__':
    main()
