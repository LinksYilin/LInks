# -*- coding: utf-8 -*-
"""fix_verifier_2.py — 修正独立验证者报告的第 3,5,6,7,8,10,11,12,16 项"""
from docx import Document

DOC = r'D:\GED_mutation\manuscript_routeA.docx'
DEL = '\u0394'
A = '\u03b1'
B = '\u03b2'
MIN = '\u2212'

RULES = [
    # 3) 两种子 → 三种子（§3.6）
    ('The capacity ladder used two seeds (42 and 123) because each rung required retraining '
     'under four graph definitions;',
     'The capacity ladder used three seeds (42, 123 and 2024) for the four stable encoders and '
     'an architecture-matched three-seed set for the equivariant encoder, because each rung '
     'required retraining under four graph definitions;'),
    ('the sequence-baseline and fusion models used three seeds (42, 123 and 2024).',
     'the sequence-baseline and fusion models used the same three seeds (42, 123 and 2024).'),

    # 3) §5.5 两种子局限 → 改为真实的局限
    ('Sixth, the capacity ladder used two seeds per rung, so its paired intervals quantify '
     'sampling uncertainty over proteins but not optimisation variance.',
     'Sixth, the capacity ladder used three seeds per rung, so its joint intervals quantify '
     'both protein sampling and inter-seed optimisation variance; the seed and protein '
     'components are resampled together rather than separately.'),
    ('We report the two-seed bounds as the resolution currently achieved rather than as a limit '
     'of the design.',
     'We report the joint seed and protein bounds as the resolution currently achieved rather '
     'than as a limit of the design.'),

    # 5) 甘氨酸率（505 对，8 Å）
    ('74.6% of glycine-involving pairs versus 6.8% of others',
     '72.4% of glycine-involving pairs versus 1.8% of others'),

    # 8) 氢修复句（538 对）
    ('for the 538 modelled pairs, this inflated the mean number of broken contacts from 2.1 to '
     '26.3.',
     'for the 538 modelled pairs, this inflated the mean number of broken contacts from 9.5 to '
     '26.3.'),

    # 12) ssym ridge 0.306 → 0.305
    ('the physicochemical baseline on ssym reached r = 0.306',
     'the physicochemical baseline on ssym reached r = 0.305'),
    ('reached r = 0.306 on ssym', 'reached r = 0.305 on ssym'),
    ('physicochemical baseline (r = 0.306)', 'physicochemical baseline (r = 0.305)'),

    # 16) 训练集蛋白来源数
    ('420 proteins (232 from MegaScale and 188 from ThermoMutDB)',
     '420 distinct proteins (232 sampled from MegaScale and 191 from ThermoMutDB, with three '
     'proteins present in both sources)'),
    ('232 from MegaScale and 188 from ThermoMutDB',
     '232 sampled from MegaScale and 191 from ThermoMutDB'),

    # 10) Figure 1 描述错误
    ('Table 1 summarises the mutation-level performance; Figure 1 shows seed-specific '
     'uncertainty and within-protein ranking.',
     'Table 1 summarises the mutation-level performance, and Figure 1 shows the study design '
     'and evaluation protocol.'),
    ('Figure 1 shows seed-specific uncertainty and within-protein ranking',
     'Figure 1 shows the study design and evaluation protocol'),
]


def main():
    d = Document(DOC)
    n = 0
    applied = []
    for p in d.paragraphs:
        if not p.runs:
            continue
        full = ''.join(r.text for r in p.runs)
        new = full
        for a, b in RULES:
            if a in new:
                new = new.replace(a, b)
                applied.append(a[:52])
        if new != full:
            p.runs[0].text = new
            for r in p.runs[1:]:
                r.text = ''
            n += 1
    d.save(DOC)
    print(f'修改 {n} 段; 命中 {len(applied)} 条规则:')
    for x in applied:
        print(f'  ✓ {x}')

    # 复核残留
    d2 = Document(DOC)
    txt = '\n'.join(p.text for p in d2.paragraphs)
    print()
    for pat, name in [('two seeds', '残留"two seeds"'),
                      ('74.6%', '残留 74.6%'),
                      ('from 2.1 to 26.3', '残留 2.1→26.3'),
                      ('r = 0.306', '残留 0.306'),
                      ('188 from ThermoMutDB', '残留 188'),
                      ('seed-specific uncertainty', '残留 Figure1 误述')]:
        print(f'  {name}: {txt.count(pat)}')


if __name__ == '__main__':
    main()
