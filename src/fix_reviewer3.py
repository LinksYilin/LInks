# -*- coding: utf-8 -*-
"""
fix_reviewer3.py — 按审稿人 3 的确认意见修正稿件
==================================================
经逐条核实，以下 3 项为真实错误：
  R3-M5  Figure 3 图注仍称 ssym 图未构建（实际已构建，Table 1 已有 ssym 列）
  R3-M2  Methods 称种子为 42/123/2024、three seeds（阶梯实际只用 42/123）
  R3-m8  SCWRL4 交叉引用指向 Section 4.3（应为 4.7）
其余指控经核实不成立（Table 1 zero-shot 为 0.317 非 0.217；配对表为 13 行非 12；
Table 1 脚注完整未截断）。
"""
from docx import Document

DOC = r'D:\GED_mutation\manuscript_routeA.docx'

FIX = [
    # --- R3-M5：Figure 3 图注过时表述 ---
    ('Encoder results are reported as the range across definitions in Table 1; the corresponding '
     'ssym graphs were not constructed for these encoders.',
     'Encoder results are reported as the range across definitions on S669 and as the '
     'side-chain-centroid value on ssym in Table 1.'),

    # --- R3-M2：种子数与实际一致（阶梯 2 个种子，序列模型 3 个） ---
    ('The GNNs were evaluated using random seeds 42, 123, and 2024.',
     'The capacity ladder used two seeds (42 and 123) because each rung required retraining under '
     'four graph definitions; the sequence-baseline and fusion models used three seeds (42, 123 and '
     '2024). Every reported value is the mean of the per-seed predictions, and every paired interval '
     'is computed on that seed-averaged prediction.'),

    ('across three random seeds. Because S669 guided the development',
     'using the seed set stated in Section 3.6. Because S669 guided the development'),

    # --- R3-m8：SCWRL4 交叉引用 ---
    ('A matched-scope SCWRL4 analysis was used to assess cross-engine reproducibility (Section 4.3).',
     'A matched-scope SCWRL4 analysis was used to assess cross-engine reproducibility (Section 4.7).'),
]


def main():
    d = Document(DOC)
    n = 0
    for p in d.paragraphs:
        if not p.runs:
            continue
        full = ''.join(r.text for r in p.runs)
        new = full
        for a, b in FIX:
            if a in new:
                new = new.replace(a, b)
                n += 1
        if new != full:
            p.runs[0].text = new
            for r in p.runs[1:]:
                r.text = ''
    d.save(DOC)
    print(f'修正 {n} 处')


if __name__ == '__main__':
    main()
