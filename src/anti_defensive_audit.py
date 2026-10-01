# -*- coding: utf-8 -*-
"""
anti_defensive_audit.py — 论文的防御性写作审计
================================================
按 anti-defensive-writing skill 的检查清单扫描稿件，列出：
  - 不必要的免责声明
  - 反复说明"本文并不是说……"
  - 叠加的保留措辞
  - 摘要/段首/结论中堆放的限制条件
  - 段落未说观点先强调不足
  - 多余的"然而/尽管如此/虽然"转折

输出：逐条位置 + 建议改法
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(__file__))
from paths import ROOT

DOCX = os.path.join(str(ROOT), 'manuscript_polished_final_20260930.docx')

PATTERNS = [
    ('不必要的免责声明', r'\b(we do not claim|we do not attempt|this paper does not|'
                         r'is not intended to|does not aim to|we make no claim)\b'),
    ('"并不是说"结构', r'\b(this is not to say|this should not be taken|'
                     r'not X but Y|rather than arguing)\b'),
    ('叠加保留措辞', r'\b(may potentially|could potentially|might possibly|'
                   r'it is possible that|to some extent may)\b'),
    ('值得注意类', r'\b(it is worth noting|it should be noted|to be clear|'
                 r'notably, it is|it must be emphasized)\b'),
    ('让步转折堆叠', r'\b(although .{0,60} however|while .{0,60} nevertheless)\b'),
    ('局限于段首', r'^(However|Nevertheless|Although|While|Despite)\b'),
    ('笼统自我贬低', r'\b(only a preliminary|limited reference value|'
                   r'far from sufficient|significant room for improvement)\b'),
]


def main():
    from docx import Document
    d = Document(DOCX)
    paras = [p.text.strip() for p in d.paragraphs]
    full = '\n'.join(paras)

    print('=' * 76)
    print('防御性写作审计')
    print('=' * 76)
    total = 0
    for name, pat in PATTERNS:
        hits = []
        for i, t in enumerate(paras):
            for m in re.finditer(pat, t, re.I):
                hits.append((i, m.group(0), t[max(0, m.start() - 40):m.end() + 60]))
        if hits:
            print(f'\n### {name}: {len(hits)} 处')
            for i, frag, ctx in hits[:8]:
                print(f'  [段{i}] 「{frag}」')
                print(f'        …{ctx}…')
        total += len(hits)

    print()
    print('=' * 76)
    print(f'合计 {total} 处可疑表述')

    # 摘要中是否堆限制
    print()
    print('=== 摘要中的限定词密度 ===')
    abs_txt = ' '.join(t for t in paras if t.startswith(('Motivation.', 'Results.', 'Conclusion.',
                                                         'Availability.')))
    n_words = len(abs_txt.split())
    hedge = len(re.findall(r'\b(may|might|could|suggest|appear|likely|possibly|'
                           r'exploratory|preliminary|not claim)\b', abs_txt, re.I))
    print(f'  摘要 {n_words} 词, 保留/弱化词 {hedge} 个, 密度 {hedge/max(n_words,1)*100:.1f}%')
    if hedge / max(n_words, 1) > 0.03:
        print('  ⚠️ 密度偏高，建议收敛（摘要只保留会影响理解的关键限定）')
    else:
        print('  ✅ 密度可接受')

    print()
    print('=== 结论段中的限制堆叠 ===')
    for i, t in enumerate(paras):
        if t.startswith('Conclusion.'):
            n = len(re.findall(r'\b(however|although|while|but|yet|limitation|'
                               r'boundary|caveat|not)\b', t, re.I))
            print(f'  段{i}: {len(t.split())} 词, 转折/否定词 {n} 个')
            print(f'    {t[:220]}…')
            break


if __name__ == '__main__':
    main()
