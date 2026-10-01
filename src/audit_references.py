# -*- coding: utf-8 -*-
"""audit_references.py — 检查参考文献的完整性与正文引用一致性"""
import re

from docx import Document

DOC = r'D:\GED_mutation\manuscript_routeA.docx'
d = Document(DOC)
paras = [p.text.strip() for p in d.paragraphs]
body_txt = '\n'.join(paras)

# 抓取参考文献列表
refs = {}
for t in paras:
    m = re.match(r'^\[(\d+)\]\s+(.+)$', t)
    if m:
        refs[int(m.group(1))] = m.group(2)

print(f'参考文献条目: {len(refs)} 条')
print()

# 正文引用（排除参考文献列表本身）
cite_region = body_txt
for n, txt in refs.items():
    cite_region = cite_region.replace(f'[{n}] {txt}', '')
cited = set()
for m in re.finditer(r'\[(\d+)(?:[,–\-](\d+))?\]', cite_region):
    a = int(m.group(1))
    cited.add(a)
    if m.group(2):
        cited.add(int(m.group(2)))

print(f'正文引用编号: {sorted(cited)}')
missing = sorted(set(refs) - cited)
extra = sorted(cited - set(refs))
print(f'  列表有但正文未引用: {missing if missing else "无"}')
print(f'  正文引用但列表缺失: {extra if extra else "无"}')
print()

# 检查每条参考文献的字段完整性
print('=== 参考文献字段检查 ===')
issues = []
for n in sorted(refs):
    t = refs[n]
    has_year = bool(re.search(r'\((19|20)\d{2}\)', t))
    has_title = len(t) > 60
    has_venue = bool(re.search(r'\.\s+(J\.|Nucleic|Nature|Science|Bioinformatics|Proc\.|Cell|'
                               r'Protein|Bioinformatics|arXiv|bioRxiv|Commun\.|IEEE|ACM)', t))
    tag = []
    if not has_year:
        tag.append('无年份')
    if not has_title:
        tag.append('标题过短')
    if not has_venue:
        tag.append('未见期刊/会议')
    if tag:
        issues.append((n, tag, t[:90]))
    print(f'  [{n}] {"✅" if not tag else "⚠ " + ",".join(tag)}  {t[:88]}')

print()
print(f'字段问题: {len(issues)} 条')
for n, tag, t in issues:
    print(f'  [{n}] {tag}: {t}')

# 检查引用顺序（Vancouver 要求按出现顺序编号）
order = []
for m in re.finditer(r'\[(\d+)\]', cite_region):
    v = int(m.group(1))
    if v not in order:
        order.append(v)
print()
print(f'首次出现顺序: {order}')
if order != sorted(order):
    print('  ⚠ 引用顺序不是升序（Vancouver 要求按首次出现编号）')
else:
    print('  ✅ 编号顺序正确')
