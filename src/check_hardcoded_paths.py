# -*- coding: utf-8 -*-
"""check_hardcoded_paths.py — 检查发布脚本中的硬编码盘符路径"""
import os
import re

REL = r'D:\GED_mutation\release'
hits = {}
for dp, _, fs in os.walk(os.path.join(REL, 'src')):
    for f in fs:
        if not f.endswith('.py'):
            continue
        p = os.path.join(dp, f)
        t = open(p, encoding='utf-8').read()
        for m in re.finditer(r'[A-Za-z]:\\', t):
            line = t[:m.start()].count('\n') + 1
            src = t.split('\n')[line - 1].strip()[:110]
            hits.setdefault(f, []).append((line, src))

print(f'含盘符路径的脚本: {len(hits)} 个 / 共 {sum(len(v) for v in hits.values())} 处')
print()
grouped = {}
for f, v in hits.items():
    for line, src in v:
        key = 'paths.py (路径集中管理)' if f == 'paths.py' else '其它脚本'
        grouped.setdefault(key, []).append((f, line, src))
for key, items in grouped.items():
    print(f'--- {key}: {len(items)} 处 ---')
    for f, line, src in items[:10]:
        print(f'  {f}:{line}  {src}')
    print()
