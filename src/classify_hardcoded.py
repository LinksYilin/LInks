# -*- coding: utf-8 -*-
"""classify_hardcoded.py — 分类硬编码路径的脚本：核心流程 vs 一次性编辑工具"""
import os
import re

REL = r'D:\GED_mutation\release'
SRC = os.path.join(REL, 'src')

EDITORIAL = re.compile(
    r'^(add_|apply_|fix_|fill_|rebuild_table|insert_figures|clean_orphan|export_current|'
    r'audit_docx|consistency_sweep|number_audit|verify_|forensic_audit|validate_release|'
    r'check_hardcoded|fix_release|audit_docx_health)')

rows = []
for f in sorted(os.listdir(SRC)):
    if not f.endswith('.py'):
        continue
    t = open(os.path.join(SRC, f), encoding='utf-8').read()
    n = len(re.findall(r'[A-Za-z]:\\', t))
    uses_paths = 'from paths import' in t or 'import paths' in t
    if n == 0:
        continue
    kind = 'editorial' if EDITORIAL.match(f) else 'PIPELINE'
    rows.append((kind, f, n, uses_paths))

core = [r for r in rows if r[0] == 'PIPELINE']
edit = [r for r in rows if r[0] == 'editorial']
print(f'硬编码脚本共 {len(rows)} 个')
print(f'  一次性编辑/核验工具: {len(edit)} 个（不影响复现）')
print(f'  核心流程脚本:       {len(core)} 个  ← 需要处理')
print()
if core:
    print('核心流程脚本明细：')
    for _, f, n, up in core:
        print(f'  {f:<42} {n:>2} 处硬编码  uses_paths={up}')
