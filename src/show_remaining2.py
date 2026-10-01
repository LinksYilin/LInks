# -*- coding: utf-8 -*-
"""show_remaining2.py — 显示剩余核心脚本的硬编码行"""
import os
import re

RELSRC = r'D:\GED_mutation\release\src'
EDITORIAL = re.compile(
    r'^(add_|apply_|fix_|fill_|rebuild_table|insert_figures|clean_orphan|export_current|'
    r'audit_docx|consistency_sweep|number_audit|verify_|forensic_audit|validate_release|'
    r'check_hardcoded|classify_hardcoded|show_hardcoded|show_remaining|fix_release|audit_docx_health)')

for f in sorted(os.listdir(RELSRC)):
    if not f.endswith('.py') or EDITORIAL.match(f):
        continue
    p = os.path.join(RELSRC, f)
    lines = open(p, encoding='utf-8').read().split('\n')
    found = [(i + 1, l.strip()) for i, l in enumerate(lines) if re.search(r'[A-Za-z]:\\', l)]
    if not found:
        continue
    print(f'=== {f} ===')
    for i, l in found:
        print(f'  L{i}: {l[:108]}')
