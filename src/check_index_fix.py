# -*- coding: utf-8 -*-
"""check_index_fix.py — 核实两文件的索引修复"""
import os

for f in ['mechanism_analysis_corrected.py', 'locality_corrected.py']:
    p = os.path.join(r'D:\GED_mutation\src', f)
    t = open(p, encoding='utf-8').read()
    print(f'== {f} ==')
    for i, l in enumerate(t.split('\n'), 1):
        if '_node_idx' in l or '_pdb_res_idx' in l:
            print(f'  L{i}: {l.strip()[:110]}')
    print()
