# -*- coding: utf-8 -*-
"""check_skip_log.py — 检查审计日志与实际分析集的一致性"""
import os

import pandas as pd

D = r'D:\GED_mutation\data'
log = pd.read_csv(os.path.join(D, 'data_flow_skip_log.csv'))
print('列:', log.columns.tolist())
print('行数:', len(log))
print()
print(log.head(20).to_string(index=False))

# 实际分析集
edits = pd.read_csv(os.path.join(D, 'edits_corrected.csv'))
e8 = edits[edits.threshold == 8.0]
set505 = set(zip(e8.pdb_id.astype(str), e8.mut_info.astype(str)))
print(f'\n505 集: {len(set505)} 对')

lad = pd.read_csv(os.path.join(D, 'ladder_predictions.csv'))
set511 = set(zip(lad[lad.benchmark == 's669'].protein_id.astype(str),
                 lad[lad.benchmark == 's669'].mutation_id.astype(str)))
print(f'511 集: {len(set511)} 对')

# 日志中标记为排除但仍在分析集中的
if 'in_final_gnn_analysis' in log.columns:
    ex = log[log.in_final_gnn_analysis.astype(str).str.lower() == 'false']
    print(f'\n日志标记为排除: {len(ex)} 行')
    cols = [c for c in log.columns if 'pdb' in c.lower() or 'mut' in c.lower()]
    print('  相关列:', cols)
    if len(cols) >= 2:
        bad = []
        for _, r in ex.iterrows():
            k = (str(r[cols[0]]), str(r[cols[1]]))
            if k in set505 or k in set511:
                bad.append(k)
        print(f'  其中仍在 505/511 集中的: {len(bad)}')
        for k in bad[:12]:
            print(f'    {k}')
