# -*- coding: utf-8 -*-
"""regen_skip_log.py — 从实际最终集合重新生成审计日志

原日志把 11 对标记为 excluded，但它们实际都在分析集中；1G3P 被标为
"无接触图" 却在 505 集中。改为按真实集合成员资格标注，并保留原因列作为
历史说明。
"""
import os

import pandas as pd

D = r'D:\GED_mutation\data'
ROOT = r'D:\GED_mutation'
log = pd.read_csv(os.path.join(D, 'data_flow_skip_log.csv'))

# 实际集合
edits = pd.read_csv(os.path.join(D, 'edits_corrected.csv'))
e8 = edits[edits.threshold == 8.0]
set505 = {(str(a), str(b)) for a, b in zip(e8.pdb_id, e8.mut_info)}

lad = pd.read_csv(os.path.join(D, 'ladder_predictions.csv'))
s669 = lad[lad.benchmark == 's669']
set511 = {(str(a), str(b)) for a, b in zip(s669.protein_id, s669.mutation_id)}

bench = pd.read_csv(os.path.join(D, 'benchmarks_s669_clean.csv'))
set543 = {(str(a), str(b)) for a, b in zip(bench.pdb_id, bench.mut_info)}

log['key'] = list(zip(log.pdb_id.astype(str), log.mutation.astype(str)))
log['in_543'] = [k in set543 for k in log.key]
log['in_511_predictive'] = [k in set511 for k in log.key]
log['in_505_contact_change'] = [k in set505 for k in log.key]
log['in_final_gnn_analysis'] = log['in_511_predictive'] | log['in_505_contact_change']
log['reason_recorded_historically'] = log['reason']

out = log[['pdb_id', 'mutation', 'reason_recorded_historically', 'in_543',
           'in_511_predictive', 'in_505_contact_change', 'in_final_gnn_analysis']]
out = out.sort_values(['in_final_gnn_analysis', 'pdb_id', 'mutation'],
                      ascending=[False, True, True])

path = os.path.join(D, 'data_flow_skip_log.csv')
out.to_csv(path, index=False)
out.to_csv(os.path.join(ROOT, 'data', 'data_flow_skip_log.csv'), index=False)

n_in = int(out.in_final_gnn_analysis.sum())
n_out = int((~out.in_final_gnn_analysis).sum())
print(f'重新生成: {len(out)} 行')
print(f'  实际在分析集中: {n_in}')
print(f'  实际排除:       {n_out}')
print()
print('被纠正的条目（原标为排除，实为在集合中）:')
chk = out[out.in_final_gnn_analysis]
sub = chk[chk.reason_recorded_historically.astype(str).str.contains('越界|无接触图', na=False)]
for _, r in sub.iterrows():
    print(f'  {r.pdb_id} {r.mutation}: {r.reason_recorded_historically}  '
          f'(511={r.in_511_predictive}, 505={r.in_505_contact_change})')
