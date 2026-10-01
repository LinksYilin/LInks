"""
diag_index_convention.py — 核查 _pdb_res_idx 的索引约定
=======================================================
问题：_pdb_res_idx 是"整条链"索引还是"切片后"索引？
若不一致，GNN 评估会把突变特征写到错误的残基上，或直接跳过。
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from contact_graph_defs import extract_residues

from pathlib import Path as _P
DATA = str(_P(__file__).resolve().parent.parent / 'data')
ss = pd.read_csv(os.path.join(DATA, 'benchmarks_s669_clean.csv'))
align = pd.read_csv(os.path.join(DATA, 'alignment_map_s669.csv'))
amap = {r['pdb_id']: r for _, r in align.iterrows()}

n_total = 0
n_in_range_sliced = 0
n_ok_minus_offset = 0
mismatch = []

for _, row in ss.iterrows():
    pid, mut = row['pdb_id'], row['mut_info']
    a = amap.get(pid)
    if a is None:
        continue
    L = len(row['wt_seq'])
    off = int(a['offset'])
    idx = int(row['_pdb_res_idx'])
    n_total += 1

    # 切片后序列
    wr_s, sw_s = extract_residues(os.path.join(DATA, 'structures', f'pdb{pid.lower()}.ent'),
                                  a['chain'], off, L)
    if sw_s is None:
        continue
    if idx < len(sw_s):
        n_in_range_sliced += 1
        aa_sliced = sw_s[idx]
    else:
        aa_sliced = None
    # 减去 offset
    j = idx - off
    aa_minus = sw_s[j] if 0 <= j < len(sw_s) else None

    wt_claim = mut[0]
    ok_sliced = (aa_sliced == wt_claim)
    ok_minus = (aa_minus == wt_claim)
    if not ok_sliced and not ok_minus:
        mismatch.append((pid, mut, idx, off, L, aa_sliced, aa_minus))

print(f'总突变: {n_total}')
print(f'_pdb_res_idx 直接作为切片索引时"在范围内"的: {n_in_range_sliced}')
print(f'两种解释都与声称 wt 不符的: {len(mismatch)}')
print()
print('=== 按 offset 分组：直接索引是否匹配声称的 wt 氨基酸 ===')
n_direct_ok = 0
n_minus_ok = 0
n_direct_inrange = 0
for _, row in ss.iterrows():
    pid, mut = row['pdb_id'], row['mut_info']
    a = amap.get(pid)
    if a is None:
        continue
    L = len(row['wt_seq'])
    off = int(a['offset'])
    idx = int(row['_pdb_res_idx'])
    wr_s, sw_s = extract_residues(os.path.join(DATA, 'structures', f'pdb{pid.lower()}.ent'),
                                  a['chain'], off, L)
    if sw_s is None:
        continue
    wt_claim = mut[0]
    if idx < len(sw_s):
        n_direct_inrange += 1
        if sw_s[idx] == wt_claim:
            n_direct_ok += 1
    j = idx - off
    if 0 <= j < len(sw_s) and sw_s[j] == wt_claim:
        n_minus_ok += 1

print(f'  _pdb_res_idx 直接用（切片后序列）匹配 wt 的: {n_direct_ok}')
print(f'  _pdb_res_idx - offset 匹配 wt 的: {n_minus_ok}')
print()
print('=== 按 offset>0 的子集 ===')
sub = ss[ss['pdb_id'].map(lambda p: int(amap[p]['offset']) if p in amap else 0) > 0]
print(f'  offset>0 的突变数: {len(sub)} / {len(ss)}')
d_ok = m_ok = 0
for _, row in sub.iterrows():
    pid, mut = row['pdb_id'], row['mut_info']
    a = amap[pid]
    L = len(row['wt_seq'])
    off = int(a['offset']); idx = int(row['_pdb_res_idx'])
    wr_s, sw_s = extract_residues(os.path.join(DATA, 'structures', f'pdb{pid.lower()}.ent'),
                                  a['chain'], off, L)
    if sw_s is None:
        continue
    if idx < len(sw_s) and sw_s[idx] == mut[0]:
        d_ok += 1
    j = idx - off
    if 0 <= j < len(sw_s) and sw_s[j] == mut[0]:
        m_ok += 1
print(f'  直接索引匹配: {d_ok}')
print(f'  减去 offset 匹配: {m_ok}')
