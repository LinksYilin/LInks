"""quantify_index_bug.py — 量化索引 bug 的实际影响范围"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from contact_graph_defs import extract_residues

DATA = r'D:\GED_mutation\data'


def check(name, label_csv, align_csv):
    ss = pd.read_csv(os.path.join(DATA, label_csv))
    align = pd.read_csv(os.path.join(DATA, align_csv))
    amap = {r['pdb_id']: r for _, r in align.iterrows()}
    rows = []
    for _, row in ss.iterrows():
        pid, mut = row['pdb_id'], row['mut_info']
        a = amap.get(pid)
        if a is None:
            continue
        L = len(row['wt_seq'])
        off = int(a['offset'])
        idx = int(row['_pdb_res_idx'])
        node_idx = idx - off
        wr, sw = extract_residues(os.path.join(DATA, 'structures', f'pdb{pid.lower()}.ent'),
                                  a['chain'], off, L)
        if sw is None:
            continue
        direct_ok = idx < len(sw)
        node_ok = 0 <= node_idx < len(sw)
        rows.append({
            'pdb_id': pid, 'mut_info': mut, 'offset': off, 'chain_idx': idx,
            'node_idx': node_idx, 'seq_len': len(sw),
            'direct_in_range': direct_ok, 'node_in_range': node_ok,
            'direct_aa': sw[idx] if direct_ok else None,
            'node_aa': sw[node_idx] if node_ok else None,
            'claim': mut[0],
        })
    d = pd.DataFrame(rows)
    print(f'=== {name} ===')
    print(f'  总突变: {len(d)}')
    print(f'  offset>0 的突变: {(d["offset"] > 0).sum()}')
    aff = d[d['offset'] > 0]
    print(f'    其中"直接用 _pdb_res_idx"越界被丢弃的: {(~aff["direct_in_range"]).sum()}')
    print(f'    其中"直接用"未越界但残基不匹配的: {((aff["direct_in_range"]) & (aff["direct_aa"] != aff["claim"])).sum()}')
    print(f'  修正后（减 offset）可用: {d["node_in_range"].sum()}')
    print(f'  直接用可用: {d["direct_in_range"].sum()}')
    return d


d669 = check('S669', 'benchmarks_s669_clean.csv', 'alignment_map_s669.csv')
print()
dssym = check('ssym', 'benchmarks_ssym_clean.csv', 'alignment_map_ssym.csv')
print()
print('=== S669 受影响的突变（offset>0）===')
aff = d669[(d669['offset'] > 0)]
print(aff[['pdb_id', 'mut_info', 'offset', 'chain_idx', 'node_idx', 'seq_len',
           'direct_in_range', 'direct_aa', 'node_aa', 'claim']].to_string(index=False))
d669.to_csv(os.path.join(DATA, 'index_bug_audit_s669.csv'), index=False)
dssym.to_csv(os.path.join(DATA, 'index_bug_audit_ssym.csv'), index=False)
print('\n已保存 index_bug_audit_s669.csv / index_bug_audit_ssym.csv')
