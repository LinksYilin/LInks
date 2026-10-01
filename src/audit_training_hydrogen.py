# -*- coding: utf-8 -*-
"""
audit_training_hydrogen.py — 全量审计训练/测试图的氢一致性
============================================================
关键问题：megascale_sc（训练）疑似含氢，而 s669_sc/ssym_sc（测试）排氢
          → 训练/测试表示不一致。
本脚本对全部图逐个判定，输出准确清单。
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(__file__))
from contact_graph_defs import extract_residues, repr_coords

DATA = r'D:\GED_mutation\data'
STRUCT = os.path.join(DATA, 'raw', 'megascale_structures', 'AlphaFold_model_PDBs.parquet')
TMP = os.path.join(DATA, 'raw', '_tmp_audit.pdb')
BACKBONE = {'N', 'CA', 'C', 'O', 'OXT'}


def centroid_edges(pdb_path, exclude_h, cutoff=8.0):
    res, _ = extract_residues(pdb_path, None, 0, None)
    cents = []
    for r in res:
        atoms = [a for a in r.get_atoms() if a.get_name() not in BACKBONE
                 and ((a.element != 'H') if exclude_h else True)]
        cents.append(np.mean([a.get_coord() for a in atoms], axis=0) if atoms
                     else r['CA'].get_coord())
    cents = np.array(cents)
    tree = cKDTree(cents)
    return {tuple(sorted((int(a), int(b))))
            for a, b in tree.query_pairs(r=cutoff, output_type='ndarray')}


print('=' * 62)
print('1) MegaScale 训练图（全量 239）')
print('=' * 62)
struct_df = pd.read_parquet(STRUCT)
struct_map = {r['name']: r['pdb'] for _, r in struct_df.iterrows()}
sc_dir = os.path.join(DATA, 'contact_graphs_megascale_sc')
files = sorted(f for f in os.listdir(sc_dir) if f.endswith('.npz'))
stat = {'含氢': 0, '排氢': 0, '其他': 0, '无结构': 0}
for f in files:
    pid = f[:-4]
    z = np.load(os.path.join(sc_dir, f), allow_pickle=True)
    sc = {tuple(sorted((int(a), int(b)))) for a, b in z['edge_index'].T}
    cand = [k for k in struct_map if k.replace('/', '_').replace('.pdb', '') == pid]
    if not cand:
        stat['无结构'] += 1
        continue
    with open(TMP, 'w') as fh:
        fh.write(struct_map[cand[0]])
    try:
        e_noH = centroid_edges(TMP, True)
        e_H = centroid_edges(TMP, False)
    finally:
        if os.path.exists(TMP):
            os.remove(TMP)
    if e_H == sc:
        stat['含氢'] += 1
    elif e_noH == sc:
        stat['排氢'] += 1
    else:
        stat['其他'] += 1
print(f'  含氢 {stat["含氢"]}, 排氢 {stat["排氢"]}, 其他 {stat["其他"]}, 无结构 {stat["无结构"]}')
print(f'  → megascale 结论: {"含氢（旧定义，需修复）" if stat["含氢"] > stat["排氢"] else "排氢"}')

print()
print('=' * 62)
print('2) 其他数据集：与各自 definitions 的一致性')
print('=' * 62)
for name in ['s669', 'ssym', 'thermomutdb']:
    sc = os.path.join(DATA, f'contact_graphs_{name}_sc')
    cent = os.path.join(DATA, f'contact_graphs_{name}_centroid')
    if not (os.path.isdir(sc) and os.path.isdir(cent)):
        print(f'  {name}: 目录缺失')
        continue
    same = diff = 0
    for f in os.listdir(sc):
        if not f.endswith('.npz'):
            continue
        a = np.load(os.path.join(sc, f), allow_pickle=True)
        b = np.load(os.path.join(cent, f), allow_pickle=True)
        ea = {tuple(sorted((int(x), int(y)))) for x, y in a['edge_index'].T}
        eb = {tuple(sorted((int(x), int(y)))) for x, y in b['edge_index'].T}
        if ea == eb:
            same += 1
        else:
            diff += 1
    verdict = '排氢（与新定义一致）' if diff == 0 else f'不一致 {diff} 个'
    print(f'  {name}_sc vs _centroid(排氢): 一致 {same}, 不同 {diff}  → {verdict}')

print()
print('=' * 62)
print('3) ca / cb 是否受氢影响')
print('=' * 62)
print('  Cα / Cβ 是单个重原子坐标，不受氢影响 → 无需检查')
