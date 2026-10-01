"""
fix_index_bug.py — 修复索引约定 bug（第 7 个问题）
====================================================
_pdb_res_idx 是"整链索引"（align_and_filter.py: idx = offset + pos），
而接触图节点是"按 offset 切片后"的索引。
正确节点索引 = _pdb_res_idx - offset。

本脚本给基准 CSV 增加 _node_idx 列（= _pdb_res_idx - offset），
并验证其残基与声称的 wt 一致。
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from contact_graph_defs import extract_residues
from paths import DATA, FIGURES, PUB_FIGURES, ROOT, SRC, TOOLS  # 路径集中管理

SRC_PATH = str(SRC)
DATA_PATH = str(DATA)
ROOT_PATH = str(ROOT)
FIG_PATH = str(FIGURES)
PUB_PATH = str(PUB_FIGURES)
TOOLS_PATH = str(TOOLS)




DATA = r'D:\GED_mutation\data'

for label_csv, align_csv, tag in [
    ('benchmarks_s669_clean.csv', 'alignment_map_s669.csv', 'S669'),
    ('benchmarks_ssym_clean.csv', 'alignment_map_ssym.csv', 'ssym'),
]:
    p = os.path.join(DATA, label_csv)
    df = pd.read_csv(p)
    align = pd.read_csv(os.path.join(DATA, align_csv))
    amap = {r['pdb_id']: r for _, r in align.iterrows()}

    node_idx, ok = [], []
    for _, row in df.iterrows():
        a = amap.get(row['pdb_id'])
        if a is None:
            node_idx.append(-1); ok.append(False); continue
        L = len(row['wt_seq'])
        off = int(a['offset'])
        ni = int(row['_pdb_res_idx']) - off
        node_idx.append(ni)
        wr, sw = extract_residues(os.path.join(DATA, 'structures', f'pdb{row["pdb_id"].lower()}.ent'),
                                  a['chain'], off, L)
        good = sw is not None and 0 <= ni < len(sw) and sw[ni] == str(row['mut_info'])[0]
        ok.append(good)

    df['_node_idx'] = node_idx
    df['_node_idx_ok'] = ok
    df.to_csv(p, index=False)
    print(f'{tag}: 加了 _node_idx 列；校验通过 {sum(ok)}/{len(df)}')
    if sum(ok) < len(df):
        bad = df[~df['_node_idx_ok']]
        print('   未通过样例（前 5）:')
        print(bad[['pdb_id', 'mut_info', '_pdb_res_idx', '_node_idx']].head(5).to_string(index=False))
