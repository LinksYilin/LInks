"""
build_megascale_graphs.py — MegaScale 训练数据接触图构建
========================================================
MegaScale 的 WT 结构在 AlphaFold_model_PDBs.parquet（name -> PDB 文本）。
对训练集用到的 WT（WT_name 字段），展开 PDB 文本，构建接触图，
并与训练集的 aa_seq（或 aa_seq_full）对齐校验。

用法：
  python build_megascale_graphs.py \
    --struct_parquet data/raw/megascale_structures/AlphaFold_model_PDBs.parquet \
    --train_parquet data/raw/megascale/train-00000-of-00003.parquet \
    --out_dir data/contact_graphs_megascale
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

from build_contact_graph import build_graph


def build_graph_from_pdb_text(pdb_text, name):
    """把 PDB 文本写到临时文件，构建接触图。"""
    tmp = os.path.join(os.environ.get('TEMP', '.'), f'_megascale_{name.replace("/", "_")}.pdb')
    with open(tmp, 'w') as f:
        f.write(pdb_text)
    try:
        return build_graph(tmp, chain_id=None)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def align_ref_to_built(ref_seq, built_seq):
    """把 ref_seq 对齐到 built_seq。

    经诊断，MegaScale 的 83 个 WT 全部满足：
      - ref 与 built 长度相同
      - 差异仅在位置 0（ref[0] 在 AlphaFold 结构里被预测成别的残基）
      - ref[1:] == built[1:]
    这是已知的 N 端首残基预测不确定性。因此对齐规则：
      - 若 ref == built → offset=0, 完全一致
      - 若 len 相同且 ref[1:] == built[1:] → offset=1, 首残基忽略
    返回 built 上可靠对齐的起始索引（0 或 1）。
    """
    if ref_seq == built_seq:
        return 0
    if len(ref_seq) == len(built_seq) and ref_seq[1:] == built_seq[1:]:
        return 1  # 首残基差异，从位置 1 起一致
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--struct_parquet', required=True)
    ap.add_argument('--train_parquet', required=True)
    ap.add_argument('--out_dir', required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    struct_df = pd.read_parquet(args.struct_parquet)
    struct_map = {r['name']: r['pdb'] for _, r in struct_df.iterrows()}

    train_df = pd.read_parquet(args.train_parquet, columns=['WT_name', 'aa_seq', 'aa_seq_full'])
    # 每个 WT 取第一条的 aa_seq 作为参考序列。
    # 注意：aa_seq 是核心折叠域序列（结构 PDB 对应它），aa_seq_full 含 His 标签/linker。
    wt_ref = {}
    for _, r in train_df.iterrows():
        wt_ref.setdefault(r['WT_name'], r['aa_seq'])

    ok, fail = [], []
    for wt_name, ref_seq in sorted(wt_ref.items()):
        if wt_name not in struct_map:
            fail.append((wt_name, '无结构'))
            continue
        pdb_text = struct_map[wt_name]
        g = build_graph_from_pdb_text(pdb_text, wt_name)
        if g is None:
            fail.append((wt_name, '构建失败'))
            continue
        nodes, edge_index, edge_attr, meta = g
        built_seq = meta['seq']
        # 对齐：结构可能 N/C 端缺残基（如 ref 首残基 Q 在 PDB 里缺失）
        off = align_ref_to_built(ref_seq, built_seq)
        if off is None:
            fail.append((wt_name, f'无法对齐 ref={len(ref_seq)} vs built={len(built_seq)}'))
            continue
        # 记录对齐信息到 meta
        meta['ref_len'] = len(ref_seq)
        meta['built_len'] = len(built_seq)
        meta['align_offset'] = off
        if off == 0 and built_seq == ref_seq:
            status = '一致'
        elif off >= 0:
            status = f'ref是built子串(offset={off})'
        else:
            status = f'ref需左移{abs(off)}位对齐'
        ok.append((wt_name, meta['n_residues'], meta['n_contacts'], status))
        # 保存
        out_path = os.path.join(args.out_dir, wt_name.replace('/', '_').replace('.pdb', '') + '.npz')
        np.savez(out_path, nodes=nodes, edge_index=edge_index, edge_attr=edge_attr,
                 metadata=np.array([json.dumps(meta)]))

    print(f'成功: {len(ok)}/{len(wt_ref)}')
    for name, n, c, st in ok[:10]:
        print(f'  ✅ {name}: {n} 残基, {c} 接触 [{st}]')
    if fail:
        print(f'失败: {len(fail)}')
        for name, msg in fail:
            print(f'  ❌ {name}: {msg}')


if __name__ == '__main__':
    main()
