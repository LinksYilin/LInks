"""
build_thermomutdb_graphs.py — ThermoMutDB 残基核对 + 接触图构建
================================================================
ThermoMutDB 对齐策略（见 ThermoMutDB对齐策略.md）：
  mutation_code 如 'E49M'：PDB 第 49 号残基（E）突变成 M。
  直接用 PDB 残基号定位，无需序列比对。

流程：
  1. 解析 mutation_code → (wt_aa, pdb_resnum, mt_aa)
  2. 下载/定位 PDB，找到 pdb_resnum 对应的残基
  3. 核对残基 == wt_aa
  4. 一致 → 建接触图（取含该残基的链）；不一致 → 剔除

用法：
  python build_thermomutdb_graphs.py \
    --label_csv data/thermomutdb_clean.csv \
    --struct_dir data/structures \
    --out_dir data/contact_graphs_thermomutdb \
    --out_clean data/thermomutdb_aligned.csv
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, Polypeptide

AA3TO1 = Polypeptide.protein_letters_3to1


def parse_mutation_code(mc):
    """'E49M' → (wt_aa='E', resnum=49, mt_aa='M')。残基号去前导零。"""
    wt_aa = mc[0]
    mt_aa = mc[-1]
    resnum_str = mc[1:-1]
    try:
        resnum = int(resnum_str)  # 去前导零：'03' -> 3
    except ValueError:
        resnum = resnum_str
    return wt_aa, resnum, mt_aa


def pdb_from_filename(fn):
    return fn.split('.')[0].replace('pdb', '').upper()


def find_pdb_file(struct_dir, pid):
    for fn in os.listdir(struct_dir):
        if (fn.lower().endswith(('.ent', '.pdb'))) and (pdb_from_filename(fn) == pid.upper()):
            return os.path.join(struct_dir, fn)
    return None


def find_residue_by_number(pdb_path, resnum):
    """按 PDB 残基号找残基，返回 (chain_id, 该链的残基列表, 该残基在链内的索引)。
    返回 None 若找不到。resnum 可为 int 或 str（统一按 int 比较，兼容前导零）。"""
    parser = PDBParser(QUIET=True)
    s = parser.get_structure('x', pdb_path)
    m = s[0]
    # 目标残基号：int 化比较
    try:
        target = int(resnum)
    except (ValueError, TypeError):
        target = None
    for chain in m.get_chains():
        residues = []
        for r in chain.get_residues():
            if Polypeptide.is_aa(r, standard=True) or r.get_resname().strip() == 'MSE':
                aa = 'M' if r.get_resname().strip() == 'MSE' else AA3TO1.get(r.get_resname().strip(), '')
                residues.append((r.id[1], aa))  # (PDB resnum, aa)
        for idx, (rn, aa) in enumerate(residues):
            if target is not None:
                try:
                    if int(rn) == target:
                        return chain.id, residues, idx
                except (ValueError, TypeError):
                    pass
            elif str(rn) == str(resnum):
                return chain.id, residues, idx
    return None


def build_graph_from_pdb(pdb_path, chain_id, residues):
    """用 BioPython 构建指定链的接触图（复用 build_contact_graph 的 build_graph）。"""
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    from build_contact_graph import build_graph
    return build_graph(pdb_path, chain_id=chain_id)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--label_csv', required=True)
    ap.add_argument('--struct_dir', required=True)
    ap.add_argument('--out_dir', required=True)
    ap.add_argument('--out_clean', required=True)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    df = pd.read_csv(args.label_csv)

    # 蛋白级缓存：每个 pid 只解析一次
    pdb_cache = {}  # pid -> (pdb_path 或 None)
    residue_cache = {}  # pid -> {chain_id: (residues, aa_list)} 或 None

    def get_pdb(pid):
        if pid not in pdb_cache:
            pdb_cache[pid] = find_pdb_file(args.struct_dir, pid)
        return pdb_cache[pid]

    def get_residues(pid, pdb_path):
        if pid not in residue_cache:
            # 解析所有链的残基
            parser = PDBParser(QUIET=True)
            s = parser.get_structure('x', pdb_path)
            m = s[0]
            chains = {}
            for chain in m.get_chains():
                res_list = []
                for r in chain.get_residues():
                    if Polypeptide.is_aa(r, standard=True) or r.get_resname().strip() == 'MSE':
                        aa = 'M' if r.get_resname().strip() == 'MSE' else AA3TO1.get(r.get_resname().strip(), '')
                        res_list.append((r.id[1], aa))
                chains[chain.id] = res_list
            residue_cache[pid] = chains
        return residue_cache[pid]

    kept, excluded = [], []
    graph_count = 0
    for _, row in df.iterrows():
        pid = row['pdb_id']
        if not isinstance(pid, str) or not pid.strip():
            excluded.append((str(pid), str(row['mutation_code']), 'pdb_id 无效'))
            continue
        pid = pid.strip()
        mc = row['mutation_code']
        wt_aa, resnum, mt_aa = parse_mutation_code(mc)
        pdb_path = get_pdb(pid)
        if pdb_path is None:
            excluded.append((pid, mc, '无结构文件'))
            continue
        chains = get_residues(pid, pdb_path)
        # 在所有链中找 resnum 匹配且 wt_aa 一致的残基
        found = None
        try:
            target = int(resnum)
        except (ValueError, TypeError):
            target = None
        for chain_id, residues in chains.items():
            for idx, (rn, aa) in enumerate(residues):
                match = False
                if target is not None:
                    try:
                        match = (int(rn) == target)
                    except (ValueError, TypeError):
                        match = (str(rn) == str(resnum))
                else:
                    match = (str(rn) == str(resnum))
                if (match) and (aa == wt_aa):
                    found = (chain_id, residues, idx, aa)
                else:
                    found = ('MISMATCH', chain_id, residues, idx, aa)
                break
            if found:
                break
        if found is None:
            excluded.append((pid, mc, f'残基号 {resnum} 未找到'))
            continue
        if found[0] == 'MISMATCH':
            chain_id, residues, idx, aa = found[1], found[2], found[3], found[4]
            excluded.append((pid, mc, f'残基不一致 PDB={aa} vs WT={wt_aa}'))
            continue
        chain_id, residues, idx, aa = found
        # 建接触图（整条链，每 pid 缓存）
        out_path = os.path.join(args.out_dir, f'{pid}.npz')
        if not os.path.exists(out_path):
            g = build_graph_from_pdb(pdb_path, chain_id, residues)
            if g is None:
                excluded.append((pid, mc, '接触图构建失败'))
                continue
            nodes, edge_index, edge_attr, meta = g
            np.savez(out_path, nodes=nodes, edge_index=edge_index, edge_attr=edge_attr,
                     metadata=np.array([json.dumps(meta)]))
            graph_count += 1
        # 记录该突变的链内索引
        new_row = row.to_dict()
        new_row['_chain'] = chain_id
        new_row['_pdb_res_idx'] = idx
        kept.append(new_row)

    kept_df = pd.DataFrame(kept)
    kept_df.to_csv(args.out_clean, index=False)
    excl_df = pd.DataFrame(excluded, columns=['pdb_id', 'mutation_code', 'reason'])
    excl_df.to_csv(args.out_clean.replace('.csv', '_excluded.csv'), index=False)

    print(f'保留: {len(kept)}, 剔除: {len(excluded)}, 接触图: {graph_count}')
    print('剔除原因分布:')
    print(excl_df['reason'].value_counts().to_string() if len(excl_df) else '无')


if __name__ == '__main__':
    main()
