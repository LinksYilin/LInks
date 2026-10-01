"""
align_and_filter.py — 对齐 + 突变级残基核对 + 剔除
====================================================
结合 alignment_map + 突变位置残基核对，产出"最终可用的评估集"。

流程：
  1. 读 benchmark CSV + alignment map。
  2. 对每个突变：用对齐的 (链, offset) 把 pos 映射到 PDB 残基，核对残基 == wt_aa。
  3. 一致 → 保留；不一致/越界/无结构 → 剔除并记录原因。
  4. 产出：可用评估集 CSV + 剔除清单 CSV。

用法：
  python align_and_filter.py \
    --benchmark data/raw/benchmarks/s669.csv \
    --alignment data/alignment_map_s669.csv \
    --struct_dir data/structures \
    --out_clean data/benchmarks_s669_clean.csv \
    --out_excluded data/benchmarks_s669_excluded.csv
"""
import argparse
import os

from Bio.PDB import PDBParser, Polypeptide

AA3TO1 = Polypeptide.protein_letters_3to1


def chain_residues(pdb_path, chain_id):
    """返回某链的标准氨基酸单字母列表。"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('x', pdb_path)
    model = structure[0]
    res_list = []
    for chain in model.get_chains():
        if chain.id == chain_id:
            for res in chain.get_residues():
                if Polypeptide.is_aa(res, standard=True):
                    aa = AA3TO1.get(res.get_resname().strip(), '')
                    if aa:
                        res_list.append(aa)
    return res_list


def pdb_from_filename(fn):
    return fn.split('.')[0].replace('pdb', '').upper()


def find_pdb_file(struct_dir, pid):
    for fn in os.listdir(struct_dir):
        if (fn.lower().endswith(('.ent', '.pdb'))) and (pdb_from_filename(fn) == pid.upper()):
                return os.path.join(struct_dir, fn)
    return None


def process(benchmark_csv, alignment_csv, struct_dir):
    import pandas as pd
    df = pd.read_csv(benchmark_csv)
    align = pd.read_csv(alignment_csv)
    align_map = {r['pdb_id']: r for _, r in align.iterrows()}

    kept = []
    excluded = []
    for _, row in df.iterrows():
        pid = row['pdb_id']
        wt_aa = row['mut_info'][0]
        pos = int(row['pos'])
        a = align_map.get(pid)
        reason = None
        if a is None or a['status'] != 'OK':
            reason = '对齐失败或无对齐记录'
        else:
            pdb_path = find_pdb_file(struct_dir, pid)
            if pdb_path is None:
                reason = '无结构文件'
            else:
                residues = chain_residues(pdb_path, a['chain'])
                offset = int(a['offset'])
                idx = offset + pos
                if idx < 0 or idx >= len(residues):
                    reason = f'位置越界(idx={idx},链长={len(residues)})'
                elif residues[idx] != wt_aa:
                    reason = f'残基不一致(PDB={residues[idx]} vs WT={wt_aa})'
                else:
                    # 保留，附上映射信息
                    new_row = row.to_dict()
                    new_row['_chain'] = a['chain']
                    new_row['_pdb_res_idx'] = idx
                    kept.append(new_row)
                    continue
        if reason is not None:
            ex = row.to_dict()
            ex['_exclusion_reason'] = reason
            excluded.append(ex)

    return pd.DataFrame(kept), pd.DataFrame(excluded)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--benchmark', required=True)
    ap.add_argument('--alignment', required=True)
    ap.add_argument('--struct_dir', required=True)
    ap.add_argument('--out_clean', required=True)
    ap.add_argument('--out_excluded', required=True)
    args = ap.parse_args()

    kept, excluded = process(args.benchmark, args.alignment, args.struct_dir)
    kept.to_csv(args.out_clean, index=False)
    excluded.to_csv(args.out_excluded, index=False)
    print(f'保留: {len(kept)} 突变, 剔除: {len(excluded)} 突变')
    print(f'已写入 {args.out_clean} 与 {args.out_excluded}')


if __name__ == '__main__':
    main()
