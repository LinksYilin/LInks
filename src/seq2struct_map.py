"""
seq2struct_map.py — 序列→结构对齐（把 wt_seq 的突变位置映射到 PDB 残基）
============================================================================
解决 S669/ssym 的 wt_seq ≠ PDB 链完整序列的问题：
  找 wt_seq 在 PDB 链序列中的偏移 offset，使 pos → PDB 残基号。

用法：
  python seq2struct_map.py --benchmark s669.csv --struct_dir structures/ \
      --out alignment_map.csv
"""
import argparse
import csv
import os

from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser, Polypeptide

AA3TO1 = Polypeptide.protein_letters_3to1


def chain_sequences(pdb_path):
    """返回 {链id: 完整单字母序列}，只取标准氨基酸。"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('x', pdb_path)
    model = structure[0]
    out = {}
    for chain in model.get_chains():
        seq = ''
        for res in chain.get_residues():
            if Polypeptide.is_aa(res, standard=True):
                aa = AA3TO1.get(res.get_resname().strip(), '')
                if aa:
                    seq += aa
        out[chain.id] = seq
    return out


def find_offset(wt_seq, chain_seq):
    """子串精确匹配，返回 offset（wt_seq 在 chain_seq 中的起始 0-based 位置）或 None。"""
    idx = chain_seq.find(wt_seq)
    return idx if idx >= 0 else None


def fuzzy_offset(wt_seq, chain_seq):
    """局部比对（PairwiseAligner），返回 (begin, end, score)。begin 是 chain_seq 上的起始。"""
    if not wt_seq or not chain_seq:
        return None
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    aln = aligner.align(chain_seq, wt_seq)
    if not aln:
        return None
    a = aln[0]
    # a.aligned[0] 是 chain_seq(target) 上的对齐块，取第一个块的起止
    target_blocks = a.aligned[0]
    if len(target_blocks) == 0:
        return None
    begin = int(target_blocks[0][0])
    end = int(target_blocks[-1][1])
    score = a.score
    return begin, end, score


def pdb_from_filename(fn):
    base = fn.split('.')[0]
    return base.replace('pdb', '').upper()


def find_pdb_file(struct_dir, pdb_id):
    for fn in os.listdir(struct_dir):
        if (fn.lower().endswith(('.ent', '.pdb'))) and (pdb_from_filename(fn) == pdb_id.upper()):
                return os.path.join(struct_dir, fn)
    return None


def map_benchmark(benchmark_csv, struct_dir):
    """对 benchmark 里每个蛋白做对齐，返回结果列表。"""
    import pandas as pd
    df = pd.read_csv(benchmark_csv)
    # 每个 pdb_id 取第一条 wt_seq
    wtseq_by_pdb = {}
    for _, row in df.iterrows():
        wtseq_by_pdb.setdefault(row['pdb_id'], row['wt_seq'])

    results = []
    for pdb_id, wt_seq in sorted(wtseq_by_pdb.items()):
        pdb_path = find_pdb_file(struct_dir, pdb_id)
        if pdb_path is None:
            results.append({'pdb_id': pdb_id, 'chain': None, 'offset': None,
                            'wt_len': len(wt_seq), 'chain_len': None,
                            'status': 'FAIL', 'note': '无结构文件'})
            continue
        try:
            chains = chain_sequences(pdb_path)
        except Exception as e:
            results.append({'pdb_id': pdb_id, 'chain': None, 'offset': None,
                            'wt_len': len(wt_seq), 'chain_len': None,
                            'status': 'FAIL', 'note': f'解析失败 {e}'})
            continue

        matched = False
        for chain_id, chain_seq in chains.items():
            off = find_offset(wt_seq, chain_seq)
            if off is not None:
                results.append({'pdb_id': pdb_id, 'chain': chain_id, 'offset': off,
                                'wt_len': len(wt_seq), 'chain_len': len(chain_seq),
                                'status': 'OK', 'note': '精确子串'})
                matched = True
                break
        if not matched:
            # 模糊匹配
            best = None
            for chain_id, chain_seq in chains.items():
                r = fuzzy_offset(wt_seq, chain_seq)
                if r is not None and (best is None or r[2] > best[1][2]):
                    best = (chain_id, r)
            if best is not None:
                chain_id, (begin, end, score) = best
                results.append({'pdb_id': pdb_id, 'chain': chain_id, 'offset': begin,
                                'wt_len': len(wt_seq), 'chain_len': len(chains[chain_id]),
                                'status': 'OK', 'note': f'模糊匹配(score={score:.0f})'})
            else:
                results.append({'pdb_id': pdb_id, 'chain': None, 'offset': None,
                                'wt_len': len(wt_seq), 'chain_len': max(len(s) for s in chains.values()),
                                'status': 'FAIL', 'note': '无法对齐'})
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--benchmark', required=True)
    ap.add_argument('--struct_dir', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    results = map_benchmark(args.benchmark, args.struct_dir)
    ok = [r for r in results if r['status'] == 'OK']
    fail = [r for r in results if r['status'] == 'FAIL']
    print(f'对齐结果: OK {len(ok)}/{len(results)}, FAIL {len(fail)}')
    for r in fail:
        print(f"  ❌ {r['pdb_id']}: {r['note']}")

    with open(args.out, 'w', newline='') as f:
        fieldnames = ['pdb_id', 'chain', 'offset', 'wt_len', 'chain_len', 'status', 'note']
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(results)
    print(f'已写入 {args.out}')


if __name__ == '__main__':
    main()
