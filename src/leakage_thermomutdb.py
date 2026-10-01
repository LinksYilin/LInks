"""
leakage_thermomutdb.py — ThermoMutDB 相对 S669/ssym 的同源性防泄漏过滤
=====================================================================
MegaScale 的 train_s669 划分已由 ThermoMPNN 做过防泄漏，但 ThermoMutDB 没有。
两者都是天然蛋白，与 S669 可能同源，必须单独过滤。

流程：
  1. 训练端序列：ThermoMutDB 各蛋白接触图的链序列（metadata['seq']）
  2. 测试端序列：S669 + ssym 的 wt_seq
  3. BLASTp（测试序列查询训练库），按 Cuturello 2024 三条件 AND 剔除训练蛋白
  4. 产出：剔除清单 + 过滤后的 thermomutdb_aligned_noleak.csv
"""
import json
import os
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from filter_leakage import parse_blast

from pathlib import Path as _P
DATA = str(_P(__file__).resolve().parent.parent / 'data')
BLAST_BIN = os.environ.get('BLAST_BIN', str(__import__('pathlib').Path(__file__).resolve().parent.parent / 'tools' / 'blast' / 'ncbi-blast-2.17.0+' / 'bin'))


def write_fasta(records, path):
    with open(path, 'w') as f:
        f.writelines(f'>{rid}\n{seq}\n' for rid, seq in records)


def main():
    # 训练端：ThermoMutDB 接触图链序列
    tm = pd.read_csv(os.path.join(DATA, 'thermomutdb_aligned.csv'))
    train_recs = []
    for pid in sorted(tm['pdb_id'].unique()):
        p = os.path.join(DATA, 'contact_graphs_thermomutdb', f'{pid}.npz')
        if not os.path.exists(p):
            continue
        meta = json.loads(str(np.load(p, allow_pickle=True)['metadata'][0]))
        if meta.get('seq'):
            train_recs.append((pid, meta['seq']))

    # 测试端：S669 + ssym 的 wt_seq（每蛋白一条）
    test_recs = []
    for name in ['benchmarks_s669_clean.csv', 'benchmarks_ssym_clean.csv']:
        df = pd.read_csv(os.path.join(DATA, name))
        for pid, seq in df.groupby('pdb_id')['wt_seq'].first().items():
            test_recs.append((f'{name[11:15]}_{pid}', seq))

    print(f'训练端 ThermoMutDB 蛋白: {len(train_recs)}，测试端蛋白: {len(test_recs)}')

    with tempfile.TemporaryDirectory() as tmp:
        train_fa = os.path.join(tmp, 'train.fasta')
        test_fa = os.path.join(tmp, 'test.fasta')
        write_fasta(train_recs, train_fa)
        write_fasta(test_recs, test_fa)
        db = os.path.join(tmp, 'train_db')
        out = os.path.join(tmp, 'blast.tsv')
        import subprocess
        subprocess.run([os.path.join(BLAST_BIN, 'makeblastdb.exe'), '-in', train_fa,
                        '-dbtype', 'prot', '-out', db], check=True, capture_output=True)
        subprocess.run([os.path.join(BLAST_BIN, 'blastp.exe'), '-db', db, '-query', test_fa,
                        '-out', out, '-evalue', '1',
                        '-outfmt', '6 qseqid sseqid pident evalue qlen slen length'],
                       check=True, capture_output=True)
        exclude = parse_blast(out)

    exclude = sorted(exclude)
    with open(os.path.join(DATA, 'thermomutdb_leakage_exclude.txt'), 'w') as f:
        f.write('\n'.join(exclude))

    kept = tm[~tm['pdb_id'].isin(exclude)]
    kept.to_csv(os.path.join(DATA, 'thermomutdb_aligned_noleak.csv'), index=False)
    print(f'需剔除 ThermoMutDB 蛋白: {len(exclude)}')
    print(f'突变: {len(tm)} → {len(kept)}（剔除 {len(tm) - len(kept)}）')
    print(f'蛋白: {tm["pdb_id"].nunique()} → {kept["pdb_id"].nunique()}')


if __name__ == '__main__':
    main()
