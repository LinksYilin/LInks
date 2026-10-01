"""
filter_leakage.py — 防泄漏过滤：训练集 vs 测试集同源性剔除
============================================================
目标：确保训练集中没有与 S669/S461/ssym 测试集蛋白同源的序列。

标准（Cuturello et al. 2024, Bioinformatics btae447）：
    剔除条件（三条件【同时满足】才剔除训练蛋白，即 AND 关系）：
      1. 序列同一性 > 25%（相对比对长度）
      2. e-value < 0.01
      3. 比对覆盖 > 50%（比对长度 / 查询长度）

实现：调用 BLASTp（需安装 blast+ 并加入 PATH）。
    提供 MMseqs2 作为更快的替代（注释中的命令）。

用法：
  python filter_leakage.py \
    --train_fasta train.fasta \
    --test_fasta S669.fasta \
    --out exclude_list.txt
"""
import argparse
import os
import subprocess
import tempfile


def run_blastp(train_fasta, test_fasta, tmpdir):
    """对 test 中每条序列在 train 库中做 blastp，返回命中列表。"""
    train_db = os.path.join(tmpdir, 'train_db')
    subprocess.run(['makeblastdb', '-in', train_fasta, '-dbtype', 'prot',
                    '-out', train_db, '-parse_seqids'],
                   check=True, capture_output=True)

    out = os.path.join(tmpdir, 'blast_out.txt')
    subprocess.run(['blastp', '-db', train_db, '-query', test_fasta,
                    '-out', out, '-outfmt', '6 qseqid sseqid pident evalue qlen slen length'],
                   check=True, capture_output=True)
    return out


# BLAST outfmt 6 固定 7 列：qseqid sseqid pident evalue qlen slen length
# （sseqid 是训练库中的序列 = 要剔除的对象；qseqid 是测试集查询序列）
BLAST_COLS = 7


def parse_blast(out_file):
    """返回需要剔除的训练序列 id 集合（三条件同时满足才剔除）。"""
    exclude = set()
    with open(out_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < BLAST_COLS:
                continue
            sseqid = parts[1]  # parts[0] 为 qseqid，本流程不使用
            pident = float(parts[2])
            evalue = float(parts[3])
            qlen = int(parts[4])
            alen = int(parts[6])
            identity_ok = pident > 25.0
            evalue_ok = evalue < 0.01
            # 覆盖 = 比对长度 / 查询长度（与 Cuturello 原文一致）
            coverage_ok = (alen / qlen) > 0.5 if qlen > 0 else False
            if identity_ok and evalue_ok and coverage_ok:
                exclude.add(sseqid)
    return exclude


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--train_fasta', required=True)
    ap.add_argument('--test_fasta', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        blast_out = run_blastp(args.train_fasta, args.test_fasta, tmpdir)
        exclude = parse_blast(blast_out)

    with open(args.out, 'w') as f:
        f.writelines(sid + '\n' for sid in sorted(exclude))
    print(f'需剔除 {len(exclude)} 条训练序列，清单已写入 {args.out}')


if __name__ == '__main__':
    main()

# ----------------------------------------------------------------------
# MMseqs2 更快替代（若训练集很大，推荐）：
#   mmseqs createdb train.fasta train_db
#   mmseqs createdb test.fasta test_db
#   mmseqs search test_db train_db result tmp --min-seq-id 0.25 -e 0.01 \
#            --cov-mode 1 -c 0.5
#   mmseqs createtsv test_db train_db result hits.tsv
# 然后按同上的三条件筛选 hits.tsv 的 sseqid。
# ----------------------------------------------------------------------
