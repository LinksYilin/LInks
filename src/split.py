"""
split.py — 训练/验证集划分（防同源泄漏）
========================================
原则：按"蛋白"划分，同一蛋白（及其同源蛋白）不能同时出现在 train 和 val。

方案 A（推荐）：CD-HIT 聚类（30% 序列同一性）→ 按聚类簇做 GroupKFold。
方案 B（无需外部工具）：按 PDB 蛋白 id 直接 GroupKFold（同源蛋白需已由
    filter_leakage.py 相对 test 集过滤；此处保证 val 与 train 无同蛋白即可）。

用法：
  python split.py --csv clean_train.csv --protein_col pdb_id --out splits/
"""
import argparse
import os

import pandas as pd

try:
    from sklearn.model_selection import GroupKFold
except ImportError:
    raise SystemExit("请先安装 scikit-learn: pip install scikit-learn")


def cdhit_cluster(fasta_path, identity=0.3, out_prefix='cdhit'):
    """用 CD-HIT 聚类（需安装 cd-hit 并加入 PATH）。返回 {序列id: 簇id}。

    CD-HIT .clstr 格式示例：
        >Cluster 0
        0\t1234aa, >seqA... *
        1\t1234aa, >seqB... at 95.00%
        >Cluster 1
        ...
    每行以数字开头（0=代表序列，>0=簇内其它成员）。所有行都必须归入当前簇，
    否则同簇成员会漏配、被 GroupKFold 拆到 train/val 两侧。
    """
    import subprocess
    clstr_out = out_prefix + '.clstr'
    subprocess.run(['cd-hit', '-i', fasta_path, '-o', out_prefix,
                    '-c', str(identity), '-n', '4', '-T', '4'],
                   check=True, capture_output=True)
    seq2cluster = {}
    cluster_id = -1
    with open(clstr_out) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>Cluster'):
                cluster_id += 1
                continue
            if not line.strip():
                continue
            # 簇成员行：形如 "0\t1234aa, >seqid... *" 或 "1\t... >seqid... at 95%"
            if '>' in line:
                # 提取 >seqid 部分（seqid 到空格或 "..." 为止）
                after_gt = line.split('>', 1)[1]
                seqid = after_gt.split('...')[0].split()[0].rstrip('*')
                if seqid:
                    seq2cluster[seqid] = cluster_id
    return seq2cluster


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', required=True, help='清洗后的训练 CSV，需含蛋白 id 列')
    ap.add_argument('--protein_col', default='pdb_id')
    ap.add_argument('--fasta', default=None, help='（可选）用于 CD-HIT 聚类的 fasta')
    ap.add_argument('--out', default='splits')
    ap.add_argument('--n_splits', type=int, default=5)
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    df = pd.read_csv(args.csv)
    os.makedirs(args.out, exist_ok=True)

    if args.fasta and os.path.exists(args.fasta):
        # 方案 A：按 CD-HIT 簇分组
        seq2cluster = cdhit_cluster(args.fasta)
        groups = df[args.protein_col].map(lambda p: seq2cluster.get(p, p)).values
        print(f'按 CD-HIT 簇分组，共 {len(set(groups))} 个簇')
    else:
        # 方案 B：按蛋白 id 分组
        groups = df[args.protein_col].values
        print(f'按蛋白 id 分组，共 {len(set(groups))} 个蛋白')

    # 留一折做 val，其余做 train；产出多折供交叉验证
    gkf = GroupKFold(n_splits=args.n_splits)
    for fold, (train_idx, val_idx) in enumerate(gkf.split(df, groups=groups)):
        df.iloc[train_idx].to_csv(os.path.join(args.out, f'train_fold{fold}.csv'), index=False)
        df.iloc[val_idx].to_csv(os.path.join(args.out, f'val_fold{fold}.csv'), index=False)
        # 报告泄漏自检：train/val 是否有共享蛋白
        shared = set(df.iloc[train_idx][args.protein_col]) & set(df.iloc[val_idx][args.protein_col])
        print(f'fold{fold}: train={len(train_idx)}, val={len(val_idx)}, 共享蛋白={len(shared)}')
        assert len(shared) == 0, f'fold{fold} 存在蛋白泄漏！'


if __name__ == '__main__':
    main()
