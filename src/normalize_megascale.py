"""
normalize_megascale.py — MegaScale 训练标签规范化
==================================================
按发现确定 MegaScale 训练数据：
  1. 标签 = ddG_ML（不是 deltaG；README 确认 ddG_ML 是标准训练标签）
  2. 符号：ddG_ML 正值=stabilizing，与标准（正值=destabilizing）相反 → 翻号
  3. 用 split_name 里的 train_s669 划分（ThermoMPNN 提供的、已对 S669 防泄漏的训练子集）
  4. 去重：同一 (WT_name, mut_type) 多条记录取一条
  5. 突变位置解析：mut_type 形如 'E1Q'（WT残基 + 位置 + MT残基，1-based）

用法：
  python normalize_megascale.py \
    --input_glob 'data/raw/megascale/train-*.parquet' \
    --out data/megascale_clean.csv
"""
import argparse
import glob

import pandas as pd


def parse_mut_type(mut_type):
    """解析 'E1Q' → (wt_aa, pos_1based, mt_aa)。"""
    wt_aa = mut_type[0]
    mt_aa = mut_type[-1]
    pos = int(mut_type[1:-1])
    return wt_aa, pos, mt_aa


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input_glob', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--split', default='train_s669', help='用哪个 split_name')
    args = ap.parse_args()

    files = sorted(glob.glob(args.input_glob))
    dfs = []
    for f in files:
        df = pd.read_parquet(f, columns=['name', 'WT_name', 'mut_type', 'ddG_ML', 'split_name'])
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True)

    # 1. 只用指定 split（train_s669 = 对 S669 防泄漏的训练子集）
    df = df[df['split_name'] == args.split].copy()
    print(f'{args.split} 划分的记录数: {len(df)}')

    # 2. 去重：同一 (WT_name, mut_type) 取第一条
    df = df.drop_duplicates(subset=['WT_name', 'mut_type'], keep='first')
    print(f'去重后: {len(df)}')

    # 3. ddG_ML 转数值 + 翻号（正值=stabilizing → 标准正值=destabilizing）
    df['ddg'] = -pd.to_numeric(df['ddG_ML'], errors='coerce')
    df = df.dropna(subset=['ddg'])

    # 4. 解析突变位置
    parsed = df['mut_type'].apply(parse_mut_type)
    df['wt_aa'] = parsed.apply(lambda x: x[0])
    df['pos_1based'] = parsed.apply(lambda x: x[1])
    df['mt_aa'] = parsed.apply(lambda x: x[2])

    # 5. 只保留单点突变（MegaScale 单点数据集本身应都是单点，但保险起见）
    df = df[df['mut_type'].str.len() >= 3]

    # 输出精简列
    out = df[['WT_name', 'mut_type', 'wt_aa', 'pos_1based', 'mt_aa', 'ddg']].copy()
    out.to_csv(args.out, index=False)
    print(f'最终训练标签: {len(out)} 条')
    print(f'ddg 范围: [{out["ddg"].min():.2f}, {out["ddg"].max():.2f}], 均值 {out["ddg"].mean():.2f}')
    print(f'翻号后正值(失稳)占比: {(out["ddg"]>0).mean()*100:.1f}%')


if __name__ == '__main__':
    main()
