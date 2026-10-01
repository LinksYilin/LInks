"""
foldx_batch_buildmodel.py — 批量 FoldX BuildModel 生成突变体结构
==================================================================
对 S669 干净集的每个突变，调用 FoldX BuildModel 生成突变体结构。

FoldX 输入（已验证）：
  - 链 ID：_chain
  - PDB 残基号：mut_info 里的数字（如 H461D → 461）
  - 野生残基：mut_info[0]
  - 突变残基：mut_info[-1]

用法：
  python foldx_batch_buildmodel.py \
    --label_csv data/benchmarks_s669_clean.csv \
    --struct_dir data/structures \
    --foldx_bin D:/GED_mutation/tools/foldx/foldx.exe \
    --out_dir data/mutant_structures_s669 \
    --limit 5   # 调试用，先跑 5 个
"""
import argparse
import os
import re
import subprocess


def parse_mut_info(mut_info):
    """'H461D' → (wt_aa='H', pdb_resnum=461, mt_aa='D')。"""
    m = re.match(r'^([A-Z])(\d+)([A-Z])$', mut_info)
    if not m:
        return None
    return m.group(1), int(m.group(2)), m.group(3)


def find_pdb_file(struct_dir, pid):
    for fn in os.listdir(struct_dir):
        if fn.lower().endswith(('.ent', '.pdb')):
            base = fn.split('.')[0].replace('pdb', '').upper()
            if base == pid.upper():
                return os.path.join(struct_dir, fn)
    return None


def run_foldx(foldx_bin, wt_pdb, chain, wt_aa, resnum, mt_aa, out_dir):
    """对单个突变调用 FoldX BuildModel。返回输出 PDB 路径。"""
    os.makedirs(out_dir, exist_ok=True)
    import shutil

    # 工作目录：把 pdb 复制成 .pdb，mutant-file 也放这里（相对路径）
    work_dir = os.path.join(out_dir, 'work')
    os.makedirs(work_dir, exist_ok=True)
    base = os.path.splitext(os.path.basename(wt_pdb))[0]
    pdb_copy = os.path.join(work_dir, f'{base}.pdb')
    shutil.copy(wt_pdb, pdb_copy)

    # FoldX individual_list 格式：野生残基+链+残基号+突变残基;
    mutant_spec = f'{wt_aa}{chain}{resnum}{mt_aa};'
    mf_path = os.path.join(work_dir, 'individual_list.txt')
    with open(mf_path, 'w') as f:
        f.write(mutant_spec + '\n')

    out_sub = os.path.join(work_dir, 'out')
    os.makedirs(out_sub, exist_ok=True)

    # 全部用相对路径，在 work_dir 下运行
    cmd = [foldx_bin, '--command=BuildModel', f'--pdb={base}.pdb',
           '--mutant-file=individual_list.txt', '--output-dir=out']
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)  # noqa: PLW1510 (返回值由调用方检查)
    if result.returncode != 0:
        return None, result.stderr[-500:]

    # FoldX 输出命名：{base}_1.pdb
    out_pdb = os.path.join(out_sub, f'{base}_1.pdb')
    if os.path.exists(out_pdb):
        return out_pdb, None
    return None, f'输出未找到: {out_pdb}'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--label_csv', required=True)
    ap.add_argument('--struct_dir', required=True)
    ap.add_argument('--foldx_bin', required=True)
    ap.add_argument('--out_dir', required=True)
    ap.add_argument('--limit', type=int, default=0, help='0=全部，否则只跑前 N 个')
    args = ap.parse_args()

    import pandas as pd
    df = pd.read_csv(args.label_csv)
    if args.limit > 0:
        df = df.head(args.limit)

    os.makedirs(args.out_dir, exist_ok=True)
    success, fail = 0, 0
    for _, row in df.iterrows():
        pid = row['pdb_id']
        parsed = parse_mut_info(row['mut_info'])
        if parsed is None:
            print(f'  [跳过] {pid} {row["mut_info"]}: 无法解析')
            fail += 1
            continue
        wt_aa, resnum, mt_aa = parsed
        chain = row['_chain']
        pdb_path = find_pdb_file(args.struct_dir, pid)
        if pdb_path is None:
            print(f'  [跳过] {pid}: 无结构文件')
            fail += 1
            continue

        sub_dir = os.path.join(args.out_dir, pid, row['mut_info'])
        out_pdb, err = run_foldx(args.foldx_bin, pdb_path, chain, wt_aa, resnum, mt_aa, sub_dir)
        if out_pdb:
            success += 1
        else:
            print(f'  [失败] {pid} {row["mut_info"]}: {err}')
            fail += 1

    print(f'\n完成: 成功 {success}, 失败 {fail}')


if __name__ == '__main__':
    main()
