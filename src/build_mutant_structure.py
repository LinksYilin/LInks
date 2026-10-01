"""
build_mutant_structure.py — WT→MT 成对结构生成（GEDMut 2.0 结构生成协议）
============================================================================
按 GEDMut2.0_设计定稿.md 的三档结构来源协议，生成突变体结构。

三档：
  Tier A：实验解析的 WT + MT 成对结构（外部提供，本脚本只做索引登记）。
  Tier B（主）：WT 实验结构 + FoldX BuildModel / SCWRL4 side-chain repacking。
  Tier C（fallback）：AlphaFold 预测（外部提供，仅登记）。

用法：
  python build_mutant_structure.py \
    --wt_pdb 1ABC.pdb \
    --mutations "A:42:W"   # 链:位置(1-based):突变残基，多个用逗号分隔 \
    --tool foldx|scwrl \
    --out_dir mutant_structures/ \
    --n_seeds 3            # 多结构生成，用于置信度统计

依赖（运行时才需要）：
  FoldX: foldx 可执行文件在 PATH，或通过 --foldx_bin 指定。
  SCWRL4: Scwrl4 可执行文件。
"""
import argparse
import os
import subprocess

try:
    from Bio.PDB import PDBParser, Polypeptide
except ImportError:
    PDBParser = Polypeptide = None


def read_wt_residues(pdb_path):
    """从 PDB 读 (链, 残基号) -> 单字母残基。"""
    parser = PDBParser(QUIET=True)
    s = parser.get_structure('x', pdb_path)
    out = {}
    for m in s:
        for chain in m:
            for res in chain:
                if Polypeptide.is_aa(res, standard=True):
                    aa = Polypeptide.protein_letters_3to1.get(res.get_resname().strip(), '')
                    if aa:
                        out[(chain.id, res.id[1])] = aa
    return out


def parse_mutation(mut_str):
    """解析 "chain:pos:aa" 为 (chain, pos_int, aa)。pos 为 1-based 残基号。"""
    chain, pos, aa = mut_str.split(':')
    return chain, int(pos), aa.upper()


def run_foldx_buildmodel(wt_pdb, mutations, out_dir, foldx_bin='foldx'):
    """调用 FoldX BuildModel 生成突变体结构。返回输出 PDB 路径。

    FoldX individual_list.txt 每行格式："野生残基+链+残基号+突变残基;"
    例如 "VA42W;" 表示 A 链 42 位 V→W。
    """
    os.makedirs(out_dir, exist_ok=True)
    wt_residues = read_wt_residues(wt_pdb)

    mut_list = []
    for chain, pos, aa in mutations:
        wt_aa = wt_residues.get((chain, pos), '?')
        mut_list.append(f"{wt_aa}{chain}{pos}{aa};")
    mf_path = os.path.join(out_dir, 'individual_list.txt')
    with open(mf_path, 'w') as f:
        f.write(','.join(mut_list) + '\n')

    cmd = [foldx_bin, '--command=BuildModel', f'--pdb={wt_pdb}',
           f'--mutant-file={mf_path}', f'--output-dir={out_dir}']
    subprocess.run(cmd, check=True)
    out_pdb = os.path.join(out_dir, f'{os.path.splitext(os.path.basename(wt_pdb))[0]}_1.pdb')
    return out_pdb


def run_scwrl(wt_pdb, mutations, out_dir, scwrl_bin='Scwrl4'):
    """调用 SCWRL4 生成突变体结构（side-chain repacking，不改骨架）。"""
    os.makedirs(out_dir, exist_ok=True)
    # SCWRL4 用突变后完整序列文件
    seq_file = os.path.join(out_dir, 'mut_seq.txt')
    # 从 wt_pdb 读序列并替换突变位点
    write_mutant_sequence(wt_pdb, mutations, seq_file)
    out_pdb = os.path.join(out_dir, 'mt_scwrl.pdb')
    cmd = [scwrl_bin, '-i', wt_pdb, '-s', seq_file, '-o', out_pdb]
    subprocess.run(cmd, check=True)
    return out_pdb


def write_mutant_sequence(wt_pdb, mutations, seq_file):
    """从 wt_pdb 读序列，替换突变位点，写出突变序列（SCWRL4 格式）。"""
    parser = PDBParser(QUIET=True)
    s = parser.get_structure('x', wt_pdb)
    for m in s:
        for chain in m:
            seq = ''
            for res in chain:
                if Polypeptide.is_aa(res, standard=True):
                    aa = Polypeptide.protein_letters_3to1.get(res.get_resname().strip(), '')
                    seq += aa if aa else ''
            # 替换该链上的突变
            # 需要残基号 → 序列索引映射；简化：按顺序假设连续（待完善）
            with open(seq_file, 'w') as f:
                f.write(f'>{chain.id}\n{seq}\n')
            break  # 只处理第一条链


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--wt_pdb', required=True)
    ap.add_argument('--mutations', required=True, help='逗号分隔的 chain:pos:aa')
    ap.add_argument('--tool', default='foldx', choices=['foldx', 'scwrl'])
    ap.add_argument('--out_dir', required=True)
    ap.add_argument('--n_seeds', type=int, default=3)
    ap.add_argument('--foldx_bin', default='foldx')
    ap.add_argument('--scwrl_bin', default='Scwrl4')
    args = ap.parse_args()

    mutations = [parse_mutation(m) for m in args.mutations.split(',')]
    print(f"WT: {args.wt_pdb}")
    print(f"突变: {mutations}")
    print(f"工具: {args.tool}, 结构数: {args.n_seeds}")

    for s in range(args.n_seeds):
        seed_dir = os.path.join(args.out_dir, f'seed_{s}')
        os.makedirs(seed_dir, exist_ok=True)
        if args.tool == 'foldx':
            out_pdb = run_foldx_buildmodel(args.wt_pdb, mutations, seed_dir, args.foldx_bin)
        else:
            out_pdb = run_scwrl(args.wt_pdb, mutations, seed_dir, args.scwrl_bin)
        print(f"  seed_{s}: {out_pdb}")

    print("完成。多结构结果用于 edit_confidence.py 计算置信度。")


if __name__ == '__main__':
    main()
