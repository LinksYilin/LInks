"""
test_scwrl4.py — 验证 SCWRL4 的 -s 突变机制（先小规模测试，再规模化）
=====================================================================
测试：用 1BFM 野生型骨架 + 改过字母的序列文件，看能否正确生成 M35W 突变体。
并检验"输出残基顺序是否与输入一致"这一关键问题。
"""
import os
import subprocess

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1

SCWRL = r'D:\GED_mutation\tools\scwrl4\Scwrl4.exe'
WORK = r'D:\GED_mutation\tools\scwrl4\test_mut'
WT = r'D:\GED_mutation\data\structures\pdb1bfm.ent'
CHAIN = 'A'
MUT_POS_0 = 34          # 0-based，对应序列第 35 位（M）
MUT_FROM, MUT_TO = 'M', 'W'

os.makedirs(WORK, exist_ok=True)
parser = PDBParser(QUIET=True)


def extract_chain(src, dst, chain_id):
    """把指定链写成单链 PDB。"""
    st = parser.get_structure('p', src)
    ch = [c for c in st[0].get_chains() if c.id == chain_id][0]
    with open(dst, 'w') as f:
        n = 0
        for res in ch:
            if not res.has_id('CA'):
                continue
            for a in res.get_atoms():
                if a.element == 'H':
                    continue
                f.write(
                    f"ATOM  {n+1:>5} {a.get_name():<4}{'':1}{res.get_resname():>3} "
                    f"{chain_id}{res.id[1]:>4}{'':1}   "
                    f"{a.get_coord()[0]:>8.3f}{a.get_coord()[1]:>8.3f}{a.get_coord()[2]:>8.3f}"
                    f"  1.00  0.00          {a.element:>2}\n")
                n += 1
            f.write('TER\n')
        f.write('END\n')
    return ch


def seq_of(chain):
    out = []
    for res in chain:
        if not res.has_id('CA'):
            continue
        nm = res.get_resname().strip()
        aa = 'M' if nm == 'MSE' else (protein_letters_3to1.get(nm) if is_aa(res, standard=True) else None)
        if aa is None:
            aa = 'G'
        out.append((aa, res.id[1]))
    return out


# 1) 抽链
single = os.path.join(WORK, 'wt_chainA.pdb')
ch = extract_chain(WT, single, CHAIN)
seq = seq_of(ch)
print(f'输入链 {CHAIN}: {len(seq)} 个残基')
print(f'  序列: {"".join(a for a, _ in seq)}')
print(f'  第 35 位(1-based) 是 {seq[MUT_POS_0][0]}（残基号 {seq[MUT_POS_0][1]}）')

# 2) 构造突变序列：突变位点大写，其余小写（保留原构象）
wt_seq = ''.join(a for a, _ in seq)
mut_seq = list(wt_seq)
mut_seq[MUT_POS_0] = MUT_TO
# 全大写 = 全部重排（与 FoldX BuildModel 的"重排突变位点+邻居"不同，但先测通）
seq_full = ''.join(mut_seq).upper()
seq_file = os.path.join(WORK, 'seq_mut.txt')
with open(seq_file, 'w') as f:
    f.write(seq_full + '\n')
print(f'  突变后序列(全大写): {seq_full}')

# 3) 运行 SCWRL4
out_pdb = os.path.join(WORK, 'mt_M35W.pdb')
cmd = [SCWRL, '-i', 'wt_chainA.pdb', '-o', 'mt_M35W.pdb', '-s', 'seq_mut.txt', '-h']
r = subprocess.run(cmd, cwd=WORK, capture_output=True, text=True)
print(f'\nSCWRL4 返回码: {r.returncode}')
if r.returncode != 0:
    print('STDERR:', r.stderr[-500:])
    raise SystemExit(1)

# 4) 检验输出
st2 = parser.get_structure('m', out_pdb)
ch2 = [c for c in st2[0].get_chains()][0]
seq2 = seq_of(ch2)
print(f'\n输出: {len(seq2)} 个残基')
print(f'  序列: {"".join(a for a, _ in seq2)}')

# 5) 关键检验：残基号顺序是否与输入一致
in_ids = [rid for _, rid in seq]
out_ids = [rid for _, rid in seq2]
print('\n=== 残基号顺序检查 ===')
print(f'  输入顺序: {in_ids[:12]} ... {in_ids[-4:]}')
print(f'  输出顺序: {out_ids[:12]} ... {out_ids[-4:]}')
print(f'  顺序是否相同: {in_ids == out_ids}')

# 6) 突变位点是否正确
target_rid = seq[MUT_POS_0][1]
mt_aa_out = next((a for a, rid in seq2 if rid == target_rid), None)
print('\n=== 突变检查 ===')
print(f'  残基号 {target_rid}: 输入 {MUT_FROM} -> 输出 {mt_aa_out}（期望 {MUT_TO}）')
print(f'  {"✅ 突变成功" if mt_aa_out == MUT_TO else "❌ 突变失败"}')

# 7) 输出是否有氢
st_h = parser.get_structure('h', out_pdb)
n_h = sum(1 for res in st_h[0].get_residues() for a in res.get_atoms() if a.element == 'H')
print(f'  输出氢原子数: {n_h}（用 -h 应为 0）')
