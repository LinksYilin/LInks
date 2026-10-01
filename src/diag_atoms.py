"""诊断 2：WT 与 FoldX MT 结构的原子组成是否一致 + 含氢质心下的编辑数"""
import numpy as np
from Bio.PDB import PDBParser
from scipy.spatial import cKDTree

WT = r'D:\GED_mutation\data\structures\pdb1bfm.ent'
MT = r'D:\GED_mutation\data\mutant_structures_s669\1BFM\M35W\work\out\pdb1bfm_1.pdb'
BACKBONE = {'N', 'CA', 'C', 'O', 'OXT'}

parser = PDBParser(QUIET=True)


def load_chain(path, chain_id='A'):
    st = parser.get_structure('p', path)
    ch = [c for c in st[0].get_chains() if c.id == chain_id]
    if not ch:
        return None
    res = [r for r in ch[0] if r.has_id('CA')]
    return res


for tag, path in [('WT', WT), ('MT(FoldX)', MT)]:
    res = load_chain(path)
    if res is None:
        print(f'{tag}: 无 A 链')
        continue
    n_all = sum(1 for r in res for a in r.get_atoms())
    n_h = sum(1 for r in res for a in r.get_atoms() if a.element == 'H')
    n_heavy = n_all - n_h
    print(f'{tag}: {len(res)} 残基, 原子 {n_all} (氢 {n_h}, 重原子 {n_heavy})')
    # 看几个残基的原子
    r0 = res[34] if len(res) > 34 else res[0]
    print(f'   残基 {r0.get_resname()}{r0.id[1]} 原子: {[a.get_name() for a in r0.get_atoms()]}')
    # 氢原子名样例
    hs = [a.get_name() for r in res for a in r.get_atoms() if a.element == 'H'][:8]
    print(f'   氢原子名样例: {hs}')

print()
print('=== 含氢质心下的 WT vs MT 编辑 ===')


def edges(res, include_h, cut=8.0):
    coords = []
    for r in res:
        a = [x for x in r.get_atoms()
             if x.get_name() not in BACKBONE and (include_h or x.element != 'H')]
        if not a:
            a = [r['CA']]
        coords.append(np.mean([x.get_coord() for x in a], axis=0))
    coords = np.array(coords)
    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=cut, output_type='ndarray')
    return {(min(int(i), int(j)), max(int(i), int(j))) for i, j in pairs}


wt_res = load_chain(WT)
mt_res = load_chain(MT)
for inc in [True, False]:
    try:
        W = edges(wt_res, inc)
        M = edges(mt_res, inc)
        tag = '含氢' if inc else '排除氢'
        print(f'  {tag:<6} WT={len(W):>4} MT={len(M):>4} 断边={len(W-M):>3} 成边={len(M-W):>3}')
    except Exception as e:
        print(f'  {"含氢" if inc else "排除氢"}: 失败 {e}')
