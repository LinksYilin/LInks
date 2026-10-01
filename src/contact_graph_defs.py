"""
contact_graph_defs.py — 多原子定义 + 多阈值的接触图构建
=========================================================
支持 4 种残基代表性定义：
  ca       : Cα 原子
  cb       : Cβ 原子（甘氨酸回退 Cα）
  centroid : 侧链质心（除主链外的重原子质心；甘氨酸回退 Cα）
  allatom  : 残基间全原子最小距离（KD-tree）

支持任意接触阈值（Å）。一次解析结构，输出各定义下"距离 < max_cut"的残基对及距离，
再由调用方按阈值切分，避免重复解析。
"""
import numpy as np
from scipy.spatial import cKDTree

from paths import DATA, FIGURES, PUB_FIGURES, ROOT, SRC, TOOLS  # 路径集中管理

SRC_PATH = str(SRC)
DATA_PATH = str(DATA)
ROOT_PATH = str(ROOT)
FIG_PATH = str(FIGURES)
PUB_PATH = str(PUB_FIGURES)
TOOLS_PATH = str(TOOLS)

try:
    from Bio.PDB import PDBParser
    from Bio.PDB.Polypeptide import is_aa, protein_letters_3to1



except ImportError:
    raise SystemExit('请先安装 biopython')

AA3TO1 = protein_letters_3to1
BACKBONE = {'N', 'CA', 'C', 'O', 'OXT'}

AA_PROPERTIES = {
    'A': [1.8, 88.6, 0], 'R': [-4.5, 173.4, 1], 'N': [-3.5, 114.1, 0],
    'D': [-3.5, 111.1, -1], 'C': [2.5, 108.5, 0], 'Q': [-3.5, 143.8, 0],
    'E': [-3.5, 138.4, -1], 'G': [-0.4, 60.1, 0], 'H': [-3.2, 153.2, 0],
    'I': [4.5, 166.7, 0], 'L': [3.8, 166.7, 0], 'K': [-3.9, 168.6, 1],
    'M': [1.9, 162.9, 0], 'F': [2.8, 189.9, 0], 'P': [-1.6, 112.7, 0],
    'S': [-0.8, 89.0, 0], 'T': [-0.7, 116.1, 0], 'W': [-0.9, 227.8, 0],
    'Y': [-1.3, 193.6, 0], 'V': [4.2, 140.0, 0],
}
AA_ORDER = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
AA_INDEX = {a: i for i, a in enumerate(AA_ORDER)}


def extract_residues(pdb_path, chain_id=None, offset=0, length=None):
    """解析结构，返回 (residue 列表, 单字母序列)。与既有流程保持一致的筛选规则。"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('p', pdb_path)
    model = structure[0]
    chains = list(model.get_chains())
    if chain_id is not None:
        chains = [c for c in chains if c.id == chain_id]
    if not chains:
        return None, None
    chain = chains[0]

    residues, seq = [], []
    for res in chain:
        name = res.get_resname().strip()
        if name == 'MSE':
            aa = 'M'
        else:
            if not is_aa(res, standard=True):
                continue
            aa = AA3TO1.get(name)
            if aa is None or aa not in AA_INDEX:
                continue
        if not res.has_id('CA'):
            continue
        residues.append(res)
        seq.append(aa)

    if offset > 0 or length is not None:
        end = (offset + length) if length is not None else None
        residues = residues[offset:end]
        seq = seq[offset:end]
    return residues, ''.join(seq)


def _centroid(res):
    atoms = [a for a in res.get_atoms() if a.get_name() not in BACKBONE and a.element != 'H']
    if atoms:
        return np.mean([a.get_coord() for a in atoms], axis=0)
    return res['CA'].get_coord()


def _atom_coord(res, name):
    if res.has_id(name):
        return res[name].get_coord()
    return None


def repr_coords(residues, atom_def):
    """返回各残基的代表性坐标 (N,3)。"""
    out = []
    for res in residues:
        if atom_def == 'ca':
            c = _atom_coord(res, 'CA')
        elif atom_def == 'cb':
            c = _atom_coord(res, 'CB')
            if c is None:
                c = _atom_coord(res, 'CA')   # 甘氨酸
        elif atom_def == 'centroid':
            c = _centroid(res)
        else:
            raise ValueError(atom_def)
        out.append(c)
    return np.array(out, dtype=np.float64)


def pairwise_within(coords, max_cut):
    """返回 {(i,j): d}，仅距离 < max_cut 的残基对（i<j）。"""
    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=max_cut, output_type='ndarray')
    best = {}
    for a, b in pairs:
        i, j = (int(a), int(b)) if a < b else (int(b), int(a))
        d = float(np.linalg.norm(coords[i] - coords[j]))
        key = (i, j)
        if key not in best or d < best[key]:
            best[key] = d
    return best


def allatom_min_pairs(residues, max_cut):
    """残基间全原子最小距离：返回 {(i,j): min_d}，仅 min_d < max_cut。"""
    coords, ridx = [], []
    for i, res in enumerate(residues):
        for a in res.get_atoms():
            if a.element == 'H':
                continue
            coords.append(a.get_coord())
            ridx.append(i)
    if not coords:
        return {}
    coords = np.array(coords, dtype=np.float64)
    ridx = np.array(ridx)
    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=max_cut, output_type='ndarray')
    best = {}
    for a, b in pairs:
        ri, rj = int(ridx[a]), int(ridx[b])
        if ri == rj:
            continue
        i, j = (ri, rj) if ri < rj else (rj, ri)
        d = float(np.linalg.norm(coords[a] - coords[b]))
        key = (i, j)
        if key not in best or d < best[key]:
            best[key] = d
    return best


def compute_all_defs(pdb_path, chain_id, offset, length, max_cut=10.0):
    """一次解析，返回 {atom_def: {(i,j): d}}（距离 < max_cut）。"""
    residues, seq = extract_residues(pdb_path, chain_id, offset, length)
    if residues is None or len(residues) < 2:
        return None, None
    out = {}
    for ad in ['ca', 'cb', 'centroid']:
        out[ad] = pairwise_within(repr_coords(residues, ad), max_cut)
    out['allatom'] = allatom_min_pairs(residues, max_cut)
    return out, seq


def edges_at(pairs_dict, threshold):
    """按阈值切出边集合。"""
    return {k for k, d in pairs_dict.items() if d < threshold}


if __name__ == '__main__':
    # 自测：1BFM M35W 的 WT/MT 在 4 种定义下的边数
    import os
    wt = os.path.join(DATA_PATH, 'structures', 'pdb1bfm.ent')
    mt = os.path.join(DATA_PATH, 'mutant_structures_s669', '1BFM', 'M35W', 'work', 'out',
                      'pdb1bfm_1.pdb')
    for tag, path in [('WT', wt), ('MT', mt)]:
        d, seq = compute_all_defs(path, 'A', 0, None, max_cut=10.0)
        if d is None:
            print(f'{tag}: 解析失败')
            continue
        print(f'{tag} ({len(seq)} 残基):')
        for ad in ['ca', 'cb', 'centroid', 'allatom']:
            print(f'   {ad:<9} 8Å 边数 = {len(edges_at(d[ad], 8.0))}')
