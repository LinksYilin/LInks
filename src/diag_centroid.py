"""诊断：既有 WT 图 243 边对应哪种质心定义？"""
import numpy as np
from Bio.PDB import PDBParser
from scipy.spatial import cKDTree

PDB = r'D:\GED_mutation\data\structures\pdb1bfm.ent'
d = np.load(r'D:\GED_mutation\data\contact_graphs_s669_sc\1BFM.npz', allow_pickle=True)
OLD = {(min(int(a), int(b)), max(int(a), int(b))) for a, b in d['edge_index'].T}
print(f'既有 WT 边数: {len(OLD)}')

parser = PDBParser(QUIET=True)
st = parser.get_structure('p', PDB)
model = st[0]
chain = [c for c in model.get_chains() if c.id == 'A'][0]

residues = []
for res in chain:
    if res.get_resname().strip() not in ('MSE',) and not res.has_id('CA'):
        continue
    if res.get_resname().strip() == 'MSE':
        pass
    residues.append(res)
print(f'残基数: {len(residues)}')
# 检查是否有氢原子
n_h = sum(1 for res in residues for a in res.get_atoms() if a.element == 'H')
n_all = sum(1 for res in residues for a in res.get_atoms())
print(f'原子总数 {n_all}，其中氢 {n_h}')

BACKBONE = {'N', 'CA', 'C', 'O', 'OXT'}


def edges_from_coords(coords, cut=8.0):
    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=cut, output_type='ndarray')
    return {(min(int(a), int(b)), max(int(a), int(b))) for a, b in pairs}


def test(name, coord_fn):
    try:
        coords = np.array([coord_fn(r) for r in residues], dtype=np.float64)
    except Exception as e:
        print(f'  {name}: 失败 {e}')
        return
    E = edges_from_coords(coords, 8.0)
    inter = len(E & OLD)
    print(f'  {name:<46} 边数={len(E):>4} 与既有交集={inter:>4} 仅既有={len(OLD-E):>3} 仅新={len(E-OLD):>3}')


def centroid_noH(r):
    a = [x for x in r.get_atoms() if x.get_name() not in BACKBONE and x.element != 'H']
    return np.mean([x.get_coord() for x in a], axis=0) if a else r['CA'].get_coord()


def centroid_all(r):
    a = [x for x in r.get_atoms() if x.get_name() not in BACKBONE]
    return np.mean([x.get_coord() for x in a], axis=0) if a else r['CA'].get_coord()


def cb(r):
    return r['CB'].get_coord() if r.has_id('CB') else r['CA'].get_coord()


def ca(r):
    return r['CA'].get_coord()


def centroid_heavy_with_gly(r):
    """含甘氨酸也用 CA（与既有实现一致）"""
    a = [x for x in r.get_atoms() if x.get_name() not in BACKBONE and x.element != 'H']
    return np.mean([x.get_coord() for x in a], axis=0) if a else r['CA'].get_coord()


def centroid_fake(r):
    return centroid_noH(r)


print('\n不同坐标定义对比 8Å 边集:')
test('Cα', ca)
test('Cβ（GLY回退Cα）', cb)
test('侧链质心（排除H，GLY回退Cα）', centroid_noH)
test('侧链质心（含H，GLY回退Cα）', centroid_all)
# 另一种：排除H与OXT但含C
def centroid_sidechain_strict(r):
    a = [x for x in r.get_atoms() if x.get_name() not in ('N', 'CA', 'C', 'O', 'OXT') and x.element != 'H']
    return np.mean([x.get_coord() for x in a], axis=0) if a else r['CA'].get_coord()
test('侧链质心（同 definition 2）', centroid_sidechain_strict)

# 直接看既有 metadata 提示的 n_contacts
print(f'\n既有 metadata n_contacts: {d["metadata"][0] if "metadata" in d else "NA"}')
