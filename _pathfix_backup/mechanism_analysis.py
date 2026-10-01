"""
mechanism_analysis.py — 机制分析：解释"结构敏感但不提升预测"
==============================================================
对应路线图第 4 项。计算：
4.1 编辑特征与 ΔΔG 的相关 + 控制理化特征后的偏相关（冗余性）
4.2 断裂接触"距突变位点距离"与 ΔΔG 的关联（距离依赖）
4.3 按体积变化幅度分层
所有统计带蛋白簇 bootstrap 95% CI（B=1000）。
"""
import os

import numpy as np
import pandas as pd
from numpy.linalg import lstsq
from scipy.stats import spearmanr

DATA = r'D:\GED_mutation\data'

KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,
      'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
VOL = {'A':88.6,'R':173.4,'N':114.1,'D':111.1,'C':108.5,'Q':143.8,'E':138.4,'G':60.1,'H':153.2,
       'I':166.7,'L':166.7,'K':168.6,'M':162.9,'F':189.9,'P':112.7,'S':89.0,'T':116.1,'W':227.8,'Y':193.6,'V':140.0}
CHG = {'D':-1,'E':-1,'K':1,'R':1}


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def partial_corr(x, y, Z):
    """控制协变量 Z 后 x 与 y 的偏相关。"""
    x, y, Z = np.asarray(x, float), np.asarray(y, float), np.asarray(Z, float)
    Z1 = np.column_stack([np.ones(len(Z)), Z])
    bx = lstsq(Z1, x, rcond=None)[0]
    by = lstsq(Z1, y, rcond=None)[0]
    rx = x - Z1 @ bx
    ry = y - Z1 @ by
    return pearson(rx, ry)


def cluster_bootstrap_ci(values, protein_ids, metric, B=1000, seed=0):
    """蛋白簇 bootstrap 95% CI。metric 为一个函数 f(idx)->float。"""
    rng = np.random.default_rng(seed)
    proteins = np.array(protein_ids)
    unique, inverse = np.unique(proteins, return_inverse=True)
    pidx = [np.where(inverse == i)[0] for i in range(len(unique))]
    stats = []
    for _ in range(B):
        sampled = rng.integers(0, len(unique), len(unique))
        mask = np.concatenate([pidx[i] for i in sampled])
        s = metric(mask)
        if s is not None and not np.isnan(s):
            stats.append(s)
    if not stats:
        return float('nan'), float('nan')
    return float(np.percentile(stats, 2.5)), float(np.percentile(stats, 97.5))


def main():
    te = pd.read_csv(os.path.join(DATA, 'true_edits_s669_sc.csv'))
    df = pd.read_csv(os.path.join(DATA, 'benchmarks_s669_clean.csv'))
    df = df[['pdb_id', 'mut_info', 'ddg']]
    m = te.merge(df, on=['pdb_id', 'mut_info'], how='inner')
    print(f'合并样本: {len(m)} 条（有真实编辑 + 有 ΔΔG）')

    m = m.reset_index(drop=True)
    m['wt'] = m['mut_info'].str[0]
    m['mt'] = m['mut_info'].str[-1]
    m['d_vol'] = m.apply(lambda r: VOL.get(r['wt'], 140) - VOL.get(r['mt'], 140), axis=1)
    m['d_kd'] = m.apply(lambda r: KD.get(r['wt'], 0) - KD.get(r['mt'], 0), axis=1)
    m['vol_mt'] = m['mt'].map(VOL).fillna(140)
    m['kd_mt'] = m['mt'].map(KD).fillna(0)
    m['chg_mt'] = m['mt'].map(CHG).fillna(0)
    m['n_edit'] = m['n_broken'] + m['n_formed']

    phys = m[['d_vol', 'd_kd', 'vol_mt', 'kd_mt', 'chg_mt']].values
    y = m['ddg'].values
    pid = m['pdb_id'].values

    print('\n=== 4.1 编辑特征 vs ΔΔG（含控制理化后的偏相关）===')
    print(f'{"feature":<12} {"Pearson":>9} {"95% CI":>18} {"Spearman":>9} {"partial r|phys":>15}')
    for name in ['n_broken', 'n_formed', 'n_edit']:
        x = m[name].values
        r = pearson(x, y)
        lo, hi = cluster_bootstrap_ci(x, pid, lambda idx: pearson(x[idx], y[idx]))
        rho = float(spearmanr(x, y)[0]) if x.std() > 0 else float('nan')
        pr = partial_corr(x, y, phys)
        print(f'{name:<12} {r:>9.3f} [{lo:>7.3f},{hi:>7.3f}] {rho:>9.3f} {pr:>15.3f}')

    # 编辑特征与理化特征的相关（冗余性证据）
    print('\n=== 编辑特征与理化特征的最大 |相关|（冗余性证据）===')
    for name in ['n_broken', 'n_formed', 'n_edit']:
        cors = [abs(pearson(m[name].values, m[c].values)) for c in ['d_vol', 'd_kd', 'vol_mt', 'kd_mt', 'chg_mt']]
        print(f'  {name}: max |r| with physicochemical = {max(cors):.3f}')

    # 4.3 按体积变化幅度分层
    print('\n=== 4.3 按体积变化幅度分层（|Δvol| 三分位）===')
    m['abs_dvol'] = m['d_vol'].abs()
    q1, q2 = m['abs_dvol'].quantile([1/3, 2/3])
    for label, sub in [('small', m[m['abs_dvol'] <= q1]),
                       ('medium', m[(m['abs_dvol'] > q1) & (m['abs_dvol'] <= q2)]),
                       ('large', m[m['abs_dvol'] > q2])]:
        if len(sub) < 10:
            continue
        r = pearson(sub['n_edit'].values, sub['ddg'].values)
        rho = float(spearmanr(sub['n_edit'], sub['ddg'])[0])
        print(f'  {label:<7} n={len(sub):>3}  |Δvol| 范围 [{sub["abs_dvol"].min():.1f},{sub["abs_dvol"].max():.1f}]  '
              f'r(n_edit,ddg)={r:>6.3f}  rho={rho:>6.3f}')

    # 4.2 距离依赖（断边到突变位点距离 与 ΔΔG）
    print('\n=== 4.2 断裂接触距突变位点距离 vs ΔΔG ===')
    import ast

    from Bio.PDB import PDBParser, Polypeptide
    parser = PDBParser(QUIET=True)
    BACKBONE = {'N', 'CA', 'C', 'O', 'OXT'}

    def sc_centroid(res):
        atoms = [a for a in res.get_atoms() if a.get_name() not in BACKBONE]
        if atoms:
            return np.mean([a.get_coord() for a in atoms], axis=0)
        if res.has_id('CA'):
            return res['CA'].get_coord()
        return None

    cache = {}
    mean_d, n_edit_list, ddg_list, pid_list = [], [], [], []
    for _, row in m.iterrows():
        pid_ = row['pdb_id']
        wt_pdb = os.path.join(DATA, 'structures', f'pdb{pid_.lower()}.ent')
        if not os.path.exists(wt_pdb):
            continue
        if pid_ not in cache:
            s = parser.get_structure('x', wt_pdb)
            coords = []
            for c in s[0].get_chains():
                for r in c.get_residues():
                    if Polypeptide.is_aa(r, standard=True):
                        sc = sc_centroid(r)
                        if sc is not None:
                            coords.append(sc)
                break
            cache[pid_] = np.array(coords)
        coords = cache[pid_]
        mi = int(row['mut_idx'])
        if mi >= len(coords):
            continue
        mc = coords[mi]
        broken = ast.literal_eval(row['broken_edges']) if isinstance(row['broken_edges'], str) and row['broken_edges'] else []
        if not broken:
            continue
        ds = []
        for (a, b) in broken:
            if a < len(coords) and b < len(coords):
                ds.append(min(np.linalg.norm(coords[a] - mc), np.linalg.norm(coords[b] - mc)))
        if ds:
            mean_d.append(np.mean(ds))
            n_edit_list.append(row['n_edit'])
            ddg_list.append(row['ddg'])
            pid_list.append(pid_)

    mean_d = np.array(mean_d)
    ddg_a = np.array(ddg_list)
    n_edit_a = np.array(n_edit_list)
    print(f'  有效样本: {len(mean_d)}')
    print(f'  平均断边距离 vs ΔΔG: r = {pearson(mean_d, ddg_a):.3f}')
    print(f'  平均断边距离 vs 编辑总数: r = {pearson(mean_d, n_edit_a):.3f}')
    # 距离分箱
    print('  距离分箱（断边到突变位点距离）:')
    for lo, hi in [(0, 5), (5, 10), (10, 15), (15, 20), (20, 100)]:
        mask = (mean_d >= lo) & (mean_d < hi)
        if mask.sum() >= 10:
            print(f'    [{lo:>2},{hi:>3}) Å: n={mask.sum():>3}, mean ΔΔG={ddg_a[mask].mean():>5.2f}, '
                  f'r(dist,ddg)={pearson(mean_d[mask], ddg_a[mask]):>6.3f}')


if __name__ == '__main__':
    main()
