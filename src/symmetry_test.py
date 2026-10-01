"""
symmetry_test.py — 突变方向对称性检验（第 3 项）
==================================================
思路：把同一个突变的"正向"与"反向"都喂给模型，检验预测是否反号。

正向输入：WT 图（突变位点节点 = 突变后氨基酸）
反向输入：MT 图（突变位点节点 = 突变前氨基酸）

理想情况：pred_reverse ≈ −pred_forward（ΔΔG(A→B) = −ΔΔG(B→A)）。
报告：
  - pred_reverse 与 pred_forward 的相关系数
  - |pred_rev + pred_fwd| 的分布（偏离理想对称的程度）
  - 各模型的对称性得分
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader

sys.path.insert(0, os.path.dirname(__file__))
from contact_graph_defs import extract_residues, pairwise_within, repr_coords
from gnn_local_baseline import GNNLocal
from paths import DATA, FIGURES, PUB_FIGURES, ROOT, SRC, TOOLS  # 路径集中管理
from train_gnn_baseline import AA_INDEX, GNNRegressor

SRC_PATH = str(SRC)
DATA_PATH = str(DATA)
ROOT_PATH = str(ROOT)
FIG_PATH = str(FIGURES)
PUB_PATH = str(PUB_FIGURES)
TOOLS_PATH = str(TOOLS)




DATA = DATA_PATH  # 来自 paths.py，可用 GED_ROOT 环境变量覆盖
GRAPH_DIRS = {'megascale': os.path.join(DATA, 'contact_graphs_megascale_sc'),
              'thermomutdb': os.path.join(DATA, 'contact_graphs_thermomutdb_sc')}
AA_ORDER = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
IDX2AA = {i: a for i, a in enumerate(AA_ORDER)}
SEEDS = [42, 123, 2024]
TH, MAXC = 8.0, 10.0
BACKBONE = {'N', 'CA', 'C', 'O', 'OXT'}


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def load_train():
    df = pd.read_csv(os.path.join(DATA, 'training_merged_noleak_sc.csv'))
    ms = df[df['source'] == 'megascale']
    tm = df[df['source'] == 'thermomutdb']
    n = min(len(tm), len(ms))
    ms = ms.sample(n=n, random_state=42)
    df = pd.concat([ms, tm], ignore_index=True)
    cache, out = {}, []
    for _, r in df.iterrows():
        p = os.path.join(GRAPH_DIRS.get(r['source']), f"{r['protein']}.npz")
        if not os.path.exists(p):
            continue
        key = f"{r['source']}:{r['protein']}"
        if key not in cache:
            d = np.load(p, allow_pickle=True)
            cache[key] = (d['nodes'], d['edge_index'], d['edge_attr'])
        nodes, ei, ea = cache[key]
        mi = int(r['mut_idx'])
        if mi < 0 or mi >= nodes.shape[0]:
            continue
        nodes = nodes.copy()
        mt = str(r['mt_aa'])
        if mt in AA_INDEX:
            nodes[mi, :20] = 0.0
            nodes[mi, AA_INDEX[mt]] = 1.0
        flag = np.zeros((nodes.shape[0], 1), dtype=np.float32)
        flag[mi, 0] = 1.0
        x = torch.tensor(np.concatenate([nodes, flag], axis=1), dtype=torch.float32)
        out.append(Data(x=x, edge_index=torch.tensor(ei, dtype=torch.long),
                        edge_attr=torch.tensor(ea, dtype=torch.float32),
                        y=torch.tensor([float(r['ddg'])], dtype=torch.float32)))
    return out


def build_from_coords(coords, pairs, nodes_aa, mut_idx, aa_at_mut, flag_val=1.0):
    """由坐标+边对构建 Data。"""
    edges = [k for k, d in pairs.items() if d < TH]
    if edges:
        ei = np.array([[i, j] for i, j in edges] + [[j, i] for i, j in edges], dtype=np.int64).T
    else:
        ei = np.zeros((2, 0), dtype=np.int64)
    ea = np.zeros((ei.shape[1], 2), dtype=np.float32)
    feats = []
    for aa in nodes_aa:
        onehot = np.zeros(20, dtype=np.float32)
        if aa in AA_INDEX:
            onehot[AA_INDEX[aa]] = 1.0
        props = np.array([KD.get(aa, 0.0), VOL.get(aa, 140.0) / 250.0, CHG.get(aa, 0)], dtype=np.float32)
        feats.append(np.concatenate([onehot, props]))
    nodes = np.stack(feats)
    nodes[mut_idx, :20] = 0.0
    if aa_at_mut in AA_INDEX:
        nodes[mut_idx, AA_INDEX[aa_at_mut]] = 1.0
    flag = np.zeros((nodes.shape[0], 1), dtype=np.float32)
    flag[mut_idx, 0] = flag_val
    x = torch.tensor(np.concatenate([nodes, flag], axis=1), dtype=torch.float32)
    return Data(x=x, edge_index=torch.tensor(ei, dtype=torch.long),
                edge_attr=torch.tensor(ea, dtype=torch.float32),
                y=torch.tensor([0.0], dtype=torch.float32))


KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,
      'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
VOL = {'A':88.6,'R':173.4,'N':114.1,'D':111.1,'C':108.5,'Q':143.8,'E':138.4,'G':60.1,'H':153.2,
       'I':166.7,'L':166.7,'K':168.6,'M':162.9,'F':189.9,'P':112.7,'S':89.0,'T':116.1,'W':227.8,'Y':193.6,'V':140.0}
CHG = {'D':-1,'E':-1,'K':1,'R':1}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--epochs', type=int, default=20)
    args = ap.parse_args()

    df = pd.read_csv(os.path.join(DATA, 'benchmarks_s669_clean.csv'))
    align = pd.read_csv(os.path.join(DATA, 'alignment_map_s669.csv'))
    amap = {r['pdb_id']: r for _, r in align.iterrows()}
    mt_dir = os.path.join(DATA, 'mutant_structures_s669')

    fwd, rev, metas = [], [], []
    for _, row in df.iterrows():
        pid, mut = row['pdb_id'], row['mut_info']
        a = amap.get(pid)
        if a is None:
            continue
        base = f'pdb{pid.lower()}'
        mt_path = os.path.join(mt_dir, pid, mut, 'work', 'out', f'{base}_1.pdb')
        wt_path = os.path.join(DATA, 'structures', f'{base}.ent')
        if not (os.path.exists(mt_path) and os.path.exists(wt_path)):
            continue
        L = len(row['wt_seq'])
        ch, off = a['chain'], int(a['offset'])
        wr, sw = extract_residues(wt_path, ch, off, L)
        mr, sm = extract_residues(mt_path, ch, off, L)
        if wr is None or mr is None or len(sw) != len(sm):
            continue
        if '_node_idx' not in row.index or pd.isna(row['_node_idx']):
            raise ValueError(
                f"缺少已验证的 _node_idx（{(row['pdb_id'], row['mut_info'])}）；"
                "拒绝回退到有问题的 _pdb_res_idx，请先运行 fix_index_bug.py")
        mi = int(row['_node_idx'])
        if mi >= len(sw):
            continue
        wt_aa, mt_aa = sw[mi], sm[mi]
        # 正向：WT 坐标 + 突变后氨基酸
        wc = repr_coords(wr, 'centroid')
        wp = pairwise_within(wc, MAXC)
        fwd.append(build_from_coords(wc, wp, list(sw), mi, mt_aa))
        # 反向：MT 坐标 + 突变前氨基酸
        mc_ = repr_coords(mr, 'centroid')
        mp = pairwise_within(mc_, MAXC)
        rev.append(build_from_coords(mc_, mp, list(sm), mi, wt_aa))
        metas.append({'pdb_id': pid, 'mut_info': mut, 'ddg': row['ddg']})

    print(f'正向/反向样本: {len(fwd)}')
    tr = load_train()
    in_dim = tr[0].x.shape[1]
    y = np.array([m['ddg'] for m in metas])

    results = {}
    for seed in SEEDS:
        for tag, Model in [('gnn_local', GNNLocal), ('gnn_global', GNNRegressor)]:
            torch.manual_seed(seed); np.random.seed(seed)
            model = Model(in_dim, hid=64)
            opt = torch.optim.Adam(model.parameters(), lr=1e-3)
            loader = DataLoader(tr, batch_size=64, shuffle=True)
            for _ in range(args.epochs):
                model.train()
                for b in loader:
                    opt.zero_grad(); F.mse_loss(model(b), b.y).backward(); opt.step()
            model.eval()
            def pred(samples, model=model):  # 绑定循环变量（B023）
                out = []
                with torch.no_grad():
                    for b in DataLoader(samples, batch_size=64):
                        out.append(model(b).numpy())
                return np.concatenate(out)
            pf, pr = pred(fwd), pred(rev)
            results[f'{tag}_s{seed}'] = (pf, pr)
            print(f'  {tag} seed {seed}: r(fwd,rev)={pearson(pf,pr):.3f}, '
                  f'r(fwd,ddg)={pearson(pf,y):.3f}, mean|fwd+rev|={np.mean(np.abs(pf+pr)):.3f}')

    print()
    print('=== 对称性汇总 ===')
    print(f'{"model":<18} {"r(fwd,rev)":>11} {"mean|fwd+rev|":>14} {"r(fwd,ddg)":>10} {"r(rev,-ddg)":>12}')
    rows = []
    for k, (pf, pr) in results.items():
        rows.append({'model': k, 'r_fwd_rev': pearson(pf, pr),
                     'mean_abs_sum': float(np.mean(np.abs(pf + pr))),
                     'r_fwd_ddg': pearson(pf, y), 'r_rev_negsym_ddg': pearson(pr, -y)})
        print(f'{k:<18} {pearson(pf,pr):>11.3f} {np.mean(np.abs(pf+pr)):>14.3f} '
              f'{pearson(pf,y):>10.3f} {pearson(pr,-y):>12.3f}')
    pd.DataFrame(rows).to_csv(os.path.join(DATA, 'symmetry_test.csv'), index=False)
    print(f'\n已保存 {os.path.join(DATA, "symmetry_test.csv")}')


if __name__ == '__main__':
    main()
