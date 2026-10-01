"""
benchmark_eval_full.py — 双基准（S669 + ssym）完整评估
========================================================
一次补齐两个缺失证据：
  A1. ssym 结果（Methods 声称做了二次评估但 Results 未报告）
  A2. edge-aware GNN 结果（§3.5 声称在 §4.4 评估但未报告）

模型（各 3 seed）：
  ridge        : 5 特征理化 Ridge
  gnn_global   : GCN 全局池化
  gnn_local    : GCN 局部（突变位点）池化
  gnn_local_blosum : 局部池化 + 固定 BLOSUM 特征
  gnn_edge     : GINEConv（使用边特征）

统计：蛋白簇 bootstrap 95% CI（B=1000，含重复簇），Spearman/MAE/RMSE，
      以及蛋白内 Pearson（Fisher-z，min_mut=5）。
逐样本 prediction 全部保存。
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from scipy.stats import spearmanr
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader

sys.path.insert(0, os.path.dirname(__file__))
from seed_utils import set_seed  # 完整确定性控制（含 cuBLAS/cudnn 标志）

SEED_DEFAULT = 42

from edit_cost import node_substitution_cost
from gnn_edge_baseline import GNNEdge
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


# ---------------- 指标 ----------------
# ---------------------------------------------------------------------------
# 逐样本排除日志（审查要求：不得静默改变样本数）
# ---------------------------------------------------------------------------
SKIP_LOG = []


def log_skip(stage, pid, mut, reason):
    SKIP_LOG.append({'stage': stage, 'pdb_id': pid, 'mut_info': mut, 'reason': reason})


def dump_skips(path):
    """把排除记录写 CSV，并打印分布；便于审稿人核对样本数变化。"""
    import pandas as pd



    if not SKIP_LOG:
        print('  [排除日志] 无样本被排除')
        return
    df = pd.DataFrame(SKIP_LOG)
    df.to_csv(path, index=False)
    print(f'  [排除日志] 共 {len(df)} 条，已写入 {path}')
    for (stage, reason), g in df.groupby(['stage', 'reason']):
        print(f'      {stage} / {reason}: {len(g)}')


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def mae(a, b):
    return float(np.mean(np.abs(np.asarray(a) - np.asarray(b))))


def rmse(a, b):
    return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def spearman(a, b):
    return float(spearmanr(a, b)[0])


def cluster_bootstrap(y_true, y_pred, protein_ids, metric, B=1000, seed=0):
    """蛋白簇 bootstrap：有放回抽蛋白，重复簇重复计数。"""
    rng = np.random.default_rng(seed)
    prot = np.array(protein_ids)
    uniq, inv = np.unique(prot, return_inverse=True)
    pidx = [np.where(inv == i)[0] for i in range(len(uniq))]
    stats = []
    for _ in range(B):
        sampled = rng.integers(0, len(uniq), len(uniq))
        mask = np.concatenate([pidx[i] for i in sampled])
        if len(mask) < 3:
            continue
        s = metric(y_true[mask], y_pred[mask])
        if not np.isnan(s):
            stats.append(s)
    if not stats:
        return float('nan'), float('nan'), float('nan')
    st = np.array(stats)
    return float(np.mean(st)), float(np.percentile(st, 2.5)), float(np.percentile(st, 97.5))


def protein_level_pearson(y_true, y_pred, protein_ids, min_mut=5):
    prot = np.array(protein_ids)
    rs = []
    for p in np.unique(prot):
        m = prot == p
        if m.sum() < min_mut:
            continue
        a, b = y_true[m], y_pred[m]
        if a.std() == 0 or b.std() == 0:
            continue
        rs.append(np.clip(np.corrcoef(a, b)[0, 1], -0.999, 0.999))
    if not rs:
        return float('nan'), 0
    return float(np.tanh(np.mean(np.arctanh(np.array(rs))))), len(rs)


# ---------------- 样本构造 ----------------
def build_sample(nodes, ei, ea, mut_idx, wt, mt, y, add_blosum):
    nodes = nodes.copy()
    if mt in AA_INDEX:
        nodes[mut_idx, :20] = 0.0
        nodes[mut_idx, AA_INDEX[mt]] = 1.0
    flag = np.zeros((nodes.shape[0], 1), dtype=np.float32)
    flag[mut_idx, 0] = 1.0
    extra = [flag]
    if add_blosum:
        bc = node_substitution_cost(wt, mt)
        col = np.zeros((nodes.shape[0], 1), dtype=np.float32)
        col[mut_idx, 0] = bc
        extra.append(col)
    x = torch.tensor(np.concatenate([nodes] + extra, axis=1), dtype=torch.float32)
    return Data(x=x, edge_index=torch.tensor(ei, dtype=torch.long),
                edge_attr=torch.tensor(ea, dtype=torch.float32),
                y=torch.tensor([float(y)], dtype=torch.float32))


def load_train(add_blosum=False):
    df = pd.read_csv(os.path.join(DATA, 'training_merged_noleak_sc.csv'))
    ms = df[df['source'] == 'megascale']
    tm = df[df['source'] == 'thermomutdb']
    n = min(len(tm), len(ms))
    ms = ms.sample(n=n, random_state=42)
    df = pd.concat([ms, tm], ignore_index=True)
    cache, out = {}, []
    for _, r in df.iterrows():
        gdir = GRAPH_DIRS.get(r['source'])
        p = os.path.join(gdir, f"{r['protein']}.npz")
        if not os.path.exists(p):
            log_skip('load', r.get('pdb_id', r.get('protein', '?')),
                     r.get('mut_info', r.get('mt_aa', '?')), 'graph_missing')
            continue
        key = f"{r['source']}:{r['protein']}"
        if key not in cache:
            d = np.load(p, allow_pickle=True)
            cache[key] = (d['nodes'], d['edge_index'], d['edge_attr'])
        nodes, ei, ea = cache[key]
        mi = int(r['mut_idx'])
        if mi < 0 or mi >= nodes.shape[0]:
            log_skip('index', r.get('pdb_id', r.get('protein', '?')),
                     r.get('mut_info', r.get('mt_aa', '?')),
                     f'index_out_of_range({mi})')
            continue
        wt = IDX2AA.get(int(np.argmax(nodes[mi, :20])), 'X')
        out.append(build_sample(nodes, ei, ea, mi, wt, str(r['mt_aa']), r['ddg'], add_blosum))
    return out


def load_benchmark(graph_dir, label_csv, add_blosum=False):
    df = pd.read_csv(label_csv)
    cache, out, meta = {}, [], []
    for _, r in df.iterrows():
        pid = r['pdb_id']
        p = os.path.join(graph_dir, f'{pid}.npz')
        if not os.path.exists(p):
            log_skip('load', r.get('pdb_id', r.get('protein', '?')),
                     r.get('mut_info', r.get('mt_aa', '?')), 'graph_missing')
            continue
        if pid not in cache:
            d = np.load(p, allow_pickle=True)
            cache[pid] = (d['nodes'], d['edge_index'], d['edge_attr'])
        nodes, ei, ea = cache[pid]
        # ★ 索引修正：_pdb_res_idx 是整链索引，图节点是切片后索引 → 用 _node_idx
        if '_node_idx' not in r.index or pd.isna(r['_node_idx']):
            raise ValueError(
                f"缺少已验证的 _node_idx（{(r['pdb_id'], r['mut_info'])}）；"
                "拒绝回退到有问题的 _pdb_res_idx，请先运行 fix_index_bug.py")
        mi = int(r['_node_idx'])
        if mi < 0 or mi >= nodes.shape[0]:
            log_skip('index', r.get('pdb_id', r.get('protein', '?')),
                     r.get('mut_info', r.get('mt_aa', '?')),
                     f'index_out_of_range({mi})')
            continue
        mi_info = str(r['mut_info'])
        out.append(build_sample(nodes, ei, ea, mi, mi_info[0], mi_info[-1], r['ddg'], add_blosum))
        meta.append({'pdb_id': pid, 'mut_info': mi_info, 'ddg': r['ddg']})
    return out, meta


def train_model(model, samples, epochs, seed):
    set_seed(seed)
    loader = DataLoader(samples, batch_size=64, shuffle=True)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    for _ in range(epochs):
        model.train()
        for b in loader:
            opt.zero_grad()
            loss = F.mse_loss(model(b), b.y)
            loss.backward()
            opt.step()
    return model


def predict(model, samples):
    model.eval()
    loader = DataLoader(samples, batch_size=64)
    out = []
    with torch.no_grad():
        for b in loader:
            out.append(model(b).numpy())
    return np.concatenate(out)


# ---------------- 主流程 ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--epochs', type=int, default=20)
    ap.add_argument('--out', default=os.path.join(DATA, 'benchmark_predictions.csv'))
    args = ap.parse_args()

    print('加载训练集...')
    tr_base = load_train(add_blosum=False)
    tr_bl = load_train(add_blosum=True)
    print(f'  训练样本 {len(tr_base)}')
    in_dim, in_dim_bl = tr_base[0].x.shape[1], tr_bl[0].x.shape[1]
    edge_dim = tr_base[0].edge_attr.shape[1]
    print(f'  节点维度 {in_dim} / {in_dim_bl}，边维度 {edge_dim}')

    benchmarks = {}
    for name, gdir, lcsv in [
        ('S669', os.path.join(DATA, 'contact_graphs_s669_sc'), os.path.join(DATA, 'benchmarks_s669_clean.csv')),
        ('ssym', os.path.join(DATA, 'contact_graphs_ssym_sc'), os.path.join(DATA, 'benchmarks_ssym_clean.csv')),
    ]:
        s, m = load_benchmark(gdir, lcsv, add_blosum=False)
        s_bl, _ = load_benchmark(gdir, lcsv, add_blosum=True)
        benchmarks[name] = {'base': s, 'blosum': s_bl, 'meta': m}
        n_prot = len({x['pdb_id'] for x in m})
        print(f'  {name}: {len(s)} 突变 / {n_prot} 蛋白')

    # Ridge 基线（每个基准各自用共同交集计算）
    KD = {'A':1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C':2.5,'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,
          'I':4.5,'L':3.8,'K':-3.9,'M':1.9,'F':2.8,'P':-1.6,'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V':4.2}
    VOL = {'A':88.6,'R':173.4,'N':114.1,'D':111.1,'C':108.5,'Q':143.8,'E':138.4,'G':60.1,'H':153.2,
           'I':166.7,'L':166.7,'K':168.6,'M':162.9,'F':189.9,'P':112.7,'S':89.0,'T':116.1,'W':227.8,'Y':193.6,'V':140.0}
    def feats(w, m):
        return [KD.get(w,0)-KD.get(m,0), VOL.get(w,0)-VOL.get(m,0), KD.get(m,0), VOL.get(m,0)]
    ms = pd.read_csv(os.path.join(DATA, 'megascale_clean.csv'))
    tm = pd.read_csv(os.path.join(DATA, 'thermomutdb_aligned_sc_noleak.csv'))
    tm = tm.assign(wt_aa=tm['mutation_code'].str[0], mt_aa=tm['mutation_code'].str[-1])
    ms = ms.sample(n=len(tm), random_state=42)
    tr_df = pd.concat([ms[['wt_aa','mt_aa','ddg']], tm[['wt_aa','mt_aa','ddg']]], ignore_index=True)
    from sklearn.linear_model import Ridge
    ridge = Ridge(1.0).fit(np.array([feats(w,m) for w,m in zip(tr_df['wt_aa'], tr_df['mt_aa'])]),
                           tr_df['ddg'].values)

    rows = []
    summary = []
    for bname, bd in benchmarks.items():
        y = np.array([m['ddg'] for m in bd['meta']])
        pid = [m['pdb_id'] for m in bd['meta']]
        mid = [m['mut_info'] for m in bd['meta']]
        preds = {}
        # Ridge
        preds['ridge'] = ridge.predict(np.array([feats(m[0], m[-1]) for m in mid]))
        # GNNs
        for seed in SEEDS:
            m1 = train_model(GNNRegressor(in_dim, hid=64), tr_base, args.epochs, seed)
            preds[f'gnn_global_s{seed}'] = predict(m1, bd['base'])
            m2 = train_model(GNNLocal(in_dim, hid=64), tr_base, args.epochs, seed)
            preds[f'gnn_local_s{seed}'] = predict(m2, bd['base'])
            m3 = train_model(GNNLocal(in_dim_bl, hid=64), tr_bl, args.epochs, seed)
            preds[f'gnn_local_blosum_s{seed}'] = predict(m3, bd['blosum'])
            m4 = train_model(GNNEdge(in_dim, edge_dim, hid=64), tr_base, args.epochs, seed)
            preds[f'gnn_edge_s{seed}'] = predict(m4, bd['base'])
            print(f'  [{bname}] seed {seed} 完成')

        for name, p in preds.items():
            r = pearson(y, p)
            _, lo, hi = cluster_bootstrap(y, p, pid, pearson)
            pl, npl = protein_level_pearson(y, p, pid, min_mut=5)
            summary.append({'benchmark': bname, 'model': name, 'n': len(y),
                            'n_proteins': len(set(pid)), 'pearson': r, 'ci_low': lo, 'ci_high': hi,
                            'spearman': spearman(y, p), 'mae': mae(y, p), 'rmse': rmse(y, p),
                            'protein_pearson': pl, 'n_proteins_pl': npl})
            for i in range(len(y)):
                rows.append({'benchmark': bname, 'y_true': y[i], 'y_pred': float(p[i]),
                             'protein_id': pid[i], 'mutation_id': mid[i], 'model': name})

    sm = pd.DataFrame(summary)
    sm.to_csv(os.path.join(DATA, 'benchmark_summary.csv'), index=False)
    pd.DataFrame(rows).to_csv(args.out, index=False)
    print(f'\n已保存 {os.path.join(DATA, "benchmark_summary.csv")} 和 {args.out}')

    # 打印主表
    for bname in benchmarks:
        print(f'\n===== {bname} =====')
        sub = sm[sm['benchmark'] == bname]
        print(f'{"model":<24} {"r":>7} {"95% CI":>18} {"rho":>6} {"MAE":>6} {"RMSE":>6} {"prot-r":>7}')
        for _, r in sub.iterrows():
            print(f'{r["model"]:<24} {r["pearson"]:>7.3f} [{r["ci_low"]:>6.3f},{r["ci_high"]:>6.3f}] '
                  f'{r["spearman"]:>6.3f} {r["mae"]:>6.2f} {r["rmse"]:>6.2f} {r["protein_pearson"]:>7.3f}')


if __name__ == '__main__':
    main()
