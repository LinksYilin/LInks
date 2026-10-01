"""
threshold_sensitivity_prediction.py — 预测层的阈值敏感性（第 5 项）
===================================================================
固定侧链质心定义，在 6/7/8/9/10 Å 上各训练一个相同 GNN，
评估 S669（511）与 ssym（342），检验结论对阈值是否稳健。
"""
import os
import sys

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader

sys.path.insert(0, os.path.dirname(__file__))
from seed_utils import set_seed  # 完整确定性控制（含 cuBLAS/cudnn 标志）

SEED_DEFAULT = 42

from gnn_local_baseline import GNNLocal
from paths import DATA, FIGURES, PUB_FIGURES, ROOT, SRC, TOOLS  # 路径集中管理
from train_gnn_baseline import AA_INDEX

SRC_PATH = str(SRC)
DATA_PATH = str(DATA)
ROOT_PATH = str(ROOT)
FIG_PATH = str(FIGURES)
PUB_PATH = str(PUB_FIGURES)
TOOLS_PATH = str(TOOLS)

DATA = DATA_PATH  # 来自 paths.py，可用 GED_ROOT 环境变量覆盖
SEEDS = [42, 123, 2024]
THRESHOLDS = [6.0, 7.0, 8.0, 9.0, 10.0]


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


def suffix(th):
    return 'sc' if abs(th - 8.0) < 1e-6 else f'centroid{th:g}A'


def pearson(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


def cluster_ci(y, p, pid, B=1000, seed=0):
    rng = np.random.default_rng(seed)
    prot = np.array(pid)
    uniq, inv = np.unique(prot, return_inverse=True)
    pidx = [np.where(inv == i)[0] for i in range(len(uniq))]
    st = []
    for _ in range(B):
        idx = np.concatenate([pidx[i] for i in rng.integers(0, len(uniq), len(uniq))])
        if len(idx) >= 3:
            v = pearson(y[idx], p[idx])
            if not np.isnan(v):
                st.append(v)
    return (float(np.percentile(st, 2.5)), float(np.percentile(st, 97.5))) if st else (np.nan, np.nan)


def make_sample(nodes, ei, ea, mi, mt, y):
    nodes = nodes.copy()
    if mt in AA_INDEX:
        nodes[mi, :20] = 0.0
        nodes[mi, AA_INDEX[mt]] = 1.0
    flag = np.zeros((nodes.shape[0], 1), dtype=np.float32)
    flag[mi, 0] = 1.0
    return Data(x=torch.tensor(np.concatenate([nodes, flag], axis=1), dtype=torch.float32),
                edge_index=torch.tensor(ei, dtype=torch.long),
                edge_attr=torch.tensor(ea, dtype=torch.float32),
                y=torch.tensor([float(y)], dtype=torch.float32))


def load_train(th):
    sfx = suffix(th)
    tr = pd.read_csv(os.path.join(DATA, 'training_merged_noleak_sc.csv'))
    ms = tr[tr['source'] == 'megascale']; tm = tr[tr['source'] == 'thermomutdb']
    n = min(len(ms), len(tm))
    tr = pd.concat([ms.sample(n=n, random_state=42), tm], ignore_index=True)
    dirs = {'megascale': os.path.join(DATA, f'contact_graphs_megascale_{sfx}'),
            'thermomutdb': os.path.join(DATA, f'contact_graphs_thermomutdb_{sfx}')}
    cache, out = {}, []
    for _, r in tr.iterrows():
        p = os.path.join(dirs[r['source']], f"{r['protein']}.npz")
        if not os.path.exists(p):
            log_skip('load', r.get('pdb_id', r.get('protein', '?')),
                     r.get('mut_info', r.get('mt_aa', '?')), 'graph_missing')
            continue
        k = f"{r['source']}:{r['protein']}"
        if k not in cache:
            d = np.load(p, allow_pickle=True)
            cache[k] = (d['nodes'], d['edge_index'], d['edge_attr'])
        nd, ei, ea = cache[k]
        mi = int(r['mut_idx'])
        if mi < 0 or mi >= nd.shape[0]:
            log_skip('index', r.get('pdb_id', r.get('protein', '?')),
                     r.get('mut_info', r.get('mt_aa', '?')),
                     f'index_out_of_range({mi})')
            continue
        out.append(make_sample(nd, ei, ea, mi, str(r['mt_aa']), r['ddg']))
    return out


def load_test(th, label_csv, prefix):
    sfx = suffix(th)
    df = pd.read_csv(os.path.join(DATA, label_csv))
    gdir = os.path.join(DATA, f'contact_graphs_{prefix}_{sfx}')
    cache, out, meta = {}, [], []
    for _, r in df.iterrows():
        if '_node_idx_ok' in df.columns and not bool(r['_node_idx_ok']):
            continue
        pid = r['pdb_id']
        p = os.path.join(gdir, f'{pid}.npz')
        if not os.path.exists(p):
            log_skip('load', r.get('pdb_id', r.get('protein', '?')),
                     r.get('mut_info', r.get('mt_aa', '?')), 'graph_missing')
            continue
        if pid not in cache:
            d = np.load(p, allow_pickle=True)
            cache[pid] = (d['nodes'], d['edge_index'], d['edge_attr'])
        nd, ei, ea = cache[pid]
        if '_node_idx' not in df.columns or pd.isna(r['_node_idx']):
            raise ValueError(
                f"缺少已验证的 _node_idx（{(r['pdb_id'], r['mut_info'])}）；"
                "拒绝回退到 _pdb_res_idx，请先运行 fix_index_bug.py")
        mi = int(r['_node_idx'])
        if mi < 0 or mi >= nd.shape[0]:
            log_skip('index', r.get('pdb_id', r.get('protein', '?')),
                     r.get('mut_info', r.get('mt_aa', '?')),
                     f'index_out_of_range({mi})')
            continue
        out.append(make_sample(nd, ei, ea, mi, str(r['mut_info'])[-1], r['ddg']))
        meta.append({'pdb_id': pid, 'ddg': r['ddg']})
    return out, meta


def main():
    rows = []
    for bench, lc, prefix in [('S669', 'benchmarks_s669_clean.csv', 's669'),
                              ('ssym', 'benchmarks_ssym_clean.csv', 'ssym')]:
        for th in THRESHOLDS:
            tr = load_train(th)
            te, meta = load_test(th, lc, prefix)
            if not tr or not te:
                print(f'{bench}/{th}: 数据缺失')
                continue
            y = np.array([m['ddg'] for m in meta])
            pid = [m['pdb_id'] for m in meta]
            in_dim = tr[0].x.shape[1]
            preds = []
            for s in SEEDS:
                set_seed(s)
                m_ = GNNLocal(in_dim, hid=64)
                opt = torch.optim.Adam(m_.parameters(), lr=1e-3)
                for _ in range(20):
                    m_.train()
                    for b in DataLoader(tr, batch_size=64, shuffle=True):
                        opt.zero_grad(); F.mse_loss(m_(b), b.y).backward(); opt.step()
                m_.eval()
                pr = []
                with torch.no_grad():
                    for b in DataLoader(te, batch_size=64):
                        pr.append(m_(b).numpy())
                preds.append(np.concatenate(pr))
            pm = np.mean(preds, axis=0)
            r = pearson(y, pm)
            lo, hi = cluster_ci(y, pm, pid)
            rows.append({'benchmark': bench, 'threshold': th, 'n': len(y),
                         'mean_edges': float(np.mean([s.edge_index.shape[1] // 2 for s in te])),
                         'r': r, 'ci_low': lo, 'ci_high': hi})
            print(f'  {bench} {th:g} Å: n={len(y)}, r={r:.3f} [{lo:.3f},{hi:.3f}]')

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(DATA, 'threshold_sensitivity_prediction.csv'), index=False)
    print()
    for bench in df['benchmark'].unique():
        sub = df[df['benchmark'] == bench]
        print(f'=== {bench} ===')
        for _, r in sub.iterrows():
            print(f'  {r["threshold"]:>4.0f} Å: r={r["r"]:.3f} [{r["ci_low"]:.3f},{r["ci_high"]:.3f}]')
    print('\n已保存 threshold_sensitivity_prediction.csv')


if __name__ == '__main__':
    main()
