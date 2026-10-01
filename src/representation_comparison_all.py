"""
representation_comparison_all.py — 4 种图定义 × 预测性能完整对比
================================================================
统一协议（相同 seed、epoch、训练规模、评估交集），只变图定义。
定义：ca / cb / centroid / allatom
基准：S669（共同交集 511）+ ssym（342）
输出：每种定义的 r 与 95% CI，以及相对质心的配对 Δr。
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
AA_ORDER = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
IDX2AA = {i: a for i, a in enumerate(AA_ORDER)}
SEEDS = [42, 123, 2024]
DEFS = ['ca', 'cb', 'centroid', 'allatom']
SUFFIX = {'ca': 'ca', 'cb': 'cb', 'centroid': 'sc', 'allatom': 'allatom'}


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


def paired_dr(y, p1, p2, pid, B=2000, seed=0):
    rng = np.random.default_rng(seed)
    prot = np.array(pid)
    uniq, inv = np.unique(prot, return_inverse=True)
    pidx = [np.where(inv == i)[0] for i in range(len(uniq))]
    st = []
    for _ in range(B):
        idx = np.concatenate([pidx[i] for i in rng.integers(0, len(uniq), len(uniq))])
        if len(idx) < 5:
            continue
        d = pearson(y[idx], p1[idx]) - pearson(y[idx], p2[idx])
        if not np.isnan(d):
            st.append(d)
    st = np.array(st)
    obs = pearson(y, p1) - pearson(y, p2)
    lo, hi = np.percentile(st, [2.5, 97.5])
    return obs, float(lo), float(hi), float(2 * min((st <= 0).mean(), (st >= 0).mean()))


def make_sample(nodes, ei, ea, mi, mt_aa, y):
    nodes = nodes.copy()
    if mt_aa in AA_INDEX:
        nodes[mi, :20] = 0.0
        nodes[mi, AA_INDEX[mt_aa]] = 1.0
    flag = np.zeros((nodes.shape[0], 1), dtype=np.float32)
    flag[mi, 0] = 1.0
    x = torch.tensor(np.concatenate([nodes, flag], axis=1), dtype=torch.float32)
    return Data(x=x, edge_index=torch.tensor(ei, dtype=torch.long),
                edge_attr=torch.tensor(ea, dtype=torch.float32),
                y=torch.tensor([float(y)], dtype=torch.float32))


def load_train(ad):
    sfx = SUFFIX[ad]
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


def load_test(ad, label_csv, set_prefix):
    sfx = SUFFIX[ad]
    df = pd.read_csv(os.path.join(DATA, label_csv))
    gdir = os.path.join(DATA, f'contact_graphs_{set_prefix}_{sfx}')
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
        meta.append({'pdb_id': pid, 'mut_info': r['mut_info'], 'ddg': r['ddg']})
    return out, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--epochs', type=int, default=20)
    args = ap.parse_args()

    all_preds, summaries, pairs = {}, [], []
    for bench, label_csv, prefix in [('S669', 'benchmarks_s669_clean.csv', 's669'),
                                     ('ssym', 'benchmarks_ssym_clean.csv', 'ssym')]:
        for ad in DEFS:
            tr = load_train(ad)
            te, meta = load_test(ad, label_csv, prefix)
            if not tr or not te:
                print(f'{bench}/{ad}: 数据缺失，跳过')
                continue
            y = np.array([m['ddg'] for m in meta])
            pid = [m['pdb_id'] for m in meta]
            in_dim = tr[0].x.shape[1]
            preds = []
            for s in SEEDS:
                set_seed(s)
                m_ = GNNLocal(in_dim, hid=64)
                opt = torch.optim.Adam(m_.parameters(), lr=1e-3)
                for _ in range(args.epochs):
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
            summaries.append({'benchmark': bench, 'atom_def': ad, 'n': len(y),
                              'n_proteins': len(set(pid)), 'r': r, 'ci_low': lo, 'ci_high': hi,
                              'seed_r': [round(pearson(y, p), 3) for p in preds]})
            all_preds[(bench, ad)] = (y, pm, pid)
            print(f'  {bench}/{ad}: n={len(y)}, r={r:.3f} [{lo:.3f},{hi:.3f}]')

        # 配对 vs 质心
        if ('centroid' in [s['atom_def'] for s in summaries if s['benchmark'] == bench]):
            y, pc, pid = all_preds[(bench, 'centroid')]
            for ad in [d for d in DEFS if d != 'centroid']:
                if (bench, ad) not in all_preds:
                    continue
                _, pa, _ = all_preds[(bench, ad)]
                obs, lo, hi, pv = paired_dr(y, pa, pc, pid)
                pairs.append({'benchmark': bench, 'comparison': f'{ad} vs centroid',
                              'delta_r': obs, 'ci_low': lo, 'ci_high': hi, 'p': pv,
                              'significant': 'yes' if (lo > 0 or hi < 0) else 'no'})
                print(f'    {ad} vs centroid: Δr={obs:+.3f} [{lo:+.3f},{hi:+.3f}] p={pv:.3f}')

    pd.DataFrame(summaries).to_csv(os.path.join(DATA, 'representation_all_defs.csv'), index=False)
    pd.DataFrame(pairs).to_csv(os.path.join(DATA, 'representation_all_pairs.csv'), index=False)
    print('\n已保存 representation_all_defs.csv / representation_all_pairs.csv')


if __name__ == '__main__':
    main()
