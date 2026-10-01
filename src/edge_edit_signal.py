"""
edge_edit_signal.py — 边编辑预测的信号诊断
============================================
回答关键问题：边编辑（断边 vs 保持边）是否有可预测的信号？

构造边级别数据集：
  每个样本 = WT 图中的一条边
  特征：图距离（突变位点到边的 hop 数）、接触类型、突变体积变化
  标签：这条边是否断（1=断，0=保持）

用逻辑回归 baseline，看 AUROC 是否 > 0.5。
若明显 > 0.5 → 边编辑预测有信号，GEDMut 边编辑归因可学习。
若 ≈ 0.5 → 信号太弱，需重新审视。
"""
import ast
import os
import sys
from collections import deque

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

from pathlib import Path as _P
DATA = str(_P(__file__).resolve().parent.parent / 'data')
VOL = {'A':88.6,'R':173.4,'N':114.1,'D':111.1,'C':108.5,'Q':143.8,'E':138.4,'G':60.1,'H':153.2,
       'I':166.7,'L':166.7,'K':168.6,'M':162.9,'F':189.9,'P':112.7,'S':89.0,'T':116.1,'W':227.8,'Y':193.6,'V':140.0}


def parse_edges(s):
    return ast.literal_eval(s) if isinstance(s, str) and s else []


def graph_distance(edge_index, n, src):
    """BFS 计算 src 到所有节点的图距离。"""
    adj = [[] for _ in range(n)]
    for a, b in edge_index.T:
        adj[int(a)].append(int(b))
        adj[int(b)].append(int(a))
    dist = [-1] * n
    dist[src] = 0
    q = deque([src])
    while q:
        u = q.popleft()
        for v in adj[u]:
            if dist[v] == -1:
                dist[v] = dist[u] + 1
                q.append(v)
    return dist


def main():
    ed = pd.read_csv(os.path.join(DATA, 'true_edits_s669_sc.csv'))
    X, y = [], []
    for _, row in ed.iterrows():
        pid = row['pdb_id']
        mi = int(row['mut_idx'])
        # 读 WT 接触图
        gpath = os.path.join(DATA, 'contact_graphs_s669_sc', f'{pid}.npz')
        if not os.path.exists(gpath):
            continue
        d = np.load(gpath, allow_pickle=True)
        nodes = d['nodes']
        ei = d['edge_index']
        n = nodes.shape[0]
        dist = graph_distance(ei, n, mi)

        broken = {tuple(sorted(e)) for e in parse_edges(row['broken_edges'])}
        # WT 边集（从 edge_index 提取，去重无向）
        wt_edges = set()
        for a, b in ei.T:
            wt_edges.add(tuple(sorted((int(a), int(b)))))

        # 突变体积变化
        wt_aa = str(row['mut_info'])[0]
        mt_aa = str(row['mut_info'])[-1]
        vol_diff = abs(VOL.get(wt_aa, 100) - VOL.get(mt_aa, 100))

        for (a, b) in wt_edges:
            # 特征：到突变位点的最小图距离、体积变化
            dmin = min(dist[a], dist[b]) if (dist[a] >= 0 and dist[b] >= 0) else 10
            X.append([dmin, vol_diff])
            y.append(1 if (a, b) in broken else 0)

    X = np.array(X)
    y = np.array(y)
    print(f'边级别样本: {len(y)} 条（断 {y.sum()}，保持 {len(y)-y.sum()}）')
    print(f'断边比例: {y.mean()*100:.1f}%')

    # 逻辑回归，5 折交叉验证
    from sklearn.model_selection import cross_val_score
    clf = LogisticRegression(max_iter=1000)
    aucs = cross_val_score(clf, X, y, cv=5, scoring='roc_auc')
    print(f'逻辑回归 AUROC: {aucs.mean():.3f} ± {aucs.std():.3f}')

    # 单特征 AUROC
    for i, name in [(0, '图距离'), (1, '体积变化')]:
        auc = roc_auc_score(y, X[:, i])
        print(f'  单特征 {name} AUROC: {auc:.3f}')


if __name__ == '__main__':
    main()
