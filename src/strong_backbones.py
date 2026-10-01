# -*- coding: utf-8 -*-
"""
strong_backbones.py — 强图编码器（WS1 能力阶梯的高能力端）
=============================================================
提供两个现代骨架，与既有 GNNLocal / GNNEdge 接口兼容
（输入 data.x / edge_index / edge_attr / batch，可选 coords）：

  1. DeepGINE  — 6 层 GINEConv（用边特征）+ 残差 + LayerNorm + 跳跃连接
                 + 注意力池化 + 突变位点多跳局部读出
  2. EGNN      — E(3) 等变图网络，直接使用 3D 坐标（现代结构模型）

设计要点（回应"你的 GNN 太弱"）：
  - 边特征真正参与（GINEConv，非 GCNConv）
  - 深层 + 残差，避免过平滑
  - 注意力池化替代均值池化
  - 局部读出取突变位点及其 k 跳邻域
  - EGNN 版本额外使用坐标几何
"""
import torch
import torch.nn as nn
from torch_geometric.nn import GINEConv, global_mean_pool


def _mlp(i, h, o, layers=2):
    mods, prev = [], i
    for _ in range(layers - 1):
        mods += [nn.Linear(prev, h), nn.ReLU(), nn.LayerNorm(h)]
        prev = h
    mods += [nn.Linear(prev, o)]
    return nn.Sequential(*mods)


# --------------------------------------------------------------------------
# 1) DeepGINE：深层边感知消息传递 + 注意力池化 + 局部读出
# --------------------------------------------------------------------------
class DeepGINE(nn.Module):
    def __init__(self, in_dim, edge_dim=2, hid=128, n_layers=6, out=1,
                 dropout=0.1, k_hop=2):
        super().__init__()
        self.k_hop = k_hop
        self.inp = nn.Linear(in_dim, hid)
        self.edge_enc = nn.Linear(edge_dim, hid)
        self.convs = nn.ModuleList()
        self.norms = nn.ModuleList()
        for _ in range(n_layers):
            self.convs.append(GINEConv(_mlp(hid, hid, hid), edge_dim=hid))
            self.norms.append(nn.LayerNorm(hid))
        self.drop = nn.Dropout(dropout)
        # 注意力池化
        self.att = nn.Sequential(nn.Linear(hid, hid // 2), nn.Tanh(), nn.Linear(hid // 2, 1))
        # 跳跃连接后的聚合
        self.jk = nn.Linear(hid * (n_layers + 1), hid)
        self.head = nn.Sequential(
            nn.Linear(hid * 2, hid), nn.ReLU(), nn.Dropout(dropout), nn.Linear(hid, out))

    def _adj(self, ei, n):
        return torch.sparse_coo_tensor(
            ei, torch.ones(ei.shape[1], device=ei.device), (n, n)).coalesce()

    def _khop_mask(self, ei, n, batch, flag, k):
        """突变位点 k 跳邻域掩码。用 edge_index 直接扩展，避免稀疏矩阵开销。"""
        reach = flag > 0.5
        frontier = reach.clone()
        for _ in range(k):
            if not bool(frontier.any()):
                break
            src_in = frontier[ei[0]]                 # 哪些边从 frontier 出发
            nxt = torch.zeros(n, dtype=torch.bool, device=ei.device)
            nxt[ei[1][src_in]] = True
            nxt &= ~reach
            reach = reach | nxt
            frontier = nxt
        return reach

    def embed(self, data):
        """返回拼接后的图级表示 (glob | loc)，供融合模型使用。"""
        x, ei, ea, batch = data.x, data.edge_index, data.edge_attr, data.batch
        h = self.inp(x)
        ae = self.edge_enc(ea)
        hs = [h]
        for conv, norm in zip(self.convs, self.norms):
            h2 = conv(h, ei, ae)
            h = norm(h + self.drop(h2.relu()))       # 残差 + LayerNorm
            hs.append(h)
        h = self.jk(torch.cat(hs, dim=1))            # 跳跃连接
        w = self.att(h).clamp(-20.0, 20.0)
        wexp = torch.exp(w)
        denom = torch.zeros(data.num_graphs, 1, device=h.device).index_add_(0, batch, wexp)
        alpha = wexp / (denom[batch] + 1e-9)
        glob = torch.zeros(data.num_graphs, h.shape[1], device=h.device).index_add_(
            0, batch, h * alpha)
        mask = self._khop_mask(ei, h.shape[0], batch, data.x[:, -1], self.k_hop)
        m = mask.float().unsqueeze(1)
        loc = torch.zeros(data.num_graphs, h.shape[1], device=h.device).index_add_(
            0, batch, h * m)
        cnt = torch.zeros(data.num_graphs, 1, device=h.device).index_add_(0, batch, m)
        loc = loc / (cnt + 1e-9)
        bad = (cnt.squeeze(1) < 1e-6)
        if bad.any():
            loc[bad] = glob[bad]
        return torch.cat([glob, loc], dim=1)

    def forward(self, data):
        return self.head(self.embed(data)).squeeze(-1)


# --------------------------------------------------------------------------
# 2) EGNN：E(3) 等变图网络，直接使用 3D 坐标
# --------------------------------------------------------------------------
class EGNNLayer(nn.Module):
    def __init__(self, hid, edge_dim, scaled=True):
        super().__init__()
        self.edge_mlp = nn.Sequential(
            nn.Linear(hid * 2 + 1 + edge_dim, hid), nn.SiLU(), nn.Linear(hid, hid))
        self.node_mlp = nn.Sequential(
            nn.Linear(hid * 2, hid), nn.SiLU(), nn.Linear(hid, hid))
        self.coord_mlp = nn.Sequential(nn.Linear(hid, hid), nn.SiLU(), nn.Linear(hid, 1))
        self.norm = nn.LayerNorm(hid)
        # ★ 坐标更新缩放：scaled=True 时初始为 0（坐标初始固定，逐步学习位移），
        #   scaled=False 复现原实现（位移无缩放），用于稳定性对照
        self.scaled = scaled
        if scaled:
            self.coord_scale = nn.Parameter(torch.zeros(1))

    def forward(self, h, x, ei, ea):
        i, j = ei[0], ei[1]
        d2 = ((x[i] - x[j]) ** 2).sum(-1, keepdim=True)
        m = self.edge_mlp(torch.cat([h[i], h[j], d2, ea], dim=1))
        w = self.coord_mlp(m)
        dx = torch.zeros_like(x).index_add_(0, i, (x[i] - x[j]) * w)
        denom = torch.zeros(x.shape[0], 1, device=x.device).index_add_(0, i, w.abs()) + 1e-9
        step = self.coord_scale * (dx / denom) if self.scaled else (dx / denom)
        x = x + step
        agg = torch.zeros_like(h).index_add_(0, i, m)
        h = self.norm(h + self.node_mlp(torch.cat([h, agg], dim=1)))
        return h, x


class EGNN(nn.Module):
    def __init__(self, in_dim, edge_dim=2, hid=128, n_layers=4, out=1,
                 dropout=0.1, k_hop=2, update_coords=True, scaled=True):
        super().__init__()
        self.k_hop = k_hop
        self.update_coords = update_coords
        self.inp = nn.Linear(in_dim, hid)
        self.layers = nn.ModuleList([EGNNLayer(hid, edge_dim, scaled=scaled)
                                     for _ in range(n_layers)])
        self.drop = nn.Dropout(dropout)
        self.att = nn.Sequential(nn.Linear(hid, hid // 2), nn.Tanh(), nn.Linear(hid // 2, 1))
        self.head = nn.Sequential(
            nn.Linear(hid * 2, hid), nn.ReLU(), nn.Dropout(dropout), nn.Linear(hid, out))

    def _khop_mask(self, ei, n, flag, k):
        """突变位点 k 跳邻域掩码。用 edge_index 直接扩展，避免稀疏矩阵开销。"""
        reach = flag > 0.5
        frontier = reach.clone()
        for _ in range(k):
            if not bool(frontier.any()):
                break
            src_in = frontier[ei[0]]
            nxt = torch.zeros(n, dtype=torch.bool, device=ei.device)
            nxt[ei[1][src_in]] = True
            nxt &= ~reach
            reach = reach | nxt
            frontier = nxt
        return reach

    def forward(self, data):
        h = self.inp(data.x)
        x = data.coords.float()
        ei, ea, batch = data.edge_index, data.edge_attr, data.batch
        for layer in self.layers:
            if self.update_coords:
                h2, x2 = layer(h, x, ei, ea)
                h = self.drop(h2)
                x = x2
            else:
                h2, _ = layer(h, x, ei, ea)
                h = self.drop(h2)
        # 注意力池化（确定性实现）
        w = self.att(h).clamp(-20.0, 20.0)
        wexp = torch.exp(w)
        denom = torch.zeros(data.num_graphs, 1, device=h.device).index_add_(0, batch, wexp)
        alpha = wexp / (denom[batch] + 1e-9)
        glob = torch.zeros(data.num_graphs, h.shape[1], device=h.device).index_add_(
            0, batch, h * alpha)
        mask = self._khop_mask(ei, h.shape[0], data.x[:, -1], self.k_hop)
        m = mask.float().unsqueeze(1)
        loc = torch.zeros(data.num_graphs, h.shape[1], device=h.device).index_add_(
            0, batch, h * m)
        cnt = torch.zeros(data.num_graphs, 1, device=h.device).index_add_(0, batch, m)
        loc = loc / (cnt + 1e-9)
        bad = (cnt.squeeze(1) < 1e-6)
        if bad.any():
            loc[bad] = glob[bad]
        return self.head(torch.cat([glob, loc], dim=1)).squeeze(-1)


BACKBONES = {
    'deep_gine': DeepGINE,
    'egnn': EGNN,
}


# --------------------------------------------------------------------------
# 3) 候选方法：按边类型分解的边缘池化（Type-Resolved Edge Pooling, TREP）
# --------------------------------------------------------------------------
class TypePoolGINE(torch.nn.Module):
    """
    动机（来自本项目实测）：
      - 边类型（疏水/静电）在**单变量**层面与 ΔΔG 显著相关（broken_hydro r=0.135）
      - 但带边特征的 GINE 模型（gnn_edge）**并不优于**无实体边特征的局部池化模型
      → 说明常规**节点级**消息传递把边类型信息稀释了
    做法：在突变位点的邻域，**按边类型分别池化边缘表示**，
          再与节点级全局/局部表示拼接。
    """

    def __init__(self, in_dim, edge_dim=2, hid=128, n_layers=4, out=1,
                 dropout=0.1, k_hop=2, n_types=3):
        super().__init__()
        self.k_hop = k_hop
        self.n_types = n_types
        self.inp = nn.Linear(in_dim, hid)
        self.edge_enc = nn.Linear(edge_dim, hid)
        self.convs = nn.ModuleList()
        self.norms = nn.ModuleList()
        for _ in range(n_layers):
            self.convs.append(GINEConv(_mlp(hid, hid, hid), edge_dim=hid))
            self.norms.append(nn.LayerNorm(hid))
        self.drop = nn.Dropout(dropout)
        self.att = nn.Sequential(nn.Linear(hid, hid // 2), nn.Tanh(), nn.Linear(hid // 2, 1))
        self.jk = nn.Linear(hid * (n_layers + 1), hid)
        # 头：全局 + 局部节点 + 按类型的边缘池化
        self.head = nn.Sequential(
            nn.Linear(hid * (2 + n_types), hid), nn.ReLU(),
            nn.Dropout(dropout), nn.Linear(hid, out))

    def _khop_mask(self, ei, n, flag, k):
        reach = flag > 0.5
        frontier = reach.clone()
        for _ in range(k):
            if not bool(frontier.any()):
                break
            src_in = frontier[ei[0]]
            nxt = torch.zeros(n, dtype=torch.bool, device=ei.device)
            nxt[ei[1][src_in]] = True
            nxt &= ~reach
            reach = reach | nxt
            frontier = nxt
        return reach

    def forward(self, data):
        x, ei, ea, batch = data.x, data.edge_index, data.edge_attr, data.batch
        n = x.shape[0]
        h = self.inp(x)
        ae = self.edge_enc(ea)
        hs = [h]
        for conv, norm in zip(self.convs, self.norms):
            h2 = conv(h, ei, ae)
            h = norm(h + self.drop(h2.relu()))
            hs.append(h)
        h = self.jk(torch.cat(hs, dim=1))
        # 节点级池化
        w = self.att(h).clamp(-20.0, 20.0)
        wexp = torch.exp(w)
        denom = torch.zeros(data.num_graphs, 1, device=h.device).index_add_(0, batch, wexp)
        alpha = wexp / (denom[batch] + 1e-9)
        glob = torch.zeros(data.num_graphs, h.shape[1], device=h.device).index_add_(
            0, batch, h * alpha)
        mask = self._khop_mask(ei, n, data.x[:, -1], self.k_hop)
        m = mask.float().unsqueeze(1)
        loc = torch.zeros(data.num_graphs, h.shape[1], device=h.device).index_add_(
            0, batch, h * m)
        cnt = torch.zeros(data.num_graphs, 1, device=h.device).index_add_(0, batch, m)
        loc = loc / (cnt + 1e-9)
        bad = (cnt.squeeze(1) < 1e-6)
        if bad.any():
            loc[bad] = glob[bad]
        # ★ 按边类型池化：只保留与突变位点（k 跳邻域）相连的边
        nb = mask[ei[0]]
        src_batch = batch[ei[0]]
        edge_h = h[ei[0]]                       # 源节点表示作为边的表示
        # 类型：hydro / elec / other（可同时为 hydro 与 elec 时归 hydro）
        is_hy = (ea[:, 0] > 0.5).float()
        is_el = ((ea[:, 1] > 0.5) & (ea[:, 0] <= 0.5)).float()
        is_ot = ((ea[:, 0] <= 0.5) & (ea[:, 1] <= 0.5)).float()
        type_pool = []
        for tmask in (is_hy, is_el, is_ot):
            wgt = (nb.float() * tmask).unsqueeze(1)
            acc = torch.zeros(data.num_graphs, h.shape[1], device=h.device).index_add_(
                0, src_batch, edge_h * wgt)
            den = torch.zeros(data.num_graphs, 1, device=h.device).index_add_(
                0, src_batch, wgt) + 1e-9
            type_pool.append(acc / den)
        return self.head(torch.cat([glob, loc] + type_pool, dim=1)).squeeze(-1)
