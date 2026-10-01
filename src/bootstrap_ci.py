"""
bootstrap_ci.py — 按蛋白分层的 paired bootstrap 置信区间
==========================================================
落地 统计与评估方案.md 的核心统计逻辑：
  - 以"蛋白"为抽样单元（而非突变），处理突变嵌套在蛋白内的非独立性。
  - paired bootstrap：每次重采样蛋白，算 GEDMut 与对照的配对指标差 Δ，得到 Δ 的分布。
  - 输出 mutation-level 与 protein-level 两个汇总层次。

纯 numpy 实现，可单元测试。

核心函数：
  - bootstrap_paired(y_true, y_a, y_b, protein_ids, metric, B) -> CI
      比较方法 A vs B：对配对指标差 Δ=metric(A)-metric(B) 做蛋白分层 bootstrap。
  - pearson, mae, rmse
  - protein_level_summary(y_true, y_pred, protein_ids, metric) -> 每蛋白先汇总再跨蛋白
"""
import numpy as np


# ----------------------------------------------------------------------
# 指标
# ----------------------------------------------------------------------
def pearson(y_true, y_pred):
    t = np.asarray(y_true, dtype=float)
    p = np.asarray(y_pred, dtype=float)
    if t.std() == 0 or p.std() == 0:
        return float('nan')
    return float(np.corrcoef(t, p)[0, 1])


def mae(y_true, y_pred):
    return float(np.mean(np.abs(np.asarray(y_true) - np.asarray(y_pred))))


def rmse(y_true, y_pred):
    return float(np.sqrt(np.mean((np.asarray(y_true) - np.asarray(y_pred)) ** 2)))


def _pearson_corrcoef(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.std() == 0 or b.std() == 0:
        return float('nan')
    return float(np.corrcoef(a, b)[0, 1])


# ----------------------------------------------------------------------
# 蛋白分层
# ----------------------------------------------------------------------
def group_by_protein(protein_ids):
    """返回 {蛋白id: 该蛋白所有样本的下标列表}。"""
    groups = {}
    for i, pid in enumerate(protein_ids):
        groups.setdefault(pid, []).append(i)
    return groups


# ----------------------------------------------------------------------
# paired bootstrap（核心）
# ----------------------------------------------------------------------
def bootstrap_paired(y_true, y_a, y_b, protein_ids, metric='pearson', B=1000, seed=0):
    """
    比较方法 A 与 B：对配对指标差 Δ = metric(A) − metric(B) 做蛋白分层 bootstrap。

    返回 dict:
      metric: 指标名
      delta_point: 全样本上的点估计 Δ
      ci_low, ci_high: Δ 的 95% 百分位置信区间
      n_proteins: 蛋白数
      n_samples: 样本数（突变数）
      B: bootstrap 次数
    """
    y_true = np.asarray(y_true, dtype=float)
    y_a = np.asarray(y_a, dtype=float)
    y_b = np.asarray(y_b, dtype=float)
    protein_ids = list(protein_ids)
    groups = group_by_protein(protein_ids)
    protein_list = list(groups.keys())
    n_proteins = len(protein_list)
    if n_proteins == 0:
        return None

    metric_fn = {'pearson': pearson, 'mae': mae, 'rmse': rmse}[metric]

    def delta_on_indices(idx):
        if len(idx) == 0:
            return float('nan')
        a = metric_fn(y_true[idx], y_a[idx])
        b = metric_fn(y_true[idx], y_b[idx])
        return a - b

    all_idx = list(range(len(y_true)))
    delta_point = delta_on_indices(all_idx)

    rng = np.random.default_rng(seed)
    deltas = []
    for _ in range(B):
        # 以蛋白为单元重采样
        sampled_proteins = rng.choice(protein_list, size=n_proteins, replace=True)
        idx = []
        for p in sampled_proteins:
            idx.extend(groups[p])
        d = delta_on_indices(idx)
        if not np.isnan(d):
            deltas.append(d)
    deltas = np.array(deltas)
    ci_low, ci_high = np.percentile(deltas, [2.5, 97.5])

    return {
        'metric': metric,
        'delta_point': delta_point,
        'ci_low': float(ci_low),
        'ci_high': float(ci_high),
        'n_proteins': n_proteins,
        'n_samples': len(y_true),
        'B': B,
    }


def protein_level_summary(y_true, y_pred, protein_ids, metric='pearson'):
    """
    每蛋白先算蛋白内指标，再跨蛋白取均值（避免大蛋白主导）。
    对 pearson 用 Fisher z 变换平均后反变换（更稳）。
    """
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    groups = group_by_protein(protein_ids)
    if metric == 'pearson':
        zs = []
        for pid, idx in groups.items():
            r = _pearson_corrcoef(y_true[idx], y_pred[idx])
            if not np.isnan(r):
                r = np.clip(r, -0.9999, 0.9999)
                zs.append(np.arctanh(r))
        if not zs:
            return float('nan')
        return float(np.tanh(np.mean(zs)))
    else:
        metric_fn = {'mae': mae, 'rmse': rmse}[metric]
        per_protein = [metric_fn(y_true[idx], y_pred[idx]) for idx in groups.values()]
        return float(np.mean(per_protein))


def report(y_true, y_a, y_b, protein_ids, metric='pearson', B=1000, seed=0):
    """生成一段可复制的报告文本。"""
    bp = bootstrap_paired(y_true, y_a, y_b, protein_ids, metric=metric, B=B, seed=seed)
    sig = "显著" if (bp['ci_low'] > 0 or bp['ci_high'] < 0) else "不显著"
    print(f"[{metric}] Δ(A−B) = {bp['delta_point']:.4f} "
          f"(95% CI: {bp['ci_low']:.4f}, {bp['ci_high']:.4f}) → {sig}")
    print(f"  蛋白数={bp['n_proteins']}, 样本数={bp['n_samples']}, B={bp['B']}")
    return bp


if __name__ == '__main__':
    print("bootstrap_ci.py 导入成功。请用 test_bootstrap_ci.py 做单元测试。")
