"""Generate publication-ready figures from the audited project outputs.

The script uses repository-relative paths, preserves source observations, and
exports editable PDF/SVG plus high-resolution TIFF/PNG and grayscale previews.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
DEFAULT_OUT = ROOT / "figures" / "publication"

# Okabe-Ito-derived colors with redundant marker encodings.
COLORS = {
    "linear": "#222222",
    "global": "#7A7A7A",
    "local": "#0072B2",
    "blosum": "#D55E00",
    "cb": "#CC79A7",
    "edge": "#009E73",
    "broken": "#D55E00",
    "formed": "#0072B2",
    "accent": "#009E73",
    "grid": "#D9D9D9",
}
MARKERS = {"linear": "D", "global": "o", "local": "s", "blosum": "^", "cb": "v", "edge": "P"}
FAMILY_ORDER = ["linear", "global", "cb", "local", "blosum", "edge"]
FAMILY_LABELS = {
    "linear": "Physicochemical\nridge",
    "global": "GNN global",
    "cb": "GNN local\n(Cβ)",
    "local": "GNN local\n(centroid)",
    "blosum": "GNN local\n+ BLOSUM",
    "edge": "GNN edge-aware\n(GINEConv)",
}
# 面板 b–d 的短标签（必须与 FAMILY_ORDER 一一对应）
TICK_LABELS = {
    "linear": "Ridge",
    "global": "Global",
    "cb": "C\u03b2",
    "local": "Cent.",
    "blosum": "+B62",
    "edge": "Edge",
}


def setup_style() -> None:
    matplotlib.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 7,
            "axes.labelsize": 7.5,
            "axes.titlesize": 8,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "legend.fontsize": 6.5,
            "axes.linewidth": 0.7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "xtick.major.size": 3,
            "ytick.major.size": 3,
            "lines.linewidth": 1.2,
            "patch.linewidth": 0.7,
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def add_panel_labels(fig: plt.Figure, axes, labels: str) -> None:
    for ax, label in zip(np.ravel(axes), labels):
        ax.annotate(
            label,
            xy=(0, 1),
            xycoords="axes fraction",
            xytext=(-20, 9),
            textcoords="offset points",
            fontsize=9,
            fontweight="bold",
            ha="left",
            va="bottom",
            annotation_clip=False,
        )


def export_figure(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=300, bbox_inches="tight", pad_inches=0.03)
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight", pad_inches=0.03)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.03)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.03)

    # A grayscale preview is a review artifact, not the primary publication file.
    from PIL import Image

    with Image.open(stem.with_suffix(".png")) as image:
        image.convert("L").save(stem.parent / f"{stem.name}_grayscale.png")


def family_for(model: str) -> str:
    if model in ("linear", "ridge"):
        return "linear"
    if model.startswith("gnn_global"):
        return "global"
    if model.startswith("gnn_local_blosum"):
        return "blosum"
    if model.startswith("gnn_local_cb"):
        return "cb"
    if model.startswith("gnn_local"):
        return "local"
    if model.startswith("gnn_edge"):
        return "edge"
    raise ValueError(f"Unrecognized model name: {model}")


def pearson(a, b) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def protein_cluster_ci(
    y_true, y_pred, protein_ids, *, n_boot: int = 1000, seed: int = 2026
) -> tuple[float, float]:
    """Percentile CI from paired cluster bootstrap with repeated clusters retained."""
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    proteins = np.asarray(protein_ids)
    unique, inverse = np.unique(proteins, return_inverse=True)
    grouped_indices = [np.flatnonzero(inverse == i) for i in range(len(unique))]
    # Randomness is used only to resample observed protein clusters for the bootstrap.
    rng = np.random.RandomState(seed)
    values = []
    for _ in range(n_boot):
        sampled_clusters = rng.randint(0, len(unique), size=len(unique))
        sampled_rows = np.concatenate([grouped_indices[i] for i in sampled_clusters])
        value = pearson(y_true[sampled_rows], y_pred[sampled_rows])
        if np.isfinite(value):
            values.append(value)
    if not values:
        return np.nan, np.nan
    return tuple(np.percentile(values, [2.5, 97.5]))


def protein_level_pearson(y_true, y_pred, protein_ids, min_mut: int = 5) -> tuple[float, int]:
    frame = pd.DataFrame({"y_true": y_true, "y_pred": y_pred, "protein_id": protein_ids})
    correlations = []
    for _, group in frame.groupby("protein_id", sort=False):
        if len(group) < min_mut:
            continue
        value = pearson(group["y_true"], group["y_pred"])
        if np.isfinite(value):
            correlations.append(np.clip(value, -0.999, 0.999))
    if not correlations:
        return np.nan, 0
    return float(np.tanh(np.mean(np.arctanh(correlations)))), len(correlations)


def model_metrics(predictions: pd.DataFrame, n_boot: int) -> pd.DataFrame:
    rows = []
    for model, group in predictions.groupby("model", sort=False):
        y_true = group["y_true"].to_numpy()
        y_pred = group["y_pred"].to_numpy()
        lo, hi = protein_cluster_ci(
            y_true, y_pred, group["protein_id"].to_numpy(), n_boot=n_boot
        )
        protein_r, n_proteins = protein_level_pearson(
            y_true, y_pred, group["protein_id"].to_numpy()
        )
        match = re.search(r"_s(\d+)$", model)
        rows.append(
            {
                "model": model,
                "family": family_for(model),
                "seed": int(match.group(1)) if match else -1,
                "pearson": pearson(y_true, y_pred),
                "pearson_ci_low": lo,
                "pearson_ci_high": hi,
                "spearman": float(spearmanr(y_true, y_pred).statistic),
                "mae": float(np.mean(np.abs(y_true - y_pred))),
                "rmse": float(np.sqrt(np.mean((y_true - y_pred) ** 2))),
                "protein_pearson": protein_r,
                "n_proteins": n_proteins,
                "n_mutations": len(group),
            }
        )
    return pd.DataFrame(rows)


def jitter_positions(center: float, count: int) -> np.ndarray:
    if count == 1:
        return np.array([center])
    return center + np.linspace(-0.12, 0.12, count)


def figure1(metrics: pd.DataFrame, out_dir: Path) -> None:
    fig = plt.figure(figsize=(7.2, 3.25), constrained_layout=True)
    grid = fig.add_gridspec(1, 4, width_ratios=[2.2, 1, 1, 1.15], wspace=0.35)
    axes = [fig.add_subplot(grid[0, i]) for i in range(4)]

    ax = axes[0]
    # 动态行序：按家族顺序 + seed
    fam_rank = {f: i for i, f in enumerate(FAMILY_ORDER)}
    met = metrics.copy()
    met["_fam_rank"] = met["family"].map(fam_rank).fillna(99)
    met = met.sort_values(["_fam_rank", "model"], kind="stable")
    row_order = met["model"].tolist()
    by_model = metrics.set_index("model").loc[row_order].reset_index()
    y = np.arange(len(by_model))[::-1]
    labels = []
    for ypos, row in zip(y, by_model.itertuples(index=False)):
        family = row.family
        ax.errorbar(
            row.pearson,
            ypos,
            xerr=[[row.pearson - row.pearson_ci_low], [row.pearson_ci_high - row.pearson]],
            fmt=MARKERS[family],
            ms=4.2,
            mfc=COLORS[family],
            mec="white",
            mew=0.35,
            ecolor=COLORS[family],
            elinewidth=0.8,
            capsize=1.7,
            zorder=3,
        )
        if row.model in ("linear", "ridge"):
            labels.append("Physicochemical ridge")
        else:
            labels.append(f"{FAMILY_LABELS[family].replace(chr(10), ' ')}  |  seed {row.seed}")
    ax.axvline(0, color=COLORS["grid"], lw=0.7, zorder=0)
    for boundary in [8.5, 5.5, 2.5]:
        ax.axhline(boundary, color="#E8E8E8", lw=0.6, zorder=0)
    ax.set_yticks(y, labels)
    ax.set_xlim(-0.05, 0.58)
    ax.set_xlabel("Mutation-level Pearson $r$")
    ax.set_title("Pearson $r$ with 95% CI", loc="left", pad=7)
    ax.tick_params(axis="y", labelsize=5.8, pad=3)
    for tick_label, family in zip(ax.get_yticklabels(), by_model["family"]):
        tick_label.set_color(COLORS[family])
        tick_label.set_fontweight("bold" if family == "linear" else "normal")

    for ax, column, title, ylim in [
        (axes[1], "spearman", "Rank correlation", (0.05, 0.48)),
        (axes[2], "rmse", "Prediction error", (1.36, 1.64)),
        (axes[3], "protein_pearson", "Within-protein ranking", (0.10, 0.52)),
    ]:
        for xpos, family in enumerate(FAMILY_ORDER):
            subset = metrics[metrics["family"] == family].sort_values("seed")
            xs = jitter_positions(xpos, len(subset))
            ax.scatter(
                xs,
                subset[column],
                s=22,
                marker=MARKERS[family],
                color=COLORS[family],
                edgecolor="white",
                linewidth=0.35,
                zorder=3,
            )
            if len(subset) > 1:
                ax.plot(
                    [xpos - 0.16, xpos + 0.16],
                    [subset[column].mean()] * 2,
                    color=COLORS[family],
                    lw=1.1,
                    zorder=2,
                )
        # ★ 修复：标签数量必须与 FAMILY_ORDER 一致；竖直排布避免 6 个标签拥挤
        ax.set_xticks(
            range(len(FAMILY_ORDER)),
            [TICK_LABELS[f] for f in FAMILY_ORDER],
            rotation=58,
            rotation_mode="anchor",
            ha="right",
            va="top",
            fontsize=5.5,
        )
        ax.tick_params(axis="x", pad=1.5)
        ax.set_ylim(*ylim)
        ax.set_title(title, loc="left", pad=7)
        ax.grid(axis="y", color="#E8E8E8", lw=0.55, zorder=0)
    axes[1].set_ylabel("Spearman $\\rho$")
    axes[2].set_ylabel("RMSE (kcal/mol)")
    axes[3].set_ylabel("Protein-level Pearson $r$")
    axes[3].text(
        0.02,
        0.02,
        "$n$ = 19 proteins",
        transform=axes[3].transAxes,
        fontsize=6.2,
        color="#555555",
        va="bottom",
    )
    add_panel_labels(fig, axes, "abcd")
    export_figure(fig, out_dir / "figure1_model_comparison")
    plt.close(fig)


def figure2(
    edits: pd.DataFrame,
    cb_edits: pd.DataFrame,
    ca_edits: pd.DataFrame,
    out_dir: Path,
) -> None:
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(7.2, 3.05),
        gridspec_kw={"width_ratios": [1.45, 1.2, 1]},
        constrained_layout=True,
    )

    ax = axes[0]
    for xpos, (column, label, color) in enumerate(
        [
            ("n_broken", "Broken", COLORS["broken"]),
            ("n_formed", "Formed", COLORS["formed"]),
        ]
    ):
        values = edits[column].to_numpy(dtype=float)
        parts = ax.violinplot(
            np.log1p(values),
            positions=[xpos],
            widths=0.72,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )
        for body in parts["bodies"]:
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(0.18)
        # A deterministic low-discrepancy sequence prevents stacked integer-valued points.
        jitter = ((np.arange(len(values)) * 0.61803398875) % 1.0 - 0.5) * 0.36
        ax.scatter(
            xpos + jitter,
            np.log1p(values),
            s=4,
            alpha=0.14,
            color=color,
            edgecolors="none",
            rasterized=True,
        )
        q1, median, q3 = np.percentile(np.log1p(values), [25, 50, 75])
        ax.plot([xpos, xpos], [q1, q3], color=color, lw=5, solid_capstyle="butt")
        ax.scatter([xpos], [median], s=16, color="white", edgecolor=color, lw=0.7, zorder=4)
    tick_values = np.array([0, 1, 3, 10, 30, 100, 300])
    ax.set_yticks(np.log1p(tick_values), [str(v) for v in tick_values])
    ax.set_xticks([0, 1], ["Broken", "Formed"])
    ax.set_ylabel("Contacts changed per mutation")
    ax.set_title("Centroid graphs expose contact rewiring", loc="left", pad=7)
    nonzero = ((edits["n_broken"] + edits["n_formed"]) > 0).mean() * 100
    ax.text(
        0.98,
        0.97,
        f"{nonzero:.1f}% with >=1 edit\n$n$ = {len(edits)} mutations",
        ha="right",
        va="top",
        transform=ax.transAxes,
        fontsize=6.3,
        color="#444444",
    )

    ax = axes[1]
    # 修正后改用 1R2Y R244E（一致氢排除下编辑数最多的案例之一）
    CASE_PID, CASE_MUT = "1R2Y", "R244E"
    case_ca = ca_edits[(ca_edits["pdb_id"] == CASE_PID) & (ca_edits["mut_info"] == CASE_MUT)]
    case_cb = cb_edits[(cb_edits["pdb_id"] == CASE_PID) & (cb_edits["mut_info"] == CASE_MUT)]
    case_sc = edits[(edits["pdb_id"] == CASE_PID) & (edits["mut_info"] == CASE_MUT)]
    if any(len(frame) != 1 for frame in (case_ca, case_cb, case_sc)):
        raise ValueError(f"Expected one {CASE_PID} {CASE_MUT} row in each contact representation")
    values = np.array(
        [
            [case_ca.iloc[0]["n_broken"], case_ca.iloc[0]["n_formed"]],
            [case_cb.iloc[0]["n_broken"], case_cb.iloc[0]["n_formed"]],
            [case_sc.iloc[0]["n_broken"], case_sc.iloc[0]["n_formed"]],
        ],
        dtype=float,
    )
    x = np.arange(values.shape[0])
    width = 0.32
    for idx, (label, color) in enumerate(
        [("Broken", COLORS["broken"]), ("Formed", COLORS["formed"])]
    ):
        bars = ax.bar(x + (idx - 0.5) * width, values[:, idx], width, label=label, color=color)
        for bar, value in zip(bars, values[:, idx]):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                value + 0.7,
                f"{int(value)}",
                ha="center",
                va="bottom",
                fontsize=6.3,
            )
    ax.set_xticks(x, [r"C$\alpha$", r"C$\beta$", "Side-chain\ncentroid"])
    ax.set_ylim(0, max(31, values.max() * 1.25))
    ax.set_ylabel("Contacts changed")
    ax.set_title(f"{CASE_PID} {CASE_MUT}: graph sensitivity", loc="left", pad=7)
    ax.legend(loc="upper left", ncol=1, handlelength=1.3)

    ax = axes[2]
    # 修正后：三案例 × 三运行的 FoldX 可复现性（一致氢排除）
    repro_path = DATA / "foldx_reproducibility_v2.csv"
    if repro_path.exists():
        rp = pd.read_csv(repro_path)
        agreement = pd.DataFrame({
            "case": rp["pdb_id"] + " " + rp["mut_info"],
            "broken_pct": rp["broken_agreement"] * 100,
            "formed_pct": rp["formed_agreement"] * 100,
        })
    else:
        raise FileNotFoundError(repro_path)
    ypos = np.arange(len(agreement))
    for j, (col, label) in enumerate([("broken_pct", "Broken"), ("formed_pct", "Formed")]):
        color = COLORS[label.lower()]
        ax.barh(ypos + (0.5 - j) * 0.34, agreement[col], height=0.32, color=color, label=label)
        for yp, v in zip(ypos + (0.5 - j) * 0.34, agreement[col]):
            ax.text(v + 1.5, yp, f"{v:.0f}%", va="center", fontsize=6.2, color=color)
    ax.set_yticks(ypos, agreement["case"], fontsize=6.5)
    ax.set_xlim(0, 118)
    ax.set_xlabel("Three-run agreement (%)")
    ax.set_title("FoldX repeatability", loc="left", pad=7)
    ax.legend(loc="lower right", fontsize=6, handlelength=1.2)
    ax.grid(axis="x", color="#E8E8E8", lw=0.55)

    add_panel_labels(fig, axes, "abc")
    export_figure(fig, out_dir / "figure2_contact_rewiring")
    agreement.to_csv(out_dir / "figure2_repeatability_source_data.csv", index=False)
    plt.close(fig)


def parse_foldx_rows() -> pd.DataFrame:
    benchmark = pd.read_csv(DATA / "benchmarks_s669_clean.csv")
    ddg_lookup = benchmark.set_index(["pdb_id", "mut_info"])["ddg"]
    rows = []
    for path in (DATA / "mutant_structures_s669").glob("*/*/work/out/Dif_*.fxout"):
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        header = None
        values = None
        for index, line in enumerate(lines[:-1]):
            if line.startswith("Pdb\t"):
                header = line.split("\t")
                values = lines[index + 1].split("\t")
                break
        if header is None or values is None:
            continue
        record = dict(zip(header, values))
        key = (path.parents[3].name, path.parents[2].name)
        if key not in ddg_lookup.index:
            continue
        try:
            predicted = float(record["total energy"])
        except (KeyError, TypeError, ValueError):
            continue
        rows.append(
            {
                "pdb_id": key[0],
                "mutation_id": key[1],
                "foldx_ddg": predicted,
                "experimental_ddg": float(ddg_lookup.loc[key]),
            }
        )
    return pd.DataFrame(rows).drop_duplicates(["pdb_id", "mutation_id"])


def figure3(predictions: pd.DataFrame, edits: pd.DataFrame, foldx: pd.DataFrame, out_dir: Path) -> None:
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(7.2, 3.05),
        gridspec_kw={"width_ratios": [0.95, 1.2, 1.35]},
        constrained_layout=True,
    )

    ax = axes[0]
    local = model_metrics(
        predictions[predictions["model"].str.startswith("gnn_local")], n_boot=250
    )
    plain = local[local["family"] == "local"].set_index("seed")
    blosum = local[local["family"] == "blosum"].set_index("seed")
    common_seeds = sorted(set(plain.index) & set(blosum.index))
    seed_markers = {42: "o", 123: "s", 2024: "^"}
    for seed in common_seeds:
        values = [plain.loc[seed, "pearson"], blosum.loc[seed, "pearson"]]
        ax.plot([0, 1], values, color="#AFAFAF", lw=0.8, zorder=1)
        ax.scatter(
            [0, 1],
            values,
            s=28,
            marker=seed_markers.get(seed, "o"),
            color=[COLORS["local"], COLORS["blosum"]],
            edgecolor="white",
            lw=0.45,
            zorder=3,
        )
    ax.set_xticks([0, 1], ["Local", "+ BLOSUM"])
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylim(0.33, 0.40)
    ax.set_ylabel("Pearson $r$")
    ax.set_title("BLOSUM adds no consistent gain", loc="left", pad=7)
    handles = [
        Line2D([], [], marker=seed_markers[s], color="#777777", lw=0, ms=4, label=f"seed {s}")
        for s in common_seeds
    ]
    ax.legend(handles=handles, loc="lower left", handletextpad=0.3)
    ax.grid(axis="y", color="#E8E8E8", lw=0.55)

    benchmark = pd.read_csv(DATA / "benchmarks_s669_clean.csv")[
        ["pdb_id", "mut_info", "ddg"]
    ]
    edit_signal = edits.merge(benchmark, on=["pdb_id", "mut_info"], how="inner")
    edit_signal["total_edits"] = edit_signal["n_broken"] + edit_signal["n_formed"]
    ax = axes[1]
    ax.scatter(
        edit_signal["total_edits"],
        edit_signal["ddg"],
        s=8,
        color=COLORS["accent"],
        alpha=0.34,
        edgecolors="none",
        rasterized=True,
    )
    if edit_signal["total_edits"].nunique() > 1:
        coefficients = np.polyfit(edit_signal["total_edits"], edit_signal["ddg"], 1)
        xline = np.linspace(0, edit_signal["total_edits"].quantile(0.99), 100)
        ax.plot(xline, np.polyval(coefficients, xline), color="#006D52", lw=1.2)
    contact_r = pearson(edit_signal["total_edits"], edit_signal["ddg"])
    contact_rho = float(spearmanr(edit_signal["total_edits"], edit_signal["ddg"]).statistic)
    ax.set_xlim(left=-5)
    ax.set_xlabel("Total contact edits")
    ax.set_ylabel("Experimental ΔΔG (kcal/mol)")
    ax.set_title("Edit counts are weakly associated", loc="left", pad=7)
    ax.text(
        0.97,
        0.96,
        f"Pearson r = {contact_r:.2f}\nSpearman rho = {contact_rho:.2f}\nn = {len(edit_signal)}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=6.2,
    )

    ax = axes[2]
    ax.scatter(
        foldx["experimental_ddg"],
        foldx["foldx_ddg"],
        s=8,
        color="#5E3C99",
        alpha=0.32,
        edgecolors="none",
        rasterized=True,
    )
    lower = min(foldx["experimental_ddg"].min(), foldx["foldx_ddg"].min())
    upper = max(foldx["experimental_ddg"].max(), foldx["foldx_ddg"].max())
    ax.plot([lower, upper], [lower, upper], color="#B5B5B5", lw=0.8, ls="--", label="Identity")
    coefficients = np.polyfit(foldx["experimental_ddg"], foldx["foldx_ddg"], 1)
    xline = np.linspace(foldx["experimental_ddg"].min(), foldx["experimental_ddg"].max(), 100)
    ax.plot(xline, np.polyval(coefficients, xline), color="#5E3C99", lw=1.2, label="Linear fit")
    foldx_r = pearson(foldx["experimental_ddg"], foldx["foldx_ddg"])
    foldx_mae = float(np.mean(np.abs(foldx["foldx_ddg"] - foldx["experimental_ddg"])))
    ax.set_xlabel("Experimental ΔΔG (kcal/mol)")
    ax.set_ylabel("FoldX ΔΔG (kcal/mol)")
    ax.set_title("FoldX energy signal remains modest", loc="left", pad=7)
    ax.text(
        0.04,
        0.96,
        f"Pearson r = {foldx_r:.2f}\nMAE = {foldx_mae:.2f} kcal/mol\nn = {len(foldx)}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.2,
    )
    ax.legend(loc="lower right", handlelength=1.4)

    add_panel_labels(fig, axes, "abc")
    export_figure(fig, out_dir / "figure3_limited_increment")
    edit_signal.to_csv(out_dir / "figure3_contact_signal_source_data.csv", index=False)
    foldx.to_csv(out_dir / "figure3_foldx_source_data.csv", index=False)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--bootstrap", type=int, default=1000)
    args = parser.parse_args()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    setup_style()
    predictions = pd.read_csv(DATA / "predictions_unified.csv")
    edits = pd.read_csv(DATA / "true_edits_s669_sc.csv")
    cb_edits = pd.read_csv(DATA / "true_edits_s669.csv")
    corrected_edits = pd.read_csv(DATA / "edits_corrected.csv")
    ca_edits = corrected_edits[
        (corrected_edits["atom_def"] == "ca") & (corrected_edits["threshold"] == 8.0)
    ]
    foldx = parse_foldx_rows()

    metrics = model_metrics(predictions, n_boot=args.bootstrap)
    metrics.to_csv(out_dir / "figure1_model_metrics_source_data.csv", index=False)
    edits[["pdb_id", "mut_info", "n_broken", "n_formed"]].to_csv(
        out_dir / "figure2_contact_edits_source_data.csv", index=False
    )

    figure1(metrics, out_dir)
    figure2(edits, cb_edits, ca_edits, out_dir)
    figure3(predictions, edits, foldx, out_dir)
    print(f"Generated publication figures in {out_dir}")
    print(f"Figure 1 models: {len(metrics)}; mutations per model: {metrics.n_mutations.min()}")
    print(f"Figure 2 centroid mutations: {len(edits)}; FoldX pairs: {len(foldx)}")


if __name__ == "__main__":
    main()
