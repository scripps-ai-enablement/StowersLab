#!/usr/bin/env python3
"""
vr_qc_figures.py -- deliverable QC figures for the vr_analysis pipeline.

    (a) marker_cpm_by_sample.png  -- marker CPM per sample, log scale, with the
        Rule 1 Trpc2 sort thresholds drawn as reference lines and target /
        nontarget / biology-suppressed samples visually distinguished.
    (b) technical_qc_panel.png    -- multi-loci vs too-many-loci rate per sample
        (the project's open question 1) and qualimap 5'-3' bias against the
        healthy ~1 reference.

Reads the joined table from vr_sample_qc.sample_qc_table so the figures can
never disagree with the gating table.

Usage
-----
    python vr_qc_figures.py --all-trials --fig-dir DIR
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vr_config import load_config, threshold, trials_of  # noqa: E402
from vr_sample_qc import sample_qc_table  # noqa: E402

__all__ = ["figure_marker_cpm", "figure_technical_qc"]

# Colour is bound to an entity once and reused across both figures (figure-style
# §4.1). Marker genes get one hue each; sample classes get their own encoding.
MARKER_COLORS = {
    "Trpc2": "#1b6ca8",    # pan-VNO-neuron, the Rule 1 gene -- focal
    "Omp": "#f2a03d",      # mature neuron
    "Gnai2": "#4a9c5d",    # V1R population
    "Gnao1": "#8c5bb0",    # V2R population
    "actin_sum": "#8a8f98",  # housekeeping
}
CLS_COLORS = {"target": "#1b6ca8", "nontarget": "#c2521a"}
ALARM = "#c0392b"      # reserved for suppressed/failed marks only
FLOOR = 1e-3           # plotting floor for zero CPM on a log axis


def _style(sizes=(9, 8, 7)):
    base, ann, tick = sizes
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": base, "axes.titlesize": base, "axes.labelsize": base,
        "legend.fontsize": ann, "xtick.labelsize": tick, "ytick.labelsize": tick,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.titlelocation": "left", "axes.titleweight": "normal",
        "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
        "legend.frameon": False, "figure.facecolor": "white",
    })


def check_overlaps(fig, label: str = "") -> List[str]:
    """
    Geometric self-check (figure-style §9.1): report visible text boxes that
    overlap each other or a spine they do not belong to. Printed as warnings so a
    layout regression surfaces in the run log rather than silently shipping.
    """
    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(matplotlib.text.Text)
             if t.get_text().strip() and t.get_visible()]
    spines = [(s, s.get_window_extent(r)) for ax in fig.axes
              for s in ax.spines.values() if s.get_visible()]
    ticks = {ax: set(ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both"))
             for ax in fig.axes}
    bad: List[str] = []
    for i, (ta, ba) in enumerate(texts):
        for tb, bb in texts[i + 1:]:
            if ba.overlaps(bb):
                bad.append(f"text/text {ta.get_text()!r} vs {tb.get_text()!r}")
    for t, bt in texts:
        for s, bs in spines:
            if bt.overlaps(bs) and t not in ticks.get(s.axes, ()):
                bad.append(f"text/spine {t.get_text()!r}")
    for b in bad:
        print(f"[vr_qc_figures] OVERLAP {label}: {b}")
    return bad


def _short(name: str) -> str:
    return (name.replace("100cells", "100c").replace("2cells", "2c")
                .replace("nontarget", "nonT").replace("target", "T")
                .replace("pool", "p"))


def _order(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values(["trial", "cell_type", "sample"],
                          ascending=[True, False, True]).reset_index(drop=True)


# ----------------------------------------------------------------------------
# (a) marker CPM by sample
# ----------------------------------------------------------------------------

def figure_marker_cpm(cfg, df: pd.DataFrame, path: str) -> str:
    _style()
    d = _order(df)
    n = len(d)
    x = np.arange(n)
    tmin = threshold(cfg, "target_trpc2_min")
    tcon = threshold(cfg, "target_trpc2_concern")
    nmax = threshold(cfg, "nontarget_trpc2_max")

    fig, ax = plt.subplots(figsize=(9.0, 4.6))

    # sample-class background bands: target vs nontarget, and suppressed samples
    for i, r in d.iterrows():
        if r["suppress_biology"]:
            ax.axvspan(i - 0.5, i + 0.5, color=ALARM, alpha=0.07, zorder=0, lw=0)
    # trial separator
    tb = np.flatnonzero(d["trial"].ne(d["trial"].shift()).values)[1:]
    for b in tb:
        ax.axvline(b - 0.5, color="0.35", lw=0.9, ls="-", zorder=1)

    # Rule 1 reference lines (never at a plotted data value)
    ax.axhline(tmin, color=MARKER_COLORS["Trpc2"], lw=1.0, ls="--", zorder=2)
    ax.axhline(nmax, color=CLS_COLORS["nontarget"], lw=1.0, ls=":", zorder=2)
    ax.text(n - 0.4, tmin * 1.35, f"target floor {tmin:g}", color=MARKER_COLORS["Trpc2"],
            fontsize=7, ha="right", va="bottom")
    ax.text(n - 0.4, nmax * 1.35, f"non-target ceiling {nmax:g}",
            color=CLS_COLORS["nontarget"], fontsize=7, ha="right", va="bottom")

    series = [("Trpc2", "Trpc2_cpm", "o", 46, 1.0),
              ("Omp", "Omp_cpm", "s", 26, 0.85),
              ("Gnai2", "Gnai2_cpm", "^", 26, 0.85),
              ("Gnao1", "Gnao1_cpm", "v", 26, 0.85),
              ("actin_sum", "actin_sum_cpm", "D", 22, 0.75)]
    for label, col, mk, size, alpha in series:
        vals = pd.to_numeric(d[col], errors="coerce").fillna(0.0).values
        y = np.where(vals <= 0, FLOOR, vals)
        zero = vals <= 0
        c = MARKER_COLORS[label]
        focal = label == "Trpc2"
        ax.scatter(x[~zero], y[~zero], marker=mk, s=size, color=c, alpha=alpha,
                   edgecolor="white" if focal else "none",
                   linewidth=0.6 if focal else 0, zorder=5 if focal else 4,
                   label="actin (sum)" if label == "actin_sum" else label)
        # zero values are marked at the floor as open glyphs, never dropped
        if zero.any():
            ax.scatter(x[zero], y[zero], marker=mk, s=size, facecolor="none",
                       edgecolor=c, linewidth=0.7, alpha=0.9, zorder=4)

    ax.set_yscale("log")
    ax.set_ylim(FLOOR / 2.5, 2e5)
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_yticks([1e-3, 1e-1, 1e1, 1e3, 1e5])
    ax.set_yticklabels(["0*", "0.1", "10", "1k", "100k"])
    ax.set_ylabel("expression (CPM)")
    # Claim-title checked against every plotted row (figure-style §1.4): in trial 2
    # three GFP+ libraries clear 1000 CPM and the GFP- library sits at 0.24. Trial 1
    # reads zero throughout because it is main olfactory epithelium, NOT because the
    # sort failed and NOT because of an annotation problem -- see the tissue panel.
    n_moe = int(pd.Series(d["tissue_verdict"]).isin(["MOE", "MOE_dominant_mixed"]).sum())
    ax.set_title("Trpc2 confirms the GFP+ sort in trial 2; trial 1 reads zero because "
                 f"{n_moe} of its libraries are main olfactory, not VNO")

    ax.set_xticks(x)
    labels = [_short(s) for s in d["sample"]]
    ax.set_xticklabels(labels, rotation=38, ha="right")
    for tick, (_, r) in zip(ax.get_xticklabels(), d.iterrows()):
        tick.set_color(CLS_COLORS.get(r["cell_type"], "black"))
        if r["suppress_biology"]:
            tick.set_fontstyle("italic")

    # trial band headers
    for trial, sub in d.groupby("trial"):
        idx = sub.index.values
        ax.text(idx.mean(), 1.15e5, trial, ha="center", va="top", fontsize=8,
                color="0.3")

    # Legend goes BELOW the axes: the zero-CPM glyphs occupy the floor row across
    # the full width, so any in-axes placement would sit on top of data.
    h, _ = ax.get_legend_handles_labels()
    extra = [plt.Line2D([], [], marker="o", ls="", mfc="none", mec="0.35",
                        label="0 CPM (at floor)"),
             plt.Line2D([], [], marker="s", ls="", color=ALARM, alpha=0.25,
                        ms=9, label="biology suppressed")]
    extra.append(plt.Line2D([], [], marker=r"$\bf{M}$", ls="", color=MOE_HUES[0],
                            ms=8, label="wrong tissue (main olfactory)"))
    ax.legend(handles=h + extra, loc="upper center", bbox_to_anchor=(0.5, -0.30),
              ncol=4, columnspacing=1.1, handletextpad=0.3, borderaxespad=0.0)
    # Mark the main-olfactory libraries directly, so this figure cannot be read
    # in isolation as "trial 1 failed its sort". Drawn just under the top spine
    # in axes coordinates so the glyph cannot land on a data point (pool100cells
    # carries a 39,881 CPM actin sum right at the top of the data range).
    for i, r in d.iterrows():
        if r["tissue_verdict"] in ("MOE", "MOE_dominant_mixed"):
            ax.annotate("M", (i, 0.965), xycoords=("data", "axes fraction"),
                        ha="center", va="top", fontsize=8, fontweight="bold",
                        color=MOE_HUES[0])
    ax.text(0.5, -0.44, "GFP+ (target) samples in blue, GFP− (non-target) in orange; "
            "italic = biology suppressed; M = main-olfactory tissue, so Trpc2 = 0 is "
            "expected there; 0* = zero CPM plotted at the axis floor.",
            transform=ax.transAxes, fontsize=7, color="0.35", ha="center")

    fig.canvas.draw()
    check_overlaps(fig, "marker_cpm")
    fig.savefig(path)
    plt.close(fig)
    return path


# ----------------------------------------------------------------------------
# (b) technical QC panel
# ----------------------------------------------------------------------------

def figure_technical_qc(cfg, df: pd.DataFrame, path: str) -> str:
    _style()
    d = _order(df)
    n = len(d)
    x = np.arange(n)
    labels = [_short(s) for s in d["sample"]]

    fig, axes = plt.subplots(2, 1, figsize=(9.0, 6.2), sharex=True,
                             gridspec_kw={"height_ratios": [1.15, 1.0], "hspace": 0.16})
    ax1, ax2 = axes
    C_MULTI, C_TOO = "#1b6ca8", "#b0862c"

    # --- top: the two multi-mapping channels ---
    multi = pd.to_numeric(d["pct_multi_loci"], errors="coerce").values
    too = pd.to_numeric(d["pct_too_many_loci"], errors="coerce").values
    for i in range(n):
        if np.isfinite(multi[i]) and np.isfinite(too[i]):
            ax1.plot([i, i], [max(too[i], 1e-3), multi[i]], color="0.75",
                     lw=0.9, zorder=2)
    ax1.scatter(x, np.where(multi <= 0, 1e-3, multi), s=46, color=C_MULTI,
                edgecolor="white", lw=0.6, zorder=4,
                label="mapped to multiple loci — retained, handed to Salmon EM")
    ax1.scatter(x, np.where(too <= 0, 1e-3, too), s=34, marker="v", color=C_TOO,
                edgecolor="white", lw=0.5, zorder=4,
                label="mapped to too many loci — discarded by STAR")
    ax1.set_yscale("log")
    ax1.set_ylim(8e-3, 200)
    ax1.set_yticks([0.01, 0.1, 1, 10, 100])
    ax1.set_yticklabels(["0.01", "0.1", "1", "10", "100"])
    ax1.set_ylabel("% of input reads")
    ax1.set_title("Multi-mapping loss is dominated by the retained channel, "
                  "so VR paralog ambiguity reaches the EM")
    ax1.legend(loc="upper left", ncol=1, borderaxespad=0.2, handletextpad=0.3)
    fin = np.isfinite(multi) & np.isfinite(too)
    if fin.any():
        ax1.text(0.995, 0.04, f"retained channel is {np.nanmedian(multi[fin]/np.maximum(too[fin],1e-9)):.0f}×"
                 " larger (median)", transform=ax1.transAxes, ha="right",
                 fontsize=7, color="0.3")

    # --- bottom: 5'-3' bias vs the healthy ~1 reference ---
    bias = pd.to_numeric(d["bias_5p3p"], errors="coerce").values
    lo = threshold(cfg, "bias_5p3p_low", section="qc_thresholds")
    hi = threshold(cfg, "bias_5p3p_high", section="qc_thresholds")
    ax2.axhspan(lo, hi, color="#4a9c5d", alpha=0.10, lw=0, zorder=1)
    ax2.axhline(1.0, color="#2f6b3d", lw=1.1, zorder=2)
    ok = (bias >= lo) & (bias <= hi)
    for i in range(n):
        if not np.isfinite(bias[i]):
            continue
        col = "#4a9c5d" if ok[i] else ALARM
        ax2.plot([i, i], [1.0, bias[i]], color=col, lw=0.9, alpha=0.6, zorder=3)
        ax2.scatter([i], [bias[i]], s=44, color=col, edgecolor="white", lw=0.6, zorder=4)
    miss = ~np.isfinite(bias)
    if miss.any():
        ax2.scatter(x[miss], np.full(miss.sum(), 1.0), s=44, marker="x",
                    color="0.45", zorder=4)
        for i in np.flatnonzero(miss):
            ax2.annotate("n.d.", (i, 1.0), textcoords="offset points",
                         xytext=(0, -13), ha="center", fontsize=7, color="0.45")
    ax2.set_yscale("log")
    ax2.set_ylim(0.4, 3000)
    ax2.set_yticks([0.5, 1, 10, 100, 1000])
    ax2.set_yticklabels(["0.5", "1", "10", "100", "1k"])
    ax2.set_ylabel("qualimap 5'–3' bias")
    # Title is derived from the data, not asserted, so it cannot drift out of
    # agreement with the plotted rows (figure-style §1.3/§1.4).
    n_out = int(np.sum(np.isfinite(bias) & ((bias < lo) | (bias > hi))))
    n_sev = int(np.sum(np.isfinite(bias) & (bias > 10 * hi)))
    ax2.set_title(f"Coverage evenness: {n_out} libraries fall outside the healthy "
                  f"band, {n_sev} of them by more than an order of magnitude")
    ax2.text(n - 0.45, 1.0 * 1.25, f"healthy band {lo:g}–{hi:g}", fontsize=7,
             color="#2f6b3d", ha="right", va="bottom")

    # headline values for the two extreme rows only (figure-style §2.5/§6.9)
    if np.isfinite(bias).any():
        for i in (int(np.nanargmax(bias)),):
            ax2.annotate(f"{bias[i]:.0f}×", (i, bias[i]), textcoords="offset points",
                         xytext=(0, 7), ha="center", fontsize=7, color=ALARM)

    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, rotation=38, ha="right")
    for tick, (_, r) in zip(ax2.get_xticklabels(), d.iterrows()):
        tick.set_color(CLS_COLORS.get(r["cell_type"], "black"))
        if r["suppress_biology"]:
            tick.set_fontstyle("italic")
    for ax in axes:
        ax.set_xlim(-0.6, n - 0.4)
        tb = np.flatnonzero(d["trial"].ne(d["trial"].shift()).values)[1:]
        for b in tb:
            ax.axvline(b - 0.5, color="0.35", lw=0.9, zorder=1)
    ax2.text(0.0, -0.46, "GFP+ (target) in blue, GFP− (non-target) in orange; "
             "italic = biology suppressed; n.d. = qualimap bias not reported.",
             transform=ax2.transAxes, fontsize=7, color="0.35")

    fig.canvas.draw()
    check_overlaps(fig, "technical_qc")
    fig.savefig(path)
    plt.close(fig)
    return path


# ----------------------------------------------------------------------------
# (c) tissue identity panel
# ----------------------------------------------------------------------------

# One hue family per tissue (figure-style §4.3): VNO greens/blues, MOE oranges.
VNO_SERIES = [("Trpc2", "Trpc2_cpm", "o"), ("Vmn1r (sum)", "Vmn1r_sum_cpm", "^"),
              ("Vmn2r (sum)", "Vmn2r_sum_cpm", "v")]
MOE_SERIES = [("Olfr (sum)", "Olfr_sum_cpm", "s"), ("Adcy3", "Adcy3_cpm", "P"),
              ("Cnga2", "Cnga2_cpm", "X"), ("Gnal", "Gnal_cpm", "D")]
VNO_HUES = ["#0b4f8a", "#2e7fbf", "#6fb3e0"]
MOE_HUES = ["#8a3d0b", "#c2651a", "#e08a3c", "#f0b678"]
SHARED_HUE = "#7a7f88"


def figure_tissue_identity(cfg, df: pd.DataFrame, path: str) -> str:
    """
    MOE vs VNO marker CPM per sample with the tissue floor drawn.

    This is the check that must be read BEFORE sort validation: Omp is plotted
    in neutral grey and labelled SHARED because it marks mature olfactory
    neurons of both tissues and cannot discriminate them.
    """
    _style()
    d = _order(df)
    n = len(d)
    x = np.arange(n)
    floor = threshold(cfg, "tissue_panel_floor_cpm", 100.0)

    fig, ax = plt.subplots(figsize=(9.4, 5.0))

    for i, r in d.iterrows():
        if r["tissue_verdict"] in ("MOE", "MOE_dominant_mixed"):
            ax.axvspan(i - 0.5, i + 0.5, color=MOE_HUES[1], alpha=0.09, zorder=0, lw=0)
    tb = np.flatnonzero(d["trial"].ne(d["trial"].shift()).values)[1:]
    for b in tb:
        ax.axvline(b - 0.5, color="0.35", lw=0.9, zorder=1)

    ax.axhline(floor, color="0.25", lw=1.1, ls="--", zorder=2)
    # label sits outside the data area, right of the last sample
    ax.text(1.005, floor, f"tissue floor\n{floor:g} CPM",
            transform=ax.get_yaxis_transform(), fontsize=7, color="0.25",
            ha="left", va="center")

    def draw(series, hues, jitter):
        for (label, col, mk), hue in zip(series, hues):
            vals = pd.to_numeric(d[col], errors="coerce").fillna(0.0).values
            y = np.where(vals <= 0, FLOOR, vals)
            zero = vals <= 0
            xp = x + jitter
            ax.scatter(xp[~zero], y[~zero], marker=mk, s=34, color=hue,
                       edgecolor="white", linewidth=0.4, zorder=5, label=label)
            if zero.any():
                ax.scatter(xp[zero], y[zero], marker=mk, s=30, facecolor="none",
                           edgecolor=hue, linewidth=0.7, zorder=4)

    draw(VNO_SERIES, VNO_HUES, -0.16)
    draw(MOE_SERIES, MOE_HUES, +0.16)
    omp = pd.to_numeric(d["Omp_cpm"], errors="coerce").fillna(0.0).values
    yo = np.where(omp <= 0, FLOOR, omp)
    ax.scatter(x, yo, marker="*", s=70, color=SHARED_HUE, edgecolor="white",
               linewidth=0.4, zorder=6, label="Omp — SHARED, not discriminating")

    ax.set_yscale("log")
    ax.set_ylim(FLOOR / 2.5, 3e5)
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_yticks([1e-3, 1e-1, 1e1, 1e3, 1e5])
    ax.set_yticklabels(["0*", "0.1", "10", "1k", "100k"])
    ax.set_ylabel("expression (CPM)")
    n_moe = int(d["tissue_verdict"].isin(["MOE", "MOE_dominant_mixed"]).sum())
    ax.set_title(f"Tissue identity: all {n_moe} flagged libraries are main olfactory "
                 "epithelium, not VNO")

    ax.set_xticks(x)
    ax.set_xticklabels([_short(s) for s in d["sample"]], rotation=38, ha="right")
    for tick, (_, r) in zip(ax.get_xticklabels(), d.iterrows()):
        tick.set_color(MOE_HUES[0] if r["tissue_verdict"] in ("MOE", "MOE_dominant_mixed")
                       else VNO_HUES[0] if r["tissue_verdict"].startswith("VNO")
                       else "0.45")
    for trial, sub in d.groupby("trial"):
        ax.text(sub.index.values.mean(), 2.2e5, trial, ha="center", va="top",
                fontsize=8, color="0.3")

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=4,
              columnspacing=1.1, handletextpad=0.3, borderaxespad=0.0)
    ax.text(0.5, -0.46, "VNO-specific markers in blue (left of each tick), main "
            "olfactory in orange (right); tissue label colours the sample name; "
            "0* = zero CPM at the axis floor.",
            transform=ax.transAxes, fontsize=7, color="0.35", ha="center")

    fig.canvas.draw()
    check_overlaps(fig, "tissue_identity")
    fig.savefig(path)
    plt.close(fig)
    return path


def _main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="QC deliverable figures")
    ap.add_argument("--config")
    ap.add_argument("--all-trials", action="store_true")
    ap.add_argument("--fig-dir", required=True)
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    frames = [sample_qc_table(cfg, t) for t in trials_of(cfg)]
    df = pd.concat(frames, ignore_index=True, sort=False)
    df.attrs.update(frames[0].attrs)
    os.makedirs(a.fig_dir, exist_ok=True)
    p1 = figure_marker_cpm(cfg, df, os.path.join(a.fig_dir, "marker_cpm_by_sample.png"))
    p2 = figure_technical_qc(cfg, df, os.path.join(a.fig_dir, "technical_qc_panel.png"))
    p3 = figure_tissue_identity(cfg, df,
                                os.path.join(a.fig_dir, "tissue_identity_panel.png"))
    for p in (p1, p2, p3):
        print(f"[vr_qc_figures] wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
