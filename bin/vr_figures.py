#!/usr/bin/env python3
"""
vr_figures.py -- deliverable figures for the tiered VNO receptor report.

Three figures, one function each. None of them duplicates the QC track's
marker/tissue/technical panels (results/figures/{marker_cpm_by_sample,
tissue_identity_panel,technical_qc_panel}.png), which already exist.

    fig_cluster_heatmap()   (a) samples x VR clusters. Samples suppressed by the
                                tier gate are drawn as hatched rows carrying the
                                failing tier, never silently dropped -- an absent
                                row would read as "no receptors found".
    fig_within_cluster_fractions()
                            (b) paralog fractions inside the dominant clusters,
                                with the even-split reference drawn in and read
                                support printed, so a 30-read even split is
                                distinguishable from a 3000-read one.
    fig_tier_overview()     (c) which tier each of the 10 libraries reached and
                                where it stopped.

Figure conventions follow the project's publication-grade rules; they are
implemented locally in `apply_figure_style()` rather than imported, so this
module runs on the cluster with only matplotlib available.

Usage
-----
    python vr_figures.py                    # all figures into results/figures/
    python vr_figures.py --fig tier
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from vr_config import load_config, threshold, trial_paths, trials_of  # noqa: E402
from vr_report import (  # noqa: E402
    NO_DATA,
    PASS,
    PASS_CAVEAT,
    SUPPRESSED,
    TIER_NAMES,
    TIERS,
    build_sample_report,
    load_qc_records,
    load_vr_tables,
)

__all__ = [
    "apply_figure_style",
    "COLORS",
    "fig_cluster_heatmap",
    "fig_within_cluster_fractions",
    "fig_tier_overview",
    "build_all",
]

# CVD-safe palette. One alarm hue (crimson) reserved for artifact/suppression
# marks and never reused as a data series colour; red/green opposition avoided.
COLORS: Dict[str, str] = {
    "focal": "#1f6fb2",        # reportable / VNO signal
    "secondary": "#e08214",    # comparator series
    "v1r": "#1f6fb2",
    "v2r": "#5e3c99",
    "alarm": "#c2185b",        # EM artifact / suppression
    "unique": "#0b7d6f",      # unique-read channel (independent evidence)
    "muted": "#9e9e9e",
    "suppressed_fill": "#e6e6e6",
    "grid": "#d9d9d9",
    "ink": "#222222",
}

# Role-mapped size ladder: base / annotation / tick.
SIZES: Tuple[int, int, int] = (9, 8, 7)


def apply_figure_style(sizes: Tuple[int, int, int] = SIZES) -> None:
    """Publication-grade rcParams: 3-size ladder, open frame, no chartjunk."""
    base, anno, tick = sizes
    mpl.rcParams.update({
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "font.size": base,
        "axes.titlesize": base,
        "axes.labelsize": base,
        "axes.titleweight": "regular",
        "axes.titlelocation": "left",
        "axes.edgecolor": COLORS["ink"],
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.labelcolor": COLORS["ink"],
        "text.color": COLORS["ink"],
        "legend.fontsize": anno,
        "legend.frameon": False,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "xtick.color": COLORS["ink"],
        "ytick.color": COLORS["ink"],
        "xtick.direction": "out",
        "ytick.direction": "out",
        "grid.color": COLORS["grid"],
        "grid.linewidth": 0.5,
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "hatch.linewidth": 0.6,
    })


def panel_letter(ax, letter: str, dx: float = -0.09, dy: float = 1.06) -> None:
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=SIZES[0] + 2,
            fontweight="bold", va="bottom", ha="left")


def _italic(name: str) -> str:
    """Italics via the text style, NOT mathtext: mathtext reinterprets the
    hyphen in a name like Vmn1r-ps5 as a minus sign."""
    return name


def _fmt_cpm(v: float) -> str:
    """Magnitude suffixes, never scientific notation."""
    if v >= 1000:
        return f"{f'{v / 1000:.1f}'.rstrip('0').rstrip('.')}k"
    if v >= 100:
        return f"{v:.0f}"
    if v >= 10:
        return f"{v:.1f}".rstrip("0").rstrip(".")
    if v >= 1:
        return f"{v:.2f}".rstrip("0").rstrip(".")
    return f"{v:.2f}"


def _verify_no_overlap(fig, label: str = "") -> List[str]:
    """§9.1 geometric check: report overlapping visible text boxes."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible()]
    findings = []
    for i, (a, ba) in enumerate(texts):
        for b, bb in texts[i + 1:]:
            if ba.overlaps(bb):
                findings.append(f"{label}: '{a.get_text()[:24]}' overlaps "
                                f"'{b.get_text()[:24]}'")
    return findings


# --------------------------------------------------------------------------
# Shared state assembly
# --------------------------------------------------------------------------


def gather(cfg, trials: Optional[Sequence[str]] = None) -> Dict[str, Any]:
    """Run the tier ladder for every sample and collect the VR tables."""
    trials = list(trials) if trials else trials_of(cfg)
    gates, recs, vr = [], [], {}
    for t in trials:
        rr = load_qc_records(cfg, trial=t)
        vt = load_vr_tables(cfg, t)
        vr[t] = vt
        recs += rr
        gates += [build_sample_report(cfg, r, vt) for r in rr]
    return {"gates": gates, "records": recs, "vr": vr, "trials": trials}


def _cluster_matrix(state: Dict[str, Any]) -> Tuple[pd.DataFrame, Dict[str, str]]:
    """
    samples x clusters matrix of cluster-level signal from the quantification
    track's vr_cluster_expression.tsv.

    The table carries BOTH the 200 kb `cluster` rows and the 800 kb
    `supercluster` rows under a `tier` column; only the 200 kb rows go into the
    matrix, or every region would appear twice. An empty matrix means the table
    has not landed yet.
    """
    frames = []
    prov: Dict[str, str] = {}
    for trial, vt in state["vr"].items():
        df = vt.get("cluster_expression")
        if df is None:
            prov[trial] = f"MISSING: {vt['_missing'].get('cluster_expression')}"
            continue
        prov[trial] = vt["_present"].get("cluster_expression", "")
        d = df.copy()
        if "tier" in d.columns:
            d = d[d["tier"].astype(str) == "cluster"]
        scol = next((c for c in ("sample", "sample_id", "Sample") if c in d.columns),
                    None)
        ccol = next((c for c in ("cluster_id", "cluster") if c in d.columns), None)
        vcol = next((c for c in ("cpm_sum", "cluster_cpm", "cluster_sum_cpm", "cpm",
                                 "value") if c in d.columns), None)
        if scol and ccol and vcol:
            d[vcol] = pd.to_numeric(d[vcol], errors="coerce")
            frames.append(d.pivot_table(index=scol, columns=ccol, values=vcol,
                                        aggfunc="sum"))
        elif ccol:  # wide: clusters x samples
            w = d.drop_duplicates(subset=[ccol]).set_index(ccol)
            frames.append(w.apply(pd.to_numeric, errors="coerce").T)
    if not frames:
        return pd.DataFrame(), prov
    out = pd.concat(frames).fillna(0.0)
    # one row per sample even if a sample appears in more than one input frame
    return out.groupby(level=0).sum(), prov


def _sample_labels(gates) -> List[Tuple[str, str, int, str]]:
    """(sample, trial, highest_tier, status_word) in report order."""
    out = []
    for g in gates:
        top = g.highest_tier_reported()
        out.append((g.sample, g.trial, top,
                    "reportable" if top >= 3 else "suppressed"))
    return out


# --------------------------------------------------------------------------
# (a) cluster-level heatmap
# --------------------------------------------------------------------------


def fig_cluster_heatmap(state: Dict[str, Any], out_path: str,
                        max_clusters: int = 24) -> Dict[str, Any]:
    """
    Samples x VR clusters, cluster-level signal on a log colour scale.

    Suppressed samples keep their row. The row is filled flat grey with diagonal
    hatching and carries, inside the row itself, the tier that stopped it.
    Dropping the row would read as "this library had no receptor signal" when
    the truth is "this library was never allowed to make a receptor statement".
    """
    apply_figure_style()
    mat, prov = _cluster_matrix(state)
    gates = state["gates"]
    labels = _sample_labels(gates)
    samples = [s for s, _, _, _ in labels]
    gate_by = {g.sample: g for g in gates}
    n_s = len(samples)

    def _stop_label(g) -> str:
        top = g.highest_tier_reported()
        stop = top + 1
        return (f"suppressed at tier {stop} ({TIER_NAMES[stop]})"
                if stop in g.outcomes else "suppressed")

    pending = mat.empty
    if not pending:
        mat = mat.reindex(samples)
        # Column selection is driven by the ELIGIBLE rows only. A cluster whose
        # only signal sits in a gated-off library must not occupy a column: the
        # figure would show an empty column and imply the eligible libraries
        # were tested for it and came back negative.
        eligible = [s for s, _, top, _ in labels if top >= 3]
        basis = mat.loc[[s for s in eligible if s in mat.index]] if eligible else mat
        keep = basis.fillna(0).sum(axis=0).sort_values(ascending=False)
        cols = list(keep[keep > 0].index[:max_clusters])
        if not cols:
            pending = True
        else:
            mat = mat[cols]

    n_reportable = sum(1 for g in gates if g.highest_tier_reported() >= 3)

    if pending:
        fig, ax = plt.subplots(figsize=(7.0, 0.36 * n_s + 1.6))
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.6, n_s - 0.4)
        ax.set_xticks([])
        for i, (s, trial, top, _) in enumerate(labels):
            y = n_s - 1 - i
            suppressed = top < 3
            ax.add_patch(Rectangle((0.02, y - 0.40), 0.96, 0.80,
                                   facecolor=COLORS["suppressed_fill"]
                                   if suppressed else "white",
                                   edgecolor=COLORS["muted"], linewidth=0.5,
                                   hatch="////" if suppressed else None))
            msg = (_stop_label(gate_by[s]) if suppressed
                   else "eligible — awaiting vr_cluster_expression.tsv")
            ax.text(0.06, y, msg, ha="left", va="center", fontsize=SIZES[2],
                    color=COLORS["alarm"] if suppressed else COLORS["secondary"],
                    style="normal" if suppressed else "italic")
        ax.set_yticks(range(n_s))
        ax.set_yticklabels([f"{s}  ({t})" for s, t, _, _ in labels][::-1])
        for sp in ("top", "right", "bottom", "left"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(length=0)
        ax.set_title(f"Cluster-level VR calls: {n_s - n_reportable} of {n_s} libraries "
                     "are gated off before any\nreceptor statement; cluster "
                     "quantification is still pending")
        ax.set_xlabel("VR clusters (200 kb definition; superclusters at 800 kb)")
        fig.savefig(out_path)
        findings = _verify_no_overlap(fig, "cluster_heatmap")
        plt.close(fig)
        return {"path": out_path, "pending": True, "provenance": prov,
                "overlaps": findings, "n_reportable": n_reportable}

    ncl = len(mat.columns)
    fig, ax = plt.subplots(figsize=(max(6.4, 0.42 * ncl + 4.6),
                                    0.40 * n_s + 2.4))
    data = mat.to_numpy(dtype=float)
    pos = data[data > 0]
    norm = LogNorm(vmin=max(pos.min(), 1e-2), vmax=pos.max()) if pos.size else None
    masked = np.ma.masked_less_equal(data, 0)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(color="#f7f7f7")
    im = ax.imshow(masked, aspect="auto", cmap=cmap, norm=norm)
    ax.set_xticks(range(ncl))
    ax.set_xticklabels(mat.columns, rotation=90)
    ax.set_yticks(range(n_s))
    ax.set_yticklabels([f"{s}  ({t})" for s, t, _, _ in labels])
    ax.tick_params(length=2)

    # value labels: viridis runs dark (low) -> bright (high), so light text goes
    # on the LOW end. Getting this backwards is what makes a small value on a
    # dark cell unreadable.
    if data.size <= 400 and norm is not None:
        for i in range(data.shape[0]):
            if labels[i][2] < 3:
                continue
            for j in range(data.shape[1]):
                v = data[i, j]
                if v > 0:
                    # contrast from the cell's actual luminance, not a guess at
                    # where the colormap turns light
                    rgb = cmap(norm(v))[:3]
                    lum = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
                    ax.text(j, i, _fmt_cpm(v), ha="center", va="center",
                            fontsize=SIZES[2] - 1,
                            color="white" if lum < 0.55 else COLORS["ink"],
                            zorder=5)

    # suppressed rows: hatch band plus the stopping tier INSIDE the row, so the
    # annotation cannot collide with the colorbar.
    for i, (s, trial, top, _) in enumerate(labels):
        if top >= 3:
            continue
        ax.add_patch(Rectangle((-0.5, i - 0.5), ncl, 1.0,
                               facecolor=COLORS["suppressed_fill"],
                               edgecolor=COLORS["muted"], linewidth=0.0,
                               hatch="////", zorder=3))
        ax.text(-0.4, i, _stop_label(gate_by[s]), va="center", ha="left",
                fontsize=SIZES[2], color=COLORS["alarm"], zorder=6,
                bbox=dict(facecolor="white", edgecolor="none", pad=1.6))

    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.030)
    cb.set_label("cluster-level VR signal (CPM, log scale)", fontsize=SIZES[1])
    cb.ax.tick_params(labelsize=SIZES[2])
    fig.suptitle("Cluster-level VR calls; suppressed libraries are marked in "
                 "place, not dropped", x=0.012, ha="left", fontsize=SIZES[0],
                 y=0.998)
    ax.set_xlabel("VR cluster (200 kb definition; superclusters at 800 kb)")
    ax.legend(handles=[
        Patch(facecolor=COLORS["suppressed_fill"], hatch="////",
              edgecolor=COLORS["muted"],
              label="gated off — no receptor statement permitted"),
        Patch(facecolor="#f7f7f7", edgecolor=COLORS["muted"],
              label="no signal in an eligible library"),
    ], loc="lower left", bbox_to_anchor=(0.0, 1.008), ncol=2, fontsize=SIZES[1],
        handlelength=1.6, handleheight=0.9)
    fig.subplots_adjust(top=0.90)
    fig.savefig(out_path)
    findings = _verify_no_overlap(fig, "cluster_heatmap")
    plt.close(fig)
    return {"path": out_path, "pending": False, "provenance": prov,
            "overlaps": findings, "n_clusters": ncl,
            "n_reportable": n_reportable}


# --------------------------------------------------------------------------
# (b) within-cluster paralog fractions
# --------------------------------------------------------------------------

# em_flag values that mean "nothing to see": no signal, or too little signal to
# say anything. These are not findings and are not annotated as such.
_EM_QUIET = {"no_signal", "insufficient_signal", "none", "nan", ""}

# Presentation for each em_flag the quantification track can emit. The `alarm`
# key decides colour: ONLY a redistribution call is alarm-coloured, because only
# that one says the data are misleading. A structured split is a real result.
_EM_STYLE: Dict[str, Dict[str, Any]] = {
    "suspected_em_redistribution": dict(
        alarm=True,
        short="EM redistribution suspected",
        gloss=("fractions are indistinguishable from an even split AND only one "
               "member has independent unique-read support"),
    ),
    "no_redistribution_signature": dict(
        alarm=False,
        short="no redistribution signature",
        gloss=("fractions are statistically distinguishable from an even split; "
               "consistent with a real per-paralog call"),
    ),
    "single_paralog_only": dict(
        alarm=False,
        short="single paralog only",
        gloss="one member carries the cluster's reads; no split to test",
    ),
}


def _em_flag_index(state: Dict[str, Any]) -> Tuple[Dict[Tuple[str, str], Dict[str, Any]],
                                                   Dict[str, str], Optional[float]]:
    """
    (sample, cluster_id) -> the quantification track's own artifact-flag row.

    This figure NEVER decides whether a split is EM redistribution. That
    judgement is a statistic (a Monte-Carlo test against uniformity, plus a
    unique-read gate) computed by bin/vr_quantify.py and recorded in
    vr_artifact_flags.tsv. Plotting code that re-derived it by eye would be free
    to disagree with the pipeline, and two artifacts in one results directory
    contradicting each other is worse than either being wrong alone. So we read
    `em_flag` and its supporting columns and annotate exactly what the pipeline
    concluded.

    Also returns the unique-read support threshold the track used, so the
    co-dominance annotation below can describe unique support using the
    pipeline's own cutoff rather than inventing one.
    """
    idx: Dict[Tuple[str, str], Dict[str, Any]] = {}
    prov: Dict[str, str] = {}
    thresholds: List[float] = []
    for trial, vt in state["vr"].items():
        df = vt.get("artifact_flags")
        if df is None:
            prov[trial] = f"MISSING: {vt['_missing'].get('artifact_flags')}"
            continue
        prov[trial] = vt["_present"].get("artifact_flags", "")
        d = df.copy()
        if "scope" in d.columns:
            d = d[d["scope"].astype(str) == "cluster"]
        if "unique_support_threshold" in d.columns:
            thresholds += [float(v) for v in
                           pd.to_numeric(d["unique_support_threshold"],
                                         errors="coerce").dropna().unique()]
        for _, r in d.iterrows():
            idx[(str(r.get("sample")), str(r.get("cluster_id")))] = r.to_dict()
    thr = max(thresholds) if thresholds else None
    return idx, prov, thr


def fig_within_cluster_fractions(state: Dict[str, Any], out_path: str,
                                 max_panels: int = 6,
                                 min_paralogs: int = 2) -> Dict[str, Any]:
    """
    Within-cluster paralog structure for the dominant clusters, one panel per
    (sample, cluster), on TWO channels per paralog:

        circle    EM fraction of cluster reads (Salmon's redistributed estimate)
        diamond   share of the cluster's UNIQUE reads (MAPQ-filtered, dedup) --
                  independent evidence that does not pass through the EM step

    The divergence between the channels is the whole point. An even EM split
    whose unique-read channel collapses onto ONE member is redistribution: one
    transcript's reads were divided among sequence-similar paralogs. An even EM
    split that the unique channel REPRODUCES is co-expression: two neurons in a
    multi-cell pool each expressing a different paralog of the same cluster.
    Monogenic choice is a per-CELL rule, so in a 100-cell pool the second case
    is expected biology, not a defect.

    The verdict text is read from vr_artifact_flags.tsv (`em_flag` and its
    supporting columns). This function does not compute a redistribution
    judgement of its own.
    """
    apply_figure_style()
    frames = []
    prov: Dict[str, str] = {}
    for trial, vt in state["vr"].items():
        df = vt.get("within_cluster_fractions")
        if df is None:
            prov[trial] = f"MISSING: {vt['_missing'].get('within_cluster_fractions')}"
            continue
        prov[trial] = vt["_present"].get("within_cluster_fractions", "")
        frames.append(df)
    reportable = {g.sample for g in state["gates"] if g.highest_tier_reported() >= 3}
    flag_idx, flag_prov, uniq_thr = _em_flag_index(state)
    prov.update({f"artifact_flags::{k}": v for k, v in flag_prov.items()})

    if not frames:
        fig, ax = plt.subplots(figsize=(7.2, 3.0))
        ax.axis("off")
        ax.text(0.0, 1.0, "Within-cluster paralog structure — awaiting "
                          "vr_within_cluster_fractions.tsv",
                fontsize=SIZES[0], va="top")
        ax.text(0.0, 0.80,
                "Each panel will show, for one dominant cluster, two channels per "
                "paralog: the EM\nfraction of cluster reads and the share of the "
                "cluster's UNIQUE (MAPQ-filtered)\nreads. Where the channels agree on "
                "an even split, two cells expressed two paralogs\n(co-expression). "
                "Where the EM is even but unique support collapses onto one member,\n"
                "the EM divided one transcript's reads. The verdict shown is the one "
                "recorded in\nvr_artifact_flags.tsv; this figure does not judge "
                "redistribution itself.",
                fontsize=SIZES[1], va="top", linespacing=1.6)
        ax.text(0.0, 0.10,
                "Libraries currently eligible: "
                + (", ".join(sorted(reportable)) if reportable
                   else "none — every library is gated off upstream of tier 3."),
                fontsize=SIZES[1], va="top", color=COLORS["secondary"])
        fig.savefig(out_path)
        findings = _verify_no_overlap(fig, "within_cluster")
        plt.close(fig)
        return {"path": out_path, "pending": True, "provenance": prov,
                "overlaps": findings}

    df = pd.concat(frames, ignore_index=True)
    scol = next((c for c in ("sample", "sample_id") if c in df.columns), None)
    ccol = next((c for c in ("cluster_id", "cluster") if c in df.columns), None)
    gcol = next((c for c in ("gene_name", "gene") if c in df.columns), None)
    fcol = next((c for c in ("frac_of_cluster", "within_cluster_fraction",
                             "fraction", "frac") if c in df.columns), None)
    rcol = next((c for c in ("counts", "reads", "n_reads") if c in df.columns), None)
    # the unique channel, in order of preference: deduplicated BAM uniques are
    # the strongest evidence the track produces
    ucol = next((c for c in ("unique_reads_bam_nodup", "unique_reads_bam",
                             "unique_reads_salmon_tx") if c in df.columns), None)
    pcol = next((c for c in ("is_pseudogene", "pseudogene") if c in df.columns), None)
    if scol:
        df = df[df[scol].isin(reportable)]
    df = df.copy()
    for c in (fcol, rcol, ucol):
        if c:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    groups = []
    for (s, cl), sub in df.groupby([scol, ccol]):
        # A paralog with no EM fraction and no unique reads is not part of this
        # cluster's story; dropping it keeps the panel readable. A paralog with
        # unique reads but no EM fraction is KEPT -- that asymmetry is evidence.
        keep = (sub[fcol].fillna(0) > 0)
        if ucol:
            keep = keep | (sub[ucol].fillna(0) >= (uniq_thr or 10.0))
        sub = sub[keep]
        if len(sub) < min_paralogs:
            continue
        row = flag_idx.get((str(s), str(cl)), {})
        flag = str(row.get("em_flag", "")) if row else ""
        tot = float(sub[rcol].sum()) if rcol else float("nan")
        # Ranking: a cluster the pipeline actually flagged must appear, however
        # few reads it carries. Ordering by read support alone hid the single
        # genuine redistribution call (8,974 reads) behind six larger clusters.
        alarm = bool(_EM_STYLE.get(flag, {}).get("alarm"))
        codom = 0
        tf, sf = row.get("top_frac"), row.get("second_frac")
        try:
            if float(tf) >= 0.2 and float(sf) >= 0.2:
                codom = 1
        except (TypeError, ValueError):
            pass
        groups.append({"key": (s, cl), "n_cells": row.get("n_cells"),
                       "sub": sub.sort_values(fcol, ascending=False),
                       "reads": tot, "flag": flag, "row": row,
                       "rank": (0 if alarm else 1 if codom else 2),
                       "alarm": alarm, "codominant": bool(codom)})
    groups.sort(key=lambda g: (g["rank"], -(g["reads"] if g["reads"] == g["reads"]
                                            else 0)))
    groups = groups[:max_panels]

    if not groups:
        fig, ax = plt.subplots(figsize=(7.0, 2.0))
        ax.axis("off")
        ax.text(0.0, 0.85, "Within-cluster paralog structure", fontsize=SIZES[0],
                va="top")
        ax.text(0.0, 0.45, "No cluster in a tier-3-eligible library carries "
                           f">= {min_paralogs} paralogs with signal.",
                fontsize=SIZES[1], va="top")
        fig.savefig(out_path)
        findings = _verify_no_overlap(fig, "within_cluster")
        plt.close(fig)
        return {"path": out_path, "pending": False, "provenance": prov,
                "overlaps": findings, "n_panels": 0,
                "n_redistribution_flagged": 0, "n_coexpression": 0}

    # ONE column. The verdict sentence under each panel is the figure's payload,
    # not a caption, and at two columns it either wraps into its neighbour or has
    # to be truncated into uselessness.
    ncol = 1
    nrow = len(groups)
    # Per-panel height. A 4-paralog cluster padded to a 9-paralog cluster's
    # height is dead space that makes the figure unreadable at page scale.
    heights = [0.245 * len(g["sub"]) + 0.72 for g in groups]
    fig = plt.figure(figsize=(9.6, sum(heights) + 1.05))
    gs = fig.add_gridspec(nrow, ncol, height_ratios=heights, hspace=0.62)
    axes = [[fig.add_subplot(gs[r, 0])] for r in range(nrow)]
    n_redist = 0
    n_coexpr = 0
    n_pseudo_codom = 0
    for gi, G in enumerate(groups):
        ax = axes[gi // ncol][gi % ncol]
        (s, cl), sub, row, flag = G["key"], G["sub"], G["row"], G["flag"]
        n = len(sub)
        y = np.arange(n)[::-1]
        vals = sub[fcol].fillna(0).to_numpy(dtype=float)
        reads = (sub[rcol].to_numpy(dtype=float) if rcol
                 else np.full(n, float("nan")))
        uq = (sub[ucol].fillna(0).to_numpy(dtype=float) if ucol
              else np.full(n, float("nan")))
        uq_tot = float(np.nansum(uq)) if ucol else 0.0
        uq_share = (uq / uq_tot if uq_tot > 0 else np.full(n, float("nan")))
        pseudo = ([str(v).strip() in ("1", "True", "true")
                   for v in sub[pcol]] if pcol else [False] * n)

        # The even-split reference is drawn ONLY where the pipeline itself
        # identified an even block (even_block_size >= 2). No block, no band --
        # that is what stops this figure from asserting a split the statistic
        # does not support.
        blk = pd.to_numeric(pd.Series([row.get("even_block_size")]),
                            errors="coerce").iloc[0]
        if blk == blk and blk >= 2:
            even = 1.0 / float(blk)
            band = 0.25 * even
            ax.axvspan(max(even - band, 0.0), even + band, color=COLORS["alarm"],
                       alpha=0.10, zorder=0)
            ax.axvline(even, color=COLORS["alarm"], lw=0.9, ls="--", zorder=1)
            ax.text(even, n - 0.28, f"1/{int(blk)} even split (pipeline)",
                    ha="center", va="top", fontsize=SIZES[2], color=COLORS["alarm"])

        for yi, v, rd, un, us, nm, ps in zip(y, vals, reads, uq, uq_share,
                                             sub[gcol].astype(str), pseudo):
            # connector spans the two channels: a long connector IS the
            # disagreement between EM and unique evidence
            if us == us:
                ax.hlines(yi, min(v, us), max(v, us), color=COLORS["muted"],
                          lw=1.4, zorder=2, alpha=0.55)
            ax.hlines(yi, 0, v, color=COLORS["grid"], lw=0.7, zorder=1)
            ms = 3.2 + 3.0 * np.log10(max(rd, 1.0)) if rd == rd else 4.0
            ax.plot([v], [yi], "o", ms=ms,
                    color=COLORS["alarm"] if ps else COLORS["focal"],
                    markeredgecolor="white", markeredgewidth=0.4, zorder=4)
            if us == us:
                ax.plot([us], [yi], "D", ms=4.6, color=COLORS["unique"],
                        markeredgecolor="white", markeredgewidth=0.4, zorder=5)
            lbl = f"{v:.2f}"
            if un == un:
                lbl += f"   {un:.0f}u"
            ax.text(1.03, yi, lbl, va="center", ha="left", fontsize=SIZES[2])
        ax.set_yticks(y)
        ax.set_yticklabels(list(sub[gcol].astype(str)), style="italic")
        ax.set_xlim(0, 1.40)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
        ax.set_ylim(-0.55, n - 0.12)

        # verdict line: the pipeline's own flag, never a judgement made here
        st = _EM_STYLE.get(flag)
        if flag in _EM_QUIET or flag == "":
            vtxt, vcol = (f"pipeline: {flag or 'no flag row'} — not annotated",
                          COLORS["muted"])
        else:
            vtxt = f"pipeline: {st['short'] if st else flag}"
            vcol = COLORS["alarm"] if (st and st["alarm"]) else COLORS["focal"]
        if st and st["alarm"]:
            n_redist += 1
            nsup = pd.to_numeric(
                pd.Series([row.get("n_block_members_unique_supported")]),
                errors="coerce").iloc[0]
            bu = str(row.get("block_unique_reads", ""))
            if nsup == nsup:
                vtxt += (f" — {int(nsup)} of {int(blk)} block members clear the "
                         f"unique-read threshold"
                         + (f" ({bu} unique)" if bu and bu != "nan" else ""))
        elif G["codominant"] and ucol and uq_tot > 0:
            # DESCRIPTIVE, not a verdict: state what the unique channel shows for
            # the two co-dominant members. The threshold is the pipeline's own.
            thr = uniq_thr if uniq_thr is not None else 10.0
            top2 = np.argsort(-vals)[:2]
            sup = int(sum(1 for k in top2 if uq[k] >= thr))
            names = [str(sub[gcol].astype(str).to_numpy()[k]) for k in top2]
            ps2 = [bool(pseudo[k]) for k in top2]
            if sup == 2 and not any(ps2):
                n_coexpr += 1
                ncells = str(G.get("n_cells") or "multi")
                vtxt += (f"; {names[0]} and {names[1]} are BOTH independently "
                         f"observed ({uq[top2[0]]:.0f} and {uq[top2[1]]:.0f} unique "
                         f"reads ≥ {thr:.0f}) — real co-expression across the "
                         f"{ncells}-cell pool, not redistribution")
            elif sup == 2 and any(ps2):
                # A PSEUDOGENE cannot be the receptor a cell chose. Two members
                # with unique support, one of them a pseudogene, is a
                # pseudogene-bleed / annotation question, not co-expression.
                n_pseudo_codom += 1
                which = names[ps2.index(True)]
                vtxt += (f"; both co-dominant members clear the unique-read "
                         f"threshold, but {which} is a PSEUDOGENE — this is not "
                         f"two-cell co-expression. Mechanism is unresolved between "
                         f"pseudogene transcription and mis-annotation")
            else:
                vtxt += (f"; only {sup} of the 2 co-dominant members clears the "
                         f"unique-read threshold ({thr:.0f})")
        # verdict below the axes in figure-relative terms, so it cannot collide
        # with the frame or the tick labels however tall the panel is
        dy = -44 if gi == len(groups) - 1 else -26
        ax.annotate(vtxt, xy=(0.0, 0.0), xycoords="axes fraction",
                    xytext=(0, dy), textcoords="offset points",
                    ha="left", va="top", fontsize=SIZES[2], color=vcol,
                    annotation_clip=False)

        title = f"{s}  ·  {cl}"
        if G["reads"] == G["reads"]:
            title += f"   ({G['reads']:.0f} EM reads"
            if uq_tot > 0:
                title += f", {uq_tot:.0f} unique"
            title += ")"
        pu = pd.to_numeric(pd.Series([row.get("p_uniform")]),
                           errors="coerce").iloc[0]
        if pu == pu:
            title += f"   p_uniform = {pu:.1e}"
        ax.set_title(title, fontsize=SIZES[1])
    for k in range(len(groups), nrow * ncol):
        axes[k // ncol][k % ncol].axis("off")
    for cix in range(ncol):
        used = [r for r in range(nrow) if r * ncol + cix < len(groups)]
        if used:
            axes[max(used)][cix].set_xlabel(
                "share of cluster reads  (circle = EM, diamond = unique)",
                labelpad=3)
    fig.suptitle("Within-cluster paralog structure: an even EM split is "
                 "redistribution only when the unique-read channel fails to "
                 "reproduce it",
                 x=0.010, ha="left", fontsize=SIZES[0], y=0.997)
    fig.legend(handles=[
        Line2D([], [], marker="o", ls="", color=COLORS["focal"], ms=5,
               label="EM fraction (marker size ∝ log EM reads)"),
        Line2D([], [], marker="D", ls="", color=COLORS["unique"], ms=5,
               label="share of cluster unique reads (MAPQ-filtered, dedup)"),
        Line2D([], [], marker="o", ls="", color=COLORS["alarm"], ms=5,
               label="pseudogene"),
        Line2D([], [], color=COLORS["alarm"], ls="--",
               label="even-split reference, drawn only where the pipeline found a "
                     "block"),
    ], loc="lower left", bbox_to_anchor=(0.010, 0.008), ncol=2, fontsize=SIZES[1],
        handletextpad=0.5, columnspacing=1.6)
    H = sum(heights) + 1.90
    fig.set_size_inches(9.6, H)
    fig.subplots_adjust(left=0.115, right=0.995,
                        top=1.0 - 0.55 / H, bottom=1.20 / H)
    fig.savefig(out_path)
    findings = _verify_no_overlap(fig, "within_cluster")
    plt.close(fig)
    return {"path": out_path, "pending": False, "provenance": prov,
            "n_panels": len(groups), "n_redistribution_flagged": n_redist,
            "n_coexpression": n_coexpr,
            "n_pseudogene_codominant": n_pseudo_codom, "unique_threshold": uniq_thr,
            "overlaps": findings}


# --------------------------------------------------------------------------
# (c) tier-status overview
# --------------------------------------------------------------------------

# Distinct glyph per status: two entries that look identical would make the key
# unresolvable, so the caveat state gets an open-centred marker of the same hue.
_STATUS_STYLE = {
    PASS: dict(marker="o", color=COLORS["focal"], fill=COLORS["focal"],
               label="reported"),
    PASS_CAVEAT: dict(marker="o", color=COLORS["focal"], fill="#a9cbe8",
                      label="reported, with a stated caveat"),
    "FAIL": dict(marker="X", color=COLORS["alarm"], fill=COLORS["alarm"],
                 label="failed at this tier"),
    SUPPRESSED: dict(marker="o", color=COLORS["muted"], fill="white",
                     label="suppressed by an upstream tier"),
    # NO_DATA means the tier RAN and produced no statement -- e.g. the
    # quantification track nominated no candidate because no cluster cleared
    # the signal threshold. It does NOT mean a table is missing or a step is
    # outstanding, and labelling it "awaiting"/"pending" misreads a completed
    # negative result as an unfinished one.
    NO_DATA: dict(marker="s", color=COLORS["secondary"], fill="white",
                  label="tier ran, nothing to report"),
}

_SHORT_STOP = {
    "tissue_identity": "wrong tissue (MOE)",
    "sort_validation": "sort validation failed",
    "population_id": "no Gnai2/Gnao1 signal",
    "cluster_vr": "no cluster above threshold",
    "individual_vr": "no candidate nominated",
}


def fig_tier_overview(state: Dict[str, Any], out_path: str) -> Dict[str, Any]:
    """
    Which tier each library reached and where it stopped: a ladder over the five
    tiers, with the stopping reason in its own right-hand column so the figure
    answers "why is there no receptor call for this sample?" on its own.
    """
    apply_figure_style()
    gates = state["gates"]
    n = len(gates)
    n_tier1 = sum(1 for g in gates if g.highest_tier_reported() >= 1)
    n_tier3 = sum(1 for g in gates if g.highest_tier_reported() >= 3)

    fig, (ax, axr) = plt.subplots(
        1, 2, figsize=(8.0, 0.42 * n + 2.9),
        gridspec_kw=dict(width_ratios=[1.0, 0.95], wspace=0.02))
    for i, g in enumerate(gates):
        y = n - 1 - i
        ax.plot([0, 4], [y, y], color=COLORS["grid"], lw=0.6, zorder=0)
        for t in TIERS:
            oc = g.outcomes.get(t)
            if oc is None:
                continue
            st = _STATUS_STYLE.get(oc.status, _STATUS_STYLE["FAIL"])
            ax.plot([t], [y], marker=st["marker"], ms=6.0,
                    markerfacecolor=st["fill"], markeredgecolor=st["color"],
                    markeredgewidth=1.0, ls="", zorder=3)
        top = g.highest_tier_reported()
        stop = top + 1
        if stop in g.outcomes:
            oc = g.outcomes[stop]
            reason = _SHORT_STOP.get(oc.name, oc.status)
            if oc.name == "tissue_identity" and \
                    str(g.record.get("tissue_verdict")) != "MOE":
                reason = "VNO identity unconfirmed"
            if oc.name == "sort_validation" and \
                    str(g.record.get("library_status")) in ("FAILED", "DEGENERATE"):
                reason = f"library {str(g.record.get('library_status')).lower()}"
            colour = (COLORS["secondary"] if oc.status == NO_DATA
                      else COLORS["alarm"])
            axr.text(0.02, y, f"stops at tier {stop} · {reason}", va="center",
                     ha="left", fontsize=SIZES[2], color=colour)
        else:
            axr.text(0.02, y, "all five tiers reported", va="center", ha="left",
                     fontsize=SIZES[2], color=COLORS["focal"])

    ax.set_yticks(range(n))
    ax.set_yticklabels([f"{g.sample}  ({g.trial})" for g in gates][::-1])
    ax.set_xticks(TIERS)
    ax.set_xticklabels([f"{t}\n{TIER_NAMES[t].replace('_', chr(10))}"
                        for t in TIERS], fontsize=SIZES[2])
    ax.set_xlim(-0.5, 4.5)
    ax.set_ylim(-0.9, n - 0.3)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", length=2)
    ax.set_xlabel("tier — reported only if every tier to its left passed",
                  labelpad=8)
    axr.set_xlim(0, 1)
    axr.set_ylim(-0.9, n - 0.3)
    axr.axis("off")

    ax.set_title(f"Tier reached per library: {n_tier1} of {n} clear tissue identity "
                 f"and sort validation,\n{n_tier3} are eligible for a cluster-level "
                 "VR call, and the gate names what stopped the rest",
                 loc="left", pad=12)
    handles = [Line2D([], [], ls="", marker=st["marker"], markerfacecolor=st["fill"],
                      markeredgecolor=st["color"], ms=6, label=st["label"])
               for st in (_STATUS_STYLE[k] for k in
                          (PASS, PASS_CAVEAT, "FAIL", SUPPRESSED, NO_DATA))]
    fig.subplots_adjust(left=0.235, right=0.995, top=0.845,
                        bottom=1.55 / (0.42 * n + 2.9))
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.02, 0.008),
               ncol=5, fontsize=SIZES[1], handletextpad=0.4, columnspacing=1.1)
    fig.savefig(out_path)
    findings = _verify_no_overlap(fig, "tier_overview")
    plt.close(fig)
    return {"path": out_path, "overlaps": findings, "n_samples": n,
            "n_reaching_tier1": n_tier1, "n_reaching_tier3": n_tier3}


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------


def build_all(cfg, which: Optional[Sequence[str]] = None,
              out_dir: Optional[str] = None) -> Dict[str, Any]:
    out_dir = out_dir or os.path.join(cfg["work"], "results", "figures")
    os.makedirs(out_dir, exist_ok=True)
    state = gather(cfg)
    want = set(which) if which else {"cluster", "fractions", "tier"}
    res: Dict[str, Any] = {}
    if "cluster" in want:
        res["cluster"] = fig_cluster_heatmap(
            state, os.path.join(out_dir, "vr_cluster_heatmap.png"))
    if "fractions" in want:
        res["fractions"] = fig_within_cluster_fractions(
            state, os.path.join(out_dir, "vr_within_cluster_fractions.png"))
    if "tier" in want:
        res["tier"] = fig_tier_overview(
            state, os.path.join(out_dir, "vr_tier_overview.png"))
    return res


def _main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--fig", action="append", default=None,
                    choices=["cluster", "fractions", "tier"])
    ap.add_argument("--out-dir", default=None)
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    res = build_all(cfg, args.fig, args.out_dir)
    for k, v in res.items():
        print(f"{k}\t{v.get('path')}\tpending={v.get('pending', False)}\t"
              f"overlaps={len(v.get('overlaps') or [])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
