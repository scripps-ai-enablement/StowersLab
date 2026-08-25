#!/usr/bin/env python3
"""
vr_markers.py -- sort validation (Rule 1) and population identification (Rule 4)
for the vr_analysis pipeline.

CPM DENOMINATOR (stated explicitly, per house rule 5)
-----------------------------------------------------
    CPM(gene, sample) = count(gene, sample) / total_count(sample) * 1e6
    total_count(sample) = column sum over ALL 45,706 genes in
                          star_salmon/salmon.merged.gene_counts.tsv

That is the *unscaled* Salmon merged gene counts (NOT gene_counts_scaled.tsv,
NOT TPM). Rationale: the previously recorded reference values for this project
(nontarget Trpc2 ~0.2 CPM; target 1029-1476 CPM) were derived this way, and an
all-gene library denominator is the only choice that keeps CPM comparable across
samples with very different VR content. `--counts-table` lets a caller
re-derive the same table from the scaled counts for a sensitivity check; the
choice is recorded in the output header and in `df.attrs['cpm_denominator']`.

Genes are matched on the `gene_name` column, case-insensitively. When several
gene_id rows share a gene_name (Salmon merged tables can carry duplicates) the
counts are SUMMED and the multiplicity is reported.

Rules implemented
-----------------
Rule 1  sort validation   target: Trpc2 > target_trpc2_min CPM = PASS;
                                  > target_trpc2_concern  = CONCERN; else FAIL
                          nontarget: Trpc2 < nontarget_trpc2_max = PASS else FAIL
Rule 3  actin sanity      near-zero actin sum => empty/failed library, not biology
Rule 4  population call   Gnai2:Gnao1 > v1r_dominant_ratio_min => V1R_dominant
                          < v2r_dominant_ratio_max             => V2R_dominant
                          otherwise                            => mixed
                          Gnao1 == 0 is reported as ratio = inf with a note
                          (preserved, never dropped and never a div-by-zero).
Failed-library gate       actin_sum_cpm < failed_lib_actin_cpm_max AND
                          trpc2_cpm < failed_lib_trpc2_cpm_max in a TARGET
                          sample => library_status = FAILED. Downstream code
                          must check `library_status` / `suppress_biology`
                          before reporting any biological conclusion.

Usage
-----
    python vr_markers.py --trial trial2
    python vr_markers.py --all-trials
    from vr_markers import marker_cpm_table
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vr_config import (  # noqa: E402
    load_config, marker_genes, samples_of, threshold, trial_paths, trials_of,
)

__all__ = [
    "load_gene_counts", "compute_cpm", "marker_cpm_table",
    "sort_verdict", "population_call", "library_status", "MARKER_COLUMNS",
]

CPM_DENOMINATOR_NOTE = (
    "CPM = count / (column sum over ALL genes in the merged Salmon gene-count "
    "table) * 1e6; source table = {table}"
)


# ----------------------------------------------------------------------------
# counts -> CPM
# ----------------------------------------------------------------------------

def load_gene_counts(path: str, samples: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Read salmon.merged.gene_counts.tsv. Returns a counts frame indexed by
    gene_name (duplicates summed) with sample columns. `gene_id` multiplicity per
    name is stashed in .attrs['name_multiplicity'].
    """
    df = pd.read_csv(path, sep="\t")
    if "gene_name" not in df.columns:
        raise ValueError(f"{path}: no gene_name column (got {list(df.columns)[:5]})")
    id_cols = [c for c in ("gene_id", "gene_name") if c in df.columns]
    sample_cols = [c for c in df.columns if c not in id_cols]
    if samples:
        missing = [s for s in samples if s not in sample_cols]
        if missing:
            raise KeyError(f"{path}: configured samples absent from table: {missing}")
        sample_cols = [s for s in samples if s in sample_cols]
    counts = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    counts.insert(0, "gene_name", df["gene_name"].astype(str))
    mult = counts.groupby("gene_name").size()
    agg = counts.groupby("gene_name", sort=False).sum(numeric_only=True)
    agg.attrs["name_multiplicity"] = mult[mult > 1].to_dict()
    agg.attrs["n_gene_rows"] = len(df)
    agg.attrs["source_table"] = path
    return agg


def compute_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    """CPM with an all-gene library-size denominator (see module docstring)."""
    totals = counts.sum(axis=0)
    if (totals <= 0).any():
        bad = list(totals[totals <= 0].index)
        raise ValueError(f"zero total counts for samples {bad}; cannot compute CPM")
    cpm = counts.divide(totals, axis=1) * 1e6
    cpm.attrs.update(counts.attrs)
    cpm.attrs["library_totals"] = totals.to_dict()
    return cpm


def _gene_cpm(cpm: pd.DataFrame, gene: str, sample: str) -> float:
    """Case-insensitive gene_name lookup. Absent gene -> 0.0 (not NaN)."""
    if gene in cpm.index:
        return float(cpm.at[gene, sample])
    lower = {str(i).lower(): i for i in cpm.index}
    key = lower.get(gene.lower())
    return float(cpm.at[key, sample]) if key is not None else 0.0


# ----------------------------------------------------------------------------
# rules
# ----------------------------------------------------------------------------

def prefix_sum_cpm(cpm: pd.DataFrame, prefix: str, sample: str) -> Any:
    """
    Sum CPM over every gene_name starting with `prefix` (case-insensitive).
    Returns (sum, n_genes_matched). Used for the receptor-family panels
    (Olfr* = main olfactory, Vmn1r*/Vmn2r* = vomeronasal), where the family
    total is the tissue signal and no individual paralog is trustworthy.
    """
    pl = prefix.lower()
    idx = [i for i in cpm.index if str(i).lower().startswith(pl)]
    if not idx:
        return 0.0, 0
    return float(cpm.loc[idx, sample].sum()), len(idx)


def tissue_identity(panel: Dict[str, float], cfg) -> Dict[str, Any]:
    """
    TISSUE IDENTITY -- runs BEFORE sort validation (Rule 1).

    Why this check exists, and why it is first: OMP-Cre drives GFP in mature
    MAIN OLFACTORY (MOE) neurons as well as VNO neurons, so Omp is a
    TISSUE-SHARED marker and cannot discriminate the two. A GFP+ sort taken from
    the wrong tissue -- or from VNO contaminated by MOE -- yields Omp-high,
    Trpc2-zero cells. Reading a high Omp as evidence of a VNO population is
    exactly the error this check prevents: Rule 1 asks "did the sort enrich VNO
    neurons?", which is only a meaningful question once the tissue is VNO.

    Panels compared on their MAXIMUM member, not their sum, so one strongly
    expressed marker establishes the tissue without being diluted by family
    members that happen to be silent in a given cell.

      MOE panel          Olfr* summed, Adcy3, Cnga2, Gnal
      VNO-specific panel Trpc2, Vmn1r* summed, Vmn2r* summed
      SHARED (excluded)  Omp -- mature olfactory neurons of BOTH tissues
      NOT tissue markers Gnai2 / Gnao1 -- these split V1R vs V2R populations
                         WITHIN the VNO and are broadly expressed elsewhere
                         (Gnai2 reads 165 CPM in a GFP- library with no VR
                         signal at all), so they carry no tissue information.

    An absolute CPM floor is required, not a bare ratio: below the floor a panel
    is noise and a ratio of two noise values is meaningless.
    """
    floor = threshold(cfg, "tissue_panel_floor_cpm", 100.0)
    mult = threshold(cfg, "tissue_dominance_ratio", 3.0)
    moe = max(panel["Olfr_sum_cpm"], panel["Adcy3_cpm"],
              panel["Cnga2_cpm"], panel["Gnal_cpm"])
    vno = max(panel["Trpc2_cpm"], panel["Vmn1r_sum_cpm"], panel["Vmn2r_sum_cpm"])

    if moe < floor and vno < floor:
        v = "no_tissue_signal"
        why = (f"both panels below the {floor:g} CPM floor (MOE max {moe:.1f}, "
               f"VNO-specific max {vno:.1f}): this library carries no tissue "
               "information either way -- not evidence of the wrong tissue")
    elif vno >= floor and moe < floor:
        v = "VNO"
        why = (f"VNO-specific panel {vno:.1f} CPM at/above floor with MOE panel "
               f"{moe:.1f} below it")
    elif moe >= floor and vno < floor:
        v = "MOE"
        why = (f"MAIN OLFACTORY panel {moe:.1f} CPM at/above floor with "
               f"VNO-specific panel {vno:.1f} below it -- wrong tissue")
    elif vno > mult * moe:
        v = "VNO_dominant_mixed"
        why = (f"both panels above floor but VNO-specific {vno:.1f} CPM exceeds "
               f"MOE {moe:.1f} CPM by more than {mult:g}x: VNO with a real minor "
               "MOE component (sort-purity information, not a tissue failure)")
    elif moe > mult * vno:
        v = "MOE_dominant_mixed"
        why = (f"both panels above floor but MOE {moe:.1f} CPM exceeds "
               f"VNO-specific {vno:.1f} CPM by more than {mult:g}x: predominantly "
               "wrong tissue with a minor VNO component")
    else:
        v = "ambiguous_mixed"
        why = (f"both panels above floor within {mult:g}x of each other "
               f"(MOE {moe:.1f}, VNO-specific {vno:.1f} CPM): tissue identity "
               "cannot be assigned")
    return {"tissue_verdict": v, "tissue_reason": why,
            "moe_panel_max_cpm": moe, "vno_panel_max_cpm": vno,
            "tissue_floor_cpm": floor}


def sort_verdict(cell_type: str, trpc2_cpm: float, cfg,
                 tissue_verdict: str = None) -> Dict[str, str]:
    """
    Rule 1, gated behind the tissue check.

    If the tissue is MOE, Rule 1 does not apply: a Trpc2 of zero in main
    olfactory epithelium is the expected value for that tissue, not a failed
    VNO sort. The verdict becomes FAIL_WRONG_TISSUE, which says the sample is
    unusable for VR biology AND that the remedy is wet-lab (dissection / sort
    gate), not a re-run of the quantification.
    """
    tmin = threshold(cfg, "target_trpc2_min")
    tcon = threshold(cfg, "target_trpc2_concern")
    nmax = threshold(cfg, "nontarget_trpc2_max")

    if tissue_verdict in ("MOE", "MOE_dominant_mixed"):
        return {"sort_verdict": "FAIL_WRONG_TISSUE",
                "sort_reason": (
                    f"tissue identity is {tissue_verdict}: the main-olfactory "
                    "marker panel is above the tissue floor while the "
                    "VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is not. Trpc2 = "
                    f"{trpc2_cpm:.2f} CPM here is the EXPECTED value for main "
                    "olfactory epithelium, so Rule 1 does not apply -- this is "
                    "not a failed VNO sort but the wrong tissue. Note Omp cannot "
                    "distinguish the two tissues (it marks mature olfactory "
                    "neurons of both). No VR biology is reportable; the fix is "
                    "wet-lab (dissection / sort gate), NOT re-quantification.")}
    if cell_type == "target":
        if trpc2_cpm > tmin:
            return {"sort_verdict": "PASS",
                    "sort_reason": f"target Trpc2 {trpc2_cpm:.1f} CPM > {tmin:g}"}
        if trpc2_cpm > tcon:
            return {"sort_verdict": "CONCERN",
                    "sort_reason": (f"target Trpc2 {trpc2_cpm:.1f} CPM above "
                                    f"{tcon:g} but below {tmin:g}")}
        return {"sort_verdict": "FAIL",
                "sort_reason": (f"target Trpc2 {trpc2_cpm:.2f} CPM <= {tcon:g}; "
                                "no VNO-neuron signal")}
    if cell_type == "nontarget":
        if trpc2_cpm < nmax:
            return {"sort_verdict": "PASS",
                    "sort_reason": f"nontarget Trpc2 {trpc2_cpm:.2f} CPM < {nmax:g}"}
        return {"sort_verdict": "FAIL",
                "sort_reason": (f"nontarget Trpc2 {trpc2_cpm:.1f} CPM >= {nmax:g}; "
                                "target contamination in the GFP- fraction")}
    return {"sort_verdict": "NA",
            "sort_reason": f"cell_type {cell_type!r} has no sort expectation"}


def population_call(gnai2: float, gnao1: float, cfg,
                    gnao1_reads: float = None) -> Dict[str, Any]:
    """
    Rule 4. Gnao1 == 0 -> ratio inf, preserved with an explicit note.

    `gnao1_reads` is the RAW count behind the Gnao1 CPM. The ratio is
    denominator-free (it is invariant to any per-sample normalisation), so its
    only real source of error is Poisson noise on a tiny read count: at 7 reads a
    one-read shift moves the ratio by ~15%. When the denominator rests on fewer
    than `ratio_min_support_reads` raw reads the ratio is flagged low-support --
    the CALL is still reliable when the ratio is orders of magnitude past the
    threshold, but the printed ratio must not be quoted as a precise figure.
    """
    rmin = threshold(cfg, "v1r_dominant_ratio_min")
    rmax = threshold(cfg, "v2r_dominant_ratio_max")
    min_support = threshold(cfg, "ratio_min_support_reads", 10,
                            section="qc_thresholds")
    note = ""
    if gnao1 == 0 and gnai2 == 0:
        return {"gnai2_gnao1_ratio": float("nan"), "population_call": "undetermined",
                "population_note": "both Gnai2 and Gnao1 are 0 CPM -- no population signal",
                "ratio_low_support": True, "ratio_support_reads": gnao1_reads}
    if gnao1 == 0:
        ratio = float("inf")
        note = (f"Gnao1 = 0 CPM; ratio reported as infinity (Gnai2 = {gnai2:.1f} CPM). "
                "Not a division error -- V2R marker genuinely absent.")
    else:
        ratio = gnai2 / gnao1
    if not math.isfinite(ratio) or ratio > rmin:
        call = "V1R_dominant"
    elif ratio < rmax:
        call = "V2R_dominant"
    else:
        call = "mixed"
    if not note:
        note = f"Gnai2:Gnao1 = {ratio:.2f} (thresholds >{rmin:g} V1R, <{rmax:g} V2R)"
    low = gnao1_reads is not None and gnao1_reads < min_support
    if low and math.isfinite(ratio):
        note += (f" -- LOW SUPPORT: only {gnao1_reads:g} raw Gnao1 reads, so the "
                 "ratio magnitude is Poisson-unstable (+/-1 read shifts it by "
                 f"~{100.0/max(gnao1_reads,1):.0f}%); the direction of the call is "
                 "robust, the exact ratio is not.")
    return {"gnai2_gnao1_ratio": ratio, "population_call": call,
            "population_note": note, "ratio_low_support": bool(low),
            "ratio_support_reads": gnao1_reads}


def library_status(cell_type: str, trpc2_cpm: float, actin_sum_cpm: float,
                   cfg, total_counts: float = None) -> Dict[str, Any]:
    """
    Failed-library gate (house rule 4) + Rule 3 actin sanity + a
    degenerate-library-size gate.

    The size gate exists because CPM is a RATIO: a library with only a few tens
    of thousands of assigned counts can show a huge actin CPM (all of a tiny
    pie) and so slip past the actin arm of the failed-library gate while
    carrying no usable information. Such a library is DEGENERATE and biology is
    suppressed regardless of its CPM values.
    """
    amax = threshold(cfg, "failed_lib_actin_cpm_max")
    tmax = threshold(cfg, "failed_lib_trpc2_cpm_max")
    min_total = threshold(cfg, "min_assigned_counts", 1_000_000,
                          section="qc_thresholds")
    if total_counts is not None and total_counts < min_total:
        return {"library_status": "DEGENERATE", "suppress_biology": True,
                "library_reason": (f"only {total_counts:,.0f} counts assigned across all "
                                   f"genes (< {min_total:,.0f}): CPM values are ratios of "
                                   "a near-empty library and carry no information, "
                                   f"however large they look (actin sum reads "
                                   f"{actin_sum_cpm:,.0f} CPM here). Biology suppressed.")}
    if cell_type == "target" and actin_sum_cpm < amax and trpc2_cpm < tmax:
        return {"library_status": "FAILED", "suppress_biology": True,
                "library_reason": (f"target sample with actin sum {actin_sum_cpm:.1f} "
                                   f"CPM < {amax:g} AND Trpc2 {trpc2_cpm:.2f} CPM "
                                   f"< {tmax:g}: empty/failed library, no biology "
                                   "may be reported")}
    if actin_sum_cpm < amax:
        return {"library_status": "SUSPECT", "suppress_biology": False,
                "library_reason": (f"actin sum {actin_sum_cpm:.1f} CPM < {amax:g} "
                                   "(Rule 3): near-zero housekeeping signal "
                                   "indicates a technical problem, not biology")}
    return {"library_status": "OK", "suppress_biology": False,
            "library_reason": f"actin sum {actin_sum_cpm:.1f} CPM >= {amax:g}"}


# ----------------------------------------------------------------------------
# assembly
# ----------------------------------------------------------------------------

MARKER_COLUMNS = [
    "trial", "sample", "cell_type", "n_cells", "prep_status",
    "library_total_counts",
    # tissue identity comes first: it gates Rule 1
    "tissue_verdict", "moe_panel_max_cpm", "vno_panel_max_cpm", "tissue_floor_cpm",
    "tissue_reason",
    "Trpc2_cpm", "Vmn1r_sum_cpm", "Vmn2r_sum_cpm",
    "Olfr_sum_cpm", "Adcy3_cpm", "Cnga2_cpm", "Gnal_cpm",
    "Omp_cpm",  # SHARED between MOE and VNO -- not tissue-discriminating
    "Gnai2_cpm", "Gnao1_cpm",
    "actin_sum_cpm",
    "Trpc2_reads", "Gnai2_reads", "Gnao1_reads", "trpc2_exact_zero",
    "gnai2_gnao1_ratio", "gnai2_gnao1_ratio_str",
    "ratio_low_support", "ratio_support_reads",
    "sort_verdict", "sort_reason",
    "population_call", "population_note",
    "library_status", "library_reason", "suppress_biology",
]


def marker_cpm_table(cfg, trial: str, counts_table: str = "gene_counts") -> pd.DataFrame:
    """Marker CPMs + Rule 1/3/4 verdicts, one row per configured sample."""
    P = trial_paths(cfg, trial)
    smeta = samples_of(cfg, trial)
    path = P[counts_table]
    counts = load_gene_counts(path, list(smeta))
    cpm = compute_cpm(counts)

    trpc2 = marker_genes(cfg, "pan_vno_neuron")
    omp = marker_genes(cfg, "mature_neuron")
    gnai2 = marker_genes(cfg, "v1r_population")
    gnao1 = marker_genes(cfg, "v2r_population")
    actins = marker_genes(cfg, "housekeeping_actin")
    # tissue panels. Prefix families are summed across the whole family; the
    # MOE singles are individual genes. Defaults keep the module runnable
    # against a config predating the tissue block.
    mk = cfg["markers"]
    vr_families = list(mk.get("vr_families", ["Vmn1r", "Vmn2r"]))
    moe_families = list(mk.get("moe_receptor_families", ["Olfr"]))
    moe_singles = list(mk.get("moe_transduction", ["Adcy3", "Cnga2", "Gnal"]))
    tissue_shared = list(mk.get("tissue_shared", ["Omp"]))

    rows = []
    for s, meta in smeta.items():
        r: Dict[str, Any] = {
            "trial": trial, "sample": s, "cell_type": meta["cell_type"],
            "n_cells": meta["n_cells"], "prep_status": meta["prep_status"],
            "library_total_counts": float(cpm.attrs["library_totals"][s]),
        }
        v_trpc2 = sum(_gene_cpm(cpm, g, s) for g in trpc2)
        v_omp = sum(_gene_cpm(cpm, g, s) for g in omp)
        v_gnai2 = sum(_gene_cpm(cpm, g, s) for g in gnai2)
        v_gnao1 = sum(_gene_cpm(cpm, g, s) for g in gnao1)
        r["Trpc2_cpm"], r["Omp_cpm"] = v_trpc2, v_omp
        r["Gnai2_cpm"], r["Gnao1_cpm"] = v_gnai2, v_gnao1
        actin_sum = 0.0
        for g in actins:
            v = _gene_cpm(cpm, g, s)
            r[f"{g}_cpm"] = v
            actin_sum += v
        r["actin_sum_cpm"] = actin_sum

        # raw counts behind the two ratio terms, for support/Poisson checks
        r["Gnai2_reads"] = sum(_gene_cpm(counts, g, s) for g in gnai2)
        r["Gnao1_reads"] = sum(_gene_cpm(counts, g, s) for g in gnao1)
        r["Trpc2_reads"] = sum(_gene_cpm(counts, g, s) for g in trpc2)

        # ---- TISSUE IDENTITY: computed and applied BEFORE Rule 1 ----
        for fam in vr_families:
            tot, n = prefix_sum_cpm(cpm, fam, s)
            r[f"{fam}_sum_cpm"] = tot
            r[f"{fam}_n_genes"] = n
        for fam in moe_families:
            tot, n = prefix_sum_cpm(cpm, fam, s)
            r[f"{fam}_sum_cpm"] = tot
            r[f"{fam}_n_genes"] = n
        for g in moe_singles:
            r[f"{g}_cpm"] = _gene_cpm(cpm, g, s)
        tpanel = {"Trpc2_cpm": v_trpc2,
                  "Vmn1r_sum_cpm": r.get("Vmn1r_sum_cpm", 0.0),
                  "Vmn2r_sum_cpm": r.get("Vmn2r_sum_cpm", 0.0),
                  "Olfr_sum_cpm": r.get("Olfr_sum_cpm", 0.0),
                  "Adcy3_cpm": r.get("Adcy3_cpm", 0.0),
                  "Cnga2_cpm": r.get("Cnga2_cpm", 0.0),
                  "Gnal_cpm": r.get("Gnal_cpm", 0.0)}
        ti = tissue_identity(tpanel, cfg)
        r.update(ti)

        r.update(sort_verdict(meta["cell_type"], v_trpc2, cfg,
                              tissue_verdict=ti["tissue_verdict"]))
        pc = population_call(v_gnai2, v_gnao1, cfg, gnao1_reads=r["Gnao1_reads"])
        r.update(pc)
        ratio = pc["gnai2_gnao1_ratio"]
        r["gnai2_gnao1_ratio_str"] = (
            "infinity" if ratio == float("inf")
            else "undetermined" if not np.isfinite(ratio)
            else f"{ratio:.1f}:1")
        r.update(library_status(meta["cell_type"], v_trpc2, actin_sum, cfg,
                               total_counts=r["library_total_counts"]))

        # Trpc2 detectability: ENSMUSG00000100254 has a SINGLE annotated
        # transcript here, so this column records an exact zero for audit.
        # It is REPORTING ONLY -- no verdict depends on it.
        # Do NOT read a trial-wide zero as an annotation artefact: that theory
        # was tested and disproved. The same annotation yields 17,918-37,230
        # Trpc2 counts in the trial-2 VNO libraries, so a zero elsewhere is a
        # property of the sample, not the reference. In this project the cause
        # was WRONG TISSUE (main olfactory epithelium), which tier 0 detects
        # from the MOE marker panel. Check tissue_verdict before anything else.
        r["trpc2_exact_zero"] = (r["Trpc2_reads"] == 0)
        rows.append(r)

    df = pd.DataFrame(rows)
    actin_cols = [f"{g}_cpm" for g in actins]
    order = [c for c in MARKER_COLUMNS if c != "actin_sum_cpm"]
    i = order.index("Gnao1_cpm") + 1
    order = order[:i] + actin_cols + ["actin_sum_cpm"] + order[i:]
    # keep the canonical order first, then any remaining columns -- never DROP a
    # computed column just because it is not in MARKER_COLUMNS (that silently
    # discarded the *_n_genes family sizes in an earlier revision)
    ordered = [c for c in order if c in df.columns]
    df = df[ordered + [c for c in df.columns if c not in ordered]]
    df = df.sort_values(["trial", "sample"]).reset_index(drop=True)
    df.attrs["cpm_denominator"] = CPM_DENOMINATOR_NOTE.format(table=path)
    df.attrs["counts_table"] = path
    df.attrs["library_totals"] = cpm.attrs["library_totals"]
    df.attrs["n_gene_rows"] = counts.attrs["n_gene_rows"]
    df.attrs["duplicate_gene_names"] = counts.attrs["name_multiplicity"]
    df.attrs["actin_columns"] = actin_cols
    return df


def _main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Sort validation + population ID")
    ap.add_argument("--trial")
    ap.add_argument("--all-trials", action="store_true")
    ap.add_argument("--config")
    ap.add_argument("--out-dir")
    ap.add_argument("--counts-table", default="gene_counts",
                    choices=["gene_counts", "gene_counts_scaled", "gene_tpm"],
                    help="config paths key for the source table (sensitivity check)")
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    trials = trials_of(cfg) if a.all_trials else [a.trial]
    if not trials or trials == [None]:
        ap.error("give --trial NAME or --all-trials")
    for t in trials:
        df = marker_cpm_table(cfg, t, a.counts_table)
        out_dir = a.out_dir or trial_paths(cfg, t)["out_dir"]
        os.makedirs(out_dir, exist_ok=True)
        path = os.path.join(out_dir, "marker_cpm.tsv")
        with open(path, "w") as fh:
            fh.write(f"# {df.attrs['cpm_denominator']}\n")
            fh.write(f"# genes in source table: {df.attrs['n_gene_rows']}\n")
            df.to_csv(fh, sep="\t", index=False)
        print(f"[vr_markers] {t}: {len(df)} samples -> {path}")
        print(f"[vr_markers] {t}: {df.attrs['cpm_denominator']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
