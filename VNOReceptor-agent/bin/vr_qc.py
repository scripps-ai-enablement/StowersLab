#!/usr/bin/env python3
"""
vr_qc.py -- library viability / technical QC module for the vr_analysis pipeline.

Parses, per sample, the three technical-QC sources produced by nf-core/rnaseq
(aligner star_salmon):

  1. star_salmon/log/<sample>.Log.final.out        (STAR alignment metrics)
  2. star_salmon/qualimap/<sample>/rnaseq_qc_results.txt  (coverage / genomic origin)
  3. multiqc/star_salmon/multiqc_report_data/multiqc_general_stats.txt
     (gap-filler: duplication, error rate, proper pairs. MultiQC column names
     change between versions, so columns are DISCOVERED at runtime by suffix
     matching and the discovered set is reported, never assumed.)

Multi-mapping is a named open question for this project, so both multi-mapping
channels are first-class outputs and are reported separately:

  * `pct_multi_loci`      -- reads mapped to multiple loci. These are RETAINED
                             in the BAM and handed to Salmon's EM, which is what
                             spreads reads across VR paralogs.
  * `pct_too_many_loci`   -- reads exceeding STAR's --outFilterMultimapNmax and
                             DISCARDED outright.

`dominant_multimap_channel` names whichever is larger per sample.

Usage
-----
    python vr_qc.py --trial trial2 [--out-dir DIR] [--config PATH]
    python vr_qc.py --all-trials

    # as a module
    from vr_qc import technical_qc_table
    df = technical_qc_table(cfg, "trial2")
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from typing import Any, Dict, List, Optional

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vr_config import (  # noqa: E402
    load_config, samples_of, threshold, trial_paths, trials_of,
)

__all__ = [
    "parse_star_log", "parse_qualimap", "discover_multiqc_columns",
    "parse_multiqc_general_stats", "technical_qc_table", "TECH_COLUMNS",
]

# ----------------------------------------------------------------------------
# small parsing helpers
# ----------------------------------------------------------------------------

def _num(text: Optional[str]) -> Optional[float]:
    """'80.81%' -> 80.81 ; '29,706,525' -> 29706525.0 ; junk -> None."""
    if text is None:
        return None
    s = str(text).strip().replace(",", "").replace("%", "")
    if s in ("", "NA", "N/A", "nan", "-"):
        return None
    try:
        return float(s)
    except ValueError:
        return None


# ----------------------------------------------------------------------------
# 1. STAR Log.final.out
# ----------------------------------------------------------------------------

# STAR label (lowercased, whitespace-collapsed) -> our column name
_STAR_FIELDS = {
    "number of input reads": "input_reads",
    "average input read length": "avg_input_read_length",
    "uniquely mapped reads number": "uniquely_mapped_reads",
    "uniquely mapped reads %": "pct_uniquely_mapped",
    "average mapped length": "avg_mapped_length",
    "mismatch rate per base, %": "pct_mismatch_rate",
    "number of reads mapped to multiple loci": "n_multi_loci",
    "% of reads mapped to multiple loci": "pct_multi_loci",
    "number of reads mapped to too many loci": "n_too_many_loci",
    "% of reads mapped to too many loci": "pct_too_many_loci",
    "% of reads unmapped: too short": "pct_unmapped_too_short",
    "% of reads unmapped: too many mismatches": "pct_unmapped_mismatches",
    "% of reads unmapped: other": "pct_unmapped_other",
    "number of splices: total": "n_splices_total",
    "% of chimeric reads": "pct_chimeric",
}


def parse_star_log(path: str) -> Dict[str, Any]:
    """Parse a STAR Log.final.out into our canonical column names."""
    out: Dict[str, Any] = {"star_log_found": False}
    if not os.path.exists(path):
        return out
    out["star_log_found"] = True
    with open(path) as fh:
        for line in fh:
            if "|" not in line:
                continue
            label, _, value = line.partition("|")
            key = re.sub(r"\s+", " ", label.strip().lower())
            col = _STAR_FIELDS.get(key)
            if col:
                out[col] = _num(value)
    return out


# ----------------------------------------------------------------------------
# 2. qualimap rnaseq_qc_results.txt
# ----------------------------------------------------------------------------

# qualimap label -> (column, wants_percent_in_parens)
_QMAP_FIELDS = {
    "5'-3' bias": ("bias_5p3p", False),
    "5' bias": ("bias_5p", False),
    "3' bias": ("bias_3p", False),
    "reads aligned (left/right)": ("reads_aligned_pair_str", False),
    "read pairs aligned": ("read_pairs_aligned", False),
    "total alignments": ("total_alignments", False),
    "secondary alignments": ("secondary_alignments", False),
    "non-unique alignments": ("non_unique_alignments", False),
    "aligned to genes": ("reads_aligned_to_genes", False),
    "ambiguous alignments": ("ambiguous_alignments", False),
    "no feature assigned": ("no_feature_assigned", False),
    "exonic": ("n_exonic", True),
    "intronic": ("n_intronic", True),
    "intergenic": ("n_intergenic", True),
    "overlapping exon": ("n_overlapping_exon", True),
    "reads at junctions": ("reads_at_junctions", False),
    "duplication rate": ("duplication_rate_pct", False),
}
# percent captured from a trailing "(83.82%)"
_QMAP_PCT_COL = {
    "n_exonic": "pct_exonic",
    "n_intronic": "pct_intronic",
    "n_intergenic": "pct_intergenic",
    "n_overlapping_exon": "pct_overlapping_exon",
}


def parse_qualimap(path: str) -> Dict[str, Any]:
    """Parse qualimap rnaseq_qc_results.txt (`label = value (pct%)` lines)."""
    out: Dict[str, Any] = {"qualimap_found": False}
    if not os.path.exists(path):
        return out
    out["qualimap_found"] = True
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if "=" not in line or line.startswith(">"):
                continue
            label, _, value = line.partition("=")
            key = re.sub(r"\s+", " ", label.strip().lower())
            hit = _QMAP_FIELDS.get(key)
            if not hit:
                continue
            col, _ = hit
            value = value.strip()
            # trailing percentage, e.g. "42,784,461 (83.82%)"
            m = re.search(r"\(([\d.]+)\s*%\)", value)
            if m and col in _QMAP_PCT_COL:
                out[_QMAP_PCT_COL[col]] = _num(m.group(1))
            base = re.sub(r"\(.*?\)", "", value).strip()
            if col == "reads_aligned_pair_str":
                # "29,706,525 / 29,704,929" -> keep the left mate count
                left = base.split("/")[0]
                out["reads_aligned_left"] = _num(left)
                continue
            out[col] = _num(base)
    return out


# ----------------------------------------------------------------------------
# 3. multiqc_general_stats.txt -- runtime column discovery
# ----------------------------------------------------------------------------

# Our column name -> list of case-insensitive column-name SUFFIXES to try, in
# priority order. MultiQC prefixes columns with the module name and the exact
# prefix varies by version, so we match on the suffix after the last '-'.
_MULTIQC_WANTED: Dict[str, List[str]] = {
    "pct_duplication": ["picard_mark_duplicates-percent_duplication",
                        "percent_duplication", "percent_duplicates"],
    "dupradar_intercept": ["dupradar_intercept"],
    "pct_proper_pairs": ["proper_pairs_percent", "reads_properly_paired_percent"],
    "samtools_error_rate": ["error_rate"],
    "pct_reads_mq0": ["reads_mq0_percent"],
    "non_primary_alignments": ["non_primary_alignments"],
    "pct_gc": ["percent_gc"],
    "pct_adapter_trimmed": ["percent_trimmed"],
    "raw_total_sequences": ["raw_total_sequences", "total_sequences"],
    "mqc_star_multimapped": ["star-multimapped", "multimapped"],
    "mqc_pct_uniquely_mapped": ["star-uniquely_mapped_percent",
                                "uniquely_mapped_percent"],
    "mqc_bias_5p3p": ["qualimap_rnaseq-5_3_bias", "5_3_bias"],
}


def discover_multiqc_columns(columns) -> Dict[str, str]:
    """
    Map our canonical names onto whatever columns this MultiQC version emitted.
    Returns {our_name: actual_column}. Missing entries are simply absent -- the
    caller reports which were found rather than assuming a fixed schema.
    """
    low = {str(c).lower(): str(c) for c in columns}
    found: Dict[str, str] = {}
    for ours, candidates in _MULTIQC_WANTED.items():
        for cand in candidates:
            cl = cand.lower()
            # exact match on the full column name first
            if cl in low:
                found[ours] = low[cl]
                break
            # then suffix match after the module prefix
            hits = [orig for lc, orig in low.items()
                    if lc == cl or lc.endswith("-" + cl)]
            if hits:
                found[ours] = sorted(hits, key=len)[0]
                break
    return found


def parse_multiqc_general_stats(path: str, samples) -> Dict[str, Any]:
    """
    Returns {'per_sample': {sample: {col: val}}, 'columns_found': {...},
             'columns_missing': [...], 'all_columns': [...], 'found': bool}

    MultiQC lists per-sample rows AND per-mate rows (`<sample>_1`, `<sample>_2`)
    for read-level modules. We take the exact-sample row for everything, and
    fall back to the mean of the mate rows for a column that is only populated
    at mate level (fastqc columns behave this way).
    """
    res: Dict[str, Any] = {"found": False, "per_sample": {}, "columns_found": {},
                           "columns_missing": [], "all_columns": []}
    if not os.path.exists(path):
        res["columns_missing"] = sorted(_MULTIQC_WANTED)
        return res
    df = pd.read_csv(path, sep="\t", dtype=str)
    res["found"] = True
    samp_col = df.columns[0]
    df = df.set_index(samp_col)
    res["all_columns"] = [str(c) for c in df.columns]
    colmap = discover_multiqc_columns(df.columns)
    res["columns_found"] = colmap
    res["columns_missing"] = sorted(set(_MULTIQC_WANTED) - set(colmap))
    for s in samples:
        row: Dict[str, Any] = {}
        mates = [i for i in df.index if str(i).startswith(s + "_")]
        for ours, actual in colmap.items():
            val = None
            if s in df.index:
                val = _num(df.at[s, actual])
            if val is None and mates:
                vals = [v for v in (_num(df.at[m, actual]) for m in mates)
                        if v is not None]
                if vals:
                    val = sum(vals) / len(vals)
            row[ours] = val
        res["per_sample"][s] = row
    return res


# ----------------------------------------------------------------------------
# assembly + flagging
# ----------------------------------------------------------------------------

TECH_COLUMNS = [
    "trial", "sample", "cell_type", "n_cells", "prep_status", "platform",
    "input_reads", "avg_input_read_length",
    "uniquely_mapped_reads", "pct_uniquely_mapped", "avg_mapped_length",
    "pct_mismatch_rate",
    "n_multi_loci", "pct_multi_loci",
    "n_too_many_loci", "pct_too_many_loci",
    "dominant_multimap_channel", "multimap_channel_ratio",
    "pct_unmapped_too_short", "pct_unmapped_other",
    "bias_5p3p", "bias_5p", "bias_3p",
    "pct_exonic", "pct_intronic", "pct_intergenic",
    "reads_aligned_to_genes", "secondary_alignments", "non_unique_alignments",
    "pct_duplication", "dupradar_intercept", "pct_proper_pairs",
    "samtools_error_rate", "pct_gc", "pct_adapter_trimmed",
    "tech_flags", "tech_verdict",
    "star_log_found", "qualimap_found", "multiqc_found",
]


def _tech_flags(row: Dict[str, Any], cfg) -> List[str]:
    """Threshold-driven technical flags. All cutoffs come from config."""
    T = lambda k: threshold(cfg, k, section="qc_thresholds")  # noqa: E731
    flags: List[str] = []

    def below(col, key, label):
        v, lim = row.get(col), T(key)
        if v is not None and v < lim:
            flags.append(f"{label}({v:g}<{lim:g})")

    def above(col, key, label):
        v, lim = row.get(col), T(key)
        if v is not None and v > lim:
            flags.append(f"{label}({v:g}>{lim:g})")

    if not row.get("star_log_found"):
        flags.append("star_log_missing")
    if not row.get("qualimap_found"):
        flags.append("qualimap_missing")
    below("input_reads", "min_input_reads", "low_input_reads")
    below("pct_uniquely_mapped", "min_uniquely_mapped_pct", "low_unique_mapping")
    above("pct_multi_loci", "max_multi_loci_pct", "high_multi_loci")
    above("pct_too_many_loci", "max_too_many_loci_pct", "high_too_many_loci")
    above("pct_unmapped_too_short", "max_unmapped_too_short_pct", "high_unmapped_short")
    above("pct_mismatch_rate", "max_mismatch_rate_pct", "high_mismatch")
    below("pct_exonic", "min_exonic_pct", "low_exonic")
    above("pct_intergenic", "max_intergenic_pct", "high_intergenic")

    bias = row.get("bias_5p3p")
    if bias is not None:
        lo, hi = T("bias_5p3p_low"), T("bias_5p3p_high")
        if bias < lo:
            flags.append(f"bias_3p_skew({bias:g}<{lo:g})")
        elif bias > hi:
            flags.append(f"bias_5p_skew({bias:g}>{hi:g})")
    return flags


def technical_qc_table(cfg, trial: str) -> pd.DataFrame:
    """Per-sample technical QC for one trial. One row per configured sample."""
    P = trial_paths(cfg, trial)
    smeta = samples_of(cfg, trial)
    mqc = parse_multiqc_general_stats(P["multiqc_data"] + "/multiqc_general_stats.txt",
                                      list(smeta))
    rows = []
    for s, meta in smeta.items():
        row: Dict[str, Any] = {
            "trial": trial, "sample": s, "platform": P["platform"],
            "cell_type": meta["cell_type"], "n_cells": meta["n_cells"],
            "prep_status": meta["prep_status"],
        }
        row.update(parse_star_log(os.path.join(P["star_logs"], f"{s}.Log.final.out")))
        row.update(parse_qualimap(os.path.join(P["qualimap"], s, "rnaseq_qc_results.txt")))
        row.update(mqc["per_sample"].get(s, {}))
        row["multiqc_found"] = mqc["found"] and s in mqc["per_sample"]

        # ---- the multi-mapping diagnostic (project open question 1) ----
        m, tm = row.get("pct_multi_loci"), row.get("pct_too_many_loci")
        if m is not None and tm is not None:
            row["dominant_multimap_channel"] = (
                "multiple_loci_retained" if m > tm else
                "too_many_loci_discarded" if tm > m else "equal")
            row["multimap_channel_ratio"] = (m / tm) if tm and tm > 0 else float("inf")
        else:
            row["dominant_multimap_channel"] = None
            row["multimap_channel_ratio"] = None

        flags = _tech_flags(row, cfg)
        row["tech_flags"] = ";".join(flags) if flags else ""
        row["tech_verdict"] = "PASS" if not flags else (
            "FAIL" if any(f.startswith(("low_input_reads", "low_unique_mapping",
                                        "star_log_missing")) for f in flags)
            else "WARN")
        rows.append(row)

    df = pd.DataFrame(rows)
    for c in TECH_COLUMNS:
        if c not in df.columns:
            df[c] = None
    extra = [c for c in df.columns if c not in TECH_COLUMNS]
    df = df[TECH_COLUMNS + extra].sort_values(["trial", "sample"]).reset_index(drop=True)
    df.attrs["multiqc_columns_found"] = mqc["columns_found"]
    df.attrs["multiqc_columns_missing"] = mqc["columns_missing"]
    df.attrs["multiqc_all_columns"] = mqc["all_columns"]
    return df


def _main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Technical QC / library viability")
    ap.add_argument("--trial")
    ap.add_argument("--all-trials", action="store_true")
    ap.add_argument("--config")
    ap.add_argument("--out-dir")
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    trials = trials_of(cfg) if a.all_trials else [a.trial]
    if not trials or trials == [None]:
        ap.error("give --trial NAME or --all-trials")
    for t in trials:
        df = technical_qc_table(cfg, t)
        out_dir = a.out_dir or trial_paths(cfg, t)["out_dir"]
        os.makedirs(out_dir, exist_ok=True)
        path = os.path.join(out_dir, "technical_qc.tsv")
        df.to_csv(path, sep="\t", index=False)
        print(f"[vr_qc] {t}: {len(df)} samples -> {path}")
        print(f"[vr_qc] {t}: multiqc columns found: "
              f"{sorted(df.attrs['multiqc_columns_found'])}")
        print(f"[vr_qc] {t}: multiqc columns missing: "
              f"{df.attrs['multiqc_columns_missing']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
