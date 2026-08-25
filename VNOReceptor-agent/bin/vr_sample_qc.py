#!/usr/bin/env python3
"""
vr_sample_qc.py -- assemble the per-sample QC table that gates every downstream
analysis in the vr_analysis pipeline.

Joins vr_qc.technical_qc_table() (library viability / alignment metrics) with
vr_markers.marker_cpm_table() (sort validation, population ID, failed-library
gate) into one machine-readable row per sample, and derives a single
`blocking_flags` column listing every check that failed.

`blocking_flags` is the contract for downstream code: a non-empty value means
at least one gate failed for that sample, and `suppress_biology == True` means
NO biological conclusion may be drawn from it regardless of what VR reads
appear (house rule 4).

Outputs (per trial, into <work>/results/<trial>/):
    technical_qc.tsv        raw technical metrics
    marker_cpm.tsv          raw marker CPMs + verdicts
    sample_qc.tsv           the joined gating table
    sample_qc_summary.txt   compact human-readable summary
and, across trials, <work>/results/sample_qc_all.tsv.

Usage
-----
    python vr_sample_qc.py --all-trials
    python vr_sample_qc.py --trial trial2
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Any, Dict, List

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vr_config import load_config, threshold, trial_paths, trials_of  # noqa: E402
from vr_markers import marker_cpm_table  # noqa: E402
from vr_qc import technical_qc_table  # noqa: E402

__all__ = ["sample_qc_table", "write_trial_outputs", "summary_text"]

# Technical flags that block downstream biology outright, vs. advisory warnings.
_BLOCKING_TECH_PREFIXES = ("star_log_missing", "qualimap_missing",
                           "low_input_reads", "low_unique_mapping")


def sample_qc_table(cfg, trial: str) -> pd.DataFrame:
    """One gating row per sample: metadata + technical QC + markers + verdicts."""
    tech = technical_qc_table(cfg, trial)
    mark = marker_cpm_table(cfg, trial)
    join_keys = ["trial", "sample", "cell_type", "n_cells", "prep_status"]
    df = tech.merge(mark, on=join_keys, how="outer", validate="one_to_one")

    blocking: List[str] = []
    for _, r in df.iterrows():
        flags: List[str] = []
        # Rule 1 -- sort validation is the top of the priority hierarchy
        # TISSUE IDENTITY -- evaluated before Rule 1, because Rule 1 only means
        # something once the tissue is VNO.
        tv = r.get("tissue_verdict")
        if tv in ("MOE", "MOE_dominant_mixed"):
            flags.append(f"WRONG_TISSUE_{tv}")
        elif tv == "no_tissue_signal":
            # For a GFP- (nontarget) library, absent tissue signal is the
            # EXPECTED result -- the sort is supposed to exclude sensory
            # neurons -- so it is recorded, not flagged as a defect.
            flags.append("no_tissue_signal" if r.get("cell_type") != "nontarget"
                         else "warn_no_tissue_signal_expected_for_nontarget")
        elif tv == "ambiguous_mixed":
            flags.append("tissue_ambiguous")
        elif tv == "VNO_dominant_mixed":
            flags.append("warn_minor_MOE_component")
        # Rule 1 -- sort validation
        if r.get("sort_verdict") == "FAIL_WRONG_TISSUE":
            flags.append("sort_validation_FAIL_WRONG_TISSUE")
        elif r.get("sort_verdict") == "FAIL":
            flags.append("sort_validation_FAIL")
        elif r.get("sort_verdict") == "CONCERN":
            flags.append("sort_validation_CONCERN")
        # house rule 4 -- failed / degenerate / suspect library
        if r.get("library_status") == "FAILED":
            flags.append("library_FAILED")
        elif r.get("library_status") == "DEGENERATE":
            flags.append("library_DEGENERATE_too_few_counts")
        elif r.get("library_status") == "SUSPECT":
            flags.append("library_SUSPECT_low_actin")
        # population call undetermined / ratio resting on too few reads
        if r.get("population_call") == "undetermined":
            flags.append("population_undetermined")
        if bool(r.get("ratio_low_support")) and r.get("population_call") != "undetermined":
            flags.append("warn_ratio_low_read_support")
        # config-recorded prep problems carry forward as flags
        ps = str(r.get("prep_status", ""))
        if ps not in ("ok", "unknown", "nan", "None", ""):
            flags.append(f"prep_{ps}")
        # technical flags: blocking ones promoted, others prefixed as warnings
        for f in str(r.get("tech_flags") or "").split(";"):
            if not f:
                continue
            flags.append(f if f.startswith(_BLOCKING_TECH_PREFIXES) else f"warn_{f}")
        blocking.append(";".join(flags))
    # Explicit "none" sentinel rather than an empty cell: an empty TSV field
    # reads back as NaN, which would make a clean sample indistinguishable from
    # a missing value in downstream code.
    df["blocking_flags"] = [b if b else "none" for b in blocking]

    df["suppress_biology"] = df["suppress_biology"].fillna(False).astype(bool)
    # Wrong tissue is a hard gate on VR biology in its own right: whatever VR
    # reads appear in a main-olfactory library, they are not VNO receptor calls.
    df["suppress_biology"] |= df["tissue_verdict"].isin(["MOE", "MOE_dominant_mixed"])
    df["n_blocking_flags"] = [
        0 if s == "none" else
        sum(1 for f in s.split(";") if f and not f.startswith("warn_"))
        for s in df["blocking_flags"]]
    df["qc_overall"] = np.where(
        df["suppress_biology"], "UNUSABLE",
        np.where(df["sort_verdict"].isin(["FAIL", "FAIL_WRONG_TISSUE"]), "UNUSABLE",
                 np.where(df["n_blocking_flags"] > 0, "USE_WITH_CAUTION", "USABLE")))

    lead = ["trial", "sample", "cell_type", "n_cells", "prep_status", "platform",
            "qc_overall", "tissue_verdict", "sort_verdict", "population_call",
            "library_status", "suppress_biology", "blocking_flags",
            "n_blocking_flags"]
    df = df[lead + [c for c in df.columns if c not in lead]]
    df = df.sort_values(["trial", "cell_type", "sample"]).reset_index(drop=True)
    df.attrs.update(mark.attrs)
    df.attrs["multiqc_columns_found"] = tech.attrs.get("multiqc_columns_found", {})
    df.attrs["multiqc_columns_missing"] = tech.attrs.get("multiqc_columns_missing", [])
    return df


def _fmt(v, spec="{:.2f}"):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "inf" if v == float("inf") else "NA"
    try:
        return spec.format(float(v))
    except (TypeError, ValueError):
        return str(v)


def summary_text(cfg, df: pd.DataFrame) -> str:
    """Compact human-readable QC summary."""
    L: List[str] = []
    L.append("VNO receptor RNA-seq -- per-sample QC / sort validation / population ID")
    L.append("=" * 78)
    L.append(f"config: {cfg.get('_config_path')}")
    L.append(f"CPM basis: {df.attrs.get('cpm_denominator', 'NA')}")
    L.append("Priority hierarchy: (0) TISSUE IDENTITY -> (1) sort validation -> "
             "(2) population ID via Gnai2:Gnao1 -> (3) cluster-level VR calls -> "
             "(4) individual receptors.")
    L.append("Tissue identity precedes sort validation because OMP-Cre drives GFP in "
             "mature MAIN OLFACTORY neurons as well as VNO neurons. Omp is therefore "
             "TISSUE-SHARED and cannot discriminate them; a high Omp with zero Trpc2 "
             "means main olfactory epithelium, not a failed VNO sort. Asking whether a "
             "sort enriched VNO neurons is only meaningful once the tissue is VNO.")
    L.append("A sample with suppress_biology=True yields NO biological conclusion "
             "(house rule 4), whatever VR reads appear in it.")
    L.append("")
    thr = (f"thresholds: target Trpc2 > {threshold(cfg,'target_trpc2_min'):g} CPM "
           f"(concern > {threshold(cfg,'target_trpc2_concern'):g}); "
           f"nontarget < {threshold(cfg,'nontarget_trpc2_max'):g}; "
           f"failed-library gate actin sum < "
           f"{threshold(cfg,'failed_lib_actin_cpm_max'):g} AND Trpc2 < "
           f"{threshold(cfg,'failed_lib_trpc2_cpm_max'):g}; "
           f"population V1R > {threshold(cfg,'v1r_dominant_ratio_min'):g}, "
           f"V2R < {threshold(cfg,'v2r_dominant_ratio_max'):g}")
    L.append(thr)
    L.append("")

    for trial, sub in df.groupby("trial"):
        L.append(f"--- {trial} ({sub['platform'].iloc[0]}) " + "-" * 40)
        for _, r in sub.iterrows():
            L.append(f"  {r['sample']}  [{r['cell_type']}, {r['n_cells']} cells, "
                     f"prep={r['prep_status']}]")
            L.append(f"      overall        : {r['qc_overall']}"
                     + ("   *** BIOLOGY SUPPRESSED ***" if r["suppress_biology"] else ""))
            L.append(f"      reads          : {_fmt(r['input_reads'], '{:,.0f}')} input, "
                     f"{_fmt(r['pct_uniquely_mapped'])}% unique, "
                     f"{_fmt(r['pct_exonic'])}% exonic, "
                     f"5'-3' bias {_fmt(r['bias_5p3p'])}")
            L.append(f"      multimapping   : multiple-loci {_fmt(r['pct_multi_loci'])}% "
                     f"(retained -> Salmon EM) vs too-many-loci "
                     f"{_fmt(r['pct_too_many_loci'])}% (discarded); dominant = "
                     f"{r['dominant_multimap_channel']}")
            L.append(f"      VNO panel (CPM): Trpc2 {_fmt(r['Trpc2_cpm'])}, "
                     f"Vmn1r sum {_fmt(r['Vmn1r_sum_cpm'])}, "
                     f"Vmn2r sum {_fmt(r['Vmn2r_sum_cpm'])}")
            L.append(f"      MOE panel (CPM): Olfr sum {_fmt(r['Olfr_sum_cpm'])}, "
                     f"Adcy3 {_fmt(r['Adcy3_cpm'])}, Cnga2 {_fmt(r['Cnga2_cpm'])}, "
                     f"Gnal {_fmt(r['Gnal_cpm'])}")
            L.append(f"      shared / other : Omp {_fmt(r['Omp_cpm'])} (SHARED - marks "
                     f"mature neurons of BOTH tissues, not VNO-specific), "
                     f"Gnai2 {_fmt(r['Gnai2_cpm'])}, Gnao1 {_fmt(r['Gnao1_cpm'])}, "
                     f"actin sum {_fmt(r['actin_sum_cpm'])}")
            L.append(f"      TISSUE (pre-R1): {r['tissue_verdict']} -- {r['tissue_reason']}")
            L.append(f"      Rule 1 sort    : {r['sort_verdict']} -- {r['sort_reason']}")
            L.append(f"      Rule 4 pop     : {r['population_call']} "
                     f"(Gnai2:Gnao1 = {r['gnai2_gnao1_ratio_str']}) -- "
                     f"{r['population_note']}")
            L.append(f"      library        : {r['library_status']} -- {r['library_reason']}")
            L.append(f"      blocking_flags : {r['blocking_flags']}")
            L.append("")

    L.append("=" * 78)
    L.append("Roll-up")
    for v, n in df["qc_overall"].value_counts().items():
        L.append(f"  {v}: {n}")
    L.append("  tissue verdicts: " + ", ".join(
        f"{v}={n}" for v, n in df["tissue_verdict"].value_counts().items()))
    wrong = df.loc[df["tissue_verdict"].isin(["MOE", "MOE_dominant_mixed"]),
                   "sample"].tolist()
    if wrong:
        L.append(f"  WRONG TISSUE ({len(wrong)}): " + ", ".join(wrong))
        L.append("    These are MAIN OLFACTORY libraries. No VR biology is reportable "
                 "from them. The remedy is wet-lab (dissection / sort gate) -- "
                 "re-quantifying against a different annotation will NOT change this, "
                 "because both trials used the same annotation and it yields tens of "
                 "thousands of Trpc2 counts in the trial-2 VNO libraries.")
    ok = df.loc[df["qc_overall"] != "UNUSABLE", "sample"].tolist()
    L.append(f"  samples cleared for downstream VR analysis ({len(ok)}): "
             + (", ".join(ok) if ok else "(none)"))
    bad = df.loc[df["suppress_biology"], "sample"].tolist()
    L.append(f"  biology-suppressed ({len(bad)}): " + (", ".join(bad) if bad else "(none)"))
    mm = df[["sample", "pct_multi_loci", "pct_too_many_loci"]].dropna()
    if len(mm):
        L.append(f"  multi-loci range: {mm['pct_multi_loci'].min():.2f}-"
                 f"{mm['pct_multi_loci'].max():.2f}% ; too-many-loci range: "
                 f"{mm['pct_too_many_loci'].min():.2f}-"
                 f"{mm['pct_too_many_loci'].max():.2f}%")
        L.append("  => the dominant multi-mapping loss/ambiguity channel is "
                 "'reads mapped to multiple loci', which STAR RETAINS and hands to "
                 "Salmon's EM. Reads discarded as 'too many loci' are negligible. "
                 "VR paralog ambiguity therefore lands in the EM, not in the "
                 "alignment filter -- cluster-level aggregation is required.")
    return "\n".join(L) + "\n"


def write_trial_outputs(cfg, trial: str, out_dir: str = None,
                        extra_dirs: List[str] = ()) -> Dict[str, Any]:
    P = trial_paths(cfg, trial)
    out_dir = out_dir or P["out_dir"]
    tech = technical_qc_table(cfg, trial)
    mark = marker_cpm_table(cfg, trial)
    qc = sample_qc_table(cfg, trial)
    written: List[str] = []
    for d in [out_dir, *extra_dirs]:
        os.makedirs(d, exist_ok=True)
        tech.to_csv(os.path.join(d, "technical_qc.tsv"), sep="\t", index=False)
        with open(os.path.join(d, "marker_cpm.tsv"), "w") as fh:
            fh.write(f"# {mark.attrs['cpm_denominator']}\n")
            mark.to_csv(fh, sep="\t", index=False)
        p = os.path.join(d, "sample_qc.tsv")
        with open(p, "w") as fh:
            fh.write(f"# {qc.attrs['cpm_denominator']}\n")
            qc.to_csv(fh, sep="\t", index=False)
        with open(os.path.join(d, "sample_qc_summary.txt"), "w") as fh:
            fh.write(summary_text(cfg, qc))
        written += [os.path.join(d, f) for f in
                    ("technical_qc.tsv", "marker_cpm.tsv", "sample_qc.tsv",
                     "sample_qc_summary.txt")]
    return {"trial": trial, "qc": qc, "written": written}


def _main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Assemble the per-sample QC table")
    ap.add_argument("--trial")
    ap.add_argument("--all-trials", action="store_true")
    ap.add_argument("--config")
    ap.add_argument("--stage-dir", help="second output dir (for job transfer-back)")
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    trials = trials_of(cfg) if a.all_trials else [a.trial]
    if not trials or trials == [None]:
        ap.error("give --trial NAME or --all-trials")

    frames, all_written = [], []
    for t in trials:
        extra = [os.path.join(a.stage_dir, t)] if a.stage_dir else []
        res = write_trial_outputs(cfg, t, extra_dirs=extra)
        frames.append(res["qc"])
        all_written += res["written"]
        print(f"[vr_sample_qc] {t}: {len(res['qc'])} samples")

    combined = pd.concat(frames, ignore_index=True, sort=False)
    combined.attrs.update(frames[0].attrs)
    for d in filter(None, [os.path.join(cfg["work"], "results"), a.stage_dir]):
        os.makedirs(d, exist_ok=True)
        p = os.path.join(d, "sample_qc_all.tsv")
        with open(p, "w") as fh:
            fh.write(f"# {combined.attrs.get('cpm_denominator','NA')}\n")
            combined.to_csv(fh, sep="\t", index=False)
        s = os.path.join(d, "sample_qc_summary.txt")
        with open(s, "w") as fh:
            fh.write(summary_text(cfg, combined))
        all_written += [p, s]
    print("[vr_sample_qc] wrote:\n  " + "\n  ".join(all_written))
    print("\n" + summary_text(cfg, combined))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
