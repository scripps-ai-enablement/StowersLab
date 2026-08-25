#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""vr_quantify.py -- cluster-level VR quantification + EM-artifact detection.

Part of the VNO receptor pipeline (Stowers Lab / Natalie Cole).
Reads every path, sample name and threshold from config/project.yaml via
vr_config; reuses the settled CPM convention from vr_markers (unscaled
salmon.merged.gene_counts.tsv, all-gene column-sum denominator).

WHY THIS MODULE EXISTS
----------------------
~250 V1R and ~120 V2R genes sit in genomic clusters of local duplicates at
85-95% nucleotide identity. 75bp reads cannot be uniquely assigned within a
cluster, so Salmon's EM distributes ambiguous reads across paralogs and
per-gene VR counts inside a cluster are not trustworthy. Cluster-level
aggregation is the reliable readout; individual-paralog calls require
evidence that does not come out of the EM.

=====================================================================
THE REDISTRIBUTION STATISTIC -- DEFINITION
=====================================================================
For one (sample, cluster) we take the member paralogs' EM counts
x_1..x_k, keep the members with x_i >= `detect_floor` (default 1.0
count) as the "detected" set (size k_det), and let N = round(sum x_i)
be the read support.

Descriptive quantities (whole cluster):
  frac_i    = x_i / sum(x)                     within-cluster fractions
  top_frac  = max frac_i
  evenness  = H / ln(k_det),  H = -sum frac_i ln frac_i   (0..1; 1 = uniform)
  hill      = exp(H)          "effective number of paralogs"

Inferential quantity. Uniformity of the fractions is tested against the
1/k expectation with a Pearson statistic
  X2 = sum_i (x_i - N/k)^2 / (N/k)
whose null distribution is obtained by MONTE CARLO from
Multinomial(N, 1/k) at the OBSERVED N and k -- not from the asymptotic
chi-square, because at N of a few tens the asymptotic approximation is
invalid and this is exactly the regime the flag has to survive.
  p_uniform = P(X2_null >= X2_observed)
A LARGE p_uniform means the observed split is indistinguishable from an
even split, i.e. consistent with one paralog having been expressed and
its reads redistributed.

READ-DEPTH ACCOUNTING (the point of the exercise)
A near-equal split of 30 reads across 3 paralogs is much weaker evidence
than a near-equal split of 3000, because at N=30 multinomial noise alone
swamps the difference between "even" and "monogenic". We therefore do NOT
report p_uniform on its own. For each (N, k) we also compute, by Monte
Carlo under the alternative
  p_alt = (dom_alt, (1-dom_alt)/(k-1), ...)   default dom_alt = 0.90
(i.e. what a genuinely monogenic cell would look like), the POWER of the
same test to reject uniformity at alpha:
  power_vs_dominant = P(reject uniform | monogenic truth, this N and k)
When power < `min_power` (default 0.80) the (sample, cluster) is graded
`indeterminate_low_depth` REGARDLESS of p_uniform: failure to reject
uniformity at low depth is absence of evidence, not evidence of an
artifact. This is the depth correction.

EVEN-BLOCK REFINEMENT
EM redistribution does not have to touch every member; it splits reads
among the paralogs that actually share sequence. So beyond the
whole-cluster test we find the largest m in 2..k_det such that the top-m
paralogs (a) are jointly consistent with uniform (p >= alpha) and (b)
carry at least `block_min_share` (default 0.80) of the cluster signal.
That m is `even_block_size` -- the number of paralogs among which a
single expressed transcript's reads appear to have been split. It
recovers the 45/45/10 and 33/33/33 patterns alike, and its own N and m
drive the power calculation for the graded flag.

UNIQUE-READ GATE (necessary vs sufficient)
An even fraction split is NECESSARY but NOT SUFFICIENT for an EM artifact.
Two neurons in a pool expressing two paralogs of the same cluster produce
exactly the same even split as one transcript whose reads were divided
between two sequence-similar paralogs. The uniformity test cannot tell
these apart, because both live in the fractions. What separates them is
EM-INDEPENDENT support: if each co-dominant member carries substantial
unique reads, each was independently observed and the split is real
co-expression; if all but one member sit at ~zero unique reads, one gene's
reads were redistributed.

So once an even block is found, the block members' unique-read counts
(`unique_reads_bam_nodup` by preference -- MAPQ-255 non-duplicate, hence
single-locus; falling back to `unique_reads_bam` then the Salmon
transcript-unique channel) are compared against
    threshold = max(unique_gate_min_reads, unique_gate_bg_ratio * floor)
where `floor` is the MEDIAN unique count of the cluster members OUTSIDE
the even block -- the cluster's own mismapping/ambient level, measured
rather than assumed. Defaults: an absolute floor of 10 reads and 3x the
background. The absolute floor stops a handful of stray reads counting as
independent observation in a quiet cluster; the background ratio stops
cluster-wide mismapping counting as it in a noisy one. A member clearing
BOTH is "unique-supported".

  >= 2 supported members -> `even_split_unique_supported` (co-expression)
  <= 1 supported member  -> `suspected_em_redistribution`
  no unique channel at all -> `even_split_unique_unavailable` (unresolvable;
      the two hypotheses cannot be separated, so no artifact is asserted)

WHY THIS MATTERS BEYOND THE PRESENT SAMPLES. Monogenic receptor choice is
a per-CELL rule, not a per-pool one. Any pooled library can legitimately
contain two paralogs of one cluster, and the planned stimulus-response
design pools responsive cells, so within-cluster co-expression will be
common. Without this gate the uniformity test would label that real
biology an artifact whenever the two paralogs happened to land near 50/50.
The uniformity statistic and its power are still reported in every case --
the gate changes the LABEL, not the measurement.

GRADED FLAG (`em_flag` / `em_flag_level`), never a bare boolean:
  insufficient_signal          N < min_reads; nothing to say
  single_paralog_only          k_det < 2; redistribution not possible
  no_redistribution_signature  no even block found; fractions are
                               structured, consistent with a real call
  indeterminate_low_depth      an even block exists but the test has no
                               power at this depth -- unresolvable
  even_split_unique_unavailable  even block but no unique-read channel;
                               redistribution vs co-expression unresolvable
  even_split_unique_supported  even block, but >=2 members carry independent
                               unique support -> co-expression, NOT an artifact
  suspected_em_redistribution  even block + adequate power + <=1 member with
                               independent unique support
                               level strong  : p >= 0.20, evenness >= 0.90,
                                 AND block read support >= strong_min_reads
                                 (default 200) -- an even split of a few tens
                                 of reads is graded moderate even when the
                                 test has formal power, because thin evidence
                                 should not read as a confident artifact call
                               level moderate: otherwise

=====================================================================
LIMITATIONS -- read before using the flag
=====================================================================
1. EM counts are not independent read draws. Salmon's output is smoother
   than a multinomial sample of the same total, so the multinomial null is
   too wide and p_uniform is ANTI-CONSERVATIVE: the test over-calls
   uniformity. Treat `suspected_em_redistribution` as a screen, not proof.
2. Near-uniform fractions are NOT by themselves diagnostic of an artifact,
   in a pooled library or anywhere else -- see UNIQUE-READ GATE above, which
   is what discriminates redistribution from co-expression. Two residual
   caveats survive the gate. (a) The gate inherits the limits of its own
   channel: MAPQ-255 counting is span-based, so a paralog pair whose
   ENTIRE length is unmappable at 75bp will show ~zero unique reads for
   BOTH members and land in `suspected_em_redistribution` whether or not
   both were truly expressed; such a pair is genuinely unresolvable with
   these data and the flag should be read as "not separable", not as proof
   of splitting. (b) Column `interpretation_context` (monogenic_expectation
   for <=2 cells vs pooled_ambiguous) still matters for how much weight a
   flag carries: in a 2-cell library an even split contradicts the
   monogenic prior and is informative, whereas in a 100-cell pool it is
   unremarkable and the flag rests entirely on the unique-read evidence.
3. Uniform is a caricature of EM behaviour. EM splits in proportion to
   unique-read anchors and effective lengths, so a genuine single-source
   split can land at 65/35 and will NOT be flagged. False negatives are
   expected; the flag has better specificity than sensitivity.
4. Fractional EM counts are rounded to integers to define N.
5. No multiple-testing correction across (sample, cluster) pairs. The
   flags are descriptive annotations, not a controlled family of tests.
6. Monte-Carlo p-values are granular at `n_mc` (default 4000) draws and
   floored at 1/n_mc.

=====================================================================
CANDIDATE EVIDENCE (Step 4) -- what is and is not EM-independent
=====================================================================
Two EM-independent read channels are used, in this order of preference:

  bam_unique_mapq255  -- alignments over the gene span in
      <sample>.markdup.sorted.bam with MAPQ 255, which for STAR means
      "placed at exactly one locus in the genome". Gene-level and
      genuinely unambiguous. Counted mate-wise (each mate of a pair is
      one alignment) and reported both including and excluding
      duplicate-flagged reads. Caveats: reads are assigned by alignment
      start position falling in the gene span, so a read straddling the
      span boundary is missed; where two annotated VR spans overlap a
      read is credited to both (column `span_overlaps_other_vr`).

  salmon_tx_unique   -- UniqueCount from <sample>/aux_info/ambig_info.tsv,
      row-aligned with quant.sf, summed to gene level via tx2gene.
      IMPORTANT limitation: "unique" there means unique to a TRANSCRIPT,
      so for a multi-transcript gene reads shared between its own
      isoforms fall into AmbigCount and gene-level unique support is
      UNDERSTATED. Column `n_transcripts` is carried so this is auditable
      (510 of 538 VR genes are single-transcript, so it mostly bites the
      other 28).

Two guards apply to the ranking itself:

  BACKGROUND FLOOR. In a large cluster most members carry a small,
  near-constant unique count from mismapping and ambient reads, so the
  median of the NON-TOP members is the empirical noise floor, not zero.
  `top_over_background_floor` is the top gene's unique count over that
  median, and a top gene that does not clear the floor by the same ratio
  required of rank1/rank2 is graded unresolvable. The top gene is
  excluded from its own floor, and a single-member cluster has no floor
  to clear (reported inf, test not applicable).

  EM/UNIQUE CONTRADICTION. The two channels should agree on which
  paralog carries the signal. Where the top unique-read gene holds
  almost no EM signal AND the EM-dominant paralog has essentially no
  unique reads, they name DIFFERENT genes and neither identifies the
  receptor; `evidence_contradiction` is populated, confidence is forced
  to unresolvable, and the text is prefixed to `notes`. This is observed
  in V1R_chr7_cl015 of target2cellsRep1_S3 (Vmn1r131 holds 99.9% of EM
  signal with 0 unique reads; Vmn1r103/Vmn1r104 have 92 each with 0 EM
  counts) and must be reported, not smoothed over.

If both channels are empty for a cluster the paralogs CANNOT be separated
with this data and the module says so (`evidence_type=none_available`,
`confidence=unresolvable`) rather than ranking by EM-distributed counts.
Every candidate row carries confirmation_status=tentative_unconfirmed
unless orthogonal evidence is supplied from outside this module.

=====================================================================
PSEUDOGENE BLEED -- framing (non-negotiable)
=====================================================================
Apparent pseudogene expression inside a functional cluster has (at least)
two unresolved candidate mechanisms:
  (a) multi-mapping / EM leakage from an expressed functional paralog, and
  (b) genuine transcription of the pseudogene locus, e.g. driven by
      regulatory elements shared across the cluster.
This module reports BOTH as unresolved and does not choose. Dietschi et
al. 2022 (Sci Adv 8(46) eabn7450) is NOT support for a quantitative mouse
expectation here: their pseudogene-expression result was significant in
RAT only (P=0.003) and NOT in mouse (W=1214, P=0.5704), and the mechanism
they propose is regulatory, not multi-mapping. Do not cite them as a
prior for how much bleed to expect in these libraries.

CLI
---
  python vr_quantify.py --trial trial2 [--outdir DIR] [--no-bam]
  python vr_quantify.py --all-trials
"""
from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import vr_config  # noqa: E402
import vr_markers  # noqa: E402

# ---------------------------------------------------------------- defaults
# Every one of these can be overridden from config under `vr_quant:`.
DEFAULTS: Dict[str, Any] = {
    "detect_floor_counts": 1.0,      # member counts >= this are "detected"
    "min_reads_for_em_test": 20,     # below this: insufficient_signal
    "alpha": 0.05,
    "dom_alt": 0.90,                 # monogenic alternative for power
    "min_power": 0.80,
    "strong_min_reads": 200,         # below this an even split grades moderate
    # unique-read gate on the redistribution flag (see module docstring)
    "unique_gate_min_reads": 10,     # absolute floor for "substantial" support
    "unique_gate_bg_ratio": 3.0,     # and this multiple of the cluster background
    "block_min_share": 0.80,
    "n_mc": 4000,
    "mc_seed": 20260820,
    "cluster_call_min_cpm": 10.0,    # a cluster is "called" if ...
    "cluster_call_min_share": 0.02,  # ... both of these hold
    "candidate_cluster_min_share": 0.10,
    "min_unique_reads_for_rank": 5,
    "unique_ratio_moderate": 3.0,
    "unique_ratio_low": 2.0,
    "unique_reads_moderate": 10,
    "max_candidates_reported": 3,
    "pseudogene_bleed_min_cpm": 1.0,
}

# Fixed candidate schema. Declared so the file always carries a header row even
# when no sample is cleared (trial1): a header-less empty TSV breaks downstream
# readers, and "no candidates" must be machine-readable as an empty table
# rather than as a parse error.
CANDIDATE_COLUMNS: List[str] = [
    "trial", "sample", "cluster_id", "supercluster_id", "family",
    "is_dominant_cluster", "cluster_share_of_vr", "cluster_cpm", "rank",
    "gene_id", "gene_name", "candidate_list", "n_candidates_reported",
    "evidence_type", "unique_reads_bam", "unique_reads_bam_nodup",
    "unique_reads_salmon_tx", "unique_share_of_cluster", "em_counts",
    "em_frac_of_cluster", "em_flag", "even_block_size",
    "unique_reads_cluster_median", "top_over_background_floor",
    "evidence_contradiction", "confidence",
    "confirmation_status", "marker_consistency", "population_call", "notes",
]


def qcfg(cfg: Mapping[str, Any], key: str) -> Any:
    """Threshold lookup: config `vr_quant:` block first, then DEFAULTS."""
    block = cfg.get("vr_quant") or {}
    if isinstance(block, Mapping) and key in block:
        return block[key]
    return DEFAULTS[key]


# ---------------------------------------------------------------- io helpers
def read_tsv_skipping_comments(path: str) -> pd.DataFrame:
    """Read a TSV whose leading lines may be '#' provenance comments.

    Cannot use pandas comment='#' -- that would also truncate any field
    containing a '#'. Only *leading* comment lines are skipped.
    """
    n_skip = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                n_skip += 1
            else:
                break
    return pd.read_csv(path, sep="\t", skiprows=n_skip, dtype=str)


def _as_bool(value: Any) -> bool:
    return str(value).strip().lower() in ("true", "1", "yes", "t")


def load_vr_reference(cfg: Mapping[str, Any]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return (gene_to_cluster, gene_annotation), primary-assembly only.

    JOIN ON gene_id, NEVER gene_name: Vmn1r-ps5 and Vmn2r118 each map to
    two distinct gene_ids on different chromosomes.
    """
    ref_dir = os.path.join(cfg["work"], "ref")
    g2c = pd.read_csv(os.path.join(ref_dir, "vr_gene_to_cluster.tsv"), sep="\t")
    ann = pd.read_csv(os.path.join(ref_dir, "vr_gene_annotation.tsv"), sep="\t")
    for frame, name in ((g2c, "vr_gene_to_cluster.tsv"), (ann, "vr_gene_annotation.tsv")):
        if "is_primary_assembly" not in frame.columns:
            raise KeyError(f"{name} lacks is_primary_assembly")
    # Columns wanted from the annotation table. Any that ALSO exist in
    # gene_to_cluster are dropped from the annotation side before merging --
    # otherwise pandas silently suffixes them (is_pseudogene_x/_y) and every
    # downstream `is_pseudogene` lookup misses, which silently disables the
    # pseudogene-bleed check rather than erroring.
    want = ["gene_id", "gene_name", "n_transcripts", "is_pseudogene",
            "pseudogene_criterion", "gene_biotype", "chrom", "start", "end"]
    want = [c for c in want if c in ann.columns]
    ann_slim = ann.loc[ann["is_primary_assembly"] == 1, want].copy()
    g2c_p = g2c.loc[g2c["is_primary_assembly"] == 1].copy()
    collide = [c for c in ann_slim.columns
               if c != "gene_id" and c in g2c_p.columns]
    merged = g2c_p.merge(ann_slim.drop(columns=collide),
                         on="gene_id", how="left", validate="one_to_one")
    for col in ("is_pseudogene", "pseudogene_criterion"):
        if col not in merged.columns:
            raise KeyError(
                f"{col} missing after reference merge -- the pseudogene-bleed "
                "check cannot run without it"
            )
    if any(c.endswith(("_x", "_y")) for c in merged.columns):
        raise AssertionError(
            "suffixed columns after reference merge: "
            f"{[c for c in merged.columns if c.endswith(('_x', '_y'))]}"
        )
    if merged["gene_id"].duplicated().any():
        raise ValueError("duplicate gene_id in primary-assembly gene_to_cluster")
    return merged, ann


def pseudogene_confidence(criterion: Any) -> str:
    """Grade pseudogene status from `pseudogene_criterion`, not `is_pseudogene`.

    The reference flags 538 genes with is_pseudogene, but the EVIDENCE differs:
    `biotype+name` means the Ensembl biotype and the -ps name agree (209
    genes); `name` alone or `biotype` alone means they disagree (5 genes
    total), and for those the pseudogene call is itself uncertain. Any bleed
    flag on a `conflicting_evidence` gene is doubly uncertain -- the apparent
    expression may simply be a functional gene that is misnamed.
    """
    text = str(criterion).strip().lower()
    if text in ("none", "nan", ""):
        return "not_pseudogene"
    if text == "biotype+name":
        return "concordant_evidence"
    if text in ("name", "biotype"):
        return "conflicting_evidence"
    return f"unrecognised_criterion({criterion})"


def load_counts_and_cpm(cfg: Mapping[str, Any], trial: str
                        ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Counts (gene_id-indexed), CPM, and the gene_id->gene_name map.

    WHY NOT vr_markers.load_gene_counts: that loader indexes by gene_name and
    SUMS duplicate names, which is right for marker panels but wrong here --
    Vmn1r-ps5 and Vmn2r118 each map to two distinct gene_ids on different
    chromosomes, so a name-collapsed table would fuse two loci in different
    clusters into one row. We therefore index by gene_id.

    The CPM convention is unchanged and that is asserted, not assumed: the
    all-gene column-sum denominator is invariant to whether rows are collapsed
    by name, so the library totals here must equal vr_markers' totals to
    floating-point tolerance. Mismatch raises.
    """
    paths = vr_config.trial_paths(cfg, trial)
    raw = pd.read_csv(paths["gene_counts"], sep="\t")
    for col in ("gene_id", "gene_name"):
        if col not in raw.columns:
            raise KeyError(f"{paths['gene_counts']}: missing {col} column")
    if raw["gene_id"].duplicated().any():
        dups = raw.loc[raw["gene_id"].duplicated(), "gene_id"].head(5).tolist()
        raise AssertionError(f"duplicate gene_id rows in quant table: {dups}")
    sample_cols = [c for c in raw.columns if c not in ("gene_id", "gene_name")]
    counts = (raw.set_index("gene_id")[sample_cols]
                 .apply(pd.to_numeric, errors="coerce").fillna(0.0))
    name_map = raw.set_index("gene_id")[["gene_name"]].astype(str)

    totals = counts.sum(axis=0)
    if (totals <= 0).any():
        bad = list(totals[totals <= 0].index)
        raise ValueError(f"zero total counts for samples {bad}; cannot compute CPM")
    ref_totals = vr_markers.load_gene_counts(paths["gene_counts"]).sum(axis=0)
    for s in totals.index:
        a, b = float(totals[s]), float(ref_totals[s])
        if not math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-6):
            raise AssertionError(
                "CPM denominator diverges from the settled convention for "
                f"sample {s}: gene_id-indexed total {a} vs vr_markers total {b}"
            )
    cpm = counts.divide(totals, axis=1) * 1e6
    cpm.attrs["library_totals"] = totals.to_dict()
    return counts, cpm, name_map


def assert_join_complete(counts: pd.DataFrame, name_map: pd.DataFrame,
                         g2c: pd.DataFrame, ann_full: pd.DataFrame) -> Dict[str, Any]:
    """Fail loudly on any incomplete VR join. Returns a reconciliation record."""
    quant_ids = set(counts.index)
    names = name_map["gene_name"].astype(str)
    name_vr_ids = set(names.index[names.str.startswith(("Vmn1r", "Vmn2r"))])
    ref_primary_ids = set(g2c["gene_id"])
    ref_all_ids = set(ann_full["gene_id"])

    missing_members = sorted(ref_primary_ids - quant_ids)
    if missing_members:
        raise AssertionError(
            "cluster member gene_ids absent from the quant table (would be "
            f"silently dropped): n={len(missing_members)} e.g. {missing_members[:10]}"
        )
    unmatched_quant = sorted(name_vr_ids - ref_primary_ids)
    non_primary = sorted(ref_all_ids - ref_primary_ids)
    unexplained = sorted(set(unmatched_quant) - set(non_primary))
    if unexplained:
        raise AssertionError(
            "VR-named quant rows with no primary-assembly cluster assignment and "
            f"no non-primary explanation: n={len(unexplained)} e.g. {unexplained[:10]}"
        )
    rec = {
        "n_vr_named_quant_rows": len(name_vr_ids),
        "n_ref_primary_genes": len(ref_primary_ids),
        "n_ref_all_genes": len(ref_all_ids),
        "n_non_primary_excluded": len(non_primary),
        "non_primary_gene_ids": non_primary,
        "n_quant_rows_matched": len(name_vr_ids & ref_primary_ids),
        "join_complete": True,
    }
    if rec["n_quant_rows_matched"] != len(ref_primary_ids):
        raise AssertionError(
            "join is not one-to-one: matched "
            f"{rec['n_quant_rows_matched']} of {len(ref_primary_ids)} reference genes"
        )
    return rec


# ---------------------------------------------------------------- step 1
def cluster_expression(counts: pd.DataFrame, cpm: pd.DataFrame, g2c: pd.DataFrame,
                       samples: Sequence[str], tier_col: str,
                       detect_floor: float) -> pd.DataFrame:
    """Aggregate VR counts/CPM to a cluster tier for every sample (long format)."""
    ids = g2c["gene_id"].tolist()
    sub_counts = counts.loc[ids, list(samples)]
    sub_cpm = cpm.loc[ids, list(samples)]
    tier = g2c.set_index("gene_id")[tier_col]
    fam = g2c.set_index("gene_id")["family"]
    chrom = g2c.set_index("gene_id")["chrom"]

    grp = tier.loc[ids]
    agg_counts = sub_counts.groupby(grp).sum()
    agg_cpm = sub_cpm.groupby(grp).sum()
    n_detected = (sub_counts >= detect_floor).groupby(grp).sum()
    n_members = grp.value_counts()
    fam_by_tier = fam.loc[ids].groupby(grp).first()
    chrom_by_tier = chrom.loc[ids].groupby(grp).first()

    total_counts = sub_counts.sum(axis=0)
    total_cpm = sub_cpm.sum(axis=0)

    rows = []
    for sample in samples:
        for cid in agg_counts.index:
            c = float(agg_counts.at[cid, sample])
            m = float(agg_cpm.at[cid, sample])
            rows.append({
                "sample": sample,
                "tier": "supercluster" if tier_col == "supercluster_id" else "cluster",
                "cluster_id": cid,
                "family": fam_by_tier.get(cid),
                "chrom": chrom_by_tier.get(cid),
                "n_member_genes": int(n_members.get(cid, 0)),
                "n_member_genes_detected": int(n_detected.at[cid, sample]),
                "counts_sum": c,
                "cpm_sum": m,
                "share_of_sample_vr": (c / float(total_counts[sample])
                                       if float(total_counts[sample]) > 0 else np.nan),
                "sample_total_vr_counts": float(total_counts[sample]),
                "sample_total_vr_cpm": float(total_cpm[sample]),
            })
    return pd.DataFrame(rows)


def within_cluster_fractions(counts: pd.DataFrame, cpm: pd.DataFrame,
                             g2c: pd.DataFrame, samples: Sequence[str]) -> pd.DataFrame:
    """Per (sample, cluster, member gene) counts, CPM and fraction of cluster."""
    ids = g2c["gene_id"].tolist()
    meta = g2c.set_index("gene_id")
    sub_counts = counts.loc[ids, list(samples)]
    sub_cpm = cpm.loc[ids, list(samples)]
    cl = meta["cluster_id"]

    frames = []
    for sample in samples:
        d = pd.DataFrame({
            "gene_id": ids,
            "counts": sub_counts[sample].to_numpy(dtype=float),
            "cpm": sub_cpm[sample].to_numpy(dtype=float),
        })
        d["cluster_id"] = cl.loc[ids].to_numpy()
        tot = d.groupby("cluster_id")["counts"].transform("sum")
        d["cluster_counts_sum"] = tot
        d["frac_of_cluster"] = np.where(tot > 0, d["counts"] / tot, np.nan)
        d["sample"] = sample
        frames.append(d)
    out = pd.concat(frames, ignore_index=True)
    carry = ["gene_name", "family", "chrom", "start", "end", "supercluster_id",
             "is_pseudogene", "pseudogene_criterion", "n_transcripts",
             "cluster_n_genes"]
    carry = [c for c in carry if c in meta.columns]
    out = out.merge(meta[carry].reset_index(), on="gene_id", how="left")
    out["rank_in_cluster"] = (out.groupby(["sample", "cluster_id"])["counts"]
                                .rank(ascending=False, method="first").astype(int))
    return out


# ---------------------------------------------------------------- step 2
class _MCCache:
    """Monte-Carlo null/alternative distributions of X2, cached per (N, k)."""

    def __init__(self, n_mc: int, seed: int, dom_alt: float, alpha: float):
        self.n_mc = int(n_mc)
        self.dom_alt = float(dom_alt)
        self.alpha = float(alpha)
        self.rng = np.random.default_rng(seed)
        self._null: Dict[Tuple[int, int], np.ndarray] = {}
        self._power: Dict[Tuple[int, int], float] = {}

    @staticmethod
    def _x2(obs: np.ndarray, n: int, k: int) -> np.ndarray:
        exp = n / k
        return ((obs - exp) ** 2 / exp).sum(axis=-1)

    def null(self, n: int, k: int) -> np.ndarray:
        key = (int(n), int(k))
        if key not in self._null:
            draws = self.rng.multinomial(n, [1.0 / k] * k, size=self.n_mc)
            self._null[key] = np.sort(self._x2(draws, n, k))
        return self._null[key]

    def p_uniform(self, obs_counts: Sequence[float]) -> Tuple[float, float, int, int]:
        obs = np.asarray(obs_counts, dtype=float)
        k = int(obs.size)
        n = int(round(float(obs.sum())))
        if k < 2 or n <= 0:
            return (np.nan, np.nan, n, k)
        x2 = float(self._x2(obs, n, k))
        null = self.null(n, k)
        p = float((null >= x2).sum() + 1) / float(self.n_mc + 1)
        return (p, x2, n, k)

    def power(self, n: int, k: int) -> float:
        """P(reject uniform at alpha | one paralog holds dom_alt of the signal)."""
        key = (int(n), int(k))
        if key in self._power:
            return self._power[key]
        if k < 2 or n <= 0:
            self._power[key] = np.nan
            return np.nan
        null = self.null(n, k)
        crit = float(np.quantile(null, 1.0 - self.alpha))
        rest = (1.0 - self.dom_alt) / (k - 1)
        p_alt = [self.dom_alt] + [rest] * (k - 1)
        draws = self.rng.multinomial(n, p_alt, size=self.n_mc)
        x2 = self._x2(draws, n, k)
        val = float((x2 > crit).mean())
        self._power[key] = val
        return val


def _entropy_stats(fracs: np.ndarray) -> Tuple[float, float]:
    f = fracs[fracs > 0]
    if f.size == 0:
        return (np.nan, np.nan)
    h = float(-(f * np.log(f)).sum())
    k = int(f.size)
    evenness = h / math.log(k) if k > 1 else 0.0
    return (evenness, float(math.exp(h)))


def em_redistribution_stats(member_counts: Sequence[float], mc: _MCCache,
                            detect_floor: float, min_reads: int,
                            min_power: float, block_min_share: float,
                            alpha: float, strong_min_reads: int = 200,
                            member_unique: Optional[Sequence[float]] = None,
                            unique_min_reads: int = 10,
                            unique_bg_ratio: float = 3.0
                            ) -> Dict[str, Any]:
    """Redistribution statistics for one (sample, cluster). See module docstring.

    `member_unique` must be positionally aligned with `member_counts` (EM-
    independent unique-read support per member, normally MAPQ-255 non-duplicate
    BAM counts). It gates the redistribution flag: an even fraction split is
    necessary but NOT sufficient for an artifact, because in a multi-cell pool
    two neurons expressing two paralogs of one cluster produce the same even
    split. See UNIQUE-READ GATE in the module docstring.
    """
    x_all = np.asarray(member_counts, dtype=float)
    u_all = (np.full(x_all.shape, np.nan) if member_unique is None
             else np.asarray(member_unique, dtype=float))
    if u_all.shape != x_all.shape:
        raise ValueError("member_unique must be positionally aligned with member_counts")
    keep = np.isfinite(x_all)
    x_all, u_all = x_all[keep], u_all[keep]
    total = float(x_all.sum())
    # Sort jointly so unique support stays paired with its own paralog.
    order = np.argsort(-x_all, kind="stable")
    x_s, u_s = x_all[order], u_all[order]
    sel = x_s >= detect_floor
    det, det_u = x_s[sel], u_s[sel]
    bg_u = u_s[~sel]                     # undetected members = background channel
    k_det = int(det.size)
    n_reads = int(round(total))
    out: Dict[str, Any] = {
        "k_members": int(x_all.size),
        "k_detected": k_det,
        "read_support": n_reads,
        "top_frac": float(det[0] / total) if k_det >= 1 and total > 0 else np.nan,
        "second_frac": float(det[1] / total) if k_det >= 2 and total > 0 else np.nan,
        "top_over_second": (float(det[0] / det[1]) if k_det >= 2 and det[1] > 0 else np.nan),
        "evenness": np.nan, "hill_effective_paralogs": np.nan,
        "chi2_uniform": np.nan, "p_uniform": np.nan,
        "power_vs_dominant": np.nan,
        "even_block_size": 0, "even_block_share": np.nan,
        "even_block_p_uniform": np.nan, "even_block_read_support": np.nan,
        "even_block_power": np.nan, "even_block_evenness": np.nan,
        "em_flag": "insufficient_signal", "em_flag_level": "none",
        # unique-read gate
        "unique_channel_available": int(bool(np.isfinite(u_all).any())),
        "unique_background_floor": np.nan,
        "unique_support_threshold": np.nan,
        "block_unique_reads": "",
        "n_block_members_unique_supported": np.nan,
        "block_unique_total": np.nan,
        "interpretation": "",
    }
    if total <= 0 or k_det == 0:
        out["em_flag"] = "no_signal"
        out["interpretation"] = ("no VR signal in this cluster for this sample; "
                                 "no redistribution question arises")
        return out
    fr = det / total
    out["evenness"], out["hill_effective_paralogs"] = _entropy_stats(fr)
    if n_reads < int(min_reads):
        out["em_flag"] = "insufficient_signal"
        out["interpretation"] = (
            f"only {n_reads} reads in this cluster (< {int(min_reads)}); too few to "
            "distinguish any within-cluster pattern from sampling noise")
        return out
    if k_det < 2:
        out["em_flag"] = "single_paralog_only"
        out["interpretation"] = (
            "a single member is detected, so redistribution among paralogs is not "
            "possible; this does NOT by itself confirm the paralog's identity")
        return out

    p_all, x2_all, _, _ = mc.p_uniform(det)
    out["chi2_uniform"] = x2_all
    out["p_uniform"] = p_all
    out["power_vs_dominant"] = mc.power(n_reads, k_det)

    best = None
    best_m = 0
    for m in range(2, k_det + 1):
        block = det[:m]
        share = float(block.sum() / total)
        if share < float(block_min_share):
            continue
        n_block = int(round(float(block.sum())))
        if n_block < int(min_reads):
            continue
        p_b, _, _, _ = mc.p_uniform(block)
        if p_b is not np.nan and p_b >= float(alpha):
            ev_b, _ = _entropy_stats(block / block.sum())
            best = {
                "even_block_size": m, "even_block_share": share,
                "even_block_p_uniform": p_b, "even_block_read_support": n_block,
                "even_block_power": mc.power(n_block, m),
                "even_block_evenness": ev_b,
            }
            best_m = m
    if best is None:
        out["em_flag"] = "no_redistribution_signature"
        out["em_flag_level"] = "none"
        out["interpretation"] = (
            "within-cluster fractions are structured (distinguishable from an even "
            "split); consistent with a real per-paralog call rather than EM splitting"
        )
        return out
    out.update(best)
    if not np.isfinite(best["even_block_power"]) or best["even_block_power"] < float(min_power):
        out["em_flag"] = "indeterminate_low_depth"
        out["em_flag_level"] = "unresolvable"
        out["interpretation"] = (
            "an even block exists but the uniformity test has no power at this read "
            "depth; absence of evidence, not evidence of an artifact"
        )
        return out

    # ---- UNIQUE-READ GATE ------------------------------------------------
    # An even fraction split is NECESSARY but NOT SUFFICIENT for an EM
    # artifact. Two neurons in a pool expressing two paralogs of one cluster
    # produce the same even split as one transcript whose reads were divided
    # in two. The discriminator is EM-INDEPENDENT support: if the co-dominant
    # members each carry substantial unique reads, they are independently
    # observed and the even split is co-expression, not redistribution.
    blk_u = det_u[:best_m]
    out["block_unique_reads"] = ",".join(
        ("NA" if not np.isfinite(v) else f"{v:.0f}") for v in blk_u)
    if not np.isfinite(blk_u).any():
        # No EM-independent channel for this cluster: the two hypotheses
        # cannot be separated. Do NOT assert an artifact on fractions alone.
        out["em_flag"] = "even_split_unique_unavailable"
        out["em_flag_level"] = "unresolvable"
        out["interpretation"] = (
            "fractions are indistinguishable from an even split, but no unique-read "
            "evidence is available for these members, so EM redistribution cannot be "
            "distinguished from genuine co-expression of two paralogs"
        )
        return out

    blk_u = np.nan_to_num(blk_u, nan=0.0)
    # Background floor from members that are NOT in the even block: their
    # unique counts are the cluster's own mismapping/ambient level.
    others = np.concatenate([det_u[best_m:], bg_u]) if (det_u.size + bg_u.size) else np.array([])
    others = others[np.isfinite(others)]
    floor_u = float(np.median(others)) if others.size else 0.0
    thresh = max(float(unique_min_reads), float(unique_bg_ratio) * floor_u)
    supported = int((blk_u >= thresh).sum())
    out["unique_background_floor"] = floor_u
    out["unique_support_threshold"] = thresh
    out["n_block_members_unique_supported"] = supported
    out["block_unique_total"] = float(blk_u.sum())

    if supported >= 2:
        out["em_flag"] = "even_split_unique_supported"
        out["em_flag_level"] = "none"
        out["interpretation"] = (
            f"co-expression, not redistribution: {supported} of {best_m} co-dominant "
            f"paralogs each carry >= {thresh:.0f} EM-independent unique reads "
            f"({out['block_unique_reads']}), so each is observed on its own evidence. "
            "In a multi-cell pool this is the expected result of different cells "
            "choosing different paralogs of one cluster; monogenic choice is a "
            "per-CELL rule, not a per-pool one."
        )
        return out

    out["em_flag"] = "suspected_em_redistribution"
    # `strong` requires depth as well as evenness: power says the test COULD
    # have seen a monogenic pattern, but an even split of a few tens of reads
    # remains thin evidence, so it is graded moderate regardless of power.
    strong = (best["even_block_p_uniform"] >= 0.20
              and np.isfinite(best["even_block_evenness"])
              and best["even_block_evenness"] >= 0.90
              and best["even_block_read_support"] >= int(strong_min_reads))
    out["em_flag_level"] = "strong" if strong else "moderate"
    out["interpretation"] = (
        f"EM redistribution suspected: fractions are indistinguishable from an even "
        f"split across {best_m} paralogs, yet only {supported} of them clears the "
        f"unique-read threshold ({thresh:.0f}; observed {out['block_unique_reads']}). "
        "One transcript's reads appear to have been divided among sequence-similar "
        "paralogs. NOTE the asymmetry: the CALL of redistribution rests on the "
        "near-absence of unique support for the other member(s), which is solid, but "
        "WHICH paralog is the true source rests on the small unique count of the one "
        "supported member and is weakly determined."
    )
    return out


# ---------------------------------------------------------------- step 4 evidence
def load_ambig_info(sample_dir: str, tx2gene: pd.DataFrame,
                    gene_ids: Sequence[str]) -> Optional[pd.DataFrame]:
    """Gene-level UniqueCount/AmbigCount from aux_info/ambig_info.tsv.

    ambig_info.tsv is row-aligned with quant.sf (no transcript ids of its
    own), so the alignment is verified by row count before use.
    """
    amb_path = os.path.join(sample_dir, "aux_info", "ambig_info.tsv")
    quant_path = os.path.join(sample_dir, "quant.sf")
    if not (os.path.exists(amb_path) and os.path.exists(quant_path)):
        return None
    amb = pd.read_csv(amb_path, sep="\t")
    quant = pd.read_csv(quant_path, sep="\t", usecols=["Name"])
    if len(amb) != len(quant):
        raise AssertionError(
            f"ambig_info.tsv ({len(amb)}) is not row-aligned with quant.sf "
            f"({len(quant)}) in {sample_dir}; cannot map to transcripts"
        )
    amb["transcript_id"] = quant["Name"].to_numpy()
    m = amb.merge(tx2gene[["transcript_id", "gene_id"]], on="transcript_id", how="left")
    g = (m.loc[m["gene_id"].isin(set(gene_ids))]
          .groupby("gene_id")[["UniqueCount", "AmbigCount"]].sum())
    g = g.rename(columns={"UniqueCount": "unique_reads_salmon_tx",
                          "AmbigCount": "ambig_reads_salmon_tx"})
    return g.reset_index()


def _require_samtools(use_bam: bool) -> None:
    # Fail fast, with a fixable message, before any heavy work.
    #
    # vr_quantify shells out to `samtools view -q 255` for the unique-read
    # evidence channel. `module load python/3.11.4` alone does NOT put
    # samtools on PATH; bin/run_pipeline.sh loads samtools/1.19 too, but a
    # manual invocation may not have. Without this check the failure surfaces
    # as a bare FileNotFoundError deep inside the per-gene loop, after the
    # cluster aggregation has already been computed.
    if not use_bam:
        return
    if shutil.which("samtools") is None:
        raise SystemExit(
            "samtools is not on PATH, and --no-bam was not passed.\n"
            "  fix (pick whichever applies to your environment):\n"
            "    conda install -c bioconda samtools     # conda/mamba\n"
            "    module load samtools                   # HPC with Environment Modules\n"
            "    apt install samtools / brew install samtools\n"
            "  or re-run with --no-bam to skip the unique-read channel\n"
            "  (DEGRADES the result: no individual-receptor calls are then possible).\n"
            "  or:   rerun with --no-bam (skips the unique-read evidence "
            "channel; NO individual-receptor call is then possible, and "
            "EM-redistribution flags degrade to fraction-only evidence)\n"
            "  or:   use bin/run_pipeline.sh, which loads both modules."
        )


def _samtools_region_counts(bam: str, region: str) -> Tuple[int, int]:
    """(unique alignments, unique non-duplicate alignments) at MAPQ 255."""
    cmd = ["samtools", "view", "-q", "255", bam, region]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"samtools failed on {region}: {proc.stderr[:400]}")
    n = 0
    n_nodup = 0
    for line in proc.stdout.splitlines():
        if not line or line[0] == "@":
            continue
        parts = line.split("\t", 2)
        if len(parts) < 2:
            continue
        try:
            flag = int(parts[1])
        except ValueError:
            continue
        n += 1
        if not (flag & 1024):
            n_nodup += 1
    return (n, n_nodup)


def bam_unique_counts(bam: str, genes: pd.DataFrame, threads: int = 8) -> pd.DataFrame:
    """MAPQ-255 alignment counts over each VR gene span (parallel index seeks)."""
    regions = [f"{r.chrom}:{int(r.start)}-{int(r.end)}" for r in genes.itertuples()]
    with ThreadPoolExecutor(max_workers=threads) as pool:
        res = list(pool.map(lambda rg: _samtools_region_counts(bam, rg), regions))
    return pd.DataFrame({
        "gene_id": genes["gene_id"].to_numpy(),
        "unique_reads_bam": [a for a, _ in res],
        "unique_reads_bam_nodup": [b for _, b in res],
    })


def flag_overlapping_spans(genes: pd.DataFrame) -> pd.DataFrame:
    """Mark VR genes whose span overlaps another VR gene span (double-count risk)."""
    g = genes.sort_values(["chrom", "start"]).copy()
    g["span_overlaps_other_vr"] = 0
    for _, idx in g.groupby("chrom").groups.items():
        block = g.loc[idx].sort_values("start")
        starts = block["start"].to_numpy()
        ends = block["end"].to_numpy()
        flags = np.zeros(len(block), dtype=int)
        for i in range(len(block)):
            for j in range(i + 1, len(block)):
                if starts[j] > ends[i]:
                    break
                flags[i] = 1
                flags[j] = 1
        g.loc[block.index, "span_overlaps_other_vr"] = flags
    return g[["gene_id", "span_overlaps_other_vr"]]


# ---------------------------------------------------------------- orchestration
def _sample_qc(cfg: Mapping[str, Any], trial: str) -> pd.DataFrame:
    """QC verdicts from the QC layer. Read its flag columns; never re-derive."""
    path = os.path.join(cfg["work"], "results", "sample_qc_all.tsv")
    qc = read_tsv_skipping_comments(path)
    qc = qc.loc[qc["trial"] == trial].copy()
    qc["suppress_biology_bool"] = qc["suppress_biology"].map(_as_bool)
    return qc


def _marker_table(cfg: Mapping[str, Any], trial: str) -> Optional[pd.DataFrame]:
    path = os.path.join(cfg["work"], "results", trial, "marker_cpm.tsv")
    if not os.path.exists(path):
        return None
    return read_tsv_skipping_comments(path)


def _qc_get(qc_by_sample: Mapping[str, Any], sample: str, col: str,
            default: Any = "NA") -> Any:
    """Scalar field from a QC row. Never use `or {}` on a Series (ambiguous truth)."""
    row = qc_by_sample.get(sample)
    if row is None:
        return default
    val = row.get(col, default)
    return default if val is None else val


def _is_suppressed(qc_by_sample: Mapping[str, Any], sample: str) -> bool:
    """Authoritative clearance verdict, read from the QC layer's own columns.

    `suppress_biology` is NOT sufficient on its own. It marks tissue/library
    failure, but a sample can pass that and still fail SORT VALIDATION, which
    is Rule 1 and outranks everything downstream: target2cellsRep3_S5 and
    pool2cellsRep1_S5 both carry suppress_biology=False with
    sort_verdict=FAIL and qc_overall=UNUSABLE. `qc_overall` is the composite
    gate the QC layer publishes, so it is the primary test here and
    suppress_biology is applied on top of it. A sample with no QC record is
    treated as suppressed -- absence of a verdict is not clearance.
    """
    row = qc_by_sample.get(sample)
    if row is None:
        return True
    if _as_bool(row.get("suppress_biology")):
        return True
    overall = row.get("qc_overall")
    if overall is None or str(overall).strip() == "":
        return False          # no composite column -> fall back to the above
    return str(overall).strip().upper() != "USABLE"


def _suppression_reason(qc_row: Optional[pd.Series]) -> str:
    if qc_row is None:
        return "no_qc_record"
    bits = []
    for col in ("qc_overall", "tissue_verdict", "sort_verdict", "library_status",
                "blocking_flags"):
        val = qc_row.get(col)
        if isinstance(val, str) and val and val.lower() not in ("nan", "none", ""):
            bits.append(f"{col}={val}")
    return "; ".join(bits) if bits else "suppressed_without_reason_recorded"


def run_trial(cfg: Mapping[str, Any], trial: str, outdir: Optional[str] = None,
              use_bam: bool = True, threads: int = 8) -> Dict[str, Any]:
    """Full Step 1-4 pass for one trial. Returns a summary dict."""
    samples_meta = vr_config.samples_of(cfg, trial)
    samples = list(samples_meta.keys())
    paths = vr_config.trial_paths(cfg, trial)
    outdir = outdir or os.path.join(cfg["work"], "results", trial)
    os.makedirs(outdir, exist_ok=True)

    g2c, ann_full = load_vr_reference(cfg)
    counts, cpm, name_map = load_counts_and_cpm(cfg, trial)
    counts = counts[[s for s in samples if s in counts.columns]]
    cpm = cpm[[s for s in samples if s in cpm.columns]]
    samples = list(counts.columns)
    recon = assert_join_complete(counts, name_map, g2c, ann_full)

    detect_floor = float(qcfg(cfg, "detect_floor_counts"))
    qc = _sample_qc(cfg, trial)
    qc_by_sample = {r["sample"]: r for _, r in qc.iterrows()}
    markers = _marker_table(cfg, trial)
    mk_by_sample = ({r["sample"]: r for _, r in markers.iterrows()}
                    if markers is not None else {})

    # ---- Step 1: cluster + supercluster expression
    cl_expr = cluster_expression(counts, cpm, g2c, samples, "cluster_id", detect_floor)
    sc_expr = cluster_expression(counts, cpm, g2c, samples, "supercluster_id", detect_floor)
    expr = pd.concat([cl_expr, sc_expr], ignore_index=True)

    min_cpm = float(qcfg(cfg, "cluster_call_min_cpm"))
    min_share = float(qcfg(cfg, "cluster_call_min_share"))
    expr["is_called"] = ((expr["cpm_sum"] >= min_cpm)
                         & (expr["share_of_sample_vr"] >= min_share)).astype(int)

    def _meta(sample: str, key: str, default: Any = None) -> Any:
        return samples_meta.get(sample, {}).get(key, default)

    expr["trial"] = trial
    expr["cell_type"] = expr["sample"].map(lambda s: _meta(s, "cell_type"))
    expr["n_cells"] = expr["sample"].map(lambda s: _meta(s, "n_cells"))
    expr["library_status"] = expr["sample"].map(
        lambda s: _qc_get(qc_by_sample, s, "library_status"))
    expr["qc_overall"] = expr["sample"].map(
        lambda s: _qc_get(qc_by_sample, s, "qc_overall"))
    expr["sort_verdict"] = expr["sample"].map(
        lambda s: _qc_get(qc_by_sample, s, "sort_verdict"))
    expr["suppress_biology"] = expr["sample"].map(
        lambda s: _is_suppressed(qc_by_sample, s))
    expr["suppression_reason"] = expr["sample"].map(
        lambda s: (_suppression_reason(qc_by_sample.get(s))
                   if _is_suppressed(qc_by_sample, s) else ""))
    expr["interpretation_context"] = expr["n_cells"].map(
        lambda n: "monogenic_expectation" if (n is not None and float(n) <= 2)
        else "pooled_ambiguous")

    # chr7 dual-tier reporting: the 60-paralog region is split by the 200kb
    # rule into two clusters reunited by supercluster V1R_chr7_sc013.
    chr7_ids = set(g2c.loc[g2c["cluster_id"].isin(["V1R_chr7_cl015", "V1R_chr7_cl016"]),
                           "supercluster_id"].unique())
    expr["chr7_dual_tier_region"] = (
        expr["cluster_id"].isin(["V1R_chr7_cl015", "V1R_chr7_cl016"])
        | expr["cluster_id"].isin(chr7_ids)).astype(int)

    # ---- within-cluster fractions
    frac = within_cluster_fractions(counts, cpm, g2c, samples)
    frac["trial"] = trial

    # ---- unique-read evidence
    gene_spans = g2c[["gene_id", "chrom", "start", "end"]].dropna().copy()
    gene_spans["start"] = gene_spans["start"].astype(int)
    gene_spans["end"] = gene_spans["end"].astype(int)
    overlaps = flag_overlapping_spans(gene_spans)
    cleared = [s for s in samples if not _is_suppressed(qc_by_sample, s)]

    ev_frames = []
    tx2gene = pd.read_csv(paths["tx2gene"], sep="\t")
    if "transcript_id" not in tx2gene.columns:
        tx2gene.columns = (["transcript_id", "gene_id", "gene_name"]
                           + list(tx2gene.columns[3:]))
    evidence_notes: List[str] = []
    for sample in samples:
        sdir = os.path.join(paths["per_sample_salmon"], sample)
        rec = pd.DataFrame({"gene_id": g2c["gene_id"]})
        try:
            amb = load_ambig_info(sdir, tx2gene, g2c["gene_id"])
        except AssertionError as exc:
            amb = None
            evidence_notes.append(f"{sample}: ambig_info unusable ({exc})")
        if amb is not None:
            rec = rec.merge(amb, on="gene_id", how="left")
        else:
            rec["unique_reads_salmon_tx"] = np.nan
            rec["ambig_reads_salmon_tx"] = np.nan
        bam = os.path.join(paths["per_sample_salmon"], f"{sample}.markdup.sorted.bam")
        if use_bam and sample in cleared and os.path.exists(bam):
            bc = bam_unique_counts(bam, gene_spans, threads=threads)
            rec = rec.merge(bc, on="gene_id", how="left")
        else:
            rec["unique_reads_bam"] = np.nan
            rec["unique_reads_bam_nodup"] = np.nan
            if use_bam and sample not in cleared:
                evidence_notes.append(
                    f"{sample}: BAM unique counting skipped (QC-suppressed sample)")
            elif use_bam:
                evidence_notes.append(f"{sample}: BAM not found at {bam}")
        rec["sample"] = sample
        ev_frames.append(rec)
    evidence = pd.concat(ev_frames, ignore_index=True)
    frac = frac.merge(evidence, on=["gene_id", "sample"], how="left")
    frac = frac.merge(overlaps, on="gene_id", how="left")

    # ---- Step 2/3: per (sample, cluster) flags
    mc = _MCCache(n_mc=int(qcfg(cfg, "n_mc")), seed=int(qcfg(cfg, "mc_seed")),
                  dom_alt=float(qcfg(cfg, "dom_alt")), alpha=float(qcfg(cfg, "alpha")))
    min_reads = int(qcfg(cfg, "min_reads_for_em_test"))
    min_power = float(qcfg(cfg, "min_power"))
    block_min_share = float(qcfg(cfg, "block_min_share"))
    alpha = float(qcfg(cfg, "alpha"))
    bleed_min_cpm = float(qcfg(cfg, "pseudogene_bleed_min_cpm"))

    flag_rows: List[Dict[str, Any]] = []
    for (sample, cid), grp in frac.groupby(["sample", "cluster_id"], sort=True):
        # Unique-read support must be positionally aligned with counts, so it
        # is taken from the SAME rows of the same frame. Prefer non-duplicate
        # MAPQ-255 BAM counts; fall back to the Salmon transcript-unique
        # channel where the BAM channel is absent for this sample.
        u_col = None
        for cand_col in ("unique_reads_bam_nodup", "unique_reads_bam",
                         "unique_reads_salmon_tx"):
            if cand_col in grp.columns and grp[cand_col].notna().any():
                u_col = cand_col
                break
        u_vals = (grp[u_col].to_numpy(dtype=float) if u_col is not None else None)
        stats = em_redistribution_stats(
            grp["counts"].to_numpy(dtype=float), mc, detect_floor, min_reads,
            min_power, block_min_share, alpha,
            strong_min_reads=int(qcfg(cfg, "strong_min_reads")),
            member_unique=u_vals,
            unique_min_reads=int(qcfg(cfg, "unique_gate_min_reads")),
            unique_bg_ratio=float(qcfg(cfg, "unique_gate_bg_ratio")))
        stats["unique_channel_used"] = (u_col or "none")
        det = grp.loc[grp["counts"] >= detect_floor].sort_values("counts", ascending=False)
        suppressed = _is_suppressed(qc_by_sample, sample)
        n_cells = _meta(sample, "n_cells")
        ctx = ("monogenic_expectation" if (n_cells is not None and float(n_cells) <= 2)
               else "pooled_ambiguous")

        # pseudogene bleed within a functional cluster
        is_ps = grp["is_pseudogene"].fillna(0).astype(float) > 0 if "is_pseudogene" in grp else None
        ps_cols: Dict[str, Any] = {
            "n_pseudogene_members": 0, "n_functional_members": 0,
            "pseudogene_counts": 0.0, "pseudogene_cpm": 0.0,
            "pseudogene_share_of_cluster": np.nan,
            "pseudogene_bleed_flag": "not_applicable",
            "pseudogene_genes_detected": "",
            "pseudogene_criteria_detected": "",
            "pseudogene_criterion_conflicting": 0,
            "pseudogene_evidence_grades": "",
            "pseudogene_mechanism": "",
        }
        if is_ps is not None:
            ps_cols["n_pseudogene_members"] = int(is_ps.sum())
            ps_cols["n_functional_members"] = int((~is_ps).sum())
            ps_det = grp.loc[is_ps & (grp["cpm"] >= bleed_min_cpm)]
            ps_counts = float(grp.loc[is_ps, "counts"].sum())
            ps_cpm = float(grp.loc[is_ps, "cpm"].sum())
            tot_c = float(grp["counts"].sum())
            ps_cols["pseudogene_counts"] = ps_counts
            ps_cols["pseudogene_cpm"] = ps_cpm
            ps_cols["pseudogene_share_of_cluster"] = (ps_counts / tot_c) if tot_c > 0 else np.nan
            if ps_cols["n_pseudogene_members"] == 0:
                ps_cols["pseudogene_bleed_flag"] = "no_pseudogene_members"
            elif ps_cols["n_functional_members"] == 0:
                ps_cols["pseudogene_bleed_flag"] = "pseudogene_only_cluster"
            elif len(ps_det) == 0:
                ps_cols["pseudogene_bleed_flag"] = "no_apparent_pseudogene_expression"
            else:
                ps_cols["pseudogene_bleed_flag"] = "apparent_pseudogene_expression"
                ps_cols["pseudogene_genes_detected"] = ",".join(
                    ps_det["gene_name"].astype(str).tolist())
                crits = ps_det.get("pseudogene_criterion")
                if crits is not None:
                    ps_cols["pseudogene_criteria_detected"] = ",".join(
                        crits.astype(str).tolist())
                    ps_cols["pseudogene_criterion_conflicting"] = int(
                        crits.map(pseudogene_confidence)
                             .eq("conflicting_evidence").any())
                ps_cols["pseudogene_evidence_grades"] = ",".join(
                        crits.map(pseudogene_confidence).tolist())
                ps_cols["pseudogene_mechanism"] = (
                    "UNRESOLVED: two candidate mechanisms, (a) multi-mapping/EM leakage "
                    "from an expressed functional paralog in the same cluster, (b) genuine "
                    "transcription of the pseudogene locus (e.g. cluster-shared regulatory "
                    "elements). This module does not adjudicate between them. Dietschi 2022 "
                    "is NOT a quantitative prior for mouse: their pseudogene result was "
                    "significant in rat only (P=0.003), not mouse (W=1214, P=0.5704), and "
                    "their proposed mechanism is regulatory, not multi-mapping."
                )

        row = {
            "trial": trial, "scope": "cluster", "sample": sample,
            "cell_type": _meta(sample, "cell_type"), "n_cells": n_cells,
            "cluster_id": cid,
            "family": grp["family"].iloc[0] if "family" in grp else None,
            "supercluster_id": (grp["supercluster_id"].iloc[0]
                                if "supercluster_id" in grp else None),
            "cluster_counts_sum": float(grp["counts"].sum()),
            "cluster_cpm_sum": float(grp["cpm"].sum()),
            "top_gene_name": (str(det["gene_name"].iloc[0]) if len(det) else ""),
            "detected_gene_names": ",".join(det["gene_name"].astype(str).tolist()),
            "detected_fracs": ",".join(
                f"{v:.3f}" for v in det["frac_of_cluster"].fillna(0).tolist()),
            "interpretation_context": ctx,
            "qc_overall": _qc_get(qc_by_sample, sample, "qc_overall"),
            "sort_verdict": _qc_get(qc_by_sample, sample, "sort_verdict"),
            "suppress_biology": suppressed,
            "suppression_reason": _suppression_reason(qc_by_sample.get(sample)) if suppressed else "",
            "biological_interpretation_permitted": (not suppressed),
        }
        row.update(stats)
        row.update(ps_cols)
        if suppressed:
            row["em_flag_note"] = ("numbers emitted for the record only; sample is "
                                   "QC-suppressed and no biology may be read from it")
        elif ctx == "pooled_ambiguous" and stats["em_flag"] == "suspected_em_redistribution":
            row["em_flag_note"] = ("near-uniform fractions in a multi-cell pool are ALSO the "
                                   "expected result of different cells choosing different "
                                   "paralogs of one cluster; this pair is called an artifact "
                                   "only because the co-dominant members lack independent "
                                   "unique-read support (see unique_support_threshold and "
                                   "block_unique_reads). Treat with more caution than the "
                                   "same flag on a 2-cell library.")
        elif stats["em_flag"] == "even_split_unique_supported":
            row["em_flag_note"] = ("even fractions with independent unique support on two or "
                                   "more members: co-expression, NOT an EM artifact. In a "
                                   "2-cell library this would mean the two cells chose "
                                   "different paralogs of one cluster; in a larger pool it "
                                   "is routine.")
        else:
            row["em_flag_note"] = ""
        flag_rows.append(row)
    flags = pd.DataFrame(flag_rows)

    # ---- Step 3: sample-scope expected-pattern checks
    exp_map = (cfg.get("thresholds", {}) or {}).get("expected_clusters", {}) or {}
    nt_max = float((cfg.get("thresholds", {}) or {})
                   .get("nontarget_total_vr_cpm_max", 100.0))
    sample_rows: List[Dict[str, Any]] = []
    cl_only = expr.loc[expr["tier"] == "cluster"]
    for sample in samples:
        s_expr = cl_only.loc[cl_only["sample"] == sample]
        n_called = int(s_expr["is_called"].sum())
        total_vr_cpm = float(s_expr["sample_total_vr_cpm"].iloc[0]) if len(s_expr) else np.nan
        cell_type = _meta(sample, "cell_type")
        n_cells = _meta(sample, "n_cells")
        key = (f"{cell_type}_{int(n_cells)}cell" if cell_type == "nontarget"
               else f"target_{int(n_cells)}cell")
        rng_exp = exp_map.get(key)
        lo, hi = (rng_exp[0], rng_exp[1]) if isinstance(rng_exp, (list, tuple)) else (np.nan, np.nan)
        suppressed = _is_suppressed(qc_by_sample, sample)
        in_range = (np.nan if (suppressed or not np.isfinite(float(lo)))
                    else int(lo <= n_called <= hi))
        checks = []
        if suppressed:
            # Never emit a bare "none" for a suppressed sample: an empty flag
            # list would read as "checked and passed". The numbers are still
            # emitted (per brief) but carry no biological verdict.
            checks.append("not_evaluated_sample_suppressed")
        if np.isfinite(float(lo)) and not (lo <= n_called <= hi):
            checks.append(f"cluster_count_out_of_expected_range({n_called} vs {lo}-{hi})")
        if cell_type == "nontarget" and np.isfinite(total_vr_cpm) and total_vr_cpm > nt_max:
            checks.append(f"nontarget_vr_cpm_above_max({total_vr_cpm:.1f} > {nt_max})")
        top = s_expr.sort_values("counts_sum", ascending=False).head(1)
        sample_rows.append({
            "trial": trial, "scope": "sample", "sample": sample,
            "cell_type": cell_type, "n_cells": n_cells,
            "cluster_id": "", "expected_key": key,
            "expected_clusters_lo": lo, "expected_clusters_hi": hi,
            "n_clusters_called": n_called,
            "cluster_count_in_expected_range": in_range,
            "sample_total_vr_cpm": total_vr_cpm,
            "nontarget_vr_cpm_max": nt_max if cell_type == "nontarget" else np.nan,
            "dominant_cluster_id": (str(top["cluster_id"].iloc[0]) if len(top) else ""),
            "dominant_cluster_share": (float(top["share_of_sample_vr"].iloc[0])
                                       if len(top) else np.nan),
            "expected_pattern_flags": ";".join(checks) if checks else "none",
            "qc_overall": _qc_get(qc_by_sample, sample, "qc_overall"),
            "sort_verdict": _qc_get(qc_by_sample, sample, "sort_verdict"),
            "suppress_biology": suppressed,
            "suppression_reason": _suppression_reason(qc_by_sample.get(sample)) if suppressed else "",
            "biological_interpretation_permitted": (not suppressed),
            "em_flag_note": ("numbers emitted for the record only; sample is QC-suppressed"
                             if suppressed else ""),
        })
    flags = pd.concat([flags, pd.DataFrame(sample_rows)], ignore_index=True)

    # ---- Step 4: candidate narrowing on non-EM evidence
    cand_min_share = float(qcfg(cfg, "candidate_cluster_min_share"))
    min_uniq = int(qcfg(cfg, "min_unique_reads_for_rank"))
    r_mod = float(qcfg(cfg, "unique_ratio_moderate"))
    r_low = float(qcfg(cfg, "unique_ratio_low"))
    u_mod = int(qcfg(cfg, "unique_reads_moderate"))
    max_cand = int(qcfg(cfg, "max_candidates_reported"))

    cand_rows: List[Dict[str, Any]] = []
    for sample in cleared:
        s_expr = cl_only.loc[(cl_only["sample"] == sample)].sort_values(
            "share_of_sample_vr", ascending=False)
        if not len(s_expr):
            continue
        pop_call = str(_qc_get(mk_by_sample, sample, "population_call"))
        picked = s_expr.loc[(s_expr["is_called"] == 1)
                            & (s_expr["share_of_sample_vr"] >= cand_min_share)]
        dom_id = str(s_expr["cluster_id"].iloc[0])
        if not len(picked):
            # No cluster carries enough signal to narrow. Record that fact
            # instead of ranking paralogs inside noise: the top "cluster" of a
            # sample with ~0 total VR CPM is a single stray read.
            top_row = s_expr.iloc[0]
            cand_rows.append({
                "trial": trial, "sample": sample, "cluster_id": dom_id,
                "supercluster_id": "", "family": str(top_row.get("family", "")),
                "is_dominant_cluster": 1,
                "cluster_share_of_vr": float(top_row["share_of_sample_vr"])
                if np.isfinite(top_row["share_of_sample_vr"]) else np.nan,
                "cluster_cpm": float(top_row["cpm_sum"]),
                "rank": 0, "gene_id": "", "gene_name": "", "candidate_list": "",
                "n_candidates_reported": 0,
                "evidence_type": "no_cluster_above_signal_threshold",
                "unique_reads_bam": 0, "unique_reads_bam_nodup": 0,
                "unique_reads_salmon_tx": 0, "unique_share_of_cluster": np.nan,
                "em_counts": np.nan, "em_frac_of_cluster": np.nan,
                "em_flag": "insufficient_signal", "even_block_size": 0,
                "confidence": "no_call",
                "confirmation_status": "tentative_unconfirmed",
                "marker_consistency": "", "population_call": pop_call,
                "notes": ("no cluster in this sample reaches the calling thresholds "
                          f"(>= {min_cpm} cluster CPM and >= {min_share:.2f} share of "
                          f"sample VR signal); sample total VR CPM = "
                          f"{float(top_row['sample_total_vr_cpm']):.3f}. No receptor "
                          "call is made. For a cleared target sample this is itself a "
                          "reportable negative: the library passed QC but carries no "
                          "interpretable VR signal."),
            })
            continue
        for _, crow in picked.iterrows():
            cid = str(crow["cluster_id"])
            members = frac.loc[(frac["sample"] == sample) & (frac["cluster_id"] == cid)].copy()
            fam = str(members["family"].iloc[0]) if len(members) else ""
            u_bam = members["unique_reads_bam"].fillna(0).to_numpy(dtype=float)
            u_tx = members["unique_reads_salmon_tx"].fillna(0).to_numpy(dtype=float)
            if np.nansum(u_bam) >= min_uniq:
                ev_type, ev_col = "bam_unique_mapq255", "unique_reads_bam"
            elif np.nansum(u_tx) >= min_uniq:
                ev_type, ev_col = "salmon_tx_unique", "unique_reads_salmon_tx"
            else:
                ev_type, ev_col = "none_available", None

            fl = flags.loc[(flags["scope"] == "cluster") & (flags["sample"] == sample)
                           & (flags["cluster_id"] == cid)]
            em_flag = str(fl["em_flag"].iloc[0]) if len(fl) else "NA"
            blk = int(fl["even_block_size"].iloc[0]) if len(fl) else 0
            marker_consistent = ""
            if fam == "V1R":
                marker_consistent = ("consistent" if pop_call.upper().startswith("V1R")
                                     else f"CONTRADICTION(cluster=V1R, population_call={pop_call})")
            elif fam == "V2R":
                marker_consistent = ("consistent" if pop_call.upper().startswith("V2R")
                                     else f"CONTRADICTION(cluster=V2R, population_call={pop_call})")

            if ev_type == "none_available":
                names = members.sort_values("counts", ascending=False)
                names = names.loc[names["counts"] >= detect_floor].head(max_cand)
                cand_rows.append({
                    "trial": trial, "sample": sample, "cluster_id": cid,
                    "supercluster_id": (str(members["supercluster_id"].iloc[0])
                                        if len(members) else ""),
                    "family": fam, "is_dominant_cluster": int(cid == dom_id),
                    "cluster_share_of_vr": float(crow["share_of_sample_vr"]),
                    "cluster_cpm": float(crow["cpm_sum"]),
                    "rank": 0, "gene_id": "", "gene_name": "",
                    "candidate_list": ",".join(names["gene_name"].astype(str).tolist()),
                    "n_candidates_reported": int(len(names)),
                    "evidence_type": "none_available",
                    "unique_reads_bam": 0, "unique_reads_bam_nodup": 0,
                    "unique_reads_salmon_tx": 0,
                    "unique_share_of_cluster": np.nan,
                    "em_counts": np.nan, "em_frac_of_cluster": np.nan,
                    "em_flag": em_flag, "even_block_size": blk,
                    "confidence": "unresolvable",
                    "confirmation_status": "tentative_unconfirmed",
                    "marker_consistency": marker_consistent,
                    "population_call": pop_call,
                    "notes": ("no EM-independent read support in this cluster "
                              f"(< {min_uniq} unique reads across all members in both the "
                              "MAPQ-255 BAM channel and Salmon transcript-unique channel); "
                              "paralogs CANNOT be separated with this data. Members listed "
                              "in candidate_list are ordered by EM counts for reference "
                              "ONLY and that order carries no evidential weight."),
                })
                continue

            members["_uniq"] = members[ev_col].fillna(0).astype(float)
            tot_u = float(members["_uniq"].sum())
            ranked = members.sort_values(["_uniq", "counts"], ascending=False)
            ranked = ranked.loc[ranked["_uniq"] > 0].head(max_cand)
            uu = ranked["_uniq"].to_numpy(dtype=float)
            ratio = (uu[0] / uu[1]) if len(uu) >= 2 and uu[1] > 0 else np.inf

            # Background floor. Most members of a large cluster carry a small,
            # near-constant unique count (mismapping / ambient reads), so the
            # median BACKGROUND member is an empirical noise floor rather than
            # zero. The top gene is EXCLUDED from that median -- including it
            # would make a single-member cluster its own background (ratio 1.0)
            # and spuriously demote a gene that is in fact the only candidate.
            # With no other member there is no background to clear and the
            # test is not applicable (inf).
            srt_bg = members.sort_values("_uniq", ascending=False)["_uniq"]
            bg = srt_bg.to_numpy(dtype=float)[1:]
            floor_u = float(np.median(bg)) if bg.size else np.nan
            top_over_floor = (uu[0] / floor_u) if (len(uu) and np.isfinite(floor_u)
                                                   and floor_u > 0) else np.inf

            # EM/unique contradiction. The two channels should agree about
            # which paralog carries the signal. When the top unique-read gene
            # holds almost none of the EM signal AND the EM-dominant gene has
            # essentially no unique reads, the channels point at DIFFERENT
            # genes and neither can be trusted to name the receptor. This is
            # surfaced, never smoothed over.
            em_dom = members.sort_values("counts", ascending=False).iloc[0]
            top_u_row = ranked.iloc[0] if len(ranked) else None
            contradiction = ""
            if top_u_row is not None:
                top_u_em_frac = (float(top_u_row["frac_of_cluster"])
                                 if np.isfinite(top_u_row["frac_of_cluster"]) else 0.0)
                em_dom_uniq = float(em_dom["_uniq"])
                if (top_u_em_frac < 0.05 and em_dom_uniq <= max(floor_u, 2.0)
                        and float(em_dom["frac_of_cluster"]) > 0.5):
                    contradiction = (
                        f"EM_UNIQUE_CONTRADICTION: unique reads rank "
                        f"{top_u_row['gene_name']} first (EM fraction "
                        f"{top_u_em_frac:.4f}) while the EM-dominant paralog "
                        f"{em_dom['gene_name']} (EM fraction "
                        f"{float(em_dom['frac_of_cluster']):.3f}) carries "
                        f"{em_dom_uniq:.0f} unique reads. The two evidence "
                        "channels name different genes; neither identifies the "
                        "expressed receptor in this cluster."
                    )
            for i, (_, g) in enumerate(ranked.iterrows(), start=1):
                if i == 1:
                    if contradiction:
                        conf = "unresolvable"
                    elif (uu[0] >= u_mod and ratio >= r_mod
                          and top_over_floor >= r_mod):
                        conf = "moderate"
                    elif (uu[0] >= min_uniq and ratio >= r_low
                          and top_over_floor >= r_low):
                        conf = "low"
                    else:
                        conf = "unresolvable"
                else:
                    conf = "alternative_candidate"
                cand_rows.append({
                    "trial": trial, "sample": sample, "cluster_id": cid,
                    "supercluster_id": str(g.get("supercluster_id", "")),
                    "family": fam, "is_dominant_cluster": int(cid == dom_id),
                    "cluster_share_of_vr": float(crow["share_of_sample_vr"]),
                    "cluster_cpm": float(crow["cpm_sum"]),
                    "rank": i, "gene_id": str(g["gene_id"]),
                    "gene_name": str(g.get("gene_name", "")),
                    "candidate_list": ",".join(ranked["gene_name"].astype(str).tolist()),
                    "n_candidates_reported": int(len(ranked)),
                    "evidence_type": ev_type,
                    "unique_reads_bam": float(g.get("unique_reads_bam", np.nan) or 0),
                    "unique_reads_bam_nodup": float(g.get("unique_reads_bam_nodup", np.nan) or 0),
                    "unique_reads_salmon_tx": float(g.get("unique_reads_salmon_tx", np.nan) or 0),
                    "unique_share_of_cluster": (float(g["_uniq"]) / tot_u) if tot_u > 0 else np.nan,
                    "em_counts": float(g["counts"]),
                    "em_frac_of_cluster": float(g["frac_of_cluster"])
                    if np.isfinite(g["frac_of_cluster"]) else np.nan,
                    "em_flag": em_flag, "even_block_size": blk,
                    "confidence": conf,
                    "confirmation_status": "tentative_unconfirmed",
                    "marker_consistency": marker_consistent,
                    "population_call": pop_call,
                    "unique_reads_cluster_median": floor_u,
                    "top_over_background_floor": (
                        float(top_over_floor) if np.isfinite(top_over_floor) else np.inf),
                    "evidence_contradiction": contradiction,
                    "notes": ((contradiction + " ") if contradiction else "")
                             + (f"ranked by {ev_type}; unique-read ratio rank1/rank2="
                              f"{'inf' if not np.isfinite(ratio) else format(ratio, '.2f')}; "
                              f"background floor (median unique reads of non-top "
                              f"members)={floor_u:.1f}, "
                              f"top/floor="
                              f"{'inf' if not np.isfinite(top_over_floor) else format(top_over_floor, '.2f')}; "
                              "EM counts shown for context only and were NOT used to rank. "
                              "n_transcripts="
                              f"{g.get('n_transcripts', 'NA')}"
                              + ("; transcript-unique counts understate gene-level unique "
                                 "support for multi-transcript genes"
                                 if ev_type == "salmon_tx_unique" else "")),
                })
    candidates = pd.DataFrame(cand_rows, columns=CANDIDATE_COLUMNS)

    # ---- write
    hdr = ("# vr_quantify.py -- trial=%s; CPM = count / (all-gene column sum of "
           "unscaled salmon.merged.gene_counts.tsv) * 1e6; cluster tier=200kb, "
           "supercluster tier=800kb; individual-paralog calls are tentative by "
           "construction\n" % trial)
    written = {}
    for name, df in (("vr_cluster_expression.tsv", expr),
                     ("vr_within_cluster_fractions.tsv", frac),
                     ("vr_artifact_flags.tsv", flags),
                     ("vr_candidates.tsv", candidates)):
        path = os.path.join(outdir, name)
        with open(path, "w") as fh:
            fh.write(hdr)
            df.to_csv(fh, sep="\t", index=False)
        written[name] = path
    return {
        "trial": trial, "samples": samples, "cleared_samples": cleared,
        "reconciliation": recon, "files": written,
        "evidence_notes": evidence_notes,
        "n_clusters": int(cl_only["cluster_id"].nunique()),
        "n_flag_rows": int(len(flags)), "n_candidate_rows": int(len(candidates)),
    }


def _main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--config", default=None)
    ap.add_argument("--trial", default="trial2")
    ap.add_argument("--all-trials", action="store_true")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--no-bam", action="store_true",
                    help="skip MAPQ-255 BAM unique counting (Salmon channel only)")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args(argv)
    _require_samtools(not args.no_bam)

    cfg = vr_config.load_config(args.config)
    trials = vr_config.trials_of(cfg) if args.all_trials else [args.trial]
    summary = {}
    for trial in trials:
        summary[trial] = run_trial(cfg, trial, outdir=args.outdir,
                                   use_bam=not args.no_bam, threads=args.threads)
    print(json.dumps(summary, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
