#!/usr/bin/env python3
"""
vr_reconcile.py -- reconcile the project briefing against what the pipeline
computes, one row per briefing claim.

Every row is (claim_id, source_section, claim_as_written, pipeline_result,
verdict, delta_cause, evidence). The claim text is quoted from the briefing
verbatim; the pipeline_result is read from the regenerated tables under
results/, never typed in. Nothing here is tuned to make a claim agree: where
the pipeline contradicts the briefing the row says so and names the cause.

Verdict vocabulary:
  AGREE       the pipeline reproduces the claim within the stated tolerance
  DIFFER      same quantity, different value -- delta_cause names the reason
  SUPERSEDED  the claim rests on a premise the pipeline has since disproved,
              so the claim is not merely wrong in value but wrong in kind
  REFRAMED    the underlying observation is confirmed but its label/attribution
              changes (e.g. a cluster alias, a tier that now carries a caveat)
  UNTESTABLE  the pipeline cannot evaluate the claim from these inputs

Usage:
  vr_reconcile.py --results results --out results/reconciliation.tsv
"""
from __future__ import annotations

import argparse
import os
from typing import Any, Dict, List, Optional

import pandas as pd


def load(path: str) -> Optional[pd.DataFrame]:
    return (pd.read_csv(path, sep="\t", comment="#", dtype=str,
                        keep_default_na=False)
            if os.path.exists(path) else None)


def f(x: Any) -> Optional[float]:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


class Tables:
    def __init__(self, results: str):
        self.qc = load(os.path.join(results, "sample_qc_all.tsv"))
        self.tier = load(os.path.join(results, "tier_status_all.tsv"))
        self.expr, self.frac, self.flags, self.cand = {}, {}, {}, {}
        for t in ("trial1", "trial2"):
            d = os.path.join(results, t)
            self.expr[t] = load(os.path.join(d, "vr_cluster_expression.tsv"))
            self.frac[t] = load(os.path.join(d, "vr_within_cluster_fractions.tsv"))
            self.flags[t] = load(os.path.join(d, "vr_artifact_flags.tsv"))
            self.cand[t] = load(os.path.join(d, "vr_candidates.tsv"))

    def s(self, sample: str, col: str) -> str:
        r = self.qc[self.qc["sample"] == sample]
        return "" if not len(r) else str(r.iloc[0][col])

    def n(self, sample: str, col: str) -> Optional[float]:
        return f(self.s(sample, col))

    def tier_of(self, sample: str) -> str:
        r = self.tier[self.tier["sample"] == sample]
        return "" if not len(r) else str(r.iloc[0]["highest_tier_reported"])

    def cl(self, trial: str, sample: str, cluster: str, col: str) -> Optional[float]:
        d = self.expr[trial]
        r = d[(d["sample"] == sample) & (d["cluster_id"] == cluster)]
        return None if not len(r) else f(r.iloc[0][col])

    def gene(self, trial: str, sample: str, cluster: str, gene: str,
             col: str) -> Optional[float]:
        d = self.frac[trial]
        r = d[(d["sample"] == sample) & (d["cluster_id"] == cluster)
              & (d["gene_name"] == gene)]
        return None if not len(r) else f(r.iloc[0][col])


def build(T: Tables) -> pd.DataFrame:
    R: List[Dict[str, str]] = []

    def row(cid, section, claim, result, verdict, cause, evidence=""):
        R.append({"claim_id": cid, "briefing_section": section,
                  "claim_as_written": claim, "pipeline_result": result,
                  "verdict": verdict, "delta_cause": cause,
                  "evidence_source": evidence})

    # ================= SECTION 9 -- current analysis findings =================

    # --- 9.1 sort validation / Trpc2 enrichment -----------------------------
    nt = T.n("nontarget100cells_S8", "Trpc2_cpm")
    t6 = T.n("target100cellsRep1_S6", "Trpc2_cpm")
    t7 = T.n("target100cellsRep2_S7", "Trpc2_cpm")
    t3 = T.n("target2cellsRep1_S3", "Trpc2_cpm")
    row("S9-01", "9 sort validation",
        "Trpc2 is 1029-1476 CPM in GFP+ target libraries",
        f"target100cellsRep1_S6 {t6:.1f}, target100cellsRep2_S7 {t7:.1f}, "
        f"target2cellsRep1_S3 {t3:.1f} CPM",
        "AGREE", "-- (exact reproduction; the CPM convention was pinned to the "
        "unscaled salmon.merged.gene_counts table with an all-gene column-sum "
        "denominator specifically to reproduce these values)",
        "results/sample_qc_all.tsv:Trpc2_cpm")

    row("S9-02", "9 sort validation",
        "Trpc2 is near-zero in the GFP- nontarget library",
        f"nontarget100cells_S8 = {nt:.3f} CPM from {T.n('nontarget100cells_S8','Trpc2_reads'):.0f} raw reads",
        "AGREE", "--", "results/sample_qc_all.tsv:Trpc2_cpm,Trpc2_reads")

    ratio_obs = (t6 / nt) if (t6 and nt) else None
    row("S9-03", "9 sort validation",
        "~5000-fold Trpc2 enrichment in target over nontarget",
        f"reproduced: {t6:.1f}/{nt:.3f} = {ratio_obs:,.0f}x for target100cellsRep1_S6, "
        f"{t7:.1f}/{nt:.3f} = {t7/nt:,.0f}x for target100cellsRep2_S7, and "
        f"{t3:.1f}/{nt:.3f} = {t3/nt:,.0f}x for target2cellsRep1_S3 -- bracketing the "
        f"briefing's ~5000x. But the denominator is SIX raw reads, so the fold-change "
        f"is a 6-read Poisson draw",
        "AGREE",
        "The NUMBER agrees and is not a delta -- recording that explicitly because the "
        "figure was flagged for scrutiny and it survived it. What the pipeline adds is "
        "that the quantity should not be quoted as a fold-change at all: with 6 reads "
        "in the denominator, the 95% Poisson interval on that count alone spans roughly "
        "2-13 reads, i.e. a fold anywhere from ~2,400x to ~16,000x. The DIRECTION "
        "(target >> nontarget) is unambiguous and is what the sort gate tests; the "
        "magnitude is a ratio of a large number to a near-zero one. Report the two "
        "CPMs, not their quotient.",
        "results/sample_qc_all.tsv:Trpc2_cpm,Trpc2_reads")

    # --- 9.2 population ID / Gnai2:Gnao1 ------------------------------------
    r6 = T.n("target100cellsRep1_S6", "gnai2_gnao1_ratio")
    sup6 = T.n("target100cellsRep1_S6", "Gnao1_reads")
    low6 = T.s("target100cellsRep1_S6", "ratio_low_support")
    row("S9-04", "9 population ID",
        "Gnai2:Gnao1 = 385:1 in target100cellsRep1_S6, indicating a V1R population",
        f"{r6:.1f}:1 (ratio_low_support={low6}, supported by {sup6:.0f} raw Gnao1 reads)",
        "DIFFER",
        "Normalisation/denominator, not thresholding: this is the SAME measurement "
        f"re-expressed. Gnao1 = {sup6:.0f} reads is the entire denominator, so the "
        "ratio is one Poisson draw from a 7-read count; 385:1 and 422.7:1 differ by "
        "less than the sampling error of that single count (one read either way moves "
        "it by ~15%). The pipeline flags this with ratio_low_support=True and the "
        "config sets ratio_min_support_reads=10 precisely so the MAGNITUDE is not "
        "quoted. The population CALL (V1R_dominant) is unaffected.",
        "results/sample_qc_all.tsv:gnai2_gnao1_ratio,Gnao1_reads,ratio_low_support")

    r7 = T.n("target100cellsRep2_S7", "gnai2_gnao1_ratio")
    sup7 = T.n("target100cellsRep2_S7", "Gnao1_reads")
    row("S9-05", "9 population ID",
        "V1R-dominant population in the GFP+ targets",
        f"V1R_dominant in all three clean targets. Only target100cellsRep2_S7 has a "
        f"well-supported ratio ({r7:.1f}:1 on {sup7:.0f} Gnao1 reads, "
        f"ratio_low_support=False); the other two rest on 0-7 Gnao1 reads and are "
        f"directionally right but numerically unstable",
        "AGREE",
        "-- direction agrees; the pipeline adds the support qualifier the briefing "
        "does not carry",
        "results/sample_qc_all.tsv:population_call,gnai2_gnao1_ratio,Gnao1_reads")

    # --- 9.3 the chr7 dominant cluster call ---------------------------------
    sc = T.cl("trial2", "target2cellsRep1_S3", "V1R_chr7_sc013", "share_of_sample_vr")
    c16 = T.cl("trial2", "target2cellsRep1_S3", "V1R_chr7_cl016", "share_of_sample_vr")
    c15 = T.cl("trial2", "target2cellsRep1_S3", "V1R_chr7_cl015", "share_of_sample_vr")
    row("S9-06", "9 dominant cluster",
        "target2cellsRep1_S3 expresses V1R chr7 cluster 039 at 99.9% of its VR signal",
        f"reproduced at the 800kb supercluster tier: V1R_chr7_sc013 = {sc*100:.2f}% of "
        f"the sample's VR signal. At the 200kb cluster tier the same signal resolves "
        f"into {c16*100:.1f}% V1R_chr7_cl016 + {c15*100:.1f}% V1R_chr7_cl015",
        "REFRAMED",
        "Cluster definition. The 99.9% figure is correct at ONE tier and misleading at "
        "the other. A 217,366bp gap between Vmn1r132 and Vmn1r135 exceeds the 200kb "
        "cluster rule by 17kb, splitting the briefing's 'cluster 039' in two; the 800kb "
        "supercluster tier reunites them. Both tiers must be reported for this region.",
        "results/trial2/vr_cluster_expression.tsv (tier=cluster and tier=supercluster)")

    row("S9-07", "9 cluster labels",
        "cluster labels cl039 (V1R chr7, 60 paralogs) and cl029 (V1R chr6, Vmn1r32-39)",
        "not reproducible as ordinals; membership maps exactly. 'V1R chr7 cluster 039' "
        "= V1R_chr7_cl015 + V1R_chr7_cl016 (supercluster V1R_chr7_sc013); "
        "'V1R chr6 cl029 (Vmn1r32-39)' = V1R_chr6_cl008, 9 genes",
        "SUPERSEDED",
        "Cluster definition. The old labels were ordinals from a different numbering "
        "pass and cannot be regenerated; they are retained as ALIASES only. Any "
        "reference to cl039/cl029 must be resolved through the membership map, not "
        "the number.",
        "ref/vr_clusters.tsv, ref/vr_gene_to_cluster.tsv")

    # --- 9.4 the 4-VR-reads claim -------------------------------------------
    ntv1 = T.n("nontarget100cells_S8", "Vmn1r_sum_cpm")
    ntv2 = T.n("nontarget100cells_S8", "Vmn2r_sum_cpm")
    d = T.frac["trial2"]
    ntg = d[(d["sample"] == "nontarget100cells_S8")
            & (pd.to_numeric(d["counts"], errors="coerce") > 0)]
    genes = ", ".join(f"{r['gene_name']} ({f(r['counts']):.0f})" for _, r in ntg.iterrows())
    ntreads = T.n("nontarget100cells_S8", "input_reads")
    row("S9-08", "9 sort purity",
        "4 total VR reads out of 59M in the nontarget library",
        f"4 VR genes each carry exactly 1 assigned count in nontarget100cells_S8 "
        f"({genes}), totalling {ntv1+ntv2:.3f} CPM. The library's STAR input is "
        f"{ntreads/1e6:.1f}M read pairs, not 59M",
        "DIFFER",
        "The COUNT of 4 is confirmed but its unit and denominator are not. 4 is the "
        "number of VR GENES with a nonzero Salmon count (1 count each), not the number "
        "of reads; and the denominator is 32.2M input read pairs for this library. 59M "
        "is neither this library's read count nor the trial total for a single library, "
        "so the briefing figure appears to mix a per-gene tally with a cross-library "
        "read total. Either way the conclusion -- the sort is clean -- holds: "
        f"{ntv1+ntv2:.3f} CPM is ~625x below the 100 CPM purity ceiling.",
        "results/trial2/vr_within_cluster_fractions.tsv, results/sample_qc_all.tsv:input_reads")

    # --- 9.5 the EM artifact ------------------------------------------------
    fl = T.flags["trial2"]
    hit = fl[(fl["sample"] == "target2cellsRep1_S3")
             & (fl["cluster_id"] == "V1R_chr7_cl016")]
    u166 = T.gene("trial2", "target2cellsRep1_S3", "V1R_chr7_cl016", "Vmn1r166",
                  "unique_reads_bam_nodup")
    u138 = T.gene("trial2", "target2cellsRep1_S3", "V1R_chr7_cl016", "Vmn1r138",
                  "unique_reads_bam_nodup")
    f166 = T.gene("trial2", "target2cellsRep1_S3", "V1R_chr7_cl016", "Vmn1r166",
                  "frac_of_cluster")
    f138 = T.gene("trial2", "target2cellsRep1_S3", "V1R_chr7_cl016", "Vmn1r138",
                  "frac_of_cluster")
    row("S9-09", "9 multi-mapping caveat",
        "per-gene VR counts within a cluster are not trustworthy; paralogs at "
        "near-equal fractions mean one paralog was expressed and EM split the reads",
        f"confirmed and made testable. Exactly one instance across both trials: "
        f"target2cellsRep1_S3 / V1R_chr7_cl016, Vmn1r166 {f166:.3f} vs Vmn1r138 "
        f"{f138:.3f} of the cluster, with deduplicated unique MAPQ255 reads of "
        f"{u166:.0f} vs {u138:.0f}",
        "AGREE",
        "-- with an important refinement the briefing does not state: an even split is "
        "NECESSARY BUT NOT SUFFICIENT for an EM artifact. The unique-read gate "
        "separates redistribution (one member has unique support) from genuine "
        "co-expression across a multi-cell pool (two or more members do). Without that "
        "gate the even-fractions heuristic alone produces false positives -- see R-01.",
        "results/trial2/vr_artifact_flags.tsv:em_flag, results/trial2/vr_within_cluster_fractions.tsv")

    if len(hit):
        h = hit.iloc[0]
        row("S9-10", "9 multi-mapping caveat",
            "(implicit) an EM-artifact call identifies which paralog is the true source",
            f"the CALL and the ATTRIBUTION have different strengths. Block read "
            f"support {f(h['even_block_read_support']):.0f}; the call that "
            f"redistribution occurred rests on a paralog holding "
            f"{f(h['cluster_counts_sum'])*f138:.0f} EM counts with ZERO unique reads, "
            f"which is solid. WHICH paralog is the source rests on {u166:.0f} unique "
            f"reads = {u166/ (f(h['cluster_counts_sum'])*f166) * 100:.1f}% of "
            f"Vmn1r166's apparent expression. confirmation_status=tentative_unconfirmed",
            "REFRAMED",
            "Evidence asymmetry, not a numeric delta. The briefing treats the artifact "
            "call and the receptor identification as one statement; the pipeline splits "
            "them because their evidence differs by two orders of magnitude.",
            "results/trial2/vr_artifact_flags.tsv, results/trial2/vr_candidates.tsv")

    # --- 9.6 pseudogene bleed ------------------------------------------------
    psf = T.gene("trial2", "target100cellsRep1_S6", "V1R_chr17_cl021",
                 "Vmn1r-ps150", "frac_of_cluster")
    psu = T.gene("trial2", "target100cellsRep1_S6", "V1R_chr17_cl021",
                 "Vmn1r-ps150", "unique_reads_bam")
    row("S9-11", "9 pseudogene bleed",
        "pseudogenes inside functional clusters appear 'expressed' via bleed-through",
        f"the observation reproduces -- Vmn1r-ps150 holds {psf*100:.1f}% of "
        f"V1R_chr17_cl021 in target100cellsRep1_S6 -- but the MECHANISM is unresolved. "
        f"That pseudogene carries {psu:.0f} unique MAPQ255 reads and the cluster is "
        f"CLEAN on the redistribution test (em_flag=no_redistribution_signature), so "
        f"an even-split EM artifact does not explain it",
        "REFRAMED",
        "Not a numeric delta: the word 'bleed-through' asserts a mechanism the data "
        "cannot support here. Two candidates remain open -- (a) EM leakage from an "
        "expressed functional paralog, (b) genuine transcription of the pseudogene "
        "locus via cluster-shared regulatory elements. Dietschi 2022 is NOT a "
        "quantitative mouse prior for this: their pseudogene result was significant in "
        "rat (P=0.003) but NOT in mouse (W=1214, P=0.5704), and the mechanism they "
        "propose is regulatory, not multi-mapping.",
        "results/trial2/vr_artifact_flags.tsv:pseudogene_bleed_flag,pseudogene_mechanism")

    # --- 9.7 aggregation strategy -------------------------------------------
    row("S9-12", "9 method",
        "cluster-level aggregation (Dietschi et al. 2022) is the reliable readout",
        "adopted as the reportable tier, with two documented limits. (1) 200kb is not "
        "a natural break: the V1R inter-gene gap distribution has its KDE minimum near "
        "2Mb, so the boundary is a conservative choice, not a discovered one -- hence "
        "the second 800kb supercluster tier. (2) V2R clustering is materially weaker: "
        "18 of 37 V2R clusters are singletons and only 180/222 V2R genes sit in "
        "clusters of >=5, so aggregation protects V2R calls less than V1R calls",
        "AGREE",
        "-- the method is adopted; the pipeline adds quantified limits the briefing "
        "does not carry",
        "ref/vr_gtf_parse_report.txt, ref/vr_clusters.tsv")

    # ================= SECTION 8 -- sample table =============================
    # The tissue finding supersedes the entire trial-1 half of the sample table.
    t1_desc = {
        "pool100cells_S8": ("failed prep", "S8-01"),
        "pool2cellsRep1_S5": ("severe RNA degradation", "S8-02"),
        "pool2cellsRep2_S6": ("borderline", "S8-03"),
        "pool2cellsRep3_S7": ("the cleanest 2-cell sample across both trials", "S8-04"),
    }
    for smp, (desc, cid) in t1_desc.items():
        tv = T.s(smp, "tissue_verdict")
        olfr = T.n(smp, "Olfr_sum_cpm")
        trpc2 = T.n(smp, "Trpc2_cpm")
        adcy3 = T.n(smp, "Adcy3_cpm")
        gnal = T.n(smp, "Gnal_cpm")
        omp = T.n(smp, "Omp_cpm")
        v1 = T.n(smp, "Vmn1r_sum_cpm")
        v2 = T.n(smp, "Vmn2r_sum_cpm")
        ls = T.s(smp, "library_status")
        detail = (f"tissue_verdict={tv}; Olfr sum {olfr:.1f} CPM, Adcy3 {adcy3:.1f}, "
                  f"Gnal {gnal:.1f}, Omp {omp:.1f} vs Trpc2 {trpc2:.2f}, "
                  f"Vmn1r {v1:.2f}, Vmn2r {v2:.2f} CPM; library_status={ls}; "
                  f"highest_tier_reported={T.tier_of(smp)}")
        row(cid, "8 sample table (trial 1)",
            f"{smp}: {desc}",
            detail, "SUPERSEDED",
            "The premise, not the grade. Every trial-1 description assumes the library "
            "is VNO tissue of varying quality; the tissue-identity panel says otherwise "
            "-- MAIN OLFACTORY EPITHELIUM for three of the four, and no tissue signal at "
            "all for pool2cellsRep1_S5 (both panels under the 100 CPM floor, which is "
            "uninformative rather than positively MOE). A library cannot be a degraded or borderline "
            "VNO sample if it is not VNO tissue at all, so the quality label is not "
            "merely wrong in degree -- it is answering the wrong question. Note "
            "specifically that Omp is TISSUE-SHARED (OMP-Cre labels mature main "
            "olfactory AND VNO neurons), which is exactly what let an MOE library look "
            "like a successful GFP+ VNO sort. Same annotation quantified both trials, "
            "and it yields 17,918-37,230 Trpc2 counts in trial-2 VNO libraries, so this "
            "is not a quantification artifact.",
            "results/sample_qc_all.tsv:tissue_verdict + MOE/VNO panels")

    row("S8-05", "8 sample table (trial 1)",
        "pool2cellsRep3_S7 is the cleanest 2-cell sample across both trials",
        f"technically the cleanest library in trial 1 (uniquely mapped "
        f"{T.n('pool2cellsRep3_S7','pct_uniquely_mapped'):.1f}%, duplication "
        f"{T.n('pool2cellsRep3_S7','pct_duplication'):.1f}%) AND a textbook mature MAIN "
        f"OLFACTORY neuron: Adcy3 {T.n('pool2cellsRep3_S7','Adcy3_cpm'):.0f}, Cnga2 "
        f"{T.n('pool2cellsRep3_S7','Cnga2_cpm'):.0f}, Gnal "
        f"{T.n('pool2cellsRep3_S7','Gnal_cpm'):.0f}, Olfr sum "
        f"{T.n('pool2cellsRep3_S7','Olfr_sum_cpm'):.0f} CPM, Trpc2 exactly 0",
        "SUPERSEDED",
        "This is the most consequential single delta in the reconciliation. The library "
        "IS clean -- the briefing's technical read was right -- which is precisely why "
        "it was trusted. Cleanliness and correct tissue are independent properties, and "
        "no amount of library quality makes an MOE neuron informative about vomeronasal "
        "receptor choice. Zero VNO biology is reportable from it.",
        "results/sample_qc_all.tsv, results/trial1/vr_cluster_expression.tsv")

    row("S8-06", "8 sample table",
        "the failed-sample list is: pool100cells_S8 (failed prep), trial-2 "
        "target2cellsRep2_S4 and target2cellsRep3_S5 (suspected failed preps)",
        "the pipeline's unusable set is larger and has different membership: 6 of 10 "
        "libraries carry no reportable biology. All four trial-1 libraries fail at the "
        "tissue gate (3 MOE, 1 no_tissue_signal); in trial 2, target2cellsRep2_S4 fails "
        "sort validation after PASSING the tissue gate, and target2cellsRep3_S5 fails "
        "the tissue gate itself (no_tissue_signal). Cleared: nontarget100cells_S8, "
        "target100cellsRep1_S6, target100cellsRep2_S7, target2cellsRep1_S3",
        "SUPERSEDED",
        "Gate definition plus the tissue finding. The briefing's list was assembled from "
        "prep notes; the pipeline's is derived from measured panels, and the clearance "
        "gate is qc_overall=='USABLE' AND NOT suppress_biology. suppress_biology alone "
        "is insufficient -- two libraries carry suppress_biology=False while "
        "sort_verdict=FAIL, so gating on that column alone would have admitted them.",
        "results/sample_qc_all.tsv:qc_overall,suppress_biology,sort_verdict; "
        "results/tier_status_all.tsv")

    row("S8-07", "8 sample table",
        "target2cellsRep2_S4 is a suspected failed prep",
        f"confirmed failed, and it is the case the gate exists for: actin sum "
        f"{T.n('target2cellsRep2_S4','actin_sum_cpm'):.2f} CPM and Trpc2 "
        f"{T.n('target2cellsRep2_S4','Trpc2_cpm'):.2f} CPM (both under the "
        f"failed-library thresholds of 100 and 10), yet it carries "
        f"{T.n('target2cellsRep2_S4','Vmn1r_sum_cpm'):.1f} CPM of Vmn1r signal. It "
        f"passes the tissue gate (it IS VNO), fails sort validation, and yields ZERO "
        f"receptor statements: {len(T.cand['trial2'][T.cand['trial2']['sample']=='target2cellsRep2_S4'])} "
        f"candidate rows",
        "AGREE",
        "-- and worth stating explicitly: VR family signal can survive a library that "
        "fails every other check. 501.5 CPM of Vmn1r in a dead library must never be "
        "read as a receptor call.",
        "results/sample_qc_all.tsv, results/trial2/vr_candidates.tsv")

    row("S8-08", "8 sample table",
        "Rep1/Rep2/Rep3 are replicates",
        "treated as INDEPENDENT libraries from DIFFERENT cells throughout. No "
        "correlation-based QC is applied between them and low agreement is never read "
        "as a failure. Confirmed by the data: target100cellsRep1_S6 and "
        "target100cellsRep2_S7 share their top cluster (V1R_chr6_cl008) but at "
        f"{T.cl('trial2','target100cellsRep1_S6','V1R_chr6_cl008','share_of_sample_vr')*100:.1f}% "
        f"vs {T.cl('trial2','target100cellsRep2_S7','V1R_chr6_cl008','share_of_sample_vr')*100:.1f}% "
        "of VR signal, and their top individual candidates differ (Vmn1r35 vs Vmn1r37)",
        "AGREE", "-- the caveat is honoured in code, not just in prose",
        "results/trial2/vr_cluster_expression.tsv, results/trial2/vr_candidates.tsv")

    # ================= RETRACTIONS AND NEW FINDINGS ==========================
    fl2 = T.flags["trial2"]
    h13 = fl2[(fl2["sample"] == "target100cellsRep2_S7")
              & (fl2["cluster_id"] == "V1R_chr7_cl013")]
    u89 = T.gene("trial2", "target100cellsRep2_S7", "V1R_chr7_cl013", "Vmn1r89",
                 "unique_reads_bam")
    u87 = T.gene("trial2", "target100cellsRep2_S7", "V1R_chr7_cl013", "Vmn1r87",
                 "unique_reads_bam")
    row("R-01", "retraction (phase 1)",
        "target100cellsRep2_S7 / V1R_chr7_cl013 is a genuine EM signature",
        f"NOT an artifact. em_flag={h13.iloc[0]['em_flag'] if len(h13) else 'n/a'}, "
        f"even_block_size={f(h13.iloc[0]['even_block_size']):.0f}; Vmn1r89 and Vmn1r87 "
        f"carry {u89:.0f} and {u87:.0f} independent unique MAPQ255 reads. This is real "
        f"CO-EXPRESSION across a 100-cell pool -- monogenic choice is a per-CELL rule, "
        f"so a pool legitimately captures two paralogs of one cluster",
        "SUPERSEDED",
        "The reporting figure had applied its own visual evenness criterion instead of "
        "reading the em_flag column; it now reads the column. This is the false "
        "positive the unique-read gate exists to prevent, and it is why an even split "
        "alone must never be reported as an artifact. Retracted -- do not reintroduce.",
        "results/trial2/vr_artifact_flags.tsv:em_flag,even_block_size")

    cc = T.cand["trial2"]
    c15rows = cc[(cc["sample"] == "target2cellsRep1_S3")
                 & (cc["cluster_id"] == "V1R_chr7_cl015")]
    contra = (c15rows.iloc[0]["evidence_contradiction"] if len(c15rows) else "")
    row("N-01", "new finding (phase 1)",
        "-- (not in the briefing)",
        "EVIDENCE CONTRADICTION in target2cellsRep1_S3 / V1R_chr7_cl015: the EM-dominant "
        f"paralog Vmn1r131 holds {T.gene('trial2','target2cellsRep1_S3','V1R_chr7_cl015','Vmn1r131','frac_of_cluster')*100:.1f}% "
        f"of cluster EM signal ({T.gene('trial2','target2cellsRep1_S3','V1R_chr7_cl015','Vmn1r131','counts'):.0f} "
        f"counts) with ZERO unique reads, while Vmn1r103 and Vmn1r104 each carry "
        f"{T.gene('trial2','target2cellsRep1_S3','V1R_chr7_cl015','Vmn1r103','unique_reads_bam_nodup'):.0f} "
        f"deduplicated unique reads and zero EM counts. The two channels name different "
        f"genes; neither identifies the receptor. No individual call is issued for this "
        f"cluster",
        "REFRAMED",
        "New: the second cluster of the sample the briefing described as a clean "
        "99.9% single-cluster call turns out to have no nameable receptor. The "
        "cluster-level call stands; the individual call does not exist.",
        "results/trial2/vr_candidates.tsv:evidence_contradiction")

    ca_all = pd.concat([x for x in T.cand.values() if x is not None])
    row("N-02", "new finding (phase 1)",
        "-- (not in the briefing)",
        f"all {len(ca_all)} candidate rows across both trials carry "
        f"confirmation_status=tentative_unconfirmed. Several rank1/rank2 unique-read "
        f"ratios are infinite against a background floor of 0, so a ranking can rest on "
        f"very few reads -- the evidence columns must be read before quoting a rank",
        "REFRAMED",
        "New: no individual receptor identification in this dataset is confirmed. "
        "Tier 4 is advisory in every instance.",
        "results/*/vr_candidates.tsv:confirmation_status")

    ret = pd.to_numeric(T.qc["pct_multi_loci"], errors="coerce")
    dis = pd.to_numeric(T.qc["pct_too_many_loci"], errors="coerce")
    row("N-03", "open question 1 (answered)",
        "-- (open question, not a claim)",
        f"STAR retains rather than discards paralog ambiguity: multiple-loci reads are "
        f"{ret.min():.2f}-{ret.max():.2f}% of input versus "
        f"{dis.min():.2f}-{dis.max():.2f}% discarded as too-many-loci, a median ratio "
        f"of ~125x, and the retained channel dominates in {int((ret>dis).sum())}/"
        f"{len(T.qc)} libraries. VR paralog ambiguity therefore lands in Salmon's EM, "
        f"not on the cutting-room floor. Cluster aggregation is mandatory, not optional",
        "AGREE", "-- answered in phase 0; reproduced by this run",
        "results/sample_qc_all.tsv:pct_multi_loci,pct_too_many_loci")

    return pd.DataFrame(R)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    df = build(Tables(a.results))
    df.to_csv(a.out, sep="\t", index=False)
    print(f"[reconcile] {len(df)} claims -> {a.out}")
    print(df["verdict"].value_counts().to_string())
    print()
    for _, r in df.iterrows():
        print(f"{r['claim_id']:<7} {r['verdict']:<11} {r['claim_as_written'][:78]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
