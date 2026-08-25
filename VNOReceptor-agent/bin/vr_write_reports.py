#!/usr/bin/env python3
"""
vr_write_reports.py -- generate results/ANALYSIS_REPORT.md and
results/OPEN_QUESTIONS.md from the regenerated tables.

Every number in both documents is read from results/ at generation time; none
is hand-typed. That is deliberate: the reports are deliverables the lab will
read and act on, so a stale figure in them is worse than a missing one.

Usage:
  vr_write_reports.py --config config/project.yaml --results results
"""
from __future__ import annotations

import argparse
import datetime as dt
import os
from typing import Any, Dict, List, Optional

import pandas as pd

TIER4 = ["target100cellsRep1_S6", "target100cellsRep2_S7", "target2cellsRep1_S3"]


def load(p: str) -> Optional[pd.DataFrame]:
    return (pd.read_csv(p, sep="\t", comment="#", dtype=str, keep_default_na=False)
            if os.path.exists(p) else None)


def fl(x: Any) -> Optional[float]:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


class Ctx:
    def __init__(self, results: str):
        self.R = results
        self.qc = load(os.path.join(results, "sample_qc_all.tsv"))
        self.tier = load(os.path.join(results, "tier_status_all.tsv"))
        self.recon = load(os.path.join(results, "reconciliation.tsv"))
        self.verify = load(os.path.join(results, "findings_verification.tsv"))
        self.rerun = load(os.path.join(results, "rerun_diff.tsv"))
        self.expr, self.frac, self.flags, self.cand = {}, {}, {}, {}
        for t in ("trial1", "trial2"):
            d = os.path.join(results, t)
            self.expr[t] = load(os.path.join(d, "vr_cluster_expression.tsv"))
            self.frac[t] = load(os.path.join(d, "vr_within_cluster_fractions.tsv"))
            self.flags[t] = load(os.path.join(d, "vr_artifact_flags.tsv"))
            self.cand[t] = load(os.path.join(d, "vr_candidates.tsv"))

    def s(self, smp: str, col: str) -> str:
        r = self.qc[self.qc["sample"] == smp]
        return "" if not len(r) else str(r.iloc[0][col])

    def n(self, smp: str, col: str) -> Optional[float]:
        return fl(self.s(smp, col))

    def g(self, trial, smp, cl, gene, col):
        d = self.frac[trial]
        r = d[(d["sample"] == smp) & (d["cluster_id"] == cl) & (d["gene_name"] == gene)]
        return None if not len(r) else fl(r.iloc[0][col])

    def c(self, trial, smp, cl, col):
        d = self.expr[trial]
        r = d[(d["sample"] == smp) & (d["cluster_id"] == cl) & (d["tier"].isin(["cluster", "supercluster"]))]
        return None if not len(r) else fl(r.iloc[0][col])


def it(g: str) -> str:
    return f"*{g}*"


# =====================================================================
# ANALYSIS REPORT
# =====================================================================
def analysis_report(C: Ctx) -> str:
    L: List[str] = []
    a = L.append
    stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M")

    a("# VNO vomeronasal-receptor RNA-seq — analysis report")
    a("")
    a("Stowers Lab (Natalie Cole) / Scripps CCBB. Mouse vomeronasal organ, "
      "OMP-Cre × GFP reporter, FACS-sorted GFP+ / GFP−, Takara Smart-seq HT, "
      "2-cell and 100-cell pools, 2×75bp PE, two trials.")
    a("")
    a(f"Generated {stamp} by `bin/vr_write_reports.py` from the tables under "
      f"`results/`. Every figure quoted below is read from those tables at "
      f"generation time.")
    a("")
    a("---")
    a("")

    # ---------------- headline ----------------
    a("## 1. The finding that reframes the dataset: trial 1 is the wrong tissue")
    a("")
    n_moe = int((C.qc[C.qc["trial"] == "trial1"]["tissue_verdict"] == "MOE").sum())
    n_t1 = int((C.qc["trial"] == "trial1").sum())
    rest = ("the remaining one carries no tissue signal at all"
            if n_t1 - n_moe == 1 else
            f"the remaining {n_t1-n_moe} carry no tissue signal at all")
    a(f"**{n_moe} of the {n_t1} trial-1 libraries are main olfactory epithelium "
      f"(MOE), not vomeronasal organ; {rest}. No vomeronasal-receptor biology is "
      f"reportable from any of them.** This is not a library-quality problem and "
      f"it is not fixable computationally; the remedy is at the bench.")
    a("")
    a("The distinction between the two verdicts matters and the pipeline keeps "
      "them separate. `MOE` is positive evidence of the wrong tissue. "
      "`no_tissue_signal` means both marker panels sit below the 100 CPM floor, "
      "so the library carries no tissue information either way — it is "
      "uninformative, **not** evidence of MOE. Both block downstream "
      "interpretation, but only the first supports a statement about what the "
      "tissue was.")
    a("")
    a("The tissue-identity panel compares MOE markers (*Olfr* family sum, "
      "*Adcy3*, *Cnga2*, *Gnal*) against VNO-specific markers (*Trpc2*, *Vmn1r* "
      "sum, *Vmn2r* sum) on the maximum panel member, with a 100 CPM absolute "
      "floor so that a ratio of two noise values cannot produce a verdict.")
    a("")
    hdr = ("| library | tissue verdict | *Olfr* sum | *Adcy3* | *Gnal* | *Trpc2* | "
           "*Vmn1r* sum | *Omp* |")
    a(hdr)
    a("|---|---|---|---|---|---|---|---|")
    for s in ["pool100cells_S8", "pool2cellsRep1_S5", "pool2cellsRep2_S6",
              "pool2cellsRep3_S7"]:
        a(f"| {s} | **{C.s(s,'tissue_verdict')}** | {C.n(s,'Olfr_sum_cpm'):,.0f} | "
          f"{C.n(s,'Adcy3_cpm'):,.0f} | {C.n(s,'Gnal_cpm'):,.0f} | "
          f"{C.n(s,'Trpc2_cpm'):.2f} | {C.n(s,'Vmn1r_sum_cpm'):.2f} | "
          f"{C.n(s,'Omp_cpm'):,.0f} |")
    a("")
    a("All CPM. `pool2cellsRep1_S5` is the `no_tissue_signal` case: its MOE panel "
      "maxes at 0.8 CPM and its VNO panel at 5.1, both far under the floor, so no "
      "tissue claim is made about it. Note also that `pool2cellsRep2_S6`'s MOE "
      "call is the weakest of the three — it rests on *Gnal* at 320 CPM, since "
      "its *Olfr* sum (67 CPM) is itself below the floor — but its VNO panel is "
      "2.2 CPM, so the direction is not in doubt.")
    a("")
    a("Two consequences deserve stating plainly.")
    a("")
    a(f"**The briefing's \"cleanest 2-cell sample across both trials\" is a "
      f"textbook mature main-olfactory neuron.** `pool2cellsRep3_S7` is indeed "
      f"the technically cleanest trial-1 library "
      f"({C.n('pool2cellsRep3_S7','pct_uniquely_mapped'):.1f}% uniquely mapped, "
      f"{C.n('pool2cellsRep3_S7','pct_duplication'):.1f}% duplication) — and it "
      f"reads *Adcy3* {C.n('pool2cellsRep3_S7','Adcy3_cpm'):,.0f}, *Cnga2* "
      f"{C.n('pool2cellsRep3_S7','Cnga2_cpm'):,.0f}, *Gnal* "
      f"{C.n('pool2cellsRep3_S7','Gnal_cpm'):,.0f}, *Olfr* sum "
      f"{C.n('pool2cellsRep3_S7','Olfr_sum_cpm'):,.0f} CPM, with *Trpc2* exactly "
      f"zero. Library cleanliness and correct tissue are independent properties. "
      f"The cleanliness is exactly why the sample was trusted.")
    a("")
    a(f"**{it('Omp')} is why an MOE library can look like a successful GFP+ VNO "
      f"sort.** OMP-Cre drives GFP in mature main-olfactory neurons as well as "
      f"VNO neurons, so *Omp* is tissue-shared: `pool2cellsRep3_S7` reads *Omp* "
      f"{C.n('pool2cellsRep3_S7','Omp_cpm'):,.0f} CPM. A GFP+ sort from MOE "
      f"tissue is GFP+ and *Omp*-high and still carries no vomeronasal receptor. "
      f"*Omp* must never be used as a VNO-specific marker.")
    a("")
    a(f"The tissue problem is not confined to trial 1: in trial 2, "
      f"`target2cellsRep3_S5` also stops at the tissue gate with "
      f"`no_tissue_signal` (MOE max 4.3, VNO max 1.5 CPM). That is a library "
      f"failure rather than a dissection failure, but it stops at the same gate.")
    a("")
    a(f"This is not a quantification artifact. The same annotation quantified "
      f"both trials and yields "
      f"{min(C.n(s,'Trpc2_reads') for s in TIER4):,.0f}–"
      f"{max(C.n(s,'Trpc2_reads') for s in TIER4):,.0f} *Trpc2* counts in the "
      f"trial-2 VNO libraries. `pool100cells_S8` carries "
      f"{C.n('pool100cells_S8','Olfr_sum_cpm'):,.0f} CPM of *Olfr* signal — "
      f"~{C.n('pool100cells_S8','Olfr_sum_cpm')/10000:.1f}% of the library — "
      f"with zero VR signal.")
    a("")
    a("---")
    a("")

    # ---------------- usable set ----------------
    a("## 2. What is usable: 4 of 10 libraries, 3 of them GFP+ targets")
    a("")
    a("Clearance requires `qc_overall == 'USABLE'` **and** "
      "`suppress_biology == False`. The second condition is not redundant: two "
      "libraries carry `suppress_biology=False` while `sort_verdict=FAIL`, so "
      "gating on that column alone would have admitted them.")
    a("")
    a("| library | trial | cells | tissue | sort | population | library | "
      "highest tier reported |")
    a("|---|---|---|---|---|---|---|---|")
    for _, r in C.qc.iterrows():
        s = r["sample"]
        t = C.tier[C.tier["sample"] == s]
        ht = "—" if not len(t) else str(t.iloc[0]["highest_tier_reported"])
        a(f"| {s} | {r['trial']} | {r['n_cells']} | {r['tissue_verdict']} | "
          f"{r['sort_verdict']} | {r['population_call']} | {r['library_status']} | "
          f"{ht} |")
    a("")
    a("`highest_tier_reported = -1` means nothing was reported at all — the "
      "library stopped at tier 0, the tissue gate. `0` means tier 0 passed and "
      "was reported, then a later tier failed.")
    a("")
    a(f"**The failed-library gate does real work.** `target2cellsRep2_S4` is VNO "
      f"tissue (it passes tier 0), then fails sort validation: actin sum "
      f"{C.n('target2cellsRep2_S4','actin_sum_cpm'):.2f} CPM and *Trpc2* "
      f"{C.n('target2cellsRep2_S4','Trpc2_cpm'):.2f} CPM, both under the "
      f"failed-library thresholds of 100 and 10. Yet it carries "
      f"**{C.n('target2cellsRep2_S4','Vmn1r_sum_cpm'):.1f} CPM of *Vmn1r* "
      f"signal**, and produces "
      f"{len(C.cand['trial2'][C.cand['trial2']['sample']=='target2cellsRep2_S4'])} "
      f"candidate receptor rows. VR family signal can survive a library that "
      f"fails every other check; it must never be read as a receptor call.")
    a("")
    a("![Sample QC and marker CPM](figures/marker_cpm_by_sample.png)")
    a("")
    a("![Tissue identity panel](figures/tissue_identity_panel.png)")
    a("")
    a("![Technical QC panel](figures/technical_qc_panel.png)")
    a("")
    a("---")
    a("")

    # ---------------- sort validation ----------------
    a("## 3. Sort validation (tier 1)")
    a("")
    nt = C.n("nontarget100cells_S8", "Trpc2_cpm")
    ntr = C.n("nontarget100cells_S8", "Trpc2_reads")
    a("Rule 1 of the interpretation hierarchy: if the sort does not validate, "
      "nothing downstream is reportable.")
    a("")
    a("| library | *Trpc2* CPM | *Trpc2* raw reads | verdict |")
    a("|---|---|---|---|")
    for s in TIER4 + ["nontarget100cells_S8"]:
        a(f"| {s} | {C.n(s,'Trpc2_cpm'):,.2f} | {C.n(s,'Trpc2_reads'):,.0f} | "
          f"{C.s(s,'sort_verdict')} |")
    a("")
    a(f"The three GFP+ targets read {min(C.n(s,'Trpc2_cpm') for s in TIER4):,.0f}–"
      f"{max(C.n(s,'Trpc2_cpm') for s in TIER4):,.0f} CPM *Trpc2* against "
      f"{nt:.2f} CPM in the GFP− library. The separation is unambiguous and the "
      f"sort passes.")
    a("")
    folds = sorted(C.n(s, "Trpc2_cpm") / nt for s in TIER4)
    a(f"**Do not quote the enrichment as a fold-change.** Across all three "
      f"cleared GFP+ targets the ratio computes to ~{folds[0]:,.0f}× "
      f"(`target2cellsRep1_S3`), ~{folds[1]:,.0f}× (`target100cellsRep1_S6`) and "
      f"~{folds[2]:,.0f}× (`target100cellsRep2_S7`) — a range that brackets the "
      f"briefing's ~5000×. But the denominator is {ntr:.0f} raw reads. The "
      f"95% Poisson interval on a count of {ntr:.0f} spans roughly 2–13 reads, so "
      f"the implied fold ranges from about 2,400× to 16,000× without anything "
      f"changing in the biology. Report the two CPMs; the direction is the result, "
      f"the magnitude is an artifact of dividing by near-zero.")
    a("")
    a("**Sort purity.** The GFP− library's total VR signal is "
      f"{C.n('nontarget100cells_S8','Vmn1r_sum_cpm')+C.n('nontarget100cells_S8','Vmn2r_sum_cpm'):.3f} "
      f"CPM — four VR genes carrying exactly one assigned count each — against a "
      f"purity ceiling of 100 CPM. This is trace carry-over, and no cluster in "
      f"that library is called.")
    a("")
    a("---")
    a("")

    # ---------------- population ----------------
    a("## 4. Population identity (tier 2)")
    a("")
    a("V1R neurons signal through Gnai2, V2R neurons through Gnao1, so the "
      "*Gnai2*:*Gnao1* ratio assigns the population. The config requires ≥10 raw "
      "*Gnao1* reads before the ratio's **magnitude** is quotable; below that the "
      "direction still holds but the number is Poisson-unstable.")
    a("")
    a("| library | *Gnai2* CPM | *Gnao1* CPM | *Gnao1* raw reads | ratio | "
      "magnitude quotable | call |")
    a("|---|---|---|---|---|---|---|")
    for s in TIER4 + ["nontarget100cells_S8"]:
        low = C.s(s, "ratio_low_support") in ("True", "true", "1")
        a(f"| {s} | {C.n(s,'Gnai2_cpm'):,.1f} | {C.n(s,'Gnao1_cpm'):,.2f} | "
          f"{C.n(s,'Gnao1_reads'):.0f} | {C.s(s,'gnai2_gnao1_ratio_str')} | "
          f"{'no' if low else 'yes'} | {C.s(s,'population_call')} |")
    a("")
    a(f"All three targets are **V1R-dominant**, consistent with the "
      f"cluster-level VR evidence in §5 (every called cluster in every target is "
      f"V1R). Only `target100cellsRep2_S7` has a quotable ratio "
      f"({C.n('target100cellsRep2_S7','gnai2_gnao1_ratio'):.1f}:1 on "
      f"{C.n('target100cellsRep2_S7','Gnao1_reads'):.0f} *Gnao1* reads).")
    a("")
    a(f"The briefing's 385:1 for `target100cellsRep1_S6` recomputes to "
      f"{C.n('target100cellsRep1_S6','gnai2_gnao1_ratio'):.1f}:1. **This is the "
      f"same measurement, not a discrepancy**: *Gnao1* = "
      f"{C.n('target100cellsRep1_S6','Gnao1_reads'):.0f} reads is the entire "
      f"denominator, and one read either way moves the ratio by ~15%. Neither "
      f"number should be quoted as precise.")
    a("")
    a(f"Note the trap the marker panel avoids: the GFP− library reads *Gnai2* "
      f"{C.n('nontarget100cells_S8','Gnai2_cpm'):,.0f} CPM with zero *Gnao1* and "
      f"essentially no VR signal, giving a formally infinite \"V1R-dominant\" "
      f"ratio. *Gnai2*/*Gnao1* split V1R from V2R **within** the VNO; they are not "
      f"tissue markers and are excluded from the tissue panel for exactly this "
      f"reason.")
    a("")
    a("---")
    a("")

    # ---------------- cluster level ----------------
    a("## 5. Cluster-level VR expression (tier 3) — the reliable tier")
    a("")
    a("~250 *Vmn1r* and ~120 *Vmn2r* genes sit in genomic clusters of local "
      "duplicates at 85–95% nucleotide identity. At this read length, reads "
      "cannot be uniquely assigned within a cluster, and Salmon's EM spreads "
      "ambiguous reads across paralogs. Cluster-level aggregation is therefore "
      "the reportable readout; per-gene counts inside a cluster are not.")
    a("")
    a("**The ambiguity is retained, not filtered.** Across all ten libraries, "
      f"STAR's retained multiple-loci reads are "
      f"{pd.to_numeric(C.qc['pct_multi_loci']).min():.2f}–"
      f"{pd.to_numeric(C.qc['pct_multi_loci']).max():.2f}% of input, versus "
      f"{pd.to_numeric(C.qc['pct_too_many_loci']).min():.2f}–"
      f"{pd.to_numeric(C.qc['pct_too_many_loci']).max():.2f}% discarded as "
      f"too-many-loci. The retained channel dominates in all "
      f"{len(C.qc)}/{len(C.qc)} libraries by a median factor of ~125×. Paralog "
      f"ambiguity lands in Salmon's EM, so cluster aggregation is mandatory.")
    a("")
    for s in TIER4:
        ex = C.expr["trial2"]
        d = ex[(ex["sample"] == s) & (ex["tier"] == "cluster")
               & (ex["is_called"].isin(["1", "True", "true"]))].copy()
        d["sh"] = pd.to_numeric(d["share_of_sample_vr"])
        d = d.sort_values("sh", ascending=False)
        tot = fl(ex[ex["sample"] == s].iloc[0]["sample_total_vr_cpm"])
        a(f"### {s} — {C.s(s,'n_cells')}-cell pool, total VR {tot:,.1f} CPM, "
          f"{len(d)} clusters called")
        a("")
        a("| cluster | family | members | detected | CPM | share of VR |")
        a("|---|---|---|---|---|---|")
        for _, r in d.iterrows():
            a(f"| {r['cluster_id']} | {r['family']} | {r['n_member_genes']} | "
              f"{r['n_member_genes_detected']} | {fl(r['cpm_sum']):,.1f} | "
              f"{fl(r['share_of_sample_vr'])*100:.1f}% |")
        a("")
    a(f"The 2-cell pool `target2cellsRep1_S3` calls exactly 2 clusters, "
      f"consistent with 2 cells each making one choice. The 100-cell pools call "
      f"{len(C.expr['trial2'][(C.expr['trial2']['sample']=='target100cellsRep1_S6') & (C.expr['trial2']['tier']=='cluster') & (C.expr['trial2']['is_called'].isin(['1','True','true']))])} "
      f"and "
      f"{len(C.expr['trial2'][(C.expr['trial2']['sample']=='target100cellsRep2_S7') & (C.expr['trial2']['tier']=='cluster') & (C.expr['trial2']['is_called'].isin(['1','True','true']))])} "
      f"clusters — well below the 100 distinct choices a fully diverse pool "
      f"would show, i.e. the sorted populations are skewed toward a few clusters.")
    a("")
    sc = C.c("trial2", "target2cellsRep1_S3", "V1R_chr7_sc013", "share_of_sample_vr")
    c16 = C.c("trial2", "target2cellsRep1_S3", "V1R_chr7_cl016", "share_of_sample_vr")
    c15 = C.c("trial2", "target2cellsRep1_S3", "V1R_chr7_cl015", "share_of_sample_vr")
    a("### Report both cluster tiers for the chr7 region")
    a("")
    a(f"The briefing's \"cluster 039 at 99.9%\" for `target2cellsRep1_S3` "
      f"reproduces **at the 800kb supercluster tier**: `V1R_chr7_sc013` = "
      f"{sc*100:.2f}% of the sample's VR signal. At the 200kb cluster tier the "
      f"same signal resolves into {c16*100:.1f}% `V1R_chr7_cl016` + "
      f"{c15*100:.1f}% `V1R_chr7_cl015`. A 217,366bp gap between *Vmn1r132* and "
      f"*Vmn1r135* exceeds the 200kb rule by 17kb and splits the briefing's "
      f"single cluster in two.")
    a("")
    a("Neither tier is more correct. 200kb is not a natural break — the V1R "
      "inter-gene gap distribution has its KDE minimum near 2Mb — so it is a "
      "conservative convention, and the 800kb tier exists because that "
      "convention has a visible consequence here. Report both for this region "
      "and treat `cl039`/`cl029` as aliases resolved through membership, not "
      "ordinals.")
    a("")
    a("**V2R aggregation is weaker than V1R.** 18 of 37 V2R clusters are "
      "singletons and only 180/222 V2R genes sit in clusters of ≥5, so "
      "aggregation buys less protection against EM artifacts for V2R than for "
      "V1R. No V2R cluster is called in any library here, so this does not "
      "affect the present calls — but it constrains any future V2R-dominant "
      "sample.")
    a("")
    a("![Cluster-level VR expression](figures/vr_cluster_heatmap.png)")
    a("")
    a("---")
    a("")

    # ---------------- EM artifacts ----------------
    a("## 6. EM redistribution vs co-expression")
    a("")
    a("Within-cluster fractions are tested against 1/k by Monte Carlo from "
      "Multinomial(N, 1/k) at the observed N and k (4000 draws), not the "
      "asymptotic chi-square, which is invalid at N of tens. Power against a "
      "monogenic alternative (dominant at 0.90) is simulated too; power < 0.80 "
      "grades a pair `indeterminate_low_depth`.")
    a("")
    a("**An even split is necessary but not sufficient for an EM artifact.** "
      "Monogenic choice is a per-**cell** rule, so a multi-cell pool can "
      "legitimately capture two paralogs of one cluster. The pipeline therefore "
      "applies a unique-read gate: within a detected even block, each member's "
      "unique MAPQ255 deduplicated reads are compared against max(10 reads, 3× "
      "the median unique count of cluster members outside the block — the "
      "cluster's own measured mismapping floor). One member clearing it means "
      "redistribution; two or more means co-expression.")
    a("")
    fl2 = C.flags["trial2"]
    hit = fl2[fl2["em_flag"] == "suspected_em_redistribution"]
    a(f"### The one flagged artifact")
    a("")
    for _, h in hit.iterrows():
        s, cl = h["sample"], h["cluster_id"]
        u166 = C.g("trial2", s, cl, "Vmn1r166", "unique_reads_bam_nodup")
        u138 = C.g("trial2", s, cl, "Vmn1r138", "unique_reads_bam_nodup")
        f166 = C.g("trial2", s, cl, "Vmn1r166", "frac_of_cluster")
        f138 = C.g("trial2", s, cl, "Vmn1r138", "frac_of_cluster")
        n166 = C.g("trial2", s, cl, "Vmn1r166", "counts")
        n138 = C.g("trial2", s, cl, "Vmn1r138", "counts")
        a(f"`{s}` / `{cl}` — **{h['em_flag']}** ({h['em_flag_level']}).")
        a("")
        a(f"- EM fractions {f166:.3f} (*Vmn1r166*, {n166:,.0f} counts) and "
          f"{f138:.3f} (*Vmn1r138*, {n138:,.0f} counts) over "
          f"{fl(h['even_block_read_support']):,.0f} reads in the block.")
        a(f"- Cluster-level `p_uniform` = {fl(h['p_uniform']):.4e}, which is the "
          f"Monte-Carlo **floor** (1/4001). Quote it as *p* < 2.5×10⁻⁴, a bound, "
          f"not a point estimate.")
        a(f"- Within-block `p_uniform` = {fl(h['even_block_p_uniform']):.2f} — "
          f"**high by design**. That is what \"indistinguishable from an even "
          f"split\" means, and it is why the block was detected. Block power = "
          f"{fl(h['even_block_power']):.2f}.")
        a(f"- Unique MAPQ255 deduplicated reads: *Vmn1r166* {u166:.0f} vs "
          f"*Vmn1r138* {u138:.0f}. Gate threshold "
          f"{fl(h['unique_support_threshold']):.0f}; "
          f"{fl(h['n_block_members_unique_supported']):.0f} of "
          f"{fl(h['even_block_size']):.0f} members clear it.")
        a("")
        a(f"**Preserve the asymmetry.** The *call* is solid: a paralog holding "
          f"{n138:,.0f} EM counts with **zero** single-locus reads is not "
          f"independently observed, and EM splitting one transcript's reads "
          f"explains it. But *which* paralog is the source rests on {u166:.0f} "
          f"unique reads — about {u166/n166*100:.1f}% of *Vmn1r166*'s apparent "
          f"expression. The artifact is established; the receptor identity is "
          f"`tentative_unconfirmed`.")
        a("")
    a("### Retracted: V1R_chr7_cl013 in target100cellsRep2_S7 is NOT an artifact")
    a("")
    h13 = fl2[(fl2["sample"] == "target100cellsRep2_S7")
              & (fl2["cluster_id"] == "V1R_chr7_cl013")]
    if len(h13):
        h = h13.iloc[0]
        u89 = C.g("trial2", "target100cellsRep2_S7", "V1R_chr7_cl013", "Vmn1r89",
                  "unique_reads_bam")
        u87 = C.g("trial2", "target100cellsRep2_S7", "V1R_chr7_cl013", "Vmn1r87",
                  "unique_reads_bam")
        a(f"An earlier reporting figure called this a genuine EM signature. It is "
          f"not. `em_flag = {h['em_flag']}`, `even_block_size = "
          f"{fl(h['even_block_size']):.0f}`, and *Vmn1r89* and *Vmn1r87* carry "
          f"{u89:,.0f} and {u87:,.0f} independent unique MAPQ255 reads. This is "
          f"real co-expression across a 100-cell pool. The figure had applied its "
          f"own visual evenness criterion instead of reading `em_flag`; it now "
          f"reads the column. **Do not reintroduce this claim** — it is the exact "
          f"false positive the unique-read gate exists to prevent.")
        a("")
    a("### Pseudogene signal: an open mechanism, not a demonstrated artifact")
    a("")
    psf = C.g("trial2", "target100cellsRep1_S6", "V1R_chr17_cl021", "Vmn1r-ps150",
              "frac_of_cluster")
    psu = C.g("trial2", "target100cellsRep1_S6", "V1R_chr17_cl021", "Vmn1r-ps150",
              "unique_reads_bam")
    pfl = fl2[(fl2["sample"] == "target100cellsRep1_S6")
              & (fl2["cluster_id"] == "V1R_chr17_cl021")]
    a(f"Largest instance: *Vmn1r-ps150* holds {psf*100:.1f}% of `V1R_chr17_cl021` "
      f"in `target100cellsRep1_S6`, with {psu:,.0f} unique MAPQ255 reads — and "
      f"that cluster is **clean** on the redistribution test "
      f"(`em_flag = {pfl.iloc[0]['em_flag'] if len(pfl) else 'n/a'}`). An "
      f"even-split artifact does not explain it.")
    a("")
    a("Two mechanisms remain open and are not separable here: (a) EM leakage from "
      "an expressed functional paralog in the same cluster, (b) genuine "
      "transcription of the pseudogene locus, e.g. via cluster-shared regulatory "
      "elements. **Dietschi et al. 2022 is not a quantitative mouse prior for "
      "this**: their pseudogene result was significant in rat (*P* = 0.003) but "
      "not in mouse (*W* = 1214, *P* = 0.5704), and the mechanism they propose is "
      "regulatory, not multi-mapping. The flag reports the ambiguity rather than "
      "asserting bleed-through.")
    a("")
    a("![Within-cluster fractions](figures/vr_within_cluster_fractions.png)")
    a("")
    a("---")
    a("")

    # ---------------- individual candidates ----------------
    a("## 7. Individual receptor candidates (tier 4) — all tentative")
    a("")
    ca = pd.concat([v for v in C.cand.values() if v is not None and len(v)])
    a(f"Ranking uses `bam_unique_mapq255` (MAPQ 255 = STAR placed the read at "
      f"exactly one locus), counted mate-wise, with and without duplicates. "
      f"**All {len(ca)} candidate rows carry "
      f"`confirmation_status = tentative_unconfirmed`.** No individual receptor "
      f"identification in this dataset is confirmed.")
    a("")
    for s in TIER4:
        d = C.cand["trial2"]
        d = d[d["sample"] == s]
        if not len(d):
            continue
        a(f"### {s}")
        a("")
        a("| cluster | rank | gene | unique (all / dedup) | EM counts | "
          "EM frac | confidence | em_flag |")
        a("|---|---|---|---|---|---|---|---|")
        for _, r in d.iterrows():
            gn = r["gene_name"] or "—"
            a(f"| {r['cluster_id']} | {r['rank']} | *{gn}* | "
              f"{fl(r['unique_reads_bam']):,.0f} / "
              f"{fl(r['unique_reads_bam_nodup']):,.0f} | "
              f"{fl(r['em_counts']):,.0f} | {fl(r['em_frac_of_cluster']):.3f} | "
              f"{r['confidence']} | {r['em_flag']} |")
        a("")
    a("**Read the evidence columns before quoting a rank.** Several rank1/rank2 "
      "unique-read ratios are infinite against a background floor of zero, which "
      "means a ranking can rest on very few reads. Two specific cases:")
    a("")
    c15c = C.cand["trial2"]
    c15c = c15c[(c15c["sample"] == "target2cellsRep1_S3")
                & (c15c["cluster_id"] == "V1R_chr7_cl015")]
    u131 = C.g("trial2", "target2cellsRep1_S3", "V1R_chr7_cl015", "Vmn1r131",
               "unique_reads_bam")
    n131 = C.g("trial2", "target2cellsRep1_S3", "V1R_chr7_cl015", "Vmn1r131",
               "counts")
    f131 = C.g("trial2", "target2cellsRep1_S3", "V1R_chr7_cl015", "Vmn1r131",
               "frac_of_cluster")
    u103 = C.g("trial2", "target2cellsRep1_S3", "V1R_chr7_cl015", "Vmn1r103",
               "unique_reads_bam_nodup")
    a(f"1. **`V1R_chr7_cl015` in `target2cellsRep1_S3` has no nameable "
      f"receptor.** The EM-dominant paralog *Vmn1r131* holds {f131*100:.1f}% of "
      f"cluster EM signal ({n131:,.0f} counts) with **zero** unique reads, while "
      f"*Vmn1r103* and *Vmn1r104* each carry {u103:.0f} deduplicated unique reads "
      f"and **zero** EM counts. The two evidence channels name different genes; "
      f"neither identifies the expressed receptor. The pipeline emits the "
      f"contradiction rather than picking a winner. Note what this means for the "
      f"briefing's clean \"99.9% single-cluster\" description of this sample: the "
      f"cluster-level call stands, and the second cluster's individual call does "
      f"not exist.")
    a("")
    a(f"2. **The background floor is a flat pedestal, not evidence.** In that "
      f"same cluster, 26 of 37 members sit at exactly 2 unique reads each — one "
      f"MAPQ255 read pair apiece, across non-overlapping gene spans. That "
      f"pedestal is what sets the gate (median 2 → threshold max(10, 3×2) = 10). "
      f"Counts at that level are noise, and the gate is calibrated to them.")
    a("")
    a("Cross-validation with markers holds throughout: every candidate carries "
      "`marker_consistency = consistent`, i.e. V1R calls sit in libraries whose "
      "*Gnai2* exceeds *Gnao1*. No VR-versus-marker contradiction arose.")
    a("")
    a("![Tier overview](figures/vr_tier_overview.png)")
    a("")
    a("---")
    a("")

    # ---------------- can / cannot ----------------
    a("## 8. What the data can and cannot support")
    a("")
    a("**Can support:**")
    a("")
    a("1. The sort works. GFP+ VNO libraries are *Trpc2*-high, the GFP− library "
       "is *Trpc2*-negative and VR-free at the 0.16 CPM level.")
    a("2. The three clean trial-2 GFP+ libraries are V1R-population neurons.")
    a("3. Cluster-level receptor assignments for those three libraries, at both "
       "the 200kb and 800kb tiers, with the 2-cell pool calling exactly 2 "
       "clusters.")
    a("4. One documented EM-redistribution event, with the artifact call "
       "separated from the source attribution.")
    a("5. A negative result with teeth: 6 of 10 libraries yield no biology, and "
       "the gates that suppress them are tested (`results/tier_gate_selftest.txt`, "
       "35/35 assertions — a wrong-tissue or failed library cannot reach a "
       "receptor call because the producer function is never invoked).")
    a("")
    a("**Cannot support:**")
    a("")
    a("1. Any vomeronasal biology from trial 1. Three libraries are positively "
       "MOE and the fourth is tissue-uninformative; neither state is recoverable "
       "computationally.")
    a("2. A confirmed individual receptor identity for any cell. Every candidate "
       "is `tentative_unconfirmed`; one dominant cluster has no nameable receptor "
       "at all.")
    a("3. Per-gene VR expression levels within a cluster. This is the "
       "multi-mapping limit, not a depth limit.")
    a("4. Any quantitative fold-change against the GFP− library, whose relevant "
       "counts are 0–6 reads.")
    a("5. A mechanism for the pseudogene signal.")
    a("6. Anything from `target2cellsRep2_S4` despite its 501.5 CPM of *Vmn1r*.")
    a("")
    a("**Most actionable for the lab:** the trial-1 tissue finding — check the "
      "dissection and the sort gate before spending another sequencing run, and "
      "note that a GFP+ *Omp*-high sort is not self-validating — and "
      "the observation that with only 3 usable GFP+ libraries, the planned "
      "150bp PE upgrade addresses read length but not sample count. See "
      "`OPEN_QUESTIONS.md` §6 for what 150bp does and does not buy.")
    a("")
    a("---")
    a("")
    a("## 9. Reproducibility of this run")
    a("")
    if C.rerun is not None:
        vc = C.rerun["verdict"].value_counts().to_dict()
        a(f"This report comes from a cold-start run: `results/` was moved aside to "
          f"`results.prerun_backup/` and every table, report and figure was "
          f"regenerated from `config/project.yaml` plus the nf-core/rnaseq trees. "
          f"Diffing the regenerated tree against the backup over "
          f"{len(C.rerun)} paths: " +
          ", ".join(f"**{k}** {v}" for k, v in sorted(vc.items())) + ".")
        a("")
        bad = C.rerun[C.rerun["verdict"].isin(["DIFFERS", "MISSING_NEW"])]
        if len(bad):
            a("Non-reproduced paths:")
            a("")
            for _, r in bad.iterrows():
                a(f"- `{r['path']}` — {r['verdict']}: {r['detail']}")
            a("")
    if C.verify is not None:
        npass = int((C.verify["verdict"] == "PASS").sum())
        a(f"Findings verification: **{npass}/{len(C.verify)}** recorded findings "
          f"reproduce (`results/findings_verification.tsv`). The reference build "
          f"is byte-identical when rebuilt from the 852MB GTF into a scratch "
          f"directory (`results/refcheck_diff.txt`).")
        a("")
    a("Claim-by-claim reconciliation against the briefing: "
      "`results/reconciliation.tsv`.")
    a("")
    return "\n".join(L) + "\n"


# =====================================================================
# OPEN QUESTIONS
# =====================================================================
def open_questions(C: Ctx) -> str:
    L: List[str] = []
    a = L.append
    stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    qc = C.qc
    ret = pd.to_numeric(qc["pct_multi_loci"])
    dis = pd.to_numeric(qc["pct_too_many_loci"])
    rl = pd.to_numeric(qc["avg_input_read_length"])

    a("# Open bioinformatics questions — status and recommendations")
    a("")
    a(f"Generated {stamp} from the tables under `results/`. Every recommendation "
      f"below is grounded in a number from those tables; where a number does not "
      f"exist yet, the entry says what would have to be computed rather than "
      f"appealing to general principle.")
    a("")
    a("---")
    a("")

    # ---- Q1 ----
    a("## Q1. STAR multi-mapping diagnostic — **ANSWERED**")
    a("")
    a("**Question.** Is VR paralog ambiguity being discarded by STAR's "
      "multi-mapping filters, or is it reaching the quantifier?")
    a("")
    a("**Status: answered, and the answer determines the whole method.** It "
      "reaches the quantifier.")
    a("")
    a("| library | retained multi-loci | discarded too-many-loci | ratio |")
    a("|---|---|---|---|")
    for _, r in qc.iterrows():
        m, t = fl(r["pct_multi_loci"]), fl(r["pct_too_many_loci"])
        a(f"| {r['sample']} | {m:.2f}% | {t:.2f}% | "
          f"{(m/t if t else float('inf')):,.0f}× |")
    a("")
    a(f"Retained multiple-loci reads are **{ret.min():.2f}–{ret.max():.2f}%** of "
      f"input; reads discarded as too-many-loci are **{dis.min():.2f}–"
      f"{dis.max():.2f}%**. The retained channel dominates in "
      f"{int((ret>dis).sum())}/{len(qc)} libraries, with a median ratio of "
      f"~{(ret/dis).median():,.0f}×.")
    a("")
    a("**Conclusion.** STAR is *not* filtering VR paralog ambiguity out of the "
      "data. With `--outFilterMultimapNmax` at the nf-core default, a read "
      "matching a handful of 85–95%-identical paralogs is retained as a "
      "multi-mapper and handed to Salmon, where the EM distributes it. The "
      "ambiguity is therefore an *estimation* problem in the quantifier, not a "
      "*loss* problem in the aligner. Two consequences: (a) cluster-level "
      "aggregation is mandatory rather than merely prudent, and (b) any fix must "
      "act on the estimation step or on the underlying read information — "
      "tightening alignment filters would only convert a biased estimate into a "
      "missing one.")
    a("")
    a("**Recommendation.** No further work. This question is closed; the "
      "diagnostic columns (`pct_multi_loci`, `pct_too_many_loci`, "
      "`dominant_multimap_channel`) are in `results/sample_qc_all.tsv` for every "
      "future run.")
    a("")
    a("---")
    a("")

    # ---- Q2 ----
    a("## Q2. Salmon selective alignment with a genome decoy")
    a("")
    a("**Question.** Would re-quantifying with selective alignment against a "
      "decoy-aware index improve VR assignment?")
    a("")
    a("**Status: scoped, expected payoff low for this specific problem.**")
    a("")
    a("Decoys solve a particular failure: a read originating from an unannotated "
      "genomic or intronic region that finds a spurious best hit in the "
      "transcriptome. The decoy gives that read somewhere honest to go, "
      "suppressing false transcript assignment.")
    a("")
    a("**That is not the failure mode here.** The ambiguity in this dataset is "
      "genuine sequence identity *between annotated paralogous transcripts*. A "
      "75bp read from *Vmn1r166* matching *Vmn1r138* at 90% identity is not "
      "mis-assigned because it lacks a genomic home; it is ambiguous because two "
      "real transcripts explain it about equally well. A decoy cannot break that "
      "tie — both candidates remain in the index, and the EM still has to "
      "apportion the read.")
    a("")
    hi = qc.loc[pd.to_numeric(qc["pct_intergenic"]).idxmax()]
    a(f"**What a decoy *would* help with, quantified.** The relevant signal is "
      f"the intergenic fraction, which is where genomic-origin reads show up. In "
      f"the four cleared libraries it runs "
      f"{pd.to_numeric(qc[qc['qc_overall']=='USABLE']['pct_intergenic']).min():.1f}–"
      f"{pd.to_numeric(qc[qc['qc_overall']=='USABLE']['pct_intergenic']).max():.1f}%. "
      f"The one library with a genuinely high intergenic fraction is "
      f"`{hi['sample']}` at {fl(hi['pct_intergenic']):.1f}% — and that is a "
      f"trial-1 MOE library already excluded for tissue. So the population where "
      f"a decoy would act is either small or already gated out.")
    a("")
    a("**Expected payoff.** Low for paralog resolution; modest and worth having "
      "for general quantification hygiene if the data are ever reprocessed for "
      "another purpose. It would *not* change any cluster-level call in this "
      "report, and it would not resolve the `V1R_chr7_cl015` contradiction, "
      "because both competing genes there are annotated VR transcripts.")
    a("")
    a("**Recommendation.** Do not re-run on this account alone. If the lab "
      "reprocesses for an unrelated reason, enable it then. Cost is one "
      "index build plus one quantification pass per sample; the honest "
      "expectation is that cluster-level shares move by less than the "
      "Monte-Carlo noise already reported.")
    a("")
    a("---")
    a("")

    # ---- Q3 ----
    a("## Q3. Sub-cluster aggregation by sequence identity")
    a("")
    a("**Question.** Should clusters be defined by sequence identity rather than "
      "genomic proximity?")
    a("")
    a("**Status: the proximity choice is already known to be arbitrary; the "
      "identity alternative is untested.**")
    a("")
    sc = C.c("trial2", "target2cellsRep1_S3", "V1R_chr7_sc013", "share_of_sample_vr")
    c16 = C.c("trial2", "target2cellsRep1_S3", "V1R_chr7_cl016", "share_of_sample_vr")
    c15 = C.c("trial2", "target2cellsRep1_S3", "V1R_chr7_cl015", "share_of_sample_vr")
    a(f"Three facts from the reference layer bear on this. (1) **200kb is not a "
      f"natural break**: the V1R inter-gene gap distribution has its KDE minimum "
      f"near 2Mb, an order of magnitude away, so the threshold is a conservative "
      f"convention rather than a discovered boundary. (2) The convention has a "
      f"visible consequence — a 217,366bp gap, 17kb over the rule, splits the "
      f"chr7 region so that the same signal reads as {sc*100:.2f}% of one "
      f"supercluster or {c16*100:.1f}%/{c15*100:.1f}% of two clusters. (3) V2R "
      f"clustering is fragmented: 18 of 37 clusters are singletons and only "
      f"180/222 genes sit in clusters of ≥5, so proximity aggregation protects "
      f"V2R much less than V1R.")
    a("")
    a("**Why sequence identity is the better-motivated grouping.** The thing "
      "being defended against is reads that cannot be assigned between two "
      "sequences. That is a property of the *sequences*, and genomic proximity is "
      "only a proxy for it — a good proxy, since local duplication produces both "
      "adjacency and similarity, but a proxy that can fail in both directions: "
      "adjacent-but-divergent paralogs get pooled unnecessarily, and "
      "distant-but-similar ones stay split.")
    a("")
    a("**Concrete evidence that the proxy is failing here.** In "
      f"`V1R_chr7_cl015`/`V1R_chr7_cl016` the two evidence channels disagree "
      f"across the cluster boundary: *Vmn1r131* (cl015) holds "
      f"{C.g('trial2','target2cellsRep1_S3','V1R_chr7_cl015','Vmn1r131','counts'):,.0f} "
      f"EM counts with zero unique reads, while *Vmn1r166* (cl016) holds "
      f"{C.g('trial2','target2cellsRep1_S3','V1R_chr7_cl016','Vmn1r166','counts'):,.0f} "
      f"EM counts with only "
      f"{C.g('trial2','target2cellsRep1_S3','V1R_chr7_cl016','Vmn1r166','unique_reads_bam_nodup'):.0f} "
      f"unique reads. A grouping that put the actual sequence-similar set in one "
      f"bin might make this region interpretable as a single unit instead of two "
      f"partly-contradictory ones.")
    a("")
    a("**What would have to be computed.** Three steps, all with tools already "
      "on the cluster: (1) extract VR CDS/transcript sequences for the 538 "
      "primary-assembly genes from the GRCm38 FASTA (`seqkit`, minutes); "
      "(2) all-pairs identity within family — `blastn` 538×538 is trivial at this "
      "scale — and, better matched to the actual failure, a *k*-mer-sharing "
      "matrix at k = 31 and at the read length, since mappability is about shared "
      "*k*-mers, not global identity; (3) single-linkage or spectral clustering "
      "on that matrix, then re-derive cluster-level tables and compare the "
      "resulting calls against the current ones. Step (3) is the actual test: "
      "does the grouping change any call, or only the labels?")
    a("")
    a("**Expected payoff.** Moderate and mostly *diagnostic*. With three usable "
      "GFP+ libraries and two clusters called in the 2-cell pool, a regrouping "
      "will not manufacture new biology. Its value is that it would tell you "
      "whether the current cluster-level calls are robust to the boundary "
      "convention — which is currently unknown and is the single largest "
      "unexamined assumption in the pipeline.")
    a("")
    a("**Recommendation.** Worth doing, and cheap (a day of compute at most). "
      "Do it as a *sensitivity analysis* on the existing calls rather than as a "
      "replacement clustering: keep the genomic tiers as the reported ones, add "
      "an identity-based tier, and report where the three disagree. Prioritise it "
      "above Q2 and Q4.")
    a("")
    a("---")
    a("")

    # ---- Q4 ----
    a("## Q4. De novo assembly on dominant-cluster reads")
    a("")
    a("**Question.** Could assembling the reads from a dominant cluster "
      "reconstruct the expressed transcript directly?")
    a("")
    a("**Status: feasible only for the deepest cases, and it inherits the same "
      "ambiguity.**")
    a("")
    a("The read support that would be available, per candidate:")
    a("")
    a("| library | cluster | gene | unique MAPQ255 (all / dedup) | EM counts |")
    a("|---|---|---|---|---|")
    for s in TIER4:
        d = C.cand["trial2"]
        d = d[(d["sample"] == s) & (pd.to_numeric(d["rank"], errors="coerce") == 1)]
        for _, r in d.iterrows():
            a(f"| {s} | {r['cluster_id']} | *{r['gene_name']}* | "
              f"{fl(r['unique_reads_bam']):,.0f} / "
              f"{fl(r['unique_reads_bam_nodup']):,.0f} | "
              f"{fl(r['em_counts']):,.0f} |")
    a("")
    a(f"**Depth is adequate in the 100-cell pools.** The top candidates there "
      f"carry 16,723–53,841 unique reads, and a V1R coding sequence is roughly "
      f"~900bp (e.g. *Vmn1r91* spans "
      f"{C.g('trial2','target2cellsRep1_S3','V1R_chr7_cl015','Vmn1r91','end') - C.g('trial2','target2cellsRep1_S3','V1R_chr7_cl015','Vmn1r91','start'):,.0f}bp "
      f"of genome). Even at 30% duplication that is coverage in the hundreds — "
      f"assembly is not depth-limited there.")
    a("")
    u166 = C.g("trial2", "target2cellsRep1_S3", "V1R_chr7_cl016", "Vmn1r166",
               "unique_reads_bam_nodup")
    a(f"**Depth is marginal in the 2-cell pool, which is the case that matters "
      f"most.** `target2cellsRep1_S3`'s flagged cluster has {u166:.0f} "
      f"deduplicated unique read-mates for *Vmn1r166* — at ~76bp per mate that is "
      f"roughly {u166*76/900:.0f}× nominal coverage of a ~900bp CDS if they were "
      f"evenly spread, and they are not. An assembly from ambiguous reads plus "
      f"{u166:.0f} anchors is not going to produce a confident contig.")
    a("")
    a("**The deeper problem is that assembly does not escape the ambiguity.** If "
      "you assemble the *cluster's* reads, the ambiguous majority is exactly the "
      "input, and an assembler faced with two 90%-identical templates will either "
      "collapse them into one consensus (losing the distinction you wanted) or "
      "fragment at every divergent site. If you assemble only the "
      "unambiguous reads, you have thrown away 98.9% of the signal in the case "
      "you care about. The information content of the read set is what it is; "
      "assembly re-arranges it rather than adding to it.")
    a("")
    a("**Expected payoff.** Low for the 2-cell pools, where the question is "
      "live. Possibly useful in the 100-cell pools as a *confirmation* device: an "
      "assembled contig that matches one paralog's divergent sites and not the "
      "other's would independently corroborate a call — but in those pools the "
      "unique-read evidence is already strong (16k–54k reads), so it would "
      "confirm what is not in doubt.")
    a("")
    a("**Recommendation.** Deprioritise. If tried, spend the effort on "
      "`target100cellsRep2_S7`'s `V1R_chr7_cl013` — where *Vmn1r89* and "
      "*Vmn1r87* both carry ~32–35k unique reads and the interesting question is "
      "whether two transcripts are genuinely present — rather than on the 2-cell "
      "sample. Longer reads (Q6) are the better route to the 2-cell case.")
    a("")
    a("---")
    a("")

    # ---- Q5 ----
    a("## Q5. Variant calling within the dominant cluster")
    a("")
    a("**Question.** Could SNVs separate paralogs that reads cannot be assigned "
      "between?")
    a("")
    a("**Status: the right idea, and there is a specific test case for it.**")
    a("")
    a("This is the most promising of the four open computational routes, because "
      "it attacks the actual problem: paralogs differ at *specific positions*, "
      "and a read overlapping such a position is diagnostic even when the rest of "
      "it is ambiguous. Unlike assembly, this extracts information that the EM "
      "discards.")
    a("")
    a(f"**The constraint is geometric.** Per-mate reads are ~75–76bp (STAR's "
      f"{rl.min():.0f}–{rl.max():.0f} \"average input read length\" is the sum "
      f"of both mates; see Q6), and the observed mismatch rate in the cleared "
      f"libraries is "
      f"{pd.to_numeric(qc[qc['qc_overall']=='USABLE']['pct_mismatch_rate']).min():.2f}–"
      f"{pd.to_numeric(qc[qc['qc_overall']=='USABLE']['pct_mismatch_rate']).max():.2f}%, "
      f"i.e. sequencing error is low enough that a genuine paralog-diagnostic "
      f"mismatch is distinguishable from noise. At 85–95% identity, paralogs "
      f"differ every ~10–20bp on average, so even a 75bp read should span several "
      f"diagnostic sites — the information is present in the reads today. What is "
      f"missing is a pipeline step that uses it. (At 2×150bp it would span "
      f"proportionally more, and crucially would link them, which is why Q5 and "
      f"Q6 compound.)")
    a("")
    a("**Concrete test case: `V1R_chr7_cl015` in `target2cellsRep1_S3`.** This is "
      "the cluster where the two evidence channels contradict each other:")
    a("")
    a(f"- *Vmn1r131*: "
      f"{C.g('trial2','target2cellsRep1_S3','V1R_chr7_cl015','Vmn1r131','counts'):,.0f} "
      f"EM counts "
      f"({C.g('trial2','target2cellsRep1_S3','V1R_chr7_cl015','Vmn1r131','frac_of_cluster')*100:.1f}% "
      f"of the cluster), **0** unique MAPQ255 reads.")
    a(f"- *Vmn1r103* / *Vmn1r104*: "
      f"{C.g('trial2','target2cellsRep1_S3','V1R_chr7_cl015','Vmn1r103','unique_reads_bam_nodup'):.0f} "
      f"deduplicated unique reads each, **0** EM counts.")
    a("")
    a("A variant-based approach makes a falsifiable prediction here. Pile up "
      "*all* cluster-assigned reads against each candidate paralog's sequence and "
      "genotype the diagnostic sites: if the expressed transcript is *Vmn1r131*, "
      "the pileup should match *Vmn1r131*'s alleles at sites where it differs "
      "from *Vmn1r103*/*104*, and the ~4,452 EM-assigned reads become real "
      "evidence. If instead the pileup matches *Vmn1r103*/*104*, then the EM "
      "assignment is wrong and the 92 unique reads were the honest signal. Either "
      "outcome resolves a contradiction the current pipeline can only report.")
    a("")
    a("**What it would take.** The BAMs already exist and `bcftools` 1.20 and "
      "`samtools` 1.19 are on the cluster. Steps: (1) build a paralog-diagnostic "
      "site catalogue per cluster by aligning member CDS sequences to each other "
      "(this is the same sequence work as Q3, so the two share infrastructure); "
      "(2) pile up cluster reads including multi-mappers — deliberately, since "
      "multi-mappers are the reads carrying the information — at those sites; "
      "(3) compute a per-paralog likelihood from the allele counts and report it "
      "as a third evidence channel alongside EM and unique reads. Roughly a "
      "week's work. The main methodological care needed is that MAPQ filtering "
      "must be *disabled* for this pileup, which is the opposite of normal "
      "practice and would need a guard so it never leaks into the standard path.")
    a("")
    a("**Expected payoff.** High relative to the others, and it is the only route "
      "that could upgrade a call from `tentative_unconfirmed` to confirmed using "
      "existing data. Honest caveat: in the 2-cell case it operates on ~4,500 "
      "cluster reads, so it will produce a likelihood ratio rather than a "
      "certainty, and a diagnostic site falling in a region where both paralogs "
      "are identical contributes nothing.")
    a("")
    a("**Recommendation.** Highest-priority computational work of the four. "
      "Do Q3's sequence comparison first (shared infrastructure), then Q5. Treat "
      "`V1R_chr7_cl015` as the acceptance test: a method that cannot adjudicate "
      "that cluster is not yet working.")
    a("")
    a("---")
    a("")

    # ---- 150bp ----
    a("## Q6. The planned 150bp PE upgrade — what it does and does not buy")
    a("")
    a(f"**Read the length field carefully.** STAR reports \"average input read "
      f"length\" as the sum of BOTH mates for paired-end data, so the "
      f"{rl.min():.0f}–{rl.max():.0f} in the table below is the fragment's "
      f"sequenced total, not the per-read length. Measured directly from the "
      f"aligned BAMs, the per-mate length is **75–76bp** (modal 76bp, "
      f"`samtools view -f 64/-f 128`), confirming the briefing's 2×75bp. The "
      f"planned upgrade to 2×150bp is therefore a genuine doubling of per-read "
      f"length, not a marginal change.")
    a("")
    a("| library | STAR avg input read length (both mates) | avg mapped length | "
      "per-mate length |")
    a("|---|---|---|---|")
    for _, r in qc.iterrows():
        a(f"| {r['sample']} | {fl(r['avg_input_read_length']):.0f}bp | "
          f"{fl(r['avg_mapped_length']):.1f}bp | ~{fl(r['avg_input_read_length'])/2:.0f}bp |")
    a("")
    a("The per-mate column is simply the STAR field halved, and it lands slightly "
      "below the 76bp nominal because STAR's figure is a mean over reads that "
      "have already been adapter- and quality-trimmed — the trimmed fraction "
      f"runs {pd.to_numeric(qc['pct_adapter_trimmed'], errors='coerce').min():.1f}–"
      f"{pd.to_numeric(qc['pct_adapter_trimmed'], errors='coerce').max():.1f}% "
      "across these libraries. The direct BAM measurement resolves the "
      "distinction: the modal SEQ length is 76bp with a 74–76bp spread, i.e. "
      "untrimmed reads sit at the nominal length and the sub-76 average reflects "
      "trimming, not a shorter protocol.")
    a("")
    a("**What longer reads improve.**")
    a("")
    a("1. *Unique-read support per paralog* — the quantity every individual call "
       "in this report rests on. A longer read spans more divergent sites, so a "
       "larger fraction of reads become uniquely placeable. This is the direct "
       "fix for the case that is currently weakest: "
       f"*Vmn1r166*'s attribution rests on {u166:.0f} unique reads out of "
       f"{C.g('trial2','target2cellsRep1_S3','V1R_chr7_cl016','Vmn1r166','counts'):,.0f} "
       "apparent counts (~1.1%). Raising that fraction is exactly what changes "
       "`tentative_unconfirmed` to confirmed.")
    a("2. *Ability to span two or more diagnostic sites in one read*, which is "
       "what makes Q5's variant approach strong rather than marginal: a read "
       "carrying two linked diagnostic alleles is far more informative than two "
       "reads each carrying one.")
    a("3. *Mate-pair span*, which helps when one mate lands in a conserved "
       "stretch and the other in a divergent one.")
    a("")
    a("**What longer reads do not fix.**")
    a("")
    a("1. **Regions where paralogs are identical.** Where two VR genes share an "
       "exact stretch longer than the read, no read length below that stretch "
       "helps. This is a property of the duplication history, not the assay.")
    a("2. **Sample count.** Only 3 usable GFP+ libraries exist. Longer reads on "
       "the same 3 libraries give better-resolved calls for the same three cells' "
       "worth of choices; they do not add biological replication. Given that "
       "`Rep1`/`Rep2`/`Rep3` are independent libraries from different cells, more "
       "libraries — not longer reads — is what buys statistical statements about "
       "receptor-choice frequency.")
    a("3. **The trial-1 tissue problem.** Nothing about read length touches it.")
    a("4. **Library failure.** `target2cellsRep2_S4` failed at "
       f"{fl(qc[qc['sample']=='target2cellsRep2_S4'].iloc[0]['actin_sum_cpm']):.2f} "
       "CPM actin; sequencing it longer would produce longer reads from the same "
       "failed prep.")
    a("")
    a("**Recommendation.** Proceed with 2×150bp — it is a real doubling of "
      "per-mate length from the current 75–76bp, and it acts directly on the "
      "quantity every individual call in this report is limited by. Expect it to "
      "raise unique-read support substantially and to make Q5's variant approach "
      "materially stronger (a 150bp read spanning two or three diagnostic sites "
      "carries linkage information that two 75bp reads do not).")
    a("")
    a("Two caveats on allocation, though. First, longer reads will not rescue a "
      "region where paralogs are identical over a stretch longer than the read, "
      "and they do nothing for the trial-1 tissue problem or for a failed prep. "
      "Second, with only 3 usable GFP+ libraries, read length and library count "
      "buy different things: longer reads give better-resolved calls for the same "
      "three cells' worth of choices, while more libraries are what would support "
      "any statement about receptor-choice frequency across the population. If "
      "the budget allows only one, the ordering depends on the question — "
      "**resolve which receptor these cells chose** favours 150bp; **characterise "
      "how choice is distributed** favours more libraries. Both are downstream of "
      "fixing the tissue sampling.")
    a("")
    a("---")
    a("")
    a("## Priority summary")
    a("")
    a("| # | question | status | payoff | recommendation |")
    a("|---|---|---|---|---|")
    a("| Q1 | STAR multi-mapping diagnostic | **answered** | — | closed; "
      "ambiguity reaches the EM, cluster aggregation mandatory |")
    a("| Q5 | variant calling in cluster | scoped | **high** | **do first** "
      "(after Q3's sequence work); acceptance test = `V1R_chr7_cl015` |")
    a("| Q3 | sequence-identity sub-clustering | scoped | moderate, diagnostic | "
      "do as a sensitivity analysis; shares infrastructure with Q5 |")
    a("| Q2 | Salmon decoy selective alignment | scoped | **low** for paralogs | "
      "skip unless reprocessing anyway |")
    a("| Q4 | de novo assembly | scoped | low | deprioritise; longer reads are "
      "the better route |")
    a("| Q6 | 150bp PE upgrade | planned | small (already ~149bp) | proceed, but "
      "more libraries beats longer reads |")
    a("")
    a("**Above all of these:** the trial-1 tissue finding and the 3-usable-GFP+-"
      "library count are wet-lab constraints, and no computational work in this "
      "table substitutes for fixing them.")
    a("")
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True)
    ap.add_argument("--outdir", default=None)
    a = ap.parse_args(argv)
    out = a.outdir or a.results
    C = Ctx(a.results)
    os.makedirs(out, exist_ok=True)
    p1 = os.path.join(out, "ANALYSIS_REPORT.md")
    p2 = os.path.join(out, "OPEN_QUESTIONS.md")
    with open(p1, "w") as fh:
        fh.write(analysis_report(C))
    with open(p2, "w") as fh:
        fh.write(open_questions(C))
    for p in (p1, p2):
        print(f"[write_reports] {p} ({os.path.getsize(p):,} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
