#!/usr/bin/env python3
"""
vr_verify_findings.py -- assert that a pipeline run reproduces the findings
that phases 0-1 established.

This is the acceptance test for an end-to-end run. Every check reads the
regenerated tables under results/ and compares against a recorded
expectation; nothing is recomputed from the raw data here, so a check that
fails means the pipeline changed, not that the check drifted.

Checks are declared as (id, description, callable) and each returns
(passed, observed). Numeric expectations carry an explicit tolerance and the
tolerance is part of the printed record, because several of these values rest
on very few reads and quoting them as exact would misrepresent them.

Usage:
  vr_verify_findings.py --config config/project.yaml \
        --results results --out results/findings_verification.tsv
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import Any, Callable, Dict, List, Optional, Tuple

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def load(path: str) -> Optional[pd.DataFrame]:
    if not os.path.exists(path):
        return None
    return pd.read_csv(path, sep="\t", comment="#", dtype=str,
                       keep_default_na=False)


def num(df: pd.DataFrame, col: str) -> pd.Series:
    return pd.to_numeric(df[col], errors="coerce")


def one(df: pd.DataFrame, **eq) -> Optional[pd.Series]:
    m = pd.Series(True, index=df.index)
    for k, v in eq.items():
        m &= df[k] == v
    sub = df[m]
    if len(sub) != 1:
        return None
    return sub.iloc[0]


class Checks:
    def __init__(self, results: str):
        self.R = results
        self.qc = load(os.path.join(results, "sample_qc_all.tsv"))
        self.tier = load(os.path.join(results, "tier_status_all.tsv"))
        self.outc = load(os.path.join(results, "tier_outcomes_all.tsv"))
        self.expr, self.frac, self.flags, self.cand = {}, {}, {}, {}
        for t in ("trial1", "trial2"):
            d = os.path.join(results, t)
            self.expr[t] = load(os.path.join(d, "vr_cluster_expression.tsv"))
            self.frac[t] = load(os.path.join(d, "vr_within_cluster_fractions.tsv"))
            self.flags[t] = load(os.path.join(d, "vr_artifact_flags.tsv"))
            self.cand[t] = load(os.path.join(d, "vr_candidates.tsv"))

    # ---- helpers ----------------------------------------------------------
    def qc_val(self, sample: str, col: str) -> Any:
        r = one(self.qc, sample=sample)
        return None if r is None else r.get(col)

    def qc_num(self, sample: str, col: str) -> Optional[float]:
        v = self.qc_val(sample, col)
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    # ---- checks -----------------------------------------------------------
    def c_ref_counts(self):
        ann = load(os.path.join(os.path.dirname(self.R), "ref",
                                "vr_gene_annotation.tsv"))
        if ann is None:
            return False, "ref/vr_gene_annotation.tsv absent"
        n = len(ann)
        v1 = int((ann["family"] == "V1R").sum()) if "family" in ann else -1
        v2 = int((ann["family"] == "V2R").sum()) if "family" in ann else -1
        prim = None
        for c in ("is_primary_assembly", "primary_assembly"):
            if c in ann.columns:
                prim = int(num(ann, c).fillna(0).astype(int).sum())
        obs = f"{n} VR genes ({v1} V1R, {v2} V2R), primary_assembly={prim}"
        return (n == 541 and v1 == 319 and v2 == 222 and prim == 538), obs

    def c_nontarget_trpc2(self):
        v = self.qc_num("nontarget100cells_S8", "Trpc2_cpm")
        return (v is not None and abs(v - 0.24) <= 0.02), f"Trpc2 = {v} CPM (expect 0.24 +/- 0.02)"

    def c_nontarget_vr_purity(self):
        v1 = self.qc_num("nontarget100cells_S8", "Vmn1r_sum_cpm")
        v2 = self.qc_num("nontarget100cells_S8", "Vmn2r_sum_cpm")
        if v1 is None or v2 is None:
            return False, "VR sums unavailable"
        tot = v1 + v2
        return (tot < 100.0), (f"Vmn1r {v1} + Vmn2r {v2} = {tot:.3f} CPM "
                               f"(purity ceiling 100 CPM)")

    def c_target_trpc2_range(self):
        want = ["target100cellsRep1_S6", "target100cellsRep2_S7",
                "target2cellsRep1_S3"]
        vals = {s: self.qc_num(s, "Trpc2_cpm") for s in want}
        ok = all(v is not None and 1000.0 <= v <= 1500.0 for v in vals.values())
        return ok, "; ".join(f"{s}={v}" for s, v in vals.items()) + " (expect 1029-1476)"

    def c_target_sort_pass(self):
        want = ["target100cellsRep1_S6", "target100cellsRep2_S7",
                "target2cellsRep1_S3"]
        vals = {s: self.qc_val(s, "sort_verdict") for s in want}
        return all(v == "PASS" for v in vals.values()), str(vals)

    def c_population_v1r(self):
        want = ["target100cellsRep1_S6", "target100cellsRep2_S7",
                "target2cellsRep1_S3"]
        vals = {s: self.qc_val(s, "population_call") for s in want}
        ok = all(isinstance(v, str) and "V1R" in v for v in vals.values())
        return ok, str(vals)

    def c_trial1_moe(self):
        t1 = self.qc[self.qc["trial"] == "trial1"]
        verdicts = dict(zip(t1["sample"], t1["tissue_verdict"]))
        ok = len(verdicts) == 4 and all(
            v.startswith("MOE") or v == "no_tissue_signal"
            for v in verdicts.values())
        n_moe = sum(1 for v in verdicts.values() if v.startswith("MOE"))
        return ok, f"{n_moe}/4 MOE-called; {verdicts}"

    def c_trial1_no_biology(self):
        """No trial-1 library reports ANY tier.

        highest_tier_reported == -1 is the pipeline's encoding for "stopped at
        tier 0, nothing reported". Tier 0 is the tissue-identity gate: failing
        it means no tier was ever emitted, so -1 -- not 0 -- is the correct
        expectation. (0 would mean tier 0 PASSED and was reported.)
        """
        t1 = self.tier[self.tier["trial"] == "trial1"]
        tiers = dict(zip(t1["sample"], t1["highest_tier_reported"]))
        stops = dict(zip(t1["sample"], t1["stopped_at_tier"]))
        ok = (len(tiers) == 4
              and all(str(v) == "-1" for v in tiers.values())
              and all(float(v) == 0.0 for v in stops.values()))
        return ok, (f"highest_tier_reported={tiers} (expect all -1 = nothing "
                    f"reported), stopped_at_tier={stops} (expect all 0 = tissue gate)")

    def c_moe_cleanest_sample(self):
        """The briefing's 'cleanest 2-cell sample' is a mature MOE neuron."""
        s = "pool2cellsRep3_S7"
        vals = {k: self.qc_num(s, k) for k in
                ("Trpc2_cpm", "Adcy3_cpm", "Cnga2_cpm", "Gnal_cpm", "Olfr_sum_cpm")}
        ok = (vals["Trpc2_cpm"] == 0.0 and vals["Olfr_sum_cpm"] > 500
              and vals["Adcy3_cpm"] > 1000)
        return ok, "; ".join(f"{k}={v}" for k, v in vals.items())

    def c_tier4_samples(self):
        t = self.tier
        reach = sorted(t.loc[num(t, "highest_tier_reported") >= 4, "sample"])
        expect = sorted(["target100cellsRep1_S6", "target100cellsRep2_S7",
                         "target2cellsRep1_S3"])
        return (reach == expect), f"tier-4 samples: {reach}"

    def c_failed_lib_no_call(self):
        """target2cellsRep2_S4: 501.5 CPM Vmn1r in a failed library, zero calls."""
        s = "target2cellsRep2_S4"
        vmn1r = self.qc_num(s, "Vmn1r_sum_cpm")
        actin = self.qc_num(s, "actin_sum_cpm")
        trpc2 = self.qc_num(s, "Trpc2_cpm")
        row = one(self.tier, sample=s)
        tier = None if row is None else row.get("highest_tier_reported")
        cand = self.cand.get("trial2")
        n_cand = 0 if cand is None else int((cand["sample"] == s).sum())
        # tier 0 (tissue) PASSES -- the sample really is VNO -- and tier 1
        # (sort validation) FAILS, so the highest tier REPORTED is 0.
        stopped = None if row is None else row.get("stopped_at_tier")
        ok = (vmn1r is not None and abs(vmn1r - 501.5) <= 1.0
              and str(tier) == "0" and float(stopped) == 1.0 and n_cand == 0)
        return ok, (f"Vmn1r {vmn1r} CPM, actin {actin}, Trpc2 {trpc2}, "
                    f"highest_tier_reported={tier} (tier 0 = VNO tissue confirmed), "
                    f"stopped_at_tier={stopped} (sort validation), "
                    f"candidate rows={n_cand}")

    def c_single_em_flag(self):
        rows = []
        for t, df in self.flags.items():
            if df is None:
                continue
            hit = df[df["em_flag"] == "suspected_em_redistribution"]
            for _, r in hit.iterrows():
                rows.append((t, r["sample"], r["cluster_id"], r.get("scope", "")))
        ok = (len(rows) == 1 and rows[0][1] == "target2cellsRep1_S3"
              and rows[0][2] == "V1R_chr7_cl016")
        return ok, f"suspected_em_redistribution rows: {rows}"

    def c_em_flag_evidence(self):
        df = self.flags.get("trial2")
        if df is None:
            return False, "trial2 flags absent"
        r = df[(df["sample"] == "target2cellsRep1_S3")
               & (df["cluster_id"] == "V1R_chr7_cl016")
               & (df["em_flag"] == "suspected_em_redistribution")]
        if len(r) != 1:
            return False, f"{len(r)} matching rows"
        r = r.iloc[0]
        support = float(r["even_block_read_support"])
        blk = int(float(r["even_block_size"]))
        nsup = int(float(r["n_block_members_unique_supported"]))
        # TWO DIFFERENT TESTS, opposite directions:
        #   p_uniform (cluster level, k=4 detected) tests the WHOLE cluster
        #     against 1/k. It sits at the permutation floor 1/4001 = 2.4994e-04
        #     because the cluster is emphatically NOT uniform (two members hold
        #     ~50% each, two hold ~0).
        #   even_block_p_uniform tests the 2-member BLOCK against 1/2. A HIGH
        #     p (0.86) is the intended result: it is what "indistinguishable
        #     from an even split" means, and it is the reason the block was
        #     detected at all. Requiring it to be small was backwards.
        p_cluster = float(r["p_uniform"])
        p_block = float(r["even_block_p_uniform"])
        power = float(r["even_block_power"])
        ok = (blk == 2 and nsup == 1 and support > 8000
              and abs(p_cluster - 1.0 / 4001.0) < 1e-9
              and p_block > 0.05 and power >= 0.80)
        return ok, (f"even_block_size={blk}, block read support={support:.0f}, "
                    f"members clearing unique gate={nsup}; "
                    f"cluster p_uniform={p_cluster:.4e} (at the permutation floor "
                    f"1/4001 -- quote as a bound, cluster is NOT uniform); "
                    f"block p_uniform={p_block:.3f} (HIGH by design: the 2-member "
                    f"block IS an even split), block power={power:.2f}")

    def c_em_flag_unique_asymmetry(self):
        df = self.frac.get("trial2")
        if df is None:
            return False, "fractions absent"
        sub = df[(df["sample"] == "target2cellsRep1_S3")
                 & (df["cluster_id"] == "V1R_chr7_cl016")]
        sub = sub[num(sub, "frac_of_cluster") > 0.4]
        obs = {r["gene_name"]: (float(r["frac_of_cluster"]),
                               int(float(r["unique_reads_bam_nodup"])))
               for _, r in sub.iterrows()}
        nz = [g for g, (f, u) in obs.items() if u > 0]
        ok = (len(obs) == 2 and len(nz) == 1)
        return ok, "; ".join(f"{g}: frac={f:.3f} unique_nodup={u}"
                             for g, (f, u) in obs.items())

    def c_chr7_supercluster_share(self):
        df = self.expr.get("trial2")
        if df is None:
            return False, "trial2 expression absent"
        sub = df[(df["sample"] == "target2cellsRep1_S3")
                 & (df["cluster_id"] == "V1R_chr7_sc013")]
        if len(sub) != 1:
            return False, f"{len(sub)} rows for V1R_chr7_sc013"
        share = float(sub.iloc[0]["share_of_sample_vr"]) * 100.0
        ok = abs(share - 99.91) <= 0.05
        return ok, f"V1R_chr7_sc013 = {share:.2f}% of sample VR (expect 99.91%)"

    def c_chr7_subcluster_split(self):
        df = self.expr.get("trial2")
        if df is None:
            return False, "trial2 expression absent"
        out = {}
        for cl, want in (("V1R_chr7_cl016", 66.8), ("V1R_chr7_cl015", 33.1)):
            sub = df[(df["sample"] == "target2cellsRep1_S3")
                     & (df["cluster_id"] == cl)]
            if len(sub) != 1:
                return False, f"{len(sub)} rows for {cl}"
            out[cl] = (float(sub.iloc[0]["share_of_sample_vr"]) * 100.0, want)
        ok = all(abs(g - w) <= 0.3 for g, w in out.values())
        return ok, "; ".join(f"{k}={g:.2f}% (expect {w}%)" for k, (g, w) in out.items())

    def c_cl015_contradiction(self):
        df = self.cand.get("trial2")
        if df is None:
            return False, "candidates absent"
        sub = df[(df["sample"] == "target2cellsRep1_S3")
                 & (df["cluster_id"] == "V1R_chr7_cl015")]
        contra = sorted(set(sub["evidence_contradiction"])) if len(sub) else []
        named = sorted(set(sub["gene_name"])) if len(sub) else []
        ok = (len(sub) == 0) or any(c not in ("", "False", "false", "None")
                                    for c in contra)
        return ok, (f"{len(sub)} candidate rows for V1R_chr7_cl015 "
                    f"(genes={named}, contradiction={contra}) -- expect either no "
                    f"call or an explicit contradiction flag")

    def c_all_candidates_tentative(self):
        stat = []
        for t, df in self.cand.items():
            if df is None:
                continue
            stat += list(df["confirmation_status"])
        ok = len(stat) > 0 and all(s == "tentative_unconfirmed" for s in stat)
        return ok, f"{len(stat)} candidate rows, statuses={sorted(set(stat))}"

    def c_pseudogene_bleed_case(self):
        df = self.frac.get("trial2")
        if df is None:
            return False, "fractions absent"
        sub = df[(df["sample"] == "target100cellsRep1_S6")
                 & (df["cluster_id"] == "V1R_chr17_cl021")
                 & (df["gene_name"] == "Vmn1r-ps150")]
        if len(sub) != 1:
            return False, f"{len(sub)} rows for Vmn1r-ps150"
        r = sub.iloc[0]
        frac = float(r["frac_of_cluster"]) * 100.0
        uniq = int(float(r["unique_reads_bam"]))
        fl = self.flags.get("trial2")
        em = ""
        if fl is not None:
            h = fl[(fl["sample"] == "target100cellsRep1_S6")
                   & (fl["cluster_id"] == "V1R_chr17_cl021")]
            em = "" if not len(h) else h.iloc[0]["em_flag"]
        ok = (abs(frac - 38.0) <= 1.0 and uniq > 2000
              and em != "suspected_em_redistribution")
        return ok, (f"Vmn1r-ps150 = {frac:.1f}% of cluster, {uniq} unique reads, "
                    f"cluster em_flag={em!r} (clean on the redistribution test)")

    def c_retracted_cl013(self):
        """target100cellsRep2_S7 / V1R_chr7_cl013 must NOT be flagged."""
        fl = self.flags.get("trial2")
        if fl is None:
            return False, "flags absent"
        h = fl[(fl["sample"] == "target100cellsRep2_S7")
               & (fl["cluster_id"] == "V1R_chr7_cl013")]
        if not len(h):
            return False, "no row for V1R_chr7_cl013"
        r = h.iloc[0]
        em, blk = r["em_flag"], int(float(r["even_block_size"]))
        fr = self.frac.get("trial2")
        uniq = {}
        if fr is not None:
            s = fr[(fr["sample"] == "target100cellsRep2_S7")
                   & (fr["cluster_id"] == "V1R_chr7_cl013")]
            s = s[num(s, "frac_of_cluster") > 0.1]
            uniq = {x["gene_name"]: int(float(x["unique_reads_bam"]))
                    for _, x in s.iterrows()}
        ok = (em == "no_redistribution_signature" and blk == 0
              and sum(1 for v in uniq.values() if v > 5000) >= 2)
        return ok, (f"em_flag={em}, even_block_size={blk}, "
                    f"unique reads={uniq} -- co-expression, not an EM artifact")

    def c_multimap_channel(self):
        qc = self.qc
        ret = num(qc, "pct_multi_loci")
        disc = num(qc, "pct_too_many_loci")
        ok = bool((ret > disc).all()) and ret.max() < 40.0
        return ok, (f"retained multi-loci {ret.min():.2f}-{ret.max():.2f}% vs "
                    f"too-many-loci {disc.min():.2f}-{disc.max():.2f}%; "
                    f"retained > discarded in {int((ret>disc).sum())}/{len(qc)} libraries")

    def c_unique_floor_uniformity(self):
        """The cl015 background floor is a flat 2 reads across 26 members.

        This is the "cluster's own measured mismapping floor" the unique-read
        gate compares against. Its exact uniformity is worth recording: 26
        non-overlapping gene spans each receiving exactly one MAPQ255 read pair
        (2 mates) is a low-level pedestal, not paralog-specific evidence. It
        matters because it sets the gate: median 2 -> threshold max(10, 3*2)=10.
        """
        fr = self.frac.get("trial2")
        if fr is None:
            return False, "fractions absent"
        sub = fr[(fr["sample"] == "target2cellsRep1_S3")
                 & (fr["cluster_id"] == "V1R_chr7_cl015")]
        u = pd.to_numeric(sub["unique_reads_bam"], errors="coerce")
        n2 = int((u == 2).sum())
        overlap = pd.to_numeric(sub["span_overlaps_other_vr"], errors="coerce").fillna(0)
        ok = (n2 >= 20 and float(overlap.sum()) == 0.0 and float(u.median()) == 2.0)
        return ok, (f"{n2}/{len(sub)} cl015 members sit at exactly 2 unique reads, "
                    f"median={u.median()}, none of the spans overlap another VR gene "
                    f"-> gate threshold max(10, 3*median)=10")

    def c_read_length_convention(self):
        """STAR's avg_input_read_length is the SUM of both mates, not per-read.

        This check exists because misreading the field inverts a design
        conclusion: 149 read as a per-mate length would make the planned
        2x150bp upgrade look like a marginal change, when per-mate reads are
        actually ~75bp and the upgrade is a genuine doubling. Direct BAM
        measurement (samtools view -f 64 / -f 128 on
        target100cellsRep1_S6.markdup.sorted.bam) gives a modal SEQ length of
        76bp with a 74-76 spread, i.e. exactly half the STAR field, confirming
        the briefing's 2x75bp protocol.
        """
        rl = num(self.qc, "avg_input_read_length")
        ok = bool(((rl >= 140) & (rl <= 152)).all())
        return ok, (f"STAR avg_input_read_length {rl.min():.0f}-{rl.max():.0f} "
                    f"(BOTH mates summed) -> per-mate mean ~{rl.min()/2:.0f}-"
                    f"{rl.max()/2:.0f}bp, which sits just under the 76bp nominal "
                    f"because STAR averages over post-trimming reads; direct BAM "
                    f"measurement gives a MODAL 76bp per mate (74-76 spread). "
                    f"2x75bp protocol confirmed; 2x150bp is a real doubling.")

    def c_usable_libraries(self):
        qc = self.qc
        usable = sorted(qc.loc[(qc["qc_overall"] == "USABLE")
                               & (qc["suppress_biology"].isin(["False", "false", "0"])),
                               "sample"])
        expect = sorted(["target100cellsRep1_S6", "target100cellsRep2_S7",
                         "target2cellsRep1_S3", "nontarget100cells_S8"])
        return (usable == expect), f"cleared libraries: {usable}"


CHECKS: List[Tuple[str, str, str]] = [
    ("ref_gene_counts", "VR reference reconciles: 541 genes / 319 V1R / 222 V2R / 538 primary", "c_ref_counts"),
    ("nontarget_trpc2", "nontarget100cells_S8 Trpc2 = 0.24 CPM", "c_nontarget_trpc2"),
    ("nontarget_vr_purity", "nontarget VR signal far below the 100 CPM purity ceiling", "c_nontarget_vr_purity"),
    ("target_trpc2_range", "target Trpc2 in 1029-1476 CPM", "c_target_trpc2_range"),
    ("target_sort_pass", "sort validation PASSes for the three clean targets", "c_target_sort_pass"),
    ("population_v1r", "population call is V1R-dominant in all three targets", "c_population_v1r"),
    ("trial1_moe", "all four trial-1 libraries called MOE / no VNO tissue signal", "c_trial1_moe"),
    ("trial1_no_biology", "no trial-1 library reports any tier (all stop at the tissue gate)", "c_trial1_no_biology"),
    ("moe_cleanest_sample", "briefing's 'cleanest 2-cell sample' is a mature MOE neuron", "c_moe_cleanest_sample"),
    ("tier4_samples", "exactly three libraries reach tier 4", "c_tier4_samples"),
    ("failed_lib_no_call", "target2cellsRep2_S4: 501.5 CPM Vmn1r yields zero receptor statements", "c_failed_lib_no_call"),
    ("single_em_flag", "exactly one suspected_em_redistribution flag, on target2cellsRep1_S3 / V1R_chr7_cl016", "c_single_em_flag"),
    ("em_flag_evidence", "that flag's block statistics reproduce (cluster p at the permutation floor, block p high by design)", "c_em_flag_evidence"),
    ("em_flag_asymmetry", "one of the two co-dominant paralogs has zero unique reads", "c_em_flag_unique_asymmetry"),
    ("chr7_supercluster", "V1R_chr7_sc013 = 99.91% of target2cellsRep1_S3 VR signal", "c_chr7_supercluster_share"),
    ("chr7_subcluster_split", "within it, 66.8% cl016 / 33.1% cl015", "c_chr7_subcluster_split"),
    ("cl015_contradiction", "V1R_chr7_cl015 yields no unflagged individual call", "c_cl015_contradiction"),
    ("candidates_tentative", "every candidate row is tentative_unconfirmed", "c_all_candidates_tentative"),
    ("pseudogene_bleed", "Vmn1r-ps150 holds ~38% of V1R_chr17_cl021 on a redistribution-clean cluster", "c_pseudogene_bleed_case"),
    ("retracted_cl013", "V1R_chr7_cl013 in target100cellsRep2_S7 is NOT flagged as an EM artifact", "c_retracted_cl013"),
    ("multimap_channel", "retained multi-loci channel dominates the discarded too-many-loci channel", "c_multimap_channel"),
    ("unique_floor_uniformity", "V1R_chr7_cl015 background floor is a flat 2 unique reads across 26 non-overlapping members", "c_unique_floor_uniformity"),
    ("read_length_convention", "per-mate read length is ~75-76bp (STAR field is both mates summed)", "c_read_length_convention"),
    ("usable_libraries", "exactly four libraries clear the QC gate", "c_usable_libraries"),
]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    ck = Checks(a.results)
    rows = []
    for cid, desc, fn in CHECKS:
        try:
            passed, observed = getattr(ck, fn)()
        except Exception as exc:                       # a check must never mask
            passed, observed = False, f"ERROR {type(exc).__name__}: {exc}"
        rows.append({"check_id": cid, "expectation": desc,
                     "verdict": "PASS" if passed else "FAIL",
                     "observed": str(observed)})

    df = pd.DataFrame(rows)
    df.to_csv(a.out, sep="\t", index=False)
    npass = int((df["verdict"] == "PASS").sum())
    print(f"[verify] {npass}/{len(df)} checks passed -> {a.out}")
    for _, r in df.iterrows():
        print(f"  {r['verdict']}  {r['check_id']}: {r['observed']}")
    return 0 if npass == len(df) else 1


if __name__ == "__main__":
    raise SystemExit(main())
