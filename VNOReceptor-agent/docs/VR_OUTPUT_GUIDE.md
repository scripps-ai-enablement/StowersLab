# Interpreting the `results/` tree — what each file tells you

Companion to `README.md` (how to run) and `docs/VR_PIPELINE_FLOWCHART.md` (what
each step does and why). This document answers a different question: **you are
looking at a file — what does it actually say about the data, and how should you
read it?**

Every column name, value vocabulary and row count below was read from the
committed tables, not recalled. Row counts are for the current two-trial run.

---

## Read in this order

Deliberately inverted from what curiosity suggests. A gene name read before its
flags is a name without error bars, and that is hard to un-remember.

| # | file | the question it answers |
|---|---|---|
| 1 | `tier_status_all.tsv` | which libraries support anything at all? |
| 2 | `<trial>/marker_cpm.tsv` | *why* did a library pass or fail? |
| 3 | `vr_report_all.md` | the narrative, already gated |
| 4 | `<trial>/vr_cluster_expression.tsv` | which receptor clusters are expressed — **the reliable answer** |
| 5 | `<trial>/vr_artifact_flags.tsv` | is the within-cluster structure trustworthy? |
| 6 | `<trial>/vr_candidates.tsv` | individual gene names — **last, always tentative** |

---

## 1. `tier_status_all.tsv` — start here

**10 rows (one per library), 15 columns.** The single most useful file in the
tree: one line per library saying how far down the claim chain its data reaches.

Key columns: `qc_overall`, `highest_tier_reported`, `highest_tier_name`,
`stopped_at_tier`, `stop_status`, `stop_reason`, and `tier0_status` …
`tier4_status`.

**How to read `highest_tier_reported`:**

| value | meaning |
|---|---|
| `-1` | nothing reportable — failed the tissue check itself |
| `0` | tissue confirmed, nothing beyond it |
| `3` | cluster-level receptor calls — **a genuinely useful result** |
| `4` | individual candidates exist (still tentative) |

**A `-1` is not a broken pipeline.** It is the pipeline refusing to answer a
question the data cannot address. `stop_reason` is a full sentence explaining
which gate closed and why; read it rather than inferring from the number.

**Do not collapse the five `tier*_status` columns to a boolean.** `PASS`,
`PASS_WITH_CAVEAT`, `FAIL`, `NO_DATA` and `SUPPRESSED` are five different
states, and the distinction that matters most is:

- `FAIL` — this tier was evaluated and the library did not pass it
- `SUPPRESSED` — this tier was **never evaluated**, because an upstream tier
  failed. The producer function was not called. It is not a negative result.
- `NO_DATA` — the tier ran and legitimately had nothing to report

A worked example from the current run: `target2cellsRep2_S4` reads
`tier0=PASS, tier1=FAIL, tier2/3/4=SUPPRESSED`. That library **is** genuine VNO
tissue (tier 0 passed) but is an empty library (actin 0.9 CPM, Trpc2 0.16 CPM).
It carries 501.5 CPM of `Vmn1r` signal and produces **zero** receptor
statements — the case the gate exists for.

> **Integrity check before you trust this file.** It is a *combined* table
> written only by the `finalize` stage. A stage-restricted run
> (`--stage report --trial trial2`) rewrites it with that trial alone — verified
> reproducibly: 10 rows → 6 rows, `trials: trial2` only. If the trial column is
> missing a trial you expect, run `bin/run_pipeline.sh --all --stage finalize`
> to rebuild it. The per-trial `<trial>/tier_status.tsv` files are never
> affected. Same applies to `tier_outcomes_all.tsv` (50 rows = 10 libraries × 5
> tiers).

`tier_outcomes_all.tsv` is the long-form version — one row per
(library, tier) with `blocked_by_tier` naming the upstream failure. Use it when
you need to know *which* tier blocked a suppression, not just that one did.

---

## 2. `<trial>/marker_cpm.tsv` and `sample_qc.tsv` — the *why*

`marker_cpm.tsv` is the evidence behind tiers 0–2. `sample_qc.tsv` (and the
combined `sample_qc_all.tsv`, **104 columns**) joins it to technical QC.

The columns that carry the verdicts:

| column | what it tells you |
|---|---|
| `tissue_verdict` + `tissue_reason` | VNO vs MOE, with a sentence of justification |
| `moe_panel_max_cpm` / `vno_panel_max_cpm` / `tissue_floor_cpm` | the actual comparison — panels are compared on their **maximum** member against a 100 CPM absolute floor |
| `sort_verdict` + `sort_reason` | `PASS`, `FAIL`, `FAIL_WRONG_TISSUE` |
| `library_status` + `library_reason` | `OK`, `FAILED`, `DEGENERATE` |
| `population_call` + `population_note` | `V1R_dominant`, `V2R_dominant`, `undetermined` |
| `suppress_biology` | **necessary but not sufficient** for clearance |
| `blocking_flags` / `n_blocking_flags` | the machine-readable reason list |

**Three traps in this file specifically.**

**`Omp_cpm` is not evidence of VNO.** It is tissue-shared — OMP-Cre labels
mature main-olfactory neurons as well as VNO neurons. A GFP+ sort from the wrong
tissue gives Omp-high, Trpc2-zero cells that look like a *successful* sort. This
is how four libraries in this project came to be described as failed preps when
they were main olfactory epithelium (`Olfr_sum_cpm` 48,729 — about 4.9% of one
library — with zero VR signal). If you read one number as tissue evidence, read
`tissue_verdict`, not `Omp_cpm`.

**`gnai2_gnao1_ratio` can be `inf`, and that is not an error.** It means
`Gnao1 = 0`, i.e. the V2R marker is genuinely absent — informative, not a
division bug. Use `gnai2_gnao1_ratio_str` for display. And check
`ratio_low_support` before quoting any magnitude: it fires when the raw
denominator is under 10 reads. "422.7:1" and "385:1" from the same library were
the *same measurement* with `Gnao1 = 7` reads. The call direction is robust; the
magnitude is Poisson noise.

**A large CPM in a `DEGENERATE` library means nothing.** CPM is a ratio: one
library shows 39,881 CPM of actin from only 60,744 assigned counts. Always read
`library_total_counts` alongside any CPM you plan to quote.

**Clearance gate**, for anything downstream:
`qc_overall == 'USABLE'` **AND** `suppress_biology` is false. The second
condition alone is insufficient — a library can carry `suppress_biology=False`
with `sort_verdict=FAIL`. Use the `vno_clearance()` helper rather than
re-deriving it.

`<trial>/technical_qc.tsv` holds the alignment-level metrics (mapping rates,
5'-3' bias, duplication, exonic fraction). Consult it when you want to know
*how* a library failed rather than *that* it did — `pct_unmapped_too_short`
above 50% means adapter/low-complexity content, not biology.

---

## 3. `<trial>/vr_cluster_expression.tsv` — the reliable answer

**23 columns, and exactly 100 rows per library** — so the file size tracks the
trial: 600 rows for trial 2 (6 libraries), 400 for trial 1 (4). Use the
per-library figure as the sanity check, not a total. This is the tier the method
actually delivers, and for most purposes it is the result.

Rows come in **two tiers** — filter on the `tier` column. Per library that is
**58 cluster rows + 42 supercluster rows** (trial 2: 348 + 252; trial 1:
232 + 168):

- `tier == 'cluster'` — 200kb max-gap clusters
- `tier == 'supercluster'` — 800kb, the coarser grouping

**They are not alternatives; they are two resolutions of the same data.** Do not
sum across both — you would double-count. 200kb is *not* a natural break (the
V1R inter-gene gap distribution has its density minimum near 2Mb); it was
retained as the conservative choice. Where `chr7_dual_tier_region == 1`, a
217,366bp gap splits one biological megacluster in two at the finer tier, and
the supercluster reunites it. **Report both tiers for that region.**

Key columns: `cluster_id`, `family`, `n_member_genes`,
`n_member_genes_detected`, `cpm_sum`, `share_of_sample_vr`, `is_called`,
`interpretation_context`, `suppression_reason`.

**`is_called == 1` is the filter you want** — 35 of trial 2's 600 rows, and
**0 of trial 1's 400**. Everything else is a cluster the pipeline scored and did
not call. A trial where nothing is called is not an empty file: it is 400 rows of
scored-and-rejected clusters, which is the correct output for four wrong-tissue
libraries.

**`interpretation_context` changes what a row means**, and this is the subtlety
most likely to be missed:

- `monogenic_expectation` — a 1–2 cell library. One dominant cluster is the
  expected result; two is already a question.
- `pooled_ambiguous` — a 100-cell pool. Multiple clusters are *expected*
  biology, not contamination. Monogenic choice is a **per-cell** rule; a pool
  does not obey it.

`share_of_sample_vr` is how you tell a dominant cluster from a trace one, and it
is a *relative* measure — always read it beside `cpm_sum`. The non-target control
in this run illustrates why: three clusters carry non-zero CPM
(`V1R_chr6_cl008` 0.080, `V1R_chr13_cl020` 0.040, `V1R_chr7_cl013` 0.040), and
their `share_of_sample_vr` values are a healthy-looking 0.50 / 0.25 / 0.25 — of a
sample total of **0.16 CPM**, about 625× below the 100 CPM purity ceiling. All
three correctly have `is_called == 0`. A large share of a negligible total is
still negligible.

---

## 4. `<trial>/vr_within_cluster_fractions.tsv` — per-paralog detail, handle with care

**24 columns, 538 rows per library** — one per primary-assembly VR gene (trial 2:
3,228 rows over 6 libraries; trial 1: 2,152 over 4). Every VR gene in every
scored cluster, with its EM counts and — critically — its independent
unique-read support.

The two evidence channels, side by side:

| column | channel |
|---|---|
| `counts`, `cpm`, `frac_of_cluster` | Salmon EM — **not trustworthy within a cluster** |
| `unique_reads_bam_nodup` | STAR MAPQ 255, deduplicated — the read was placed at exactly one locus |
| `unique_reads_bam` | same, with duplicates |
| `unique_reads_salmon_tx` | cross-check only; "unique" here means unique to a *transcript*, so gene-level support is understated for the 30 of 538 VR genes carrying multiple transcripts |

**The whole point of this file is that the two channels can disagree.** When
`frac_of_cluster` is high and `unique_reads_bam_nodup` is zero, the EM assigned
reads to a gene that was never independently observed. Do not read
`frac_of_cluster` on its own.

Also here: `is_pseudogene` + `pseudogene_criterion` (which rule fired, so the
call is auditable), `rank_in_cluster`, and `span_overlaps_other_vr` — genes whose
annotated spans overlap, where assignment is inherently ambiguous.

---

## 5. `<trial>/vr_artifact_flags.tsv` — read this before any gene name

**59 rows per library**, one per (library, cluster, tier) — trial 2: 354 rows,
trial 1: 236. The widest table in the tree because it exposes every input to its
own verdicts, so you can audit rather than trust.

> **The column count differs by trial, and that is informative.** Trial 2 has
> **68** columns, trial 1 only **60**. The eight absent from trial 1 are exactly
> the unique-read gate's inputs and outputs — `unique_channel_available`,
> `unique_background_floor`, `unique_support_threshold`, `block_unique_reads`,
> `n_block_members_unique_supported`, `block_unique_total`,
> `unique_channel_used`, `interpretation`. They are missing because no trial-1
> library was ever cleared for BAM reading, so that channel never ran. If you
> write code against this table, select columns by name with a default rather
> than by position, and treat their absence as "the gate did not run here" —
> not as a missing value to impute.

### `em_flag` — is the within-cluster structure real?

Exactly these five values. Trial-2 distribution:

| value | n | meaning |
|---|---|---|
| `no_signal` | 293 | nothing detected in this cluster |
| `insufficient_signal` | 29 | too few reads to test |
| `no_redistribution_signature` | 17 | tested, and the structure looks real |
| `single_paralog_only` | 8 | one gene, no ambiguity to resolve |
| `suspected_em_redistribution` | **1** | **the EM split one gene's reads across paralogs** |

**`no_redistribution_signature` is a clean result, not a warning.** It appears
17 times and means the test ran and found nothing wrong. Only
`suspected_em_redistribution` is an alarm — once, in this run.

### How that verdict is reached, and why it needs two steps

`p_uniform` tests the observed fractions against a uniform 1/k expectation by
Monte Carlo from `Multinomial(N, 1/k)` at the observed N and k. **It saturates at
2.4994e-04 = 1/4001** — the permutation floor. Quote it as `p < 2.5e-4`, never
as a point estimate. `power_vs_dominant` below 0.80 means the test could not have
detected a monogenic pattern at this depth, whatever the p-value says.

**An even split is necessary but not sufficient for an artifact**, and this is
the single most important thing in the file. Two neurons in a pool expressing two
paralogs produce *exactly the same even split* as one transcript whose reads were
divided between two similar paralogs. The fractions cannot distinguish them. The
discriminator is `n_block_members_unique_supported` against
`unique_support_threshold` (itself `max(10, 3 × unique_background_floor)`, the
cluster's own measured mismapping floor):

| unique-supported members | reading |
|---|---|
| exactly 1 | EM redistribution — one gene expressed, reads split |
| 2 or more | **genuine co-expression** across a multi-cell pool |

The worked contrast from this run:

| case | EM split | unique reads | verdict |
|---|---|---|---|
| 2-cell / `V1R_chr7_cl016` | 0.500 / 0.498 | **51 vs 0** | redistribution |
| 100-cell / `V1R_chr7_cl013` | 0.511 / 0.486 | **18,468 vs 18,571** | co-expression |

The second was initially mislabelled a "genuine EM signature" by a figure that
judged evenness itself instead of reading this column. Calling it an artifact
would have told the lab their sequencing was broken while it was working.

### `pseudogene_bleed_flag` — a separate, open question

Distribution: `no_apparent_pseudogene_expression` 177,
`pseudogene_only_cluster` 90, `no_pseudogene_members` 78,
`apparent_pseudogene_expression` **3**.

**This is not the same question as `em_flag`, and a cluster can be clean on one
and flagged on the other.** The largest instance —
`V1R_chr17_cl021`, where `Vmn1r-ps150` holds 38.0% of cluster signal with 2,824
unique reads — is **clean** on the redistribution test, so an even-split artifact
does not explain it. `pseudogene_mechanism` states plainly that EM leakage and
genuine cluster-regulatory transcription of the pseudogene locus are not
separable from this data. Treat it as an open mechanism question, not a defect.

`biological_interpretation_permitted` is False on rows belonging to suppressed
libraries — 118 of trial 2's 354, and **all 236** of trial 1's. Filter on it
before aggregating anything; in trial 1 that filter empties the table, which is
the correct result rather than a bug.

---

## 6. `<trial>/vr_candidates.tsv` — names, and their error bars

**26 rows × 30 columns** for trial 2; **0 rows** for trial 1 (no library got
that far). Read this last.

**Every row carries `confirmation_status = tentative_unconfirmed`.** Not a
hedge — a structural fact. Confirmation needs evidence that does not pass
through the EM step, which this pipeline cannot generate, so it never claims it.

`confidence` vocabulary: `moderate` (the rank-1 candidate),
`alternative_candidate` (ranks 2–3, reported so the shortlist is visible),
`unresolvable`, `no_call`.

Read these columns alongside every name:

- `rank` and `n_candidates_reported` — 2–3 candidates per cluster is the
  designed output, not a failure to decide
- `unique_reads_bam_nodup` — the evidence the rank rests on
- `top_over_background_floor` — **can be `inf`**, meaning the background was
  zero. That is a rank resting on very few reads, not infinite confidence.
- `evidence_contradiction` — when populated, the EM and unique channels name
  **different genes** and the honest output is *no call for that cluster*. In
  this run one cluster has an EM-dominant paralog with 4,452 counts and **zero**
  unique reads, while two others have 92 unique reads each and zero EM counts.
  Neither channel can name that receptor.
- `marker_consistency` — `consistent` means the family call agrees with the
  `Gnai2:Gnao1` population call. A contradiction is surfaced, never smoothed.

**Preserve asymmetries when you quote a call.** In this run, the *call* that one
cluster held a single expressed paralog is solid — 4,473 EM counts with zero
single-locus reads is not an independent observation. But *which* paralog was the
source rests on 51 unique reads, ~1.1% of that gene's apparent expression. Both
halves belong in the same sentence.

---

## 7. The figures

### `tissue_identity_panel.png`
The most consequential figure in the set. VNO-specific markers (blue) against
main-olfactory markers (orange) per library, log CPM, with the 100 CPM floor
drawn and `Omp` explicitly labelled as shared/non-discriminating. **What to look
for:** which colour sits above the floor. Trial-1 libraries show the orange
panel high and the blue panel at zero — wrong tissue, unmistakably, and the
remedy is wet-lab rather than bioinformatic.

### `marker_cpm_by_sample.png`
Tier 1–2 evidence: Trpc2, actin, `Gnai2`/`Gnao1` per library. **What to look
for:** the target/non-target Trpc2 separation (4,292–6,155× in the clean targets
of this run against a non-target at 0.2398 CPM). A target sitting near the
non-target level is a failed sort.

### `technical_qc_panel.png`
Mapping rates, 5'-3' bias, duplication, exonic fraction. **What to look for:**
libraries with `pct_unmapped_too_short` above ~50% (adapter/low-complexity
content) and extreme 5'-3' bias. Diagnoses *how* a library failed, and it is the
figure that separates a degraded prep from a wrong-tissue one.

### `vr_tier_overview.png`
All libraries × 5 tiers, with each stopping point and reason. **What to look
for:** the shape of the dataset in one glance — how many libraries support
receptor claims at all. Show this one first when presenting.

### `vr_cluster_heatmap.png`
Libraries × VR clusters, cluster-level CPM on a log scale. Gated-off libraries
keep their rows as hatched bands carrying the tier that stopped them — **absent
data is drawn, not silently dropped**, so you cannot mistake a suppressed
library for a negative result. **What to look for:** which clusters carry signal,
and whether a 2-cell library shows the expected single dominant cluster.

### `vr_within_cluster_fractions.png`
The most information-dense panel. Per paralog: a circle for the EM fraction and
a diamond for its share of cluster unique reads, connected — **the gap between
the two symbols is the signal.** The verdict under each panel is read from
`em_flag`, not recomputed. **What to look for:** a circle far right with its
diamond at zero means the EM assigned reads to a gene never independently
observed. Two paralogs each with circle *and* diamond high is real co-expression.

---

## 8. Reports and provenance

- **`vr_report_all.md`** — the cross-trial narrative, tier-gated. The document
  to hand someone who wants the answer rather than the tables.
- **`<trial>/vr_report.md`** — per-trial equivalents.
- **`ANALYSIS_REPORT.md`** — the scientific write-up, leading with the tissue
  finding.
- **`OPEN_QUESTIONS.md`** — status of the five open bioinformatics questions,
  each grounded in a number from these tables.
- **`reconciliation.tsv`** (24 rows) — every claim in the original project
  briefing against what the pipeline computes: agree / differ / superseded, with
  the cause of each delta. Read this if you have the old numbers in your head.
- **`findings_verification.tsv`** (24 rows) — independent re-derivation of the
  headline findings.
- **`tier_gate_selftest.txt`** — 35 assertions proving the gate cannot be
  bypassed, including spy producers confirming zero invocations for a failed
  library. If this file does not say `35 passed, 0 failed`, distrust everything
  else in the tree.
- **`rerun_diff.tsv`** (36 rows) — the cold-start reproducibility check.
- **`PROVENANCE_NOTES.md`** — **read this.** It records that
  `results/vr_cluster_map.png` has no producer in `bin/` and is therefore not
  regenerable from config + data.
- **`refcheck_diff.txt`** and `refcheck_*/` — the reference layer rebuilt into a
  scratch directory and diffed, so `ref/` can be verified without disturbing it.

---

## Quick reference — five ways to misread this tree

1. **Reading `Omp` as a VNO marker.** It is tissue-shared. This one cost four
   libraries.
2. **Treating `SUPPRESSED` as a negative result.** It means never evaluated,
   because an upstream gate closed.
3. **Reading `frac_of_cluster` without `unique_reads_bam_nodup`.** The EM
   fraction alone cannot distinguish an artifact from co-expression.
4. **Quoting a CPM from a `DEGENERATE` library**, or a ratio with
   `ratio_low_support` set.
5. **Expecting one cluster in a 100-cell pool.** Monogenic choice is per-cell;
   check `interpretation_context` before calling multiplicity a problem.
