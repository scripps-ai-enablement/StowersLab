# VNO receptor analysis — results

Derived results for two sequencing trials. Produced by the `vno-receptor-rnaseq`
pipeline, **version 1.0.0**, bundled 2026-08-25.

Every table here is regenerable from `bin/vr_analyze.py` plus the raw data. These
files exist to be *read*, not as the archival copy.

## Read these first

| file | what it answers |
|---|---|
| `tier_status_all.tsv` | one row per library: how far it got, and what stopped it |
| `vr_report_all.md` | the narrative version of everything below |
| `ANALYSIS_REPORT.md` | the written analysis, with the reasoning |
| `OPEN_QUESTIONS.md` | what is still unresolved |
| `figures/vr_tier_overview.png` | the same census as a picture |

## What the 10 libraries produced

| outcome | n | libraries |
|---|---|---|
| stopped at tier 0 — wrong tissue or VNO identity unconfirmed | 5 | all four trial-1 libraries, plus `target2cellsRep3_S5` |
| stopped at tier 1 — library failed viability | 1 | `target2cellsRep2_S4` |
| reached tier 3 — cluster calls, no individual candidate | 1 | `nontarget100cells_S8` (expected: it is the negative control) |
| all five tiers reported | 3 | `target100cellsRep1_S6`, `target100cellsRep2_S7`, `target2cellsRep1_S3` |

**Trial 1 produced no receptor calls at all.** All four libraries are main
olfactory epithelium, not VNO — a tissue-identity failure, not a sequencing or
analysis failure. `figures/tissue_identity_panel.png` is the evidence. The remedy
is wet-lab; re-quantifying against a different annotation will not change it.

`target2cellsRep2_S4` is the case the QC gate exists for: it passes tissue
identity as genuine VNO, carries substantial receptor signal, and still yields
**zero** receptor statements because the library itself failed. Do not quote
receptor names from it.

Trial 2 nominated 26 individual-receptor candidate rows;
trial 1 nominated 0.

## The one rule for reading these tables

**Every individual-receptor call is `tentative_unconfirmed`.** No evidence
channel in this pipeline can confirm a single receptor identity at these read
lengths — V1R/V2R paralogs sit at 85–95% nucleotide identity, so short reads are
ambiguous within a cluster and Salmon's EM distributes them. Cluster-level calls
are the reliable tier; gene-level names below that are nominations awaiting
independent evidence.

An even within-cluster split is **not** by itself evidence of an EM artifact. In
a multi-cell pool, two paralogs each carrying independent unique-read support is
real co-expression — monogenic choice is a per-*cell* rule. See
`em_artifact_adjudication` context in `ANALYSIS_REPORT.md` and
`figures/vr_within_cluster_fractions.png`.

## Layout

```
tier_status_all.tsv        per-library census, both trials
tier_outcomes_all.tsv      every tier decision with its evidence
sample_qc_all.tsv          QC verdicts, both trials
reconciliation.tsv         cross-check against the prior project's clusters
findings_verification.tsv  each claim with the check that supports it
tier_gate_selftest.txt     the gate's own assertions, all passing
PROVENANCE_NOTES.md        how these were produced
OPEN_QUESTIONS.md          unresolved items
figures/                   6 panels
trial1/ trial2/            per-trial tables (QC, markers, clusters, flags, candidates, report)
```

Absolute paths in provenance headers were replaced with `<project>`,
`<references>` and `<user>` placeholders when this bundle was assembled; the data
columns are untouched.

## Not included

Raw FASTQ, BAMs, and the merged Salmon count matrices — too large for git and
unchanged from the nf-core run. `vr_within_cluster_fractions.tsv` (~800 KB
across both trials) is also omitted; ask if you need per-paralog detail.
