# VNO receptor tiered analysis report -- trial2

*Generated 2026-08-21 12:56 by `bin/vr_report.py`. Project: VNO receptor RNA-seq (Stowers Lab / Natalie Cole).*

## Provenance

### Reference

* Genome key: **mouse-ensembl-grcm38-r91**
* GTF: `<references>/genome_references/Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf`
* FASTA: `<references>/genome_references/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa`
* Chromosome naming is Ensembl (`1`, `7`, `X`), not `chr1`. This GTF has no `gene` feature rows; gene spans in the VR reference were derived by aggregating exon/transcript rows per `gene_id`.
* VR reference build: `bin/build_vr_reference.py` -> `ref/vr_gene_annotation.tsv`, `ref/vr_clusters.tsv`, `ref/vr_gene_to_cluster.tsv`; parse report in `ref/vr_gtf_parse_report.txt`.
* 541 VR genes were annotated (319 V1R, 222 V2R); 538 are on the primary assembly. Filtering `is_primary_assembly == 1` reconciles exactly to the 318 Vmn1r / 220 Vmn2r rows present in the merged quantification table.
* **Joins are on `gene_id`, never `gene_name`**: *Vmn1r-ps5* and *Vmn2r118* each map to two distinct gene_ids on different chromosomes.

### Input tables

* QC verdicts: `<work>/results/sample_qc_all.tsv`
* trial2 (Element AVITI) quantification: `<project>/results/trail2/star_salmon/salmon.merged.gene_counts.tsv`
  * VR table `cluster_expression`: `<work>/results/trial2/vr_cluster_expression.tsv`
  * VR table `within_cluster_fractions`: `<work>/results/trial2/vr_within_cluster_fractions.tsv`
  * VR table `artifact_flags`: `<work>/results/trial2/vr_artifact_flags.tsv`
  * VR table `candidates`: `<work>/results/trial2/vr_candidates.tsv`

### CPM convention

* CPM = gene count / (column sum over **all** genes in the **unscaled** `salmon.merged.gene_counts.tsv`) x 1e6.
* The scaled and TPM tables give values 2.3-2.5x lower and must NOT be used for CPM in this project.
* As recorded in the QC table header: CPM = count / (column sum over ALL genes in the merged Salmon gene-count table) * 1e6; source table = <project>/results/trail1/star_salmon/salmon.merged.gene_counts.tsv

### Thresholds in force

| threshold | value | config section | gates |
|---|---:|---|---|
| `target_trpc2_min` | 1000 | `thresholds` | Rule 1 sort validation: target Trpc2 CPM must exceed this |
| `target_trpc2_concern` | 100 | `thresholds` | Rule 1: target Trpc2 below this is a failed sort |
| `nontarget_trpc2_max` | 10 | `thresholds` | Rule 1: GFP- nontarget Trpc2 CPM must stay below this |
| `failed_lib_actin_cpm_max` | 100 | `thresholds` | Rule 4 failed-library gate (with Trpc2 below its own max) |
| `failed_lib_trpc2_cpm_max` | 10 | `thresholds` | Rule 4 failed-library gate |
| `v1r_dominant_ratio_min` | 2.0 | `thresholds` | Rule 3 population call: Gnai2:Gnao1 above this = V1R-dominant |
| `v2r_dominant_ratio_max` | 0.5 | `thresholds` | Rule 3 population call: Gnai2:Gnao1 below this = V2R-dominant |
| `nontarget_total_vr_cpm_max` | 100 | `thresholds` | sort-purity ceiling on VR signal in the GFP- library |
| `cluster_max_gap_bp` | 200000 | `thresholds` | genomic cluster definition (see cluster caveat below) |
| `tissue_panel_floor_cpm` | 100.0 | `thresholds` | tier 0: absolute CPM floor below which a marker panel is noise |
| `tissue_dominance_ratio` | 3.0 | `thresholds` | tier 0: factor by which one tissue panel must exceed the other |
| `min_assigned_counts` | 1000000 | `qc_thresholds` | library_status=DEGENERATE below this all-gene assigned-count total |
| `ratio_min_support_reads` | 10 | `qc_thresholds` | raw Gnao1 reads required before the ratio MAGNITUDE is quotable |
| `min_input_reads` | 1000000 | `qc_thresholds` | library viability floor |
| `min_uniquely_mapped_pct` | 50.0 | `qc_thresholds` | technical QC: unique mapping rate floor |

All values are read from `<work>/config/project.yaml`; none are hardcoded in the report module.

### Cluster definition and its limits

**Cluster-threshold caveat.** Clusters are defined by a maximum inter-gene gap of
200 kb. This is NOT a natural break in the data: the kernel density of V1R
inter-gene gaps has its minimum near 2 Mb, so 200 kb is a deliberately
conservative choice that splits rather than merges. Both tiers are therefore
carried in the reference tables and both should be read:

* `cluster_id` (200 kb) -- 24 V1R and 37 V2R clusters.
* `supercluster_id` (800 kb) -- merges neighbouring clusters that the 200 kb rule
  separates.

The one region where this matters most: the historical "V1R chr7 cluster 039
(60 paralogs)" is split by the 200 kb rule into `V1R_chr7_cl015` +
`V1R_chr7_cl016`, because the gap between *Vmn1r132* and *Vmn1r135* is 217,366 bp
-- 17 kb over the rule. Supercluster `V1R_chr7_sc013` reunites them. Report both
tiers for that region. Historical label "V1R chr6 cl029 (Vmn1r32-39)" maps
exactly onto `V1R_chr6_cl008` (9 genes). The prior cl0NN ordinals are aliases,
not reproducible identifiers.

**V2R aggregation is weaker protection than V1R.** 18 of 37 V2R clusters are
singletons and only 180 of 222 V2R genes sit in a cluster of >= 5, so
cluster-level aggregation absorbs less EM ambiguity for V2R than for V1R.

**Pseudogene bleed-through is mechanistically ambiguous.** A pseudogene inside a
functional cluster showing signal may reflect (a) EM redistribution of reads from
an expressed paralog, or (b) genuine shared-regulatory-element transcription.
Dietschi et al. 2022 (Sci Adv 8(46) eabn7450) report the pseudogene effect as
significant in rat (P = 0.003) but NOT in mouse (W = 1214, P = 0.5704), and their
proposed mechanism is regulatory, not multi-mapping. This report therefore flags
the ambiguity and does not claim literature support for a quantitative mouse
expectation.

### Multi-mapping context

* STAR retained 1.66-36.35% of input reads as multi-locus alignments versus 0.02-0.30% discarded as too-many-loci (median ratio 125x). STAR is **not** filtering VR paralog ambiguity; it lands in Salmon's EM. Cluster-level aggregation is therefore mandatory, not optional.
* **Even within-cluster splits in a MULTI-cell pool are expected biology, not automatically an artifact.** Monogenic receptor choice is a per-CELL rule. In a 100-cell pool, two neurons each expressing a different paralog of the same cluster produce a ~50/50 within-cluster split that is REAL co-expression across the pool. Only the unique-read channel separates that from EM redistribution: redistribution divides ONE transcript's reads, so unique (MAPQ-filtered, deduplicated) support collapses onto a single member, whereas genuine co-expression reproduces the split in unique reads too. The quantification track therefore gates its redistribution call on both a Monte-Carlo test against uniformity AND the requirement that only one co-dominant member clears the unique-read threshold. Neither this report nor the figures re-derive that judgement; both read `em_flag` from `vr_artifact_flags.tsv`. This matters for the planned stimulus-response experiment, which pools responsive cells and will make co-expression the common case.
* Equal-fraction signature: paralogs within one cluster at near-equal fractions (e.g. 33/33/33) are a CANDIDATE redistribution artifact -- one expressed paralog whose reads the EM split. It is confirmed as an artifact only when the unique-read channel fails to reproduce the split (see the preceding caveat); in a single-cell or 2-cell library, where monogenic choice makes co-expression unlikely a priori, the artifact reading is the more probable one.
* For calibration, Dietschi et al. 2022 report that only 57% of V1rD reads were unambiguously assignable under STAR `--outFilterMultimapNmax 4` with featureCounts `-M --fraction`.

### Tier ladder

A tier is reported only when **every** upstream tier passed for that sample. This is enforced by `TierGate.emit()`, which does not invoke a tier's content producer when an upstream tier failed and substitutes an explicit suppression line instead.

0. **tissue_identity** -- tissue identity (VNO vs main olfactory epithelium)
1. **sort_validation** -- sort validation and library viability
2. **population_id** -- population identification (V1R vs V2R)
3. **cluster_vr** -- cluster-level VR calls (reliable tier)
4. **individual_vr** -- individual receptor candidates (tentative)

## Per-sample QC verdicts (read from the QC layer, not re-derived)

| trial | sample | type | cells | QC | tissue | sort | population | library | suppress_biology | Trpc2 CPM | actin CPM |
|---|---|---|---:|---|---|---|---|---|---|---:|---:|
| trial2 | `nontarget100cells_S8` | nontarget | 100 | USABLE | no_tissue_signal | PASS | V1R_dominant | OK | False | 0.24 | 3899.6 |
| trial2 | `target100cellsRep1_S6` | target | 100 | USABLE | VNO_dominant_mixed | PASS | V1R_dominant | OK | False | 1293.01 | 2591.8 |
| trial2 | `target100cellsRep2_S7` | target | 100 | USABLE | VNO_dominant_mixed | PASS | V1R_dominant | OK | False | 1476.07 | 2212.0 |
| trial2 | `target2cellsRep1_S3` | target | 2 | USABLE | VNO | PASS | V1R_dominant | OK | False | 1029.19 | 936.4 |
| trial2 | `target2cellsRep2_S4` | target | 2 | UNUSABLE | VNO | FAIL | V1R_dominant | FAILED | True | 0.16 | 0.9 |
| trial2 | `target2cellsRep3_S5` | target | 2 | UNUSABLE | no_tissue_signal | FAIL | V1R_dominant | OK | False | 0.37 | 131.2 |

## Tier outcomes by sample

### `nontarget100cells_S8` (nontarget, 100 cells, prep ok)

**Highest tier reported: 3 (cluster_vr).** stopped at tier 4 (individual_vr) [NO_DATA]: the quantification track nominated no candidate for nontarget100cells_S8 (evidence_type = no_cluster_above_signal_threshold). No individual receptor is named; the cluster-level result is the statement that stands.

**Tier 0 -- tissue identity (VNO vs main olfactory epithelium)** — `PASS_WITH_CAVEAT`

tissue_verdict=no_tissue_signal -- expected for a GFP- nontarget library: absence of the VNO-specific panel is the designed outcome, so tissue identity is not a gate for this sample. It serves as the sort-specificity control, not as VNO tissue evidence.

**Tier 1 -- sort validation and library viability** — `PASS`

nontarget Trpc2 0.24 CPM < 10

**Tier 2 -- population identification (V1R vs V2R)** — `PASS_WITH_CAVEAT`

population_call=V1R_dominant; Gnai2:Gnao1 = infinity but the denominator has too few raw reads for the ratio MAGNITUDE to be quotable (Poisson-unstable). The direction of the call stands; the number does not. Gnao1 = 0 CPM; ratio reported as infinity (Gnai2 = 165.4 CPM). Not a division error -- V2R marker genuinely absent.

**Tier 3 -- cluster-level VR calls (reliable tier)** — `PASS`

0 cluster(s) called; cluster-level aggregation is the reliable readout at 75bp

Clusters called: **0** (is_called flag from the quantification track); expected 0-0 for a nontarget 100cell
* **Sort-purity check:** GFP- control, 0 cluster(s) called, total 0.00 CPM against a 100 CPM ceiling. Within the ceiling. No receptor is called from a nontarget library either way.

**Tier 4 -- individual receptor candidates (tentative)** — `NO_DATA`

> NO DATA. the quantification track nominated no candidate for nontarget100cells_S8 (evidence_type = no_cluster_above_signal_threshold). No individual receptor is named; the cluster-level result is the statement that stands.

### `target100cellsRep1_S6` (target, 100 cells, prep ok)

**Highest tier reported: 4 (individual_vr).** all tiers reported

**Tier 0 -- tissue identity (VNO vs main olfactory epithelium)** — `PASS_WITH_CAVEAT`

tissue_verdict=VNO_dominant_mixed -- VNO dominant with a real minor MOE component; treat as a sort-purity note, not a tissue failure

**Tier 1 -- sort validation and library viability** — `PASS`

target Trpc2 1293.0 CPM > 1000

**Tier 2 -- population identification (V1R vs V2R)** — `PASS_WITH_CAVEAT`

population_call=V1R_dominant; Gnai2:Gnao1 = 422.7:1 but the denominator has too few raw reads for the ratio MAGNITUDE to be quotable (Poisson-unstable). The direction of the call stands; the number does not. Gnai2:Gnao1 = 422.71 (thresholds >2 V1R, <0.5 V2R) -- LOW SUPPORT: only 7 raw Gnao1 reads, so the ratio magnitude is Poisson-unstable (+/-1 read shifts it by ~14%); the direction of the call is robust, the exact ratio is not.

**Tier 3 -- cluster-level VR calls (reliable tier)** — `PASS`

7 cluster(s) called; cluster-level aggregation is the reliable readout at 75bp

Clusters called: **7** (is_called flag from the quantification track); expected 5-25 for a target 100cell
* **Dual-tier note:** `V1R_chr7_cl019`, `V1R_chr7_cl018` are all called and all belong to supercluster `V1R_chr7_sc014` (800 kb tier). Under the 800 kb definition this is ONE genomic region, not 2 independent calls. The 200 kb split is the conservative choice, not a natural break -- both tiers stand.

| cluster | family | CPM | share of sample VR | paralogs detected / in cluster |
|---|---|---:|---:|---:|
| `V1R_chr6_cl008` | V1R | 1177.4 | 0.249 | 3 / 9 |
| `V1R_chr7_cl013` | V1R | 962.1 | 0.203 | 1 / 6 |
| `V1R_chr7_cl019` | V1R | 850.3 | 0.180 | 1 / 1 |
| `V1R_chr7_cl018` | V1R | 706.7 | 0.149 | 1 / 2 |
| `V1R_chr7_cl017` | V1R | 477.2 | 0.101 | 3 / 18 |
| `V1R_chr6_cl007` | V1R | 377.0 | 0.080 | 4 / 21 |
| `V1R_chr17_cl021` | V1R | 176.2 | 0.037 | 3 / 24 |

**No EM-redistribution artifact was found in any cluster called for this sample.**

Pseudogene signal inside functional clusters — an OPEN mechanism question, not a redistribution verdict (each of these passed the redistribution test):
- `V1R_chr17_cl021` — pseudogene signal present (Vmn1r-ps150, 38.0% of cluster signal). redistribution test: no_redistribution_signature (within-cluster fractions are structured (distinguishable from an even split); consistent with a real per-paralog call rather than EM splitting) — so this is NOT a redistribution call. **Mechanism unresolved**: EM leakage from an expressed paralog in the same cluster, or genuine transcription of the pseudogene locus (cluster-shared regulatory elements). This module does not adjudicate. Dietschi et al. 2022 is NOT a quantitative prior for mouse: their pseudogene result was significant in rat only (P=0.003), not mouse (W=1214, P=0.5704), and their proposed mechanism is regulatory, not multi-mapping.

Clusters checked and clean on both tests (reportable results, not warnings):
- `V1R_chr6_cl007` — redistribution test: no_redistribution_signature (within-cluster fractions are structured (distinguishable from an even split); consistent with a real per-paralog call rather than EM splitting)
- `V1R_chr6_cl008` — redistribution test: no_redistribution_signature (within-cluster fractions are structured (distinguishable from an even split); consistent with a real per-paralog call rather than EM splitting)
- `V1R_chr6_cl009` — redistribution test: single_paralog_only (a single member is detected, so redistribution among paralogs is not possible; this does NOT by itself confirm the paralog's identity)
- `V1R_chr7_cl013` — redistribution test: single_paralog_only (a single member is detected, so redistribution among paralogs is not possible; this does NOT by itself confirm the paralog's identity)
- `V1R_chr7_cl017` — redistribution test: no_redistribution_signature (within-cluster fractions are structured (distinguishable from an even split); consistent with a real per-paralog call rather than EM splitting)
- `V1R_chr7_cl018` — redistribution test: single_paralog_only (a single member is detected, so redistribution among paralogs is not possible; this does NOT by itself confirm the paralog's identity)
- `V1R_chr7_cl019` — redistribution test: single_paralog_only (a single member is detected, so redistribution among paralogs is not possible; this does NOT by itself confirm the paralog's identity)

**Tier 4 -- individual receptor candidates (tentative)** — `PASS_WITH_CAVEAT`

11 tentative candidate(s) across 5 cluster(s); 0 marker contradiction(s); no identity confirmed

Every row is a TENTATIVE within-cluster assignment. At 75bp, paralogs inside one cluster share 85-95% nucleotide identity and reads are not uniquely assignable; Salmon's EM distributes them. **The cluster column, not the gene column, is the defensible unit.** The confirmation status and read support below are what a reader should weigh, not the gene name alone.

| candidate | family | cluster (200 kb) | supercluster (800 kb) | EM fraction | unique reads | unique share | EM flag | confidence | confirmation | marker cross-check |
|---|---|---|---|---:|---:|---:|---|---|---|---|
| *Vmn1r35*<br>`ENSMUSG00000060699` | V1R | `V1R_chr6_cl008` (dominant) | `V1R_chr6_sc007` | 0.946 | 33670 | 0.966 | no_redistribution_signature | moderate | tentative_unconfirmed | consistent (track check) |
| *Vmn1r36*<br>`ENSMUSG00000093764` | V1R | `V1R_chr6_cl008` (dominant) | `V1R_chr6_sc007` | 0.054 | 1183 | 0.034 | no_redistribution_signature | alternative_candidate | tentative_unconfirmed | consistent (track check) |
| *Vmn1r37*<br>`ENSMUSG00000057612` | V1R | `V1R_chr6_cl008` (dominant) | `V1R_chr6_sc007` | 0.000 | 5 | 0.000 | no_redistribution_signature | alternative_candidate | tentative_unconfirmed | consistent (track check) |
| *Vmn1r89*<br>`ENSMUSG00000095629` | V1R | `V1R_chr7_cl013` | `V1R_chr7_sc011` | 1.000 | 33139 | 1.000 | single_paralog_only | moderate | tentative_unconfirmed | consistent (track check) |
| *Vmn1r86*<br>`ENSMUSG00000070816` | V1R | `V1R_chr7_cl013` | `V1R_chr7_sc011` | 0.000 | 1 | 0.000 | single_paralog_only | alternative_candidate | tentative_unconfirmed | consistent (track check) |
| *Vmn1r87*<br>`ENSMUSG00000070815` | V1R | `V1R_chr7_cl013` | `V1R_chr7_sc011` | 0.000 | 1 | 0.000 | single_paralog_only | alternative_candidate | tentative_unconfirmed | consistent (track check) |
| *Vmn1r185*<br>`ENSMUSG00000091924` | V1R | `V1R_chr7_cl019` | `V1R_chr7_sc014` | 1.000 | 29712 | 1.000 | single_paralog_only | moderate | tentative_unconfirmed | consistent (track check) |
| *Vmn1r184*<br>`ENSMUSG00000046130` | V1R | `V1R_chr7_cl018` | `V1R_chr7_sc014` | 1.000 | 24549 | 1.000 | single_paralog_only | moderate | tentative_unconfirmed | consistent (track check) |
| *Vmn1r178*<br>`ENSMUSG00000062598` | V1R | `V1R_chr7_cl017` | `V1R_chr7_sc013` | 1.000 | 25992 | 1.000 | no_redistribution_signature | moderate | tentative_unconfirmed | consistent (track check) |
| *Vmn1r183*<br>`ENSMUSG00000066723` | V1R | `V1R_chr7_cl017` | `V1R_chr7_sc013` | 0.000 | 2 | 0.000 | no_redistribution_signature | alternative_candidate | tentative_unconfirmed | consistent (track check) |
| *Vmn1r168*<br>`ENSMUSG00000074291` | V1R | `V1R_chr7_cl017` | `V1R_chr7_sc013` | 0.000 | 2 | 0.000 | no_redistribution_signature | alternative_candidate | tentative_unconfirmed | consistent (track check) |

Candidates nominated per called cluster — `V1R_chr6_cl008`: 3, `V1R_chr7_cl013`: 3, `V1R_chr7_cl017`: 3, `V1R_chr7_cl018`: 1, `V1R_chr7_cl019`: 1. The project's stated tolerance is 2-3 candidates per cell, so a cluster reporting more than one paralog is the expected outcome of 75bp reads, not a failure of the assay.

Evidence basis: `bam_unique_mapq255`. Unique-read evidence (MAPQ-filtered BAM) is stronger than an EM fraction, because it does not depend on the redistribution step that makes within-cluster calls unsafe.

Track notes: ranked by bam_unique_mapq255; unique-read ratio rank1/rank2=12996.00; background floor (median unique reads of non-top members)=0.0, top/floor=inf; EM counts shown for context only and were NOT used to rank. n_transcripts=1 / ranked by bam_unique_mapq255; unique-read ratio rank1/rank2=28.46; background floor (median unique reads of non-top members)=1.0, top/floor=33670.00; EM counts shown for context only and were NOT used to rank. n_transcripts=1 / ranked by bam_unique_mapq255; unique-read ratio rank1/rank2=33139.00; background floor (median unique reads of non-top members)=1.0, top/floor=33139.00; EM counts shown for context only and were NOT used to rank. n_transcripts=1 / ranked by bam_unique_mapq255; unique-read ratio rank1/rank2=inf; background floor (median unique reads of non-top members)=0.0, top/floor=inf; EM counts shown for context only and were NOT used to rank. n_transcripts=1 / ranked by bam_unique_mapq255; unique-read ratio rank1/rank2=inf; background floor (median unique reads of non-top members)=nan, top/floor=inf; EM counts shown for context only and were NOT used to rank. n_transcripts=1

Every candidate in this table carries confirmation status `tentative_unconfirmed`. **No individual receptor identity is confirmed by this pipeline.** Confirmation requires evidence that does not go through the EM step — longer reads, targeted amplicon sequencing of the cluster, or in-situ/immunostaining.

### `target100cellsRep2_S7` (target, 100 cells, prep ok)

**Highest tier reported: 4 (individual_vr).** all tiers reported

**Tier 0 -- tissue identity (VNO vs main olfactory epithelium)** — `PASS_WITH_CAVEAT`

tissue_verdict=VNO_dominant_mixed -- VNO dominant with a real minor MOE component; treat as a sort-purity note, not a tissue failure

**Tier 1 -- sort validation and library viability** — `PASS`

target Trpc2 1476.1 CPM > 1000

**Tier 2 -- population identification (V1R vs V2R)** — `PASS`

population_call=V1R_dominant; Gnai2:Gnao1 = 11.3:1. Gnai2:Gnao1 = 11.34 (thresholds >2 V1R, <0.5 V2R)

**Tier 3 -- cluster-level VR calls (reliable tier)** — `PASS`

9 cluster(s) called; cluster-level aggregation is the reliable readout at 75bp

Clusters called: **9** (is_called flag from the quantification track); expected 5-25 for a target 100cell
* **Dual-tier note:** `V1R_chr7_cl013`, `V1R_chr7_cl012` are all called and all belong to supercluster `V1R_chr7_sc011` (800 kb tier). Under the 800 kb definition this is ONE genomic region, not 2 independent calls. The 200 kb split is the conservative choice, not a natural break -- both tiers stand.

| cluster | family | CPM | share of sample VR | paralogs detected / in cluster |
|---|---|---:|---:|---:|
| `V1R_chr6_cl008` | V1R | 2196.6 | 0.363 | 4 / 9 |
| `V1R_chr7_cl013` | V1R | 1238.1 | 0.204 | 3 / 6 |
| `V1R_chr17_cl021` | V1R | 643.4 | 0.106 | 4 / 24 |
| `V1R_chr7_cl017` | V1R | 481.9 | 0.080 | 2 / 18 |
| `V1R_chr7_cl010` | V1R | 304.2 | 0.050 | 1 / 14 |
| `V1R_chr13_cl020` | V1R | 272.7 | 0.045 | 3 / 64 |
| `V1R_chr7_cl011` | V1R | 257.2 | 0.042 | 3 / 12 |
| `V1R_chr6_cl007` | V1R | 205.4 | 0.034 | 1 / 21 |
| `V1R_chr7_cl012` | V1R | 196.4 | 0.032 | 4 / 29 |

**No EM-redistribution artifact was found in any cluster called for this sample.**

Pseudogene signal inside functional clusters — an OPEN mechanism question, not a redistribution verdict (each of these passed the redistribution test):
- `V1R_chr7_cl011` — pseudogene signal present (Vmn1r-ps47, 22.1% of cluster signal). redistribution test: no_redistribution_signature (within-cluster fractions are structured (distinguishable from an even split); consistent with a real per-paralog call rather than EM splitting) — so this is NOT a redistribution call. **Mechanism unresolved**: EM leakage from an expressed paralog in the same cluster, or genuine transcription of the pseudogene locus (cluster-shared regulatory elements). This module does not adjudicate. Dietschi et al. 2022 is NOT a quantitative prior for mouse: their pseudogene result was significant in rat only (P=0.003), not mouse (W=1214, P=0.5704), and their proposed mechanism is regulatory, not multi-mapping.
- `V1R_chr7_cl012` — pseudogene signal present (Vmn1r-ps64, 13.1% of cluster signal). redistribution test: no_redistribution_signature (within-cluster fractions are structured (distinguishable from an even split); consistent with a real per-paralog call rather than EM splitting) — so this is NOT a redistribution call. **Mechanism unresolved**: EM leakage from an expressed paralog in the same cluster, or genuine transcription of the pseudogene locus (cluster-shared regulatory elements). This module does not adjudicate. Dietschi et al. 2022 is NOT a quantitative prior for mouse: their pseudogene result was significant in rat only (P=0.003), not mouse (W=1214, P=0.5704), and their proposed mechanism is regulatory, not multi-mapping.

Clusters checked and clean on both tests (reportable results, not warnings):
- `V1R_chr13_cl020` — redistribution test: no_redistribution_signature (within-cluster fractions are structured (distinguishable from an even split); consistent with a real per-paralog call rather than EM splitting)
- `V1R_chr17_cl021` — redistribution test: no_redistribution_signature (within-cluster fractions are structured (distinguishable from an even split); consistent with a real per-paralog call rather than EM splitting)
- `V1R_chr6_cl006` — redistribution test: no_redistribution_signature (within-cluster fractions are structured (distinguishable from an even split); consistent with a real per-paralog call rather than EM splitting)
- `V1R_chr6_cl007` — redistribution test: single_paralog_only (a single member is detected, so redistribution among paralogs is not possible; this does NOT by itself confirm the paralog's identity)
- `V1R_chr6_cl008` — redistribution test: no_redistribution_signature (within-cluster fractions are structured (distinguishable from an even split); consistent with a real per-paralog call rather than EM splitting)
- `V1R_chr7_cl010` — redistribution test: single_paralog_only (a single member is detected, so redistribution among paralogs is not possible; this does NOT by itself confirm the paralog's identity)
- `V1R_chr7_cl013` — redistribution test: no_redistribution_signature (within-cluster fractions are structured (distinguishable from an even split); consistent with a real per-paralog call rather than EM splitting)
- `V1R_chr7_cl015` — redistribution test: no_redistribution_signature (within-cluster fractions are structured (distinguishable from an even split); consistent with a real per-paralog call rather than EM splitting)
- `V1R_chr7_cl016` — redistribution test: no_redistribution_signature (within-cluster fractions are structured (distinguishable from an even split); consistent with a real per-paralog call rather than EM splitting)
- `V1R_chr7_cl017` — redistribution test: no_redistribution_signature (within-cluster fractions are structured (distinguishable from an even split); consistent with a real per-paralog call rather than EM splitting)
- `V1R_chr7_cl019` — redistribution test: single_paralog_only (a single member is detected, so redistribution among paralogs is not possible; this does NOT by itself confirm the paralog's identity)
- `V2R_chr3_cl002` — redistribution test: no_redistribution_signature (within-cluster fractions are structured (distinguishable from an even split); consistent with a real per-paralog call rather than EM splitting)
- `V2R_chr7_cl013` — redistribution test: single_paralog_only (a single member is detected, so redistribution among paralogs is not possible; this does NOT by itself confirm the paralog's identity)

**Tier 4 -- individual receptor candidates (tentative)** — `PASS_WITH_CAVEAT`

9 tentative candidate(s) across 3 cluster(s); 0 marker contradiction(s); no identity confirmed

Every row is a TENTATIVE within-cluster assignment. At 75bp, paralogs inside one cluster share 85-95% nucleotide identity and reads are not uniquely assignable; Salmon's EM distributes them. **The cluster column, not the gene column, is the defensible unit.** The confirmation status and read support below are what a reader should weigh, not the gene name alone.

| candidate | family | cluster (200 kb) | supercluster (800 kb) | EM fraction | unique reads | unique share | EM flag | confidence | confirmation | marker cross-check |
|---|---|---|---|---:|---:|---:|---|---|---|---|
| *Vmn1r37*<br>`ENSMUSG00000057612` | V1R | `V1R_chr6_cl008` (dominant) | `V1R_chr6_sc007` | 0.743 | 53841 | 0.648 | no_redistribution_signature | unresolvable | tentative_unconfirmed | consistent (track check) |
| *Vmn1r33*<br>`ENSMUSG00000059375` | V1R | `V1R_chr6_cl008` (dominant) | `V1R_chr6_sc007` | 0.255 | 29099 | 0.350 | no_redistribution_signature | alternative_candidate | tentative_unconfirmed | consistent (track check) |
| *Vmn1r35*<br>`ENSMUSG00000060699` | V1R | `V1R_chr6_cl008` (dominant) | `V1R_chr6_sc007` | 0.001 | 113 | 0.001 | no_redistribution_signature | alternative_candidate | tentative_unconfirmed | consistent (track check) |
| *Vmn1r89*<br>`ENSMUSG00000095629` | V1R | `V1R_chr7_cl013` | `V1R_chr7_sc011` | 0.511 | 35153 | 0.520 | no_redistribution_signature | unresolvable | tentative_unconfirmed | consistent (track check) |
| *Vmn1r87*<br>`ENSMUSG00000070815` | V1R | `V1R_chr7_cl013` | `V1R_chr7_sc011` | 0.486 | 32301 | 0.478 | no_redistribution_signature | alternative_candidate | tentative_unconfirmed | consistent (track check) |
| *Vmn1r85*<br>`ENSMUSG00000070817` | V1R | `V1R_chr7_cl013` | `V1R_chr7_sc011` | 0.002 | 169 | 0.002 | no_redistribution_signature | alternative_candidate | tentative_unconfirmed | consistent (track check) |
| *Vmn1r236*<br>`ENSMUSG00000054142` | V1R | `V1R_chr17_cl021` | `V1R_chr17_sc016` | 0.352 | 16723 | 0.386 | no_redistribution_signature | unresolvable | tentative_unconfirmed | consistent (track check) |
| *Vmn1r225*<br>`ENSMUSG00000043537` | V1R | `V1R_chr17_cl021` | `V1R_chr17_sc016` | 0.365 | 12671 | 0.293 | no_redistribution_signature | alternative_candidate | tentative_unconfirmed | consistent (track check) |
| *Vmn1r224*<br>`ENSMUSG00000091151` | V1R | `V1R_chr17_cl021` | `V1R_chr17_sc016` | 0.167 | 8571 | 0.198 | no_redistribution_signature | alternative_candidate | tentative_unconfirmed | consistent (track check) |

Candidates nominated per called cluster — `V1R_chr17_cl021`: 3, `V1R_chr6_cl008`: 3, `V1R_chr7_cl013`: 3. The project's stated tolerance is 2-3 candidates per cell, so a cluster reporting more than one paralog is the expected outcome of 75bp reads, not a failure of the assay.

Evidence basis: `bam_unique_mapq255`. Unique-read evidence (MAPQ-filtered BAM) is stronger than an EM fraction, because it does not depend on the redistribution step that makes within-cluster calls unsafe.

Track notes: ranked by bam_unique_mapq255; unique-read ratio rank1/rank2=1.09; background floor (median unique reads of non-top members)=2.0, top/floor=17576.50; EM counts shown for context only and were NOT used to rank. n_transcripts=1 / ranked by bam_unique_mapq255; unique-read ratio rank1/rank2=1.32; background floor (median unique reads of non-top members)=0.0, top/floor=inf; EM counts shown for context only and were NOT used to rank. n_transcripts=1 / ranked by bam_unique_mapq255; unique-read ratio rank1/rank2=1.85; background floor (median unique reads of non-top members)=9.0, top/floor=5982.33; EM counts shown for context only and were NOT used to rank. n_transcripts=1 / ranked by bam_unique_mapq255; unique-read ratio rank1/rank2=1.85; background floor (median unique reads of non-top members)=9.0, top/floor=5982.33; EM counts shown for context only and were NOT used to rank. n_transcripts=2

Every candidate in this table carries confirmation status `tentative_unconfirmed`. **No individual receptor identity is confirmed by this pipeline.** Confirmation requires evidence that does not go through the EM step — longer reads, targeted amplicon sequencing of the cluster, or in-situ/immunostaining.

### `target2cellsRep1_S3` (target, 2 cells, prep ok)

**Highest tier reported: 4 (individual_vr).** all tiers reported

**Tier 0 -- tissue identity (VNO vs main olfactory epithelium)** — `PASS`

tissue_verdict=VNO -- VNO-specific panel establishes vomeronasal identity

**Tier 1 -- sort validation and library viability** — `PASS`

target Trpc2 1029.2 CPM > 1000

**Tier 2 -- population identification (V1R vs V2R)** — `PASS_WITH_CAVEAT`

population_call=V1R_dominant; Gnai2:Gnao1 = infinity but the denominator has too few raw reads for the ratio MAGNITUDE to be quotable (Poisson-unstable). The direction of the call stands; the number does not. Gnao1 = 0 CPM; ratio reported as infinity (Gnai2 = 3602.2 CPM). Not a division error -- V2R marker genuinely absent.

**Tier 3 -- cluster-level VR calls (reliable tier)** — `PASS`

2 cluster(s) called; cluster-level aggregation is the reliable readout at 75bp

Clusters called: **2** (is_called flag from the quantification track); expected 1-2 for a target 2cell
* **Dual-tier note:** `V1R_chr7_cl016`, `V1R_chr7_cl015` are all called and all belong to supercluster `V1R_chr7_sc013` (800 kb tier). Under the 800 kb definition this is ONE genomic region, not 2 independent calls. The 200 kb split is the conservative choice, not a natural break -- both tiers stand.
* The quantification track marked `V1R_chr7_cl016`, `V1R_chr7_cl015` as lying in the chr7 dual-tier region (historical "cluster 039", 60 paralogs, split by the 200 kb rule at the 217,366 bp Vmn1r132-Vmn1r135 gap and reunited by supercluster `V1R_chr7_sc013`).

| cluster | family | CPM | share of sample VR | paralogs detected / in cluster |
|---|---|---:|---:|---:|
| `V1R_chr7_cl016` | V1R | 515.4 | 0.668 | 4 / 23 |
| `V1R_chr7_cl015` | V1R | 255.9 | 0.331 | 2 / 37 |

**EM-redistribution findings** — the within-cluster split here is an artifact of read assignment, not per-paralog biology:
- `V1R_chr7_cl016` — **suspected_em_redistribution** — level strong — even block of 2 paralogs — 8974 reads of support

Clusters checked and clean on both tests (reportable results, not warnings):
- `V1R_chr7_cl015` — redistribution test: no_redistribution_signature (within-cluster fractions are structured (distinguishable from an even split); consistent with a real per-paralog call rather than EM splitting)

**Tier 4 -- individual receptor candidates (tentative)** — `PASS_WITH_CAVEAT`

5 tentative candidate(s) across 2 cluster(s); 0 marker contradiction(s); no identity confirmed

Every row is a TENTATIVE within-cluster assignment. At 75bp, paralogs inside one cluster share 85-95% nucleotide identity and reads are not uniquely assignable; Salmon's EM distributes them. **The cluster column, not the gene column, is the defensible unit.** The confirmation status and read support below are what a reader should weigh, not the gene name alone.

| candidate | family | cluster (200 kb) | supercluster (800 kb) | EM fraction | unique reads | unique share | EM flag | confidence | confirmation | marker cross-check |
|---|---|---|---|---:|---:|---:|---|---|---|---|
| *Vmn1r166*<br>`ENSMUSG00000096073` | V1R | `V1R_chr7_cl016` (dominant) | `V1R_chr7_sc013` | 0.500 | 70 | 0.972 | suspected_em_redistribution | moderate | tentative_unconfirmed | consistent (track check) |
| *Vmn1r143*<br>`ENSMUSG00000096071` | V1R | `V1R_chr7_cl016` (dominant) | `V1R_chr7_sc013` | 0.000 | 2 | 0.028 | suspected_em_redistribution | alternative_candidate | tentative_unconfirmed | consistent (track check) |
| *Vmn1r103*<br>`ENSMUSG00000096663` | V1R | `V1R_chr7_cl015` | `V1R_chr7_sc013` | 0.000 | 186 | 0.429 | no_redistribution_signature | unresolvable | tentative_unconfirmed | consistent (track check) |
| *Vmn1r104*<br>`ENSMUSG00000096903` | V1R | `V1R_chr7_cl015` | `V1R_chr7_sc013` | 0.000 | 186 | 0.429 | no_redistribution_signature | alternative_candidate | tentative_unconfirmed | consistent (track check) |
| *Vmn1r-ps75*<br>`ENSMUSG00000095185` | V1R | `V1R_chr7_cl015` | `V1R_chr7_sc013` | 0.000 | 10 | 0.023 | no_redistribution_signature | alternative_candidate | tentative_unconfirmed | consistent (track check) |

Candidates nominated per called cluster — `V1R_chr7_cl015`: 3, `V1R_chr7_cl016`: 2. The project's stated tolerance is 2-3 candidates per cell, so a cluster reporting more than one paralog is the expected outcome of 75bp reads, not a failure of the assay.

Evidence basis: `bam_unique_mapq255`. Unique-read evidence (MAPQ-filtered BAM) is stronger than an EM fraction, because it does not depend on the redistribution step that makes within-cluster calls unsafe.

Track notes: EM_UNIQUE_CONTRADICTION: unique reads rank Vmn1r103 first (EM fraction 0.0000) while the EM-dominant paralog Vmn1r131 (EM fraction 0.999) carries 0 unique reads. The two evidence channels name different genes; neither identifies the expressed receptor in this cluster. ranked by bam_unique_mapq255; unique-read ratio rank1/rank2=1.00; background floor (median unique reads of non-top members)=2.0, top/floor=93.00; EM counts shown for context only and were NOT used to rank. n_transcripts=1 / ranked by bam_unique_mapq255; unique-read ratio rank1/rank2=35.00; background floor (median unique reads of non-top members)=0.0, top/floor=inf; EM counts shown for context only and were NOT used to rank. n_transcripts=1

Every candidate in this table carries confirmation status `tentative_unconfirmed`. **No individual receptor identity is confirmed by this pipeline.** Confirmation requires evidence that does not go through the EM step — longer reads, targeted amplicon sequencing of the cluster, or in-situ/immunostaining.

### `target2cellsRep2_S4` (target, 2 cells, prep suspect_failed_prep)

**Highest tier reported: 0 (tissue_identity).** stopped at tier 1 (sort_validation) [FAIL]: library_status=FAILED -- target sample with actin sum 0.9 CPM < 100 AND Trpc2 0.16 CPM < 10: empty/failed library, no biology may be reported Every CPM in this sample is a ratio taken against a near-empty or failed library, so no downstream number is interpretable however large it looks.

**Tier 0 -- tissue identity (VNO vs main olfactory epithelium)** — `PASS`

tissue_verdict=VNO -- VNO-specific panel establishes vomeronasal identity

**Tier 1 -- sort validation and library viability** — `FAIL`

> FAIL. library_status=FAILED -- target sample with actin sum 0.9 CPM < 100 AND Trpc2 0.16 CPM < 10: empty/failed library, no biology may be reported Every CPM in this sample is a ratio taken against a near-empty or failed library, so no downstream number is interpretable however large it looks.

**Tier 2 -- population identification (V1R vs V2R)** — `SUPPRESSED`

> SUPPRESSED. not evaluated -- tier 1 (sort_validation) FAIL: library_status=FAILED -- target sample with actin sum 0.9 CPM < 100 AND Trpc2 0.16 CPM < 10: empty/failed library, no biology may be reported Every CPM in this sample is a ratio taken against a near-empty or failed library, so no downstream number is interpretable however large it looks.

**Tier 3 -- cluster-level VR calls (reliable tier)** — `SUPPRESSED`

> SUPPRESSED. suppressed by tier 1 (sort_validation) FAIL: library_status=FAILED -- target sample with actin sum 0.9 CPM < 100 AND Trpc2 0.16 CPM < 10: empty/failed library, no biology may be reported Every CPM in this sample is a ratio taken against a near-empty or failed library, so no downstream number is interpretable however large it looks.

**Tier 4 -- individual receptor candidates (tentative)** — `SUPPRESSED`

> SUPPRESSED. suppressed by tier 1 (sort_validation) FAIL: library_status=FAILED -- target sample with actin sum 0.9 CPM < 100 AND Trpc2 0.16 CPM < 10: empty/failed library, no biology may be reported Every CPM in this sample is a ratio taken against a near-empty or failed library, so no downstream number is interpretable however large it looks.

### `target2cellsRep3_S5` (target, 2 cells, prep suspect_failed_prep)

**Highest tier reported: none.** stopped at tier 0 (tissue_identity) [FAIL]: tissue_verdict=no_tissue_signal -- both marker panels sit below the absolute CPM floor, so this target library carries no tissue information either way. VNO identity is UNCONFIRMED (this is not positive evidence of the wrong tissue), and an unconfirmed tissue cannot support a receptor call.

**Tier 0 -- tissue identity (VNO vs main olfactory epithelium)** — `FAIL`

> FAIL. tissue_verdict=no_tissue_signal -- both marker panels sit below the absolute CPM floor, so this target library carries no tissue information either way. VNO identity is UNCONFIRMED (this is not positive evidence of the wrong tissue), and an unconfirmed tissue cannot support a receptor call.

**Tier 1 -- sort validation and library viability** — `SUPPRESSED`

> SUPPRESSED. not evaluated -- tier 0 (tissue_identity) FAIL: tissue_verdict=no_tissue_signal -- both marker panels sit below the absolute CPM floor, so this target library carries no tissue information either way. VNO identity is UNCONFIRMED (this is not positive evidence of the wrong tissue), and an unconfirmed tissue cannot support a receptor call.

**Tier 2 -- population identification (V1R vs V2R)** — `SUPPRESSED`

> SUPPRESSED. not evaluated -- tier 0 (tissue_identity) FAIL: tissue_verdict=no_tissue_signal -- both marker panels sit below the absolute CPM floor, so this target library carries no tissue information either way. VNO identity is UNCONFIRMED (this is not positive evidence of the wrong tissue), and an unconfirmed tissue cannot support a receptor call.

**Tier 3 -- cluster-level VR calls (reliable tier)** — `SUPPRESSED`

> SUPPRESSED. suppressed by tier 0 (tissue_identity) FAIL: tissue_verdict=no_tissue_signal -- both marker panels sit below the absolute CPM floor, so this target library carries no tissue information either way. VNO identity is UNCONFIRMED (this is not positive evidence of the wrong tissue), and an unconfirmed tissue cannot support a receptor call.

**Tier 4 -- individual receptor candidates (tentative)** — `SUPPRESSED`

> SUPPRESSED. suppressed by tier 0 (tissue_identity) FAIL: tissue_verdict=no_tissue_signal -- both marker panels sit below the absolute CPM floor, so this target library carries no tissue information either way. VNO identity is UNCONFIRMED (this is not positive evidence of the wrong tissue), and an unconfirmed tissue cannot support a receptor call.

