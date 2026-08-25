# VNO receptor tiered analysis report -- trial1

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
* trial1 (Illumina) quantification: `<project>/results/trail1/star_salmon/salmon.merged.gene_counts.tsv`
  * VR table `cluster_expression`: `<work>/results/trial1/vr_cluster_expression.tsv`
  * VR table `within_cluster_fractions`: `<work>/results/trial1/vr_within_cluster_fractions.tsv`
  * VR table `artifact_flags`: `<work>/results/trial1/vr_artifact_flags.tsv`
  * VR table `candidates`: `<work>/results/trial1/vr_candidates.tsv`

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
| trial1 | `pool100cells_S8` | target | 100 | UNUSABLE | MOE | FAIL_WRONG_TISSUE | undetermined | DEGENERATE | True | 0.00 | 39880.7 |
| trial1 | `pool2cellsRep1_S5` | target | 2 | UNUSABLE | no_tissue_signal | FAIL | undetermined | OK | False | 0.00 | 917.6 |
| trial1 | `pool2cellsRep2_S6` | target | 2 | UNUSABLE | MOE | FAIL_WRONG_TISSUE | V2R_dominant | FAILED | True | 0.00 | 8.3 |
| trial1 | `pool2cellsRep3_S7` | target | 2 | UNUSABLE | MOE | FAIL_WRONG_TISSUE | V2R_dominant | OK | True | 0.00 | 670.8 |

## Tier outcomes by sample

### `pool100cells_S8` (target, 100 cells, prep suspect_failed_prep)

**Highest tier reported: none.** stopped at tier 0 (tissue_identity) [FAIL]: tissue_verdict=MOE -- the main-olfactory marker panel is above the tissue floor while the VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is not. This is the WRONG TISSUE, not a failed VNO sort; Trpc2=0 is the expected value for main olfactory epithelium. No VR biology is reportable and the remedy is wet-lab (dissection / sort gate), not re-quantification.

**Tier 0 -- tissue identity (VNO vs main olfactory epithelium)** — `FAIL`

> FAIL. tissue_verdict=MOE -- the main-olfactory marker panel is above the tissue floor while the VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is not. This is the WRONG TISSUE, not a failed VNO sort; Trpc2=0 is the expected value for main olfactory epithelium. No VR biology is reportable and the remedy is wet-lab (dissection / sort gate), not re-quantification.

**Tier 1 -- sort validation and library viability** — `SUPPRESSED`

> SUPPRESSED. not evaluated -- tier 0 (tissue_identity) FAIL: tissue_verdict=MOE -- the main-olfactory marker panel is above the tissue floor while the VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is not. This is the WRONG TISSUE, not a failed VNO sort; Trpc2=0 is the expected value for main olfactory epithelium. No VR biology is reportable and the remedy is wet-lab (dissection / sort gate), not re-quantification.

**Tier 2 -- population identification (V1R vs V2R)** — `SUPPRESSED`

> SUPPRESSED. not evaluated -- tier 0 (tissue_identity) FAIL: tissue_verdict=MOE -- the main-olfactory marker panel is above the tissue floor while the VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is not. This is the WRONG TISSUE, not a failed VNO sort; Trpc2=0 is the expected value for main olfactory epithelium. No VR biology is reportable and the remedy is wet-lab (dissection / sort gate), not re-quantification.

**Tier 3 -- cluster-level VR calls (reliable tier)** — `SUPPRESSED`

> SUPPRESSED. suppressed by tier 0 (tissue_identity) FAIL: tissue_verdict=MOE -- the main-olfactory marker panel is above the tissue floor while the VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is not. This is the WRONG TISSUE, not a failed VNO sort; Trpc2=0 is the expected value for main olfactory epithelium. No VR biology is reportable and the remedy is wet-lab (dissection / sort gate), not re-quantification.

**Tier 4 -- individual receptor candidates (tentative)** — `SUPPRESSED`

> SUPPRESSED. suppressed by tier 0 (tissue_identity) FAIL: tissue_verdict=MOE -- the main-olfactory marker panel is above the tissue floor while the VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is not. This is the WRONG TISSUE, not a failed VNO sort; Trpc2=0 is the expected value for main olfactory epithelium. No VR biology is reportable and the remedy is wet-lab (dissection / sort gate), not re-quantification.

### `pool2cellsRep1_S5` (target, 2 cells, prep suspect_degraded)

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

### `pool2cellsRep2_S6` (target, 2 cells, prep borderline)

**Highest tier reported: none.** stopped at tier 0 (tissue_identity) [FAIL]: tissue_verdict=MOE -- the main-olfactory marker panel is above the tissue floor while the VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is not. This is the WRONG TISSUE, not a failed VNO sort; Trpc2=0 is the expected value for main olfactory epithelium. No VR biology is reportable and the remedy is wet-lab (dissection / sort gate), not re-quantification.

**Tier 0 -- tissue identity (VNO vs main olfactory epithelium)** — `FAIL`

> FAIL. tissue_verdict=MOE -- the main-olfactory marker panel is above the tissue floor while the VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is not. This is the WRONG TISSUE, not a failed VNO sort; Trpc2=0 is the expected value for main olfactory epithelium. No VR biology is reportable and the remedy is wet-lab (dissection / sort gate), not re-quantification.

**Tier 1 -- sort validation and library viability** — `SUPPRESSED`

> SUPPRESSED. not evaluated -- tier 0 (tissue_identity) FAIL: tissue_verdict=MOE -- the main-olfactory marker panel is above the tissue floor while the VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is not. This is the WRONG TISSUE, not a failed VNO sort; Trpc2=0 is the expected value for main olfactory epithelium. No VR biology is reportable and the remedy is wet-lab (dissection / sort gate), not re-quantification.

**Tier 2 -- population identification (V1R vs V2R)** — `SUPPRESSED`

> SUPPRESSED. not evaluated -- tier 0 (tissue_identity) FAIL: tissue_verdict=MOE -- the main-olfactory marker panel is above the tissue floor while the VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is not. This is the WRONG TISSUE, not a failed VNO sort; Trpc2=0 is the expected value for main olfactory epithelium. No VR biology is reportable and the remedy is wet-lab (dissection / sort gate), not re-quantification.

**Tier 3 -- cluster-level VR calls (reliable tier)** — `SUPPRESSED`

> SUPPRESSED. suppressed by tier 0 (tissue_identity) FAIL: tissue_verdict=MOE -- the main-olfactory marker panel is above the tissue floor while the VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is not. This is the WRONG TISSUE, not a failed VNO sort; Trpc2=0 is the expected value for main olfactory epithelium. No VR biology is reportable and the remedy is wet-lab (dissection / sort gate), not re-quantification.

**Tier 4 -- individual receptor candidates (tentative)** — `SUPPRESSED`

> SUPPRESSED. suppressed by tier 0 (tissue_identity) FAIL: tissue_verdict=MOE -- the main-olfactory marker panel is above the tissue floor while the VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is not. This is the WRONG TISSUE, not a failed VNO sort; Trpc2=0 is the expected value for main olfactory epithelium. No VR biology is reportable and the remedy is wet-lab (dissection / sort gate), not re-quantification.

### `pool2cellsRep3_S7` (target, 2 cells, prep ok)

**Highest tier reported: none.** stopped at tier 0 (tissue_identity) [FAIL]: tissue_verdict=MOE -- the main-olfactory marker panel is above the tissue floor while the VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is not. This is the WRONG TISSUE, not a failed VNO sort; Trpc2=0 is the expected value for main olfactory epithelium. No VR biology is reportable and the remedy is wet-lab (dissection / sort gate), not re-quantification.

**Tier 0 -- tissue identity (VNO vs main olfactory epithelium)** — `FAIL`

> FAIL. tissue_verdict=MOE -- the main-olfactory marker panel is above the tissue floor while the VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is not. This is the WRONG TISSUE, not a failed VNO sort; Trpc2=0 is the expected value for main olfactory epithelium. No VR biology is reportable and the remedy is wet-lab (dissection / sort gate), not re-quantification.

**Tier 1 -- sort validation and library viability** — `SUPPRESSED`

> SUPPRESSED. not evaluated -- tier 0 (tissue_identity) FAIL: tissue_verdict=MOE -- the main-olfactory marker panel is above the tissue floor while the VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is not. This is the WRONG TISSUE, not a failed VNO sort; Trpc2=0 is the expected value for main olfactory epithelium. No VR biology is reportable and the remedy is wet-lab (dissection / sort gate), not re-quantification.

**Tier 2 -- population identification (V1R vs V2R)** — `SUPPRESSED`

> SUPPRESSED. not evaluated -- tier 0 (tissue_identity) FAIL: tissue_verdict=MOE -- the main-olfactory marker panel is above the tissue floor while the VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is not. This is the WRONG TISSUE, not a failed VNO sort; Trpc2=0 is the expected value for main olfactory epithelium. No VR biology is reportable and the remedy is wet-lab (dissection / sort gate), not re-quantification.

**Tier 3 -- cluster-level VR calls (reliable tier)** — `SUPPRESSED`

> SUPPRESSED. suppressed by tier 0 (tissue_identity) FAIL: tissue_verdict=MOE -- the main-olfactory marker panel is above the tissue floor while the VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is not. This is the WRONG TISSUE, not a failed VNO sort; Trpc2=0 is the expected value for main olfactory epithelium. No VR biology is reportable and the remedy is wet-lab (dissection / sort gate), not re-quantification.

**Tier 4 -- individual receptor candidates (tentative)** — `SUPPRESSED`

> SUPPRESSED. suppressed by tier 0 (tissue_identity) FAIL: tissue_verdict=MOE -- the main-olfactory marker panel is above the tissue floor while the VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is not. This is the WRONG TISSUE, not a failed VNO sort; Trpc2=0 is the expected value for main olfactory epithelium. No VR biology is reportable and the remedy is wet-lab (dissection / sort gate), not re-quantification.

