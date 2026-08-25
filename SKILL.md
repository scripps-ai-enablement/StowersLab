---
name: vno-receptor-rnaseq
description: Identify which vomeronasal receptor (V1R/Vmn1r, V2R/Vmn2r) each mouse VNO sensory neuron or FACS-sorted cell pool expresses, from an nf-core/rnaseq star_salmon results directory. Use for VNO / vomeronasal / Trpc2 / OMP-Cre GFP sort RNA-seq, monogenic receptor choice, receptor-cluster (paralog) aggregation, and any short-read quantification where Salmon EM spreads reads across near-identical gene-family paralogs. Covers tissue-identity gating (MOE vs VNO), sort validation, cluster-level aggregation at two gap tiers, and the unique-read gate that separates EM redistribution from genuine co-expression.
---

# VNO receptor RNA-seq: tiered receptor identification

Each mouse vomeronasal sensory neuron expresses **exactly one** vomeronasal
receptor (monogenic choice). This skill runs the validated pipeline that names
that receptor — or, more often, states honestly how far the data supports a
name.

The hard part is not quantification. ~250 V1R (`Vmn1r*`) and ~120 V2R
(`Vmn2r*`) genes sit in genomic clusters of local duplicates sharing 85-95%
nucleotide identity. 75bp reads cannot be uniquely assigned within a cluster,
so Salmon's EM spreads ambiguous reads across paralogs. **Per-gene VR counts
inside a cluster are not trustworthy**; cluster-level aggregation is the
reliable readout (Dietschi et al. 2022, Sci Adv 8(46) eabn7450,
doi:10.1126/sciadv.abn7450). Pseudogenes inside functional clusters can appear
"expressed" through the same channel.

Mathematical warning signature: paralogs within one cluster at near-equal
fractions (33/33/33) usually means ONE paralog was expressed and EM split the
reads. **Usually — not always.** See the unique-read gate below; an even split
across a multi-cell pool can be real co-expression.

## Version and distribution

This skill is distributed as **per-user copies** (`skill/install_skill.py` in
the pipeline package), not from a shared registry entry. Two people can
therefore be on different versions with no warning. Call
`vno_skill_version()` and compare before trusting that two runs used the same
conventions; re-run the installer from a newer package to update.

If your organization has published this skill at org scope instead, it appears
in `host.skills.list()` with `origin: "organization"` and needs no install step
— that is the better arrangement for more than a couple of users, because one
edit reaches everyone.

## Setup — runs anywhere, nothing to build

The pipeline is a set of standalone Python modules plus a bash driver. No
package to install, no service, no site-specific paths. It needs:

| requirement | why |
|---|---|
| Python 3.9+ with pandas, numpy, matplotlib, pyyaml, scipy | the modules import only these plus the standard library |
| `samtools` on PATH | the unique-read evidence channel (`samtools view -q 255`) |
| an nf-core/rnaseq `star_salmon` output directory | the input |
| the reference GTF used for that quantification | to build the VR cluster tables |

Install whichever way suits the machine — `vno_install_check()` reports what is
missing and prints the matching remedy:

```bash
conda create -n vr python=3.11 pandas numpy matplotlib pyyaml scipy -c conda-forge
conda activate vr && conda install -c bioconda samtools
# or: python3 -m venv ~/vr-env && . ~/vr-env/bin/activate && pip install pandas numpy matplotlib pyyaml scipy
# or on an HPC: module load python3 && module load samtools
```

`samtools` is separate from the python stack and is easy to miss. Without it
`vr_quantify.py` stops at a preflight guard with the fix; `--no-bam` runs anyway
but **degrades the result** — no individual-receptor calls are possible and
EM-redistribution flags fall back to fraction-only evidence (22 candidate rows
instead of 26 on the reference dataset). Cluster-level results are unaffected.

### The one command

`bin/vr_analyze.py` is the whole interface — it scaffolds a work folder,
discovers samples, infers `target`/`nontarget` from their names, writes and
validates the config, runs every stage in order, and prints both a human report
and JSON:

```bash
python3 bin/vr_analyze.py --results /path/to/nfcore/results --gtf /path/to/genes.gtf
```

Useful flags: `--out DIR` (default `./vr_out_<trial>`), `--trial NAME`,
`--threads N`, `--no-bam`, `--force-ref`, `--json-only` (JSON on stdout,
progress on stderr — the mode to use when driving it programmatically),
`--dry-run` (scaffold and validate, run nothing), `--selftest`.

Exit codes are meaningful: `0` ran and something reached a cluster call, `3` ran
but nothing survived QC (**a result, not a failure**), `4` sample roles could not
be inferred, `5` missing dependency, `6` bad input. Never treat `3` as an error
to retry around — read the stop reasons.

**It refuses rather than guessing.** If a sample name matches neither the target
nor the nontarget vocabulary, it exits `4` and asks for `--target`/`--nontarget`
regexes, because inverting those two inverts the entire sort-validation tier.

The JSON on stdout carries `libraries[]` (with `highest_tier` and
`stopped_because` per library), `clusters_called`, `clusters_suppressed`,
`called_clusters[]`, `candidates`, `top_candidates[]`, `flags{}`, and
`unique_read_channel`. **`clusters_called` already excludes QC-suppressed
libraries**; the count that was withheld is reported separately as
`clusters_suppressed` so nothing is silently dropped. `--selftest` asserts that
gate against synthetic rows (13 assertions) — the first version of the summary
did leak a failed library's called cluster.

### Lower-level: scaffold only

`bin/vr_init.py` writes a config without running anything, if you want to edit it
before the stages:

```bash
python3 bin/vr_init.py --dest ~/vr_run --trial trial1 \
    --results /path/to/nfcore/results \
    --gtf /path/to/genes.gtf --genome-key my-grcm38
```

It discovers sample names from the quant table's header, guesses `cell_type`
from them, and validates the result through the pipeline's own `vr_config`
before exiting. **Check the guessed `cell_type` values** — `target` vs
`nontarget` cannot be inferred reliably from a filename, and the
sort-validation tier cannot interpret `UNKNOWN`. Add later trials with
`--add-trial NAME --results DIR`, which inserts inside the existing blocks
rather than appending a duplicate top-level key.

Threshold, marker and path blocks are copied verbatim from the shipped
`config/project.template.yaml`, so a generated config produces identical
verdicts to the reference run. Verified on a fresh work folder outside the
original project tree, compared field-by-field: `tissue_verdict`,
`sort_verdict`, `library_status`, `population_call`, `qc_overall`,
`suppress_biology`, `em_flag`, `pseudogene_bleed_flag`, `interpretation`, and
every candidate's `rank`/`gene_name`/`confidence`/`confirmation_status` — all
identical. (A separate clean-tarball-unpack test confirmed the package runs to
completion and agrees on tissue verdicts and candidate count; it did not repeat
the field-by-field comparison.)
Edit the template to change defaults for future work folders.

### Work-folder resolution

Helpers take `work_dir=` explicitly; otherwise `$VR_WORK`, then `$PWD` if it
contains `config/project.yaml`. There is no built-in default path.

```
<work folder>
  bin/     vr_init.py vr_config.py build_vr_reference.py vr_qc.py vr_markers.py
           vr_sample_qc.py vr_quantify.py vr_report.py vr_figures.py
           vr_qc_figures.py run_pipeline.sh
  config/  project.yaml            <- READ THIS FIRST; every path and threshold
           project.template.yaml   <- the analysis contract, copied into new configs
  ref/     vr_gene_annotation.tsv vr_clusters.tsv vr_gene_to_cluster.tsv
  results/<trial>/                 <- per-trial outputs
  docs/    VR_OUTPUT_GUIDE.md VR_PIPELINE_FLOWCHART.md QC_MODULE_README.md
```

This skill does **not** duplicate the modules — they are ~5000 lines and the
work folder is the source of truth. `kernel.py` ships the thin helpers a session
needs to locate a run, load the tables, apply the gate, and read the flag
vocabulary without re-deriving any of it.

### If the data lives on a cluster

Optional. Set `$VR_SSH_TARGET` and reach it with
`host.compute.create(vno_ssh_target())` from the `repl` tool. `call_command` has
a short timeout on such hosts, so anything that parses the GTF or reads BAMs
must go through `submit_job`. On a local run none of this applies.


## Input contract

Two things, nothing else:

1. An **nf-core/rnaseq (>=3.18) `star_salmon` results directory** containing
   `salmon.merged.gene_counts.tsv`, `salmon.merged.gene_tpm.tsv`, `tx2gene.tsv`,
   `log/<sample>.Log.final.out`, `qualimap/<sample>/rnaseq_qc_results.txt`,
   `<sample>/quant.sf`, and (for individual-receptor evidence) the per-sample
   coordinate-sorted BAM.
2. **`config/project.yaml`** — the path contract, sample table, marker panels
   and every numeric threshold.

Reference genome must be GRCm38 with Ensembl chromosome naming (`1`, `7`, `X`,
not `chr1`). The GTF used here has **no `gene` feature rows**; gene
spans are aggregated from transcript rows with an independent exon-derived
cross-check.

`vno_validate_run_dir(results_root)` checks all of this and returns what is
missing before you spend a job on it.

## Running it on a new trial

```bash
python3 bin/vr_analyze.py --results <nfcore run dir> --gtf <genes.gtf> --trial <name>
```

That is it. Everything below is the manual path, for debugging one stage or when
you need to hand-edit the config first.

<details>
<summary>Manual stage-by-stage</summary>

1. **Add the trial to `config/project.yaml`** under `trials:` (`results:`,
   `fastq:`, `platform:`) and under `samples:` with `cell_type`
   (`target` = GFP+, `nontarget` = GFP-), `n_cells`, `prep_status`.
   `vno_run_plan(work_dir, trial)` prints the exact commands.
2. **Edit the config INSIDE the existing blocks.** A naive append that creates
   a second top-level `markers:` or `thresholds:` key silently shadows the
   first — YAML resolves last-wins with no error. Run
   `vno_check_yaml_duplicate_keys(path)` after every config edit. This has
   already broken this project once. (`vr_analyze.py` and `vr_init.py` both
   check for it themselves.)
3. **Reference layer** — only if the genome build changed:
   `python bin/build_vr_reference.py --config config/project.yaml --outdir ref`
4. **QC + gating:** `python bin/vr_sample_qc.py --all-trials`
5. **Quantification:** `python bin/vr_quantify.py --trial <trial> --threads 8`
   (add `--no-bam` to skip unique-read evidence; then no individual-receptor
   call is possible)
6. **Report:** `python bin/vr_report.py` (`--selftest` first — it asserts the
   tier gate cannot be bypassed)
7. **Figures:** `python bin/vr_figures.py` and
   `python bin/vr_qc_figures.py --all-trials --fig-dir results/figures`

</details>

### What to check first, in this order

Do not look at receptor names first. Look at these:

1. `results/tier_status_all.tsv` → `tier0_status`. This, not `tissue_verdict`
   read on its own, is the verdict on tissue. `tissue_verdict` has six values
   (`VNO`, `VNO_dominant_mixed`, `MOE`, `MOE_dominant_mixed`,
   `ambiguous_mixed`, `no_tissue_signal`) and they do not map one-to-one to
   pass/fail: a GFP− nontarget library legitimately reads `no_tissue_signal`
   and still passes tier 0, because a control with no VNO neurons in it is the
   expected result, whereas the same value in a GFP+ target library is a stop.
   `MOE` / `MOE_dominant_mixed` is the wrong tissue and the remedy is wet-lab.
2. `results/<trial>/marker_cpm.tsv` → `library_status`. `FAILED` or
   `DEGENERATE` ends that sample.
3. `results/tier_status_all.tsv` → `highest_tier_reported` per sample. This is
   the one-line answer to "what does this library support?"
4. `results/<trial>/vr_cluster_expression.tsv` where `is_called == 1` — the
   reliable tier.
5. `results/<trial>/vr_artifact_flags.tsv` → `em_flag` **before** reading any
   individual gene name.
6. `results/<trial>/vr_candidates.tsv` — every row is
   `confirmation_status = tentative_unconfirmed` by construction.

## The CPM convention (settled — do not re-derive)

**CPM = count / (column sum over ALL genes of the unscaled
`salmon.merged.gene_counts.tsv`) × 1e6.**

The scaled (`gene_counts_scaled`) and TPM tables give values 2.3-2.5× lower.
Every threshold in `config/project.yaml` is calibrated against the unscaled
all-gene-denominator convention; using another table silently invalidates all
of them. `vno_cpm()` implements exactly this.

## The five-tier priority hierarchy

The hierarchy is enforced **structurally** in `vr_report.py`, not by
convention: each tier's content comes from a callable that `TierGate.emit()`
refuses to invoke unless every upstream tier passed. A wrong-tissue or failed
library cannot reach a receptor call because the producer function is never
called; the gate substitutes an explicit suppression line naming the failing
tier. `vr_report.py --selftest` asserts this (35 assertions).

| tier | name | question | source |
|---|---|---|---|
| 0 | `tissue_identity` | is this VNO at all? | MOE vs VNO marker panels |
| 1 | `sort_validation` | viable GFP+/GFP- library? | Trpc2 CPM + actin/failed-library gate |
| 2 | `population_id` | V1R or V2R? | Gnai2:Gnao1 |
| 3 | `cluster_vr` | which receptor cluster? | **the reliable tier** |
| 4 | `individual_vr` | which paralog? | always tentative |

`status` vocabulary: `PASS`, `PASS_WITH_CAVEAT`, `FAIL`, `NO_DATA`,
`SUPPRESSED`.

**Clearance gate for reporting biology:** `qc_overall == 'USABLE'` **AND**
`suppress_biology` is false. `suppress_biology` alone is insufficient — a
library can carry `suppress_biology=False` with `sort_verdict=FAIL`. Use
`vno_clearance(sample_qc_df)`, don't re-derive the condition.

## Marker and tissue panels

**Tier 0 runs before sort validation.** Compared on the *maximum* member of
each panel against a 100 CPM absolute floor (a ratio of two sub-floor noise
values carries no tissue information), then a 3× dominance ratio:

| panel | genes |
|---|---|
| VNO-specific | `Trpc2`, `Vmn1r*` sum, `Vmn2r*` sum |
| MOE | `Olfr*` sum, `Adcy3`, `Cnga2`, `Gnal` |
| tissue-shared — **never a VNO marker** | `Omp` |

**`Omp` is TISSUE-SHARED.** OMP-Cre labels mature main-olfactory neurons *and*
VNO neurons. A GFP+ sort from the wrong tissue yields Omp-high, Trpc2-zero
cells that look like a successful sort. Reading high `Omp` as evidence of a VNO
population is the exact error that makes a main-olfactory library pass for a
successful GFP+ VNO sort.

`tissue_verdict` vocabulary: `VNO`, `VNO_dominant_mixed` (VNO panel dominant
but the MOE panel is also above the floor — the common outcome for a real
100-cell VNO pool), `MOE`, `MOE_dominant_mixed`, `ambiguous_mixed`,
`no_tissue_signal` (both panels below the floor). `sort_verdict` is `PASS`,
`FAIL`, or `FAIL_WRONG_TISSUE`. Do not collapse these to a boolean — read
`tier0_status` for the pass/fail decision, which accounts for `cell_type`.

`Gnai2`/`Gnao1` are **excluded** from the tissue panel — they split V1R vs V2R
*within* the VNO and are broadly expressed elsewhere (`Gnai2` reads 165 CPM in a
GFP- library with no VR signal at all). They are tier-2 markers only.

Cross-validation rule: a **V1R call requires `Gnai2` > `Gnao1`**. A
contradiction between VR evidence and marker evidence must be surfaced, never
smoothed over.

`ratio_low_support` is set when the raw denominator is < 10 reads. **Never quote
such a ratio as precise** — "422.7:1" and "385:1" from the same library were the
same measurement with `Gnao1 = 7 reads`.

## Cluster aggregation: report both tiers

Clusters are contiguous same-family paralog runs. Two tiers, both in
`ref/vr_clusters.tsv` and `ref/vr_gene_to_cluster.tsv`:

- `cluster_id` at **200kb** max gap → 24 V1R, 37 V2R clusters
- `supercluster_id` at **800kb** → the coarser tier

**200kb is not a natural break.** The V1R inter-gene gap KDE has its minimum
near 2Mb; 200kb was retained as the conservative choice, and the threshold
sweep is recorded in `ref/vr_gtf_parse_report.txt`. The consequence is real: a
217,366bp gap splits the chr7 V1R megacluster into `V1R_chr7_cl015` +
`V1R_chr7_cl016`, which `V1R_chr7_sc013` reunites. **Report both tiers for that
region** (rows carry `chr7_dual_tier_region = 1`).

**V2R aggregation is weaker protection than V1R**: 18 of 37 V2R clusters are
singletons and only 180/222 V2R genes sit in clusters of >=5. Aggregation
cannot absorb EM ambiguity that has nowhere to go.

Join rules, both non-negotiable:

- **Join on `gene_id`, never `gene_name`.** `Vmn1r-ps5` and `Vmn2r118` each map
  to two distinct `gene_id`s on different chromosomes.
- **Filter `is_primary_assembly == 1`** (538 of 541 genes) when joining to
  quant tables. This is what reconciles exactly to the 318 `Vmn1r` / 220
  `Vmn2r` quant rows.

## The EM-artifact rule and the unique-read gate

This is the part that is easy to get wrong in both directions.

**Step 1 — is the within-cluster split even?** Uniformity of observed fractions
is tested against 1/k by **Monte Carlo** from `Multinomial(N, 1/k)` at the
observed N and k (4000 draws). The asymptotic chi-square is invalid at N of
tens. Power against a monogenic alternative (dominant at 0.90) is simulated too;
power < 0.80 grades the pair `indeterminate_low_depth`.

`p_uniform` **saturates at 2.4994e-04 = 1/4001** — that is the permutation
floor, a bound. Quote it as `p < 2.5e-4`, never as a point estimate.

**Step 2 — an even split is NECESSARY BUT NOT SUFFICIENT.** Inside a detected
even block, each member's unique MAPQ255 deduplicated reads are compared
against `max(10 reads, 3 x median unique count of cluster members OUTSIDE the
block)` — the cluster's own measured mismapping floor.

- exactly **one** co-dominant member clears it → **EM redistribution**
- **two or more** clear it → **genuine CO-EXPRESSION**, because monogenic
  choice is a *per-cell* rule and a multi-cell pool legitimately captures two
  paralogs of one cluster

That distinction matters for any pooled design. Use `vno_em_verdict(row)`; it
reads the columns rather than re-judging evenness.

`em_flag` vocabulary — **these exact values, do not invent others**:
`no_signal`, `insufficient_signal`, `no_redistribution_signature`,
`single_paralog_only`, `suspected_em_redistribution`.

**Individual-receptor evidence** is `bam_unique_mapq255` — STAR MAPQ 255 means
the read was placed at exactly one locus — counted mate-wise, with and without
duplicates. Salmon's `aux_info/ambig_info.tsv` `UniqueCount` is kept as a
cross-check only: its "unique" means unique to a *transcript*, so gene-level
support is understated for the 30 of 538 VR genes with multiple transcripts.
Read the evidence columns before quoting a rank — `top_over_background_floor`
can be `inf` against a zero background, meaning a rank resting on very few
reads.

**Pseudogene bleed is an open mechanism question.** `pseudogene_bleed_flag`
reports apparent expression; `pseudogene_mechanism` states that EM leakage and
genuine cluster-regulatory transcription of the pseudogene locus are not
separable from this data. **Dietschi is not a quantitative mouse prior**: their
pseudogene-in-cluster result is a Wilcoxon rank sum test significant in rat only
(W = 1113, P = 0.002987) and **not** significant in mouse (W = 1214,
P = 0.5704), and their proposed mechanism is regulatory (cluster-associated
transcription-stabilizing elements), not multi-mapping. Do not claim literature
support for a quantitative mouse expectation.

What that paper *does* support is the cluster tier itself: cluster identity on
expression level is ANOVA P < 2 × 10⁻¹⁶, and cluster membership on VSN
abundance (choice probability) is ANOVA P = 0.000232. They also estimated only
57% of V1rD reads assignable without ambiguity, and handled it differently from
us — STAR `--outFilterMultimapNmax 4` plus `featureCounts -M --fraction`
(fractional counting of multimappers) rather than cluster aggregation with a
MAPQ255 unique-read requirement. Verify these against the paper (open access
via PMC, doi:10.1126/sciadv.abn7450) before quoting them onward; they are
recorded here with the test that produced each so a reader can check the right
one.

## Mandatory uncertainty reporting

Every individual-receptor statement carries (a) its cluster context, (b) its
`confirmation_status`, (c) the evidence channel and read count it rests on.
Preserve asymmetries rather than averaging them: in this project the *call*
that a cluster held only one expressed paralog was solid (4,473 EM counts with
zero single-locus reads is not an independent observation), while *which*
paralog was the source rested on 51 unique reads — ~1.1% of that gene's
apparent expression. Both halves must appear in the same sentence.

When the two evidence channels name different genes,
`evidence_contradiction` is populated and the honest output is **no call for
that cluster**. Do not pick the channel you like better.

## Failure modes this project actually hit

These cost weeks. Check for them explicitly.

1. **A wrong-tissue dataset that looks like a failed sort.** Four libraries were
   briefed as "failed / degraded preps". They were main olfactory epithelium:
   `Olfr` sum 48,729 CPM (~4.9% of the library) with zero VR signal, and the
   "cleanest" sample was a textbook mature MOE neuron (Adcy3 2374, Cnga2 585,
   Gnal 4277, Trpc2 0). Both trials used the same annotation, which yields
   17,918-37,230 Trpc2 counts in genuine VNO libraries — so this is not a
   quantification artifact. **The remedy is wet-lab (dissection, sort gate), not
   re-quantification.** Verify a diagnosis against the actual annotation before
   recommending a re-run.
2. **CPM inflation in a near-empty library.** CPM is a ratio: one library showed
   39,881 CPM actin from only 60,744 assigned counts. `min_assigned_counts`
   marks such a library `DEGENERATE` and suppresses biology. A separate case
   carried 501.5 CPM of `Vmn1r` signal in a library failing every other check —
   VR family signal can survive a failed library, and it is not a receptor call.
3. **An even split that is co-expression, not an artifact.** A cluster whose two
   paralogs each carried ~18.5k *independent* unique reads was initially
   labelled a "genuine EM signature". It is real co-expression across a 100-cell
   pool; its `em_flag` is `no_redistribution_signature` with
   `even_block_size = 0`.
4. **A figure that recomputed a verdict instead of reading the flag column.**
   That mislabel came from the figure applying its own visual evenness
   criterion. Figures must read `em_flag`; a downstream consumer that re-derives
   an upstream verdict will eventually disagree with it.
5. **Duplicate YAML keys silently shadowing config.** Covered above — run
   `vno_check_yaml_duplicate_keys()` after every edit and re-read the file.

## kernel.py helpers

Loaded automatically with this skill.

| function | purpose |
|---|---|
| `vno_skill_version()` | version of this skill copy — compare with a collaborator's before trusting two runs match |
| `vno_install_check()` | which dependencies are missing, with the platform-matched remedy |
| `vno_default_work_dir()` | resolve the work folder: `$VR_WORK`, else `$PWD` if it holds `config/project.yaml` |
| `vno_validate_work_dir(work_dir=None)` | report what a work folder is missing and how to fix it |
| `vno_ssh_target()` | `$VR_SSH_TARGET` if the data lives on a cluster, else `None` |
| `vno_paths(work_dir=None)` | resolved bin/ref/config/results paths + install hints |
| `vno_validate_run_dir(results_root)` | check an nf-core star_salmon dir before spending a job; lists missing files and discovered samples |
| `vno_read_table(path)` | read a pipeline TSV, returning `(df, provenance_comment)` |
| `vno_cpm(counts, ...)` | the settled CPM convention, nothing else |
| `vno_load_results(work_dir=None, trial=...)` | load every output table for a trial into one dict |
| `vno_clearance(sample_qc)` | apply the reporting clearance gate |
| `vno_em_verdict(row)` | read `em_flag` + unique-read columns into a verdict dict; never re-judges evenness |
| `vno_flag_vocabulary()` | the exact allowed values for every status/flag column |
| `vno_check_yaml_duplicate_keys(path)` | catch the shadowing bug |
| `vno_run_plan(work_dir=None, trial=...)` | ordered shell commands for a new trial |

Deliverable figures: load the `figure-style` skill and call
`apply_figure_style()` before rendering.
