# `vr_analysis` — VNO receptor RNA-seq pipeline

Vomeronasal receptor (V1R/V2R) identification from FACS-sorted mouse VNO
sensory neurons. Stowers Lab (Natalie Cole) / Scripps CCBB.
OMP-Cre × GFP reporter, GFP+ (target) vs GFP− (nontarget), Takara Smart-seq HT,
2-cell and 100-cell pools, 2×75bp PE. Upstream quantification is
nf-core/rnaseq 3.18, aligner `star_salmon`, genome
`mouse-ensembl-grcm38-r91` (GRCm38, Ensembl chromosome naming).

Each VNO sensory neuron expresses **exactly one** vomeronasal receptor
(monogenic choice). The pipeline's job is to name it — or to state precisely
how far the data falls short of naming it.

---

## READ THIS FIRST

Two findings determine how the whole `results/` tree should be read. A
newcomer who skips them will misread it.

**1. Trial 1 is the wrong tissue — main olfactory epithelium (MOE), not VNO.**
All four trial-1 libraries were originally briefed as "failed or degraded
preps". They are not: they are MOE. `pool100cells_S8` carries an `Olfr` family
sum of 48,729 CPM (~4.9% of the library) with zero VR signal, and
`pool2cellsRep3_S7` — briefed as the *cleanest* 2-cell sample — is a textbook
mature MOE neuron (Adcy3 2374, Cnga2 585, Gnal 4277, Olfr 1113, Trpc2 0 CPM).
This is not a quantification artifact: both trials used the same annotation,
which yields 17,918–37,230 Trpc2 counts in the trial-2 VNO libraries.
**No trial-1 VR biology is reportable, and the remedy is wet-lab (dissection /
sort gate), not re-quantification.** `results/trial1/` exists for the record;
every row in it is QC-suppressed.

The error that hid this for so long: **`Omp` is TISSUE-SHARED.** OMP-Cre labels
mature main-olfactory neurons as well as VNO neurons, so a GFP+ sort from the
wrong tissue is Omp-high and looks like a successful sort. `Omp` must never be
used as a VNO-specific marker.

**2. Only 4 of 10 libraries are usable, and only 3 support receptor
candidates.** All usable libraries are trial 2:

| library | status |
|---|---|
| `target100cellsRep1_S6` | usable, reaches tier 4 (receptor candidates) |
| `target100cellsRep2_S7` | usable, reaches tier 4 |
| `target2cellsRep1_S3` | usable, reaches tier 4 |
| `nontarget100cells_S8` | usable GFP− control; tier 3 then `NO_DATA` at tier 4 |
| `target2cellsRep2_S4` | passes tier 0 (`tissue_verdict=VNO`), **fails** tier 1 sort validation — zero receptor statements despite 501.5 CPM of Vmn1r signal |
| `target2cellsRep3_S5` | stops at tier 0 (`no_tissue_signal` in a GFP+ target) |
| all four trial-1 libraries | stop at tier 0 (three `MOE`, one `no_tissue_signal`) |

Read `tier0_status` in `results/tier_status_all.tsv`, not `tissue_verdict` on
its own — the two are not interchangeable. `tissue_verdict` takes six values
(`VNO`, `VNO_dominant_mixed`, `MOE`, `MOE_dominant_mixed`, `ambiguous_mixed`,
`no_tissue_signal`) and their meaning depends on `cell_type`: the GFP−
control `nontarget100cells_S8` reads `no_tissue_signal` and **passes** tier 0,
because a control containing no VNO neurons is the expected result, while the
identical value in a GFP+ target library (`target2cellsRep3_S5`,
`pool2cellsRep1_S5`) is a stop. Both usable 100-cell target libraries read
`VNO_dominant_mixed` — VNO panel dominant with the MOE panel also above the
floor — and pass with `tier0_status=PASS_WITH_CAVEAT`.

`target2cellsRep2_S4` is the case the gate exists for: VR family signal can
survive a library that fails every other check, and it is not a receptor call.

**3. "Rep1/Rep2/Rep3" are independent libraries from different cells**, not
biological replicates. Low correlation between them is expected and is not a
QC failure.

---

## For a new user — start here

Two ways to run this. Pick one; they produce identical results.

### Option A — plain terminal (no Claude Science needed)

Works on a laptop, a login node, or inside a SLURM job.

```bash
# 1. unpack
tar -xzf vr_pipeline_standalone.tar.gz && cd vr_pipeline

# 2. dependencies (pick whichever fits your machine)
conda create -n vr python=3.11 pandas numpy matplotlib pyyaml scipy -c conda-forge
conda activate vr
conda install -c bioconda samtools
#   or: module load python3 && module load samtools        # HPC
#   or: python3 -m venv ~/vr-env && . ~/vr-env/bin/activate
#       pip install pandas numpy matplotlib pyyaml scipy   # + samtools via your pkg manager

# 3. check the install before spending a run
python3 bin/vr_analyze.py --selftest      # expect: 13 passed, 0 failed

# 4. run it — one command
python3 bin/vr_analyze.py \
    --results /path/to/your/nfcore/rnaseq/results \
    --gtf     /path/to/the/genes.gtf \
    --trial   myrun
```

`--results` is the nf-core/rnaseq output directory (the one containing
`star_salmon/`). `--gtf` is the annotation that run used — find it in
`<results>/pipeline_info/params_*.json` under `"gtf"`. Outputs land in
`./vr_out_myrun/results/`.

On a SLURM cluster, submit it rather than running on the login node:

```bash
sbatch --partition=<your partition> --time=60 --mem=16G --cpus-per-task=8 \
       --wrap="python3 bin/vr_analyze.py --results <dir> --gtf <gtf> --trial myrun"
```

### Option B — Claude Science

Install the skill once, then ask in plain language. The agent runs the pipeline,
reads the outputs, and applies the interpretation rules without you learning the
flag vocabulary.

**Install (once per user).** In a Claude Science session, in the `repl` tool:

```python
exec(open("/abs/path/to/VNOReceptor-agent/install.py").read())
```

That reads the shipped `SKILL.md` and `kernel.py`, writes them into your own
skill catalog, and publishes. It is idempotent — re-run it to pick up an update.
If you would rather not run a script, the two files are plain text under
the repository root; paste them in with `host.skills.edit(...)` and
`host.skills.publish("vno-receptor-rnaseq")`.

**Verify:**

```python
skill({skill: "vno-receptor-rnaseq"})
```

The result should list the `vno_*` helpers now loaded in your python kernel.
Then, in a `python` cell:

```python
vno_install_check()     # what is missing, with the fix for your platform
```

**Use it.** Once installed, plain requests work:

> Run the VNO receptor pipeline on `/data/run7/results`, GTF at
> `/refs/GRCm38/genes.gtf`. Tell me which receptor clusters each library
> expresses and flag anything that failed QC.

The agent will validate the inputs, run the stages, and report by tier — it will
not hand you a receptor name from a library that failed QC, because the pipeline
refuses to produce one.

Useful things to ask for by name: `tier_status_all.tsv` (which libraries support
anything), `vr_report_all.md` (the narrative), `docs/VR_OUTPUT_GUIDE.md` (what
every output column means).

### What a new user needs to supply

| input | where it comes from | notes |
|---|---|---|
| nf-core/rnaseq results dir | your own run, version >= 3.18, `--aligner star_salmon` | must contain `star_salmon/salmon.merged.gene_counts.tsv` and the per-sample BAMs |
| the reference GTF | whatever that run used | read it out of `pipeline_info/params_*.json`; a *different* GTF than the run used will give wrong coordinates |
| which samples are which | you | `target` = the sorted population you want receptor calls from (e.g. GFP+); `nontarget` = the negative control |

Sample roles are inferred from names (`target*`, `nontarget*`, `pool*`,
`GFPpos`, `GFPneg`, `*ctrl*`, `sorted*`). If a name matches neither vocabulary
the run **stops with exit 4** and asks you to say which is which:

```bash
python3 bin/vr_analyze.py --results <dir> --gtf <gtf> \
    --target 'sampleA|sampleB' --nontarget 'input|ctrl'
```

It refuses rather than guessing because reversing those two inverts the entire
sort-validation tier — every downstream verdict would be confidently wrong.

### Exit codes (they carry meaning)

| code | meaning | what to do |
|---|---|---|
| `0` | ran; at least one library reached a cluster-level call | read the report |
| `3` | ran; **nothing survived QC** | this is a *result*, not a crash — read the stop reasons |
| `4` | sample roles could not be inferred | rerun with `--target`/`--nontarget` |
| `5` | missing dependency | install what it names |
| `6` | bad input | not an nf-core `star_salmon` tree, or the GTF is unreadable |

Do not wrap this in a retry loop on any non-zero exit. `3` will never succeed on
a retry, and treating it as a failure hides a real finding.

### Species and genome scope

The shipped `ref/` tables are **mouse GRCm38** with Ensembl chromosome naming
(`1`, `7`, `X` — not `chr1`). For a different build, add `--force-ref` and the
VR reference is rebuilt from your GTF (one streaming pass, a few minutes). For a
different species, the marker panels in `config/project.template.yaml`
(`Trpc2`, `Omp`, `Gnai2`, `Gnao1`, `Adcy3`, `Cnga2`, `Gnal`, `Olfr*`) and the
`Vmn1r`/`Vmn2r` family patterns need editing first — the analysis logic
generalizes, the gene names do not.

### Two flags worth understanding before you use them

**`--no-bam`** runs without `samtools` and **degrades the result**: the
unique-read evidence channel disappears, so no individual-receptor call is
possible and EM-redistribution flags fall back to fraction-only evidence
(measured: 20 named candidates instead of 25 on the reference dataset).
Cluster-level results are unaffected. Use it only when BAMs are genuinely
unavailable.

**`--force-ref`** rebuilds the VR gene and cluster tables from your GTF. Needed
when the genome build differs from the shipped tables; otherwise it just costs
you the parse time.

## Install — reference detail

**Nothing to build.** The pipeline is standalone Python modules plus a bash
driver, importing only pandas, numpy, matplotlib, pyyaml and scipy from the
standard library outward, and shelling out to `samtools`. It runs on a laptop,
in a conda env, or on an HPC — no site-specific paths.

```bash
# conda (recommended if you have it)
conda create -n vr python=3.11 pandas numpy matplotlib pyyaml scipy -c conda-forge
conda activate vr
conda install -c bioconda samtools

# or venv
python3 -m venv ~/vr-env && . ~/vr-env/bin/activate
pip install pandas numpy matplotlib pyyaml scipy
# samtools separately: apt install samtools / brew install samtools

# or an HPC with Environment Modules
module load python3 && module load samtools
```

`bin/run_pipeline.sh` detects which of these you have and, if the interpreter
cannot import what it needs, prints the missing package names with all three
remedies rather than assuming a module system. To force a specific interpreter:
`VR_PYTHON=/path/to/python bin/run_pipeline.sh ...`

**`samtools` is separate from the python stack and easy to miss.** It powers the
unique-read evidence channel (`samtools view -q 255`). Without it
`vr_quantify.py` stops at a preflight guard naming the fix. `--no-bam` runs
anyway but **degrades the result**: no individual-receptor call is possible and
EM-redistribution flags fall back to fraction-only evidence (verified: 22
candidate rows instead of 26 on trial 2). Cluster-level results are unaffected.

If your site provides these as environment modules, load them instead of
creating an environment; `deploy/DEPLOY.md` covers installing the pipeline
itself as a shared module so users need neither step.

## Quick start

### On an existing work folder

```bash
cd <work folder>                    # the directory holding bin/ and config/
bin/run_pipeline.sh --all --finalize
```

### Creating a work folder anywhere

`bin/vr_init.py` scaffolds one and writes a valid config from your data:

```bash
python3 bin/vr_init.py --dest ~/vr_run --trial trial1 \
    --results /path/to/nfcore/results \
    --gtf /path/to/genes.gtf --genome-key my-grcm38
cd ~/vr_run && bin/run_pipeline.sh --trial trial1 --finalize
```

It discovers sample names from the quant table header, guesses `cell_type` from
them, and validates the generated config through the pipeline's own `vr_config`
before exiting — so a bad config fails at init, not five stages later.

**Check the guessed `cell_type` values.** `target` (GFP+) vs `nontarget` (GFP−)
cannot be inferred reliably from a filename, and the sort-validation tier cannot
interpret `UNKNOWN`. Add later trials with `--add-trial NAME --results DIR`,
which inserts inside the existing blocks rather than appending a duplicate
top-level key.

The `paths`, `markers`, `thresholds` and `qc_thresholds` blocks are copied
**verbatim** from `config/project.template.yaml`, never retyped, so a generated
config produces identical verdicts to this reference run. That is not a
stylistic choice: an earlier hand-typed threshold block completed without error
and produced *different verdicts* on identical input. Edit the template to
change defaults for future work folders.

Verified by two separate tests, which checked different things:

| test | what was compared |
|---|---|
| fresh work folder outside this project tree (`vr_init.py` + local `bin/`) | field-by-field against the canonical run: `tissue_verdict`, `sort_verdict`, `library_status`, `population_call`, `qc_overall`, `suppress_biology`, `em_flag`, `pseudogene_bleed_flag`, `interpretation`, and every candidate's `rank`/`gene_name`/`confidence`/`confirmation_status` — all identical |
| clean tarball unpack, arbitrary trial name | ran to completion (exit 0) and produced 26 candidate rows all `tentative_unconfirmed`, with tissue verdicts matching (2 `VNO`, 2 `VNO_dominant_mixed`, 2 `no_tissue_signal`). This test did **not** re-do the field-by-field comparison. |

So: the generated-config path is verified equal on every verdict column, and the
packaged tarball is verified to run and agree on the columns that were checked.

### Work-folder resolution

The driver derives its root from its own location, so `bin/run_pipeline.sh`
works from anywhere. Python modules and kernel helpers resolve in this order:
`$VR_WORK`, then a directory containing `config/project.yaml`, then `$PWD`.
There is no built-in default path.

### The driver

`run_pipeline.sh` runs every stage in dependency order and is the recommended
entry point — the step-by-step commands further down exist for debugging a
single stage. Its own `--help` is authoritative:

```bash
bin/run_pipeline.sh --help          # prints the stage rationale and all options
bin/run_pipeline.sh --all           # every registered trial, then finalize
bin/run_pipeline.sh --trial trial1  # one trial, no cross-trial finalize
bin/run_pipeline.sh --all --stage quant --stage report   # re-run two stages only
bin/run_pipeline.sh --refcheck      # rebuild ref/ into a scratch dir and diff
```

An unregistered trial name stops immediately with the config block to paste —
no stages run and no traceback.

Stages are `ref qc quant report finalize`. Two notes that matter:

- **`qc` always runs `--all-trials`, by design.** It writes the combined
  `results/sample_qc_all.tsv` that quantification and reporting both read;
  running it per trial would leave that table holding one trial and make the
  result depend on invocation order.
- **`ref` is skipped when `ref/vr_gene_annotation.tsv` exists.** It is
  trial-independent — rebuild it only if the genome build changed
  (`--force-ref`), or verify it non-destructively with `--refcheck`.

A cold start (empty `results/`) and a re-run over existing outputs must produce
byte-identical tables. That property has been tested over 36 paths, with this
exact result:

| outcome | n | detail |
|---|---|---|
| `IDENTICAL` | 31 | every TSV and all 6 figures, byte-for-byte |
| `TIMESTAMP_ONLY` | 3 | the three markdown reports, differing only in a generation-timestamp line |
| `EXTRA_NEW` | 1 | `refcheck_diff.txt`, produced by the verification run itself |
| `MISSING_NEW` | 1 | **a known bug — see below** |

Zero `NUMERIC_ONLY` rows: every float matched exactly, not merely within
tolerance. The reference layer was additionally rebuilt from the 852MB GTF into
a scratch directory and all three `ref/` tables came back identical.

> **Known bug: one deliverable is not regenerable.** `results/vr_cluster_map.png`
> has **no producer anywhere in `bin/`** — it was produced out-of-band during
> development, so a cold start does not recreate it. The committed `results/`
> tree therefore contains one artifact the pipeline cannot reproduce from
> `config/project.yaml` + data. Resolution options are recorded in
> `results/PROVENANCE_NOTES.md`; until one is taken, do not treat that figure as
> a pipeline output. Everything else under `results/` is reproducible.

### Running it from an agent session

The procedure is also packaged so a Claude Science session can drive it without
re-deriving any of the conventions:

```python
skill({skill: "vno-receptor-rnaseq"})   # loads guidance + kernel.py helpers
```

Then, in a `python` cell:

```python
vno_install_check()                       # what is missing, with the remedy
vno_validate_work_dir("~/vr_run")         # is this folder usable?
vno_validate_run_dir("<nfcore>/star_salmon")   # check inputs before spending a job
vno_run_plan(work_dir="~/vr_run", trial="trial1")   # exact commands
res = vno_load_results(work_dir="~/vr_run", trial="trial1")  # every output table
vno_clearance(res["sample_qc"])           # the reporting gate
```

If the data lives on a cluster, set `$VR_SSH_TARGET` and reach it with
`host.compute.create(vno_ssh_target())` from the `repl` tool; `call_command`
has a short timeout there, so any stage that parses the GTF or reads BAMs must
go through `submit_job`. For a local run none of this applies.

There is also a **VNO Receptor Specialist** agent profile with the
interpretation rules built in (see *Reusable packaging* below). Switch to it
from the session config selector when the work is VNO receptor analysis rather
than general RNA-seq.

---

## Layout

```
vr_analysis/
  bin/         pipeline modules (versioned, re-runnable, config-driven)
  ref/         VR gene annotation + genomic cluster reference — see ref/README.md
  config/      project.yaml  <- the path/threshold contract; READ IT FIRST
  results/     per-trial outputs + figures/ + combined tables
  docs/        QC_MODULE_README.md         (tissue / sort / population detail)
               VR_PIPELINE_FLOWCHART.md   (what each step does and why)
               VR_OUTPUT_GUIDE.md         (how to read every output)
  logs/        build logs
```

### `bin/` modules

| module | role |
|---|---|
| `vr_init.py` | scaffold a work folder anywhere; writes and validates config/project.yaml from your data |
| `run_pipeline.sh` | **the entry point** — loads both modules, runs every stage in dependency order, `--help` documents each stage and why it is ordered where it is |
| `vr_config.py` | loads `config/project.yaml`; every path, sample name and threshold comes through here. Nothing downstream hardcodes any of them. |
| `build_vr_reference.py` | streams the 852MB GTF once → `ref/` tables |
| `vr_qc.py` | STAR / qualimap / MultiQC parsing → `technical_qc.tsv` |
| `vr_markers.py` | CPM, tissue identity, sort validation, actin sanity, population call → `marker_cpm.tsv` |
| `vr_sample_qc.py` | joins the two, derives `blocking_flags` → `sample_qc.tsv` |
| `vr_quantify.py` | cluster-level VR expression, within-cluster fractions, EM-artifact detection, candidate nomination |
| `vr_report.py` | the tiered reporting engine; enforces the priority hierarchy structurally |
| `vr_figures.py` | the three receptor-analysis deliverable figures |
| `vr_qc_figures.py` | the three QC deliverable figures |

Cross-references, not duplicated here: **`docs/VR_OUTPUT_GUIDE.md`** (what every output file and figure tells you, and how to read it), **`docs/VR_PIPELINE_FLOWCHART.md`** (flowcharts of every decision layer with the rationale for each gate), **`ref/README.md`** (reference-layer
provenance, cluster construction, threshold sweep) and
**`docs/QC_MODULE_README.md`** (tissue-identity, sort-validation and
population-ID detail).

---

## Running the pipeline stage by stage

Use this when debugging one stage; otherwise prefer `bin/run_pipeline.sh` above.

### Environment

Activate whatever environment you installed into (see **Install** above),
or let `bin/run_pipeline.sh` detect it. `samtools` must be on PATH for the
unique-read channel; `vr_quantify.py` stops with a fixable message if it is
not, and `--no-bam` degrades the result (22 candidate rows instead of 26 on
trial 2 — cluster-level results unaffected).

Check what is missing without running anything:

```python
vno_install_check()   # kernel helper from the vno-receptor-rnaseq skill
```

**Step 1 — register the trial in `config/project.yaml`.** Add it under
`trials:` (`results:`, `fastq:`, `platform:`) and under `samples:` with
`cell_type` (`target` = GFP+, `nontarget` = GFP−), `n_cells`, `prep_status`.

> **Edit inside the existing blocks.** A naive append that creates a *second*
> top-level `markers:` or `thresholds:` key silently shadows the first — YAML
> resolves last-wins with no error and no warning. This already broke the
> config once; the pre-repair backup is `config/project.yaml.bak_pre_tissue`.
> After any edit, re-read the file and confirm each top-level key appears
> exactly once (`grep -n '^[a-z_]*:' config/project.yaml`).

Note that the upstream nf-core results directories are spelled `trail1` /
`trail2` (sic). That typo is confined to `config/project.yaml`; no module
constructs it.

**Step 2 — reference layer**, only if the genome build changed:

```bash
python bin/build_vr_reference.py --config config/project.yaml --outdir ref
```

**Step 3 — QC and gating** (must run before anything else reads biology):

```bash
python bin/vr_sample_qc.py --all-trials
```

**Step 4 — quantification:**

```bash
python bin/vr_quantify.py --trial <trial> --threads 8
```

`--no-bam` skips the MAPQ255 unique-read evidence. It runs much faster, but
**no individual-receptor call is possible without it** — the candidate ranking
and the EM unique-read gate both depend on that channel.

**Step 5 — report** (self-test first; it asserts the gate cannot be bypassed):

```bash
python bin/vr_report.py --selftest     # 35 assertions
python bin/vr_report.py
```

**Step 6 — figures:**

```bash
python bin/vr_figures.py
python bin/vr_qc_figures.py --all-trials --fig-dir results/figures
```

On SLURM: `--partition=shared`, account `it_mgmt`; use `highmem` only above
64G. Anything that parses the GTF or reads BAMs belongs in a batch job, not an
interactive login-node command.

### What to check first, in this order

Do **not** open the candidate table first.

1. `results/tier_status_all.tsv` → `tier0_status`, then
   `highest_tier_reported`. One line per library answering "did the tissue
   check pass, and what does this library actually support?" Use this rather
   than reading `tissue_verdict` yourself: the gate accounts for `cell_type`,
   which a bare verdict string does not (see the note above).
2. `results/<trial>/marker_cpm.tsv` → `tissue_verdict` and `tissue_reason` for
   *why*, and `library_status`: `FAILED` or `DEGENERATE` ends that sample.
3. `results/<trial>/vr_cluster_expression.tsv` where `is_called == 1` — the
   reliable tier.
4. `results/<trial>/vr_artifact_flags.tsv` → `em_flag`, **before** reading any
   individual gene name.
5. `results/<trial>/vr_candidates.tsv` last, and read the evidence columns
   alongside every rank.

---

## The priority hierarchy

Enforced **structurally**, not by convention. Each tier's content is produced
by a callable that `TierGate.emit()` refuses to invoke unless every upstream
tier passed for that sample. There is no code path from a wrong-tissue or
failed library to a receptor call — the producer function is never called and
the gate substitutes an explicit suppression line naming the failing tier.
`results/tier_gate_selftest.txt` records the 35 assertions that prove it.

| tier | name | question |
|---|---|---|
| 0 | `tissue_identity` | is this VNO at all? |
| 1 | `sort_validation` | viable GFP+/GFP− library? |
| 2 | `population_id` | V1R or V2R (Gnai2:Gnao1)? |
| 3 | `cluster_vr` | which receptor cluster? **the reliable tier** |
| 4 | `individual_vr` | which paralog? **always tentative** |

Status vocabulary: `PASS`, `PASS_WITH_CAVEAT`, `FAIL`, `NO_DATA`, `SUPPRESSED`.

**Clearance gate for reporting biology:** `qc_overall == 'USABLE'` **AND**
`suppress_biology` false. `suppress_biology` alone is *not* sufficient — two
libraries here carry `suppress_biology=False` with `sort_verdict=FAIL`.

---

## Output files

### `results/` (combined, both trials)

| file | contents |
|---|---|
| `sample_qc_all.tsv` | one row per library: technical metrics + marker CPMs + `tissue_verdict`, `sort_verdict`, `population_call`, `library_status`, `qc_overall`, `suppress_biology`, `blocking_flags` |
| `sample_qc_summary.txt` | human-readable version of the above |
| `tier_status_all.tsv` | one row per library: `highest_tier_reported`, `stopped_at_tier`, `stop_reason`, and the per-tier status columns |
| `tier_outcomes_all.tsv` | one row per (library × tier): `status`, `blocked_by_tier`, `reason` |
| `tier_gate_selftest.txt` | the gate self-test transcript |
| `vr_report_all.md` | the combined tiered report — **the document to read** |
| `vr_cluster_map.png` | genomic map of the VR clusters **(NOT regenerable — no producer in `bin/`; see the known-bug note in Quick start and `results/PROVENANCE_NOTES.md`)** |

### `results/<trial>/`

| file | contents |
|---|---|
| `technical_qc.tsv` | STAR/qualimap/MultiQC metrics, `tech_flags`, `tech_verdict` |
| `marker_cpm.tsv` | marker CPMs, tissue panels, sort/population/library verdicts, `ratio_low_support` |
| `sample_qc.tsv`, `sample_qc_summary.txt` | the joined gating table for this trial |
| `vr_cluster_expression.tsv` | **the reliable tier.** One row per (sample × cluster) at both the 200kb `cluster` and 800kb `supercluster` tier: `counts_sum`, `cpm_sum`, `share_of_sample_vr`, `is_called`, plus the QC context that governs it |
| `vr_within_cluster_fractions.tsv` | one row per (sample × VR gene): EM counts, `frac_of_cluster`, and the unique-read channels (`unique_reads_bam`, `unique_reads_bam_nodup`, `unique_reads_salmon_tx`) |
| `vr_artifact_flags.tsv` | one row per (sample × cluster) plus sample-scope rows: evenness statistics, the Monte Carlo `p_uniform` and power, the even-block/unique-read gate columns, `em_flag`, pseudogene-bleed columns, expected-pattern checks |
| `vr_candidates.tsv` | nominated individual receptors with rank, evidence channel, read support, `evidence_contradiction`, `confidence`, `confirmation_status`, `marker_consistency` |
| `tier_status.tsv`, `tier_outcomes.tsv` | per-trial slices of the combined tier tables |
| `vr_report.md` | the per-trial tiered report |

Every pipeline TSV carries a leading `#` provenance line recording the CPM
convention and source table. Read with `comment='#'` — or use
`vno_read_table()` from the `vno-receptor-rnaseq` skill, which returns the
provenance string alongside the frame.

### `results/figures/`

`marker_cpm_by_sample.png`, `tissue_identity_panel.png`,
`technical_qc_panel.png`, `vr_cluster_heatmap.png`,
`vr_within_cluster_fractions.png`, `vr_tier_overview.png`.

QC-suppressed samples are drawn as hatched rows carrying their failing tier,
never silently dropped — an absent row would read as "no receptors found".

---

## Verdict vocabularies (exhaustive)

Read these columns; never recompute the verdict, and never introduce a value
outside the list.

| column | allowed values |
|---|---|
| `tissue_verdict` | `VNO`, `VNO_dominant_mixed`, `MOE`, `MOE_dominant_mixed`, `ambiguous_mixed`, `no_tissue_signal` |
| `sort_verdict` | `PASS`, `FAIL`, `FAIL_WRONG_TISSUE` |
| `population_call` | `V1R_dominant`, `V2R_dominant`, `undetermined` |
| `library_status` | `OK`, `FAILED`, `DEGENERATE` |
| `tech_verdict` | `PASS`, `WARN`, `FAIL` |
| `qc_overall` | `USABLE`, `UNUSABLE` |
| tier `status` | `PASS`, `PASS_WITH_CAVEAT`, `FAIL`, `NO_DATA`, `SUPPRESSED` |
| `em_flag` | `no_signal`, `insufficient_signal`, `no_redistribution_signature`, `single_paralog_only`, `suspected_em_redistribution` |
| `em_flag_level` | `none`, `strong` |
| `pseudogene_bleed_flag` | `no_pseudogene_members`, `no_apparent_pseudogene_expression`, `apparent_pseudogene_expression`, `pseudogene_only_cluster` |
| `evidence_type` | `bam_unique_mapq255`, `no_cluster_above_signal_threshold` |
| `confidence` | `moderate`, `alternative_candidate`, `unresolvable`, `no_call` |
| `confirmation_status` | `tentative_unconfirmed` (the only value, by construction) |
| `interpretation_context` | `monogenic_expectation` (2-cell), `pooled_ambiguous` (100-cell) |

Across both trials the observed `em_flag` distribution is: trial 1 —
`no_signal` 216, `insufficient_signal` 15, `single_paralog_only` 1; trial 2 —
`no_signal` 293, `insufficient_signal` 29, `no_redistribution_signature` 17,
`single_paralog_only` 8, `suspected_em_redistribution` **1**. That single
flagged artifact is described at the end of this document.

## The CPM convention (settled — do not re-derive)

**CPM = count / (column sum over ALL genes of the *unscaled*
`salmon.merged.gene_counts.tsv`) × 1e6.**

The `gene_counts_scaled` and `gene_tpm` tables give values 2.3–2.5× lower.
Every threshold below is calibrated against the unscaled all-gene-denominator
convention; substituting another table silently invalidates all of them.

---

## Thresholds and where to change them

All of these live in `config/project.yaml`. Nothing is hardcoded in `bin/`.
`vr_config.threshold(cfg, key, section=...)` is the accessor; the QC block has
fallback defaults in `vr_config.TECH_THRESHOLD_DEFAULTS` so the modules still
run against an older config.

### `thresholds:` — biology gates

| key | value | meaning |
|---|---|---|
| `target_trpc2_min` | 1000 | Trpc2 CPM a GFP+ library should exceed |
| `target_trpc2_concern` | 100 | below this, sort validation fails |
| `nontarget_trpc2_max` | 10 | Trpc2 CPM ceiling for a GFP− library |
| `failed_lib_actin_cpm_max` | 100 | actin sum below this **and** Trpc2 below the next → FAILED |
| `failed_lib_trpc2_cpm_max` | 10 | " |
| `v1r_dominant_ratio_min` | 2.0 | Gnai2:Gnao1 above this → V1R population |
| `v2r_dominant_ratio_max` | 0.5 | below this → V2R population |
| `nontarget_total_vr_cpm_max` | 100 | sort-purity ceiling on total VR signal in a GFP− library |
| `expected_clusters` | 0/0, 1–2, 2-cell; 5–25, 100-cell | expected called-cluster count by sample type |
| `cluster_max_gap_bp` | 200000 | the 200kb cluster tier (see caveat below) |
| `tissue_panel_floor_cpm` | 100.0 | absolute CPM floor: below it a marker panel is noise and a ratio of two noise values carries no tissue information |
| `tissue_dominance_ratio` | 3.0 | with both panels above the floor, one must exceed the other by this factor |

### `qc_thresholds:` — library viability and alignment

| key | value | meaning |
|---|---|---|
| `min_input_reads` | 1e6 | minimum sequenced reads |
| `min_uniquely_mapped_pct` | 50.0 | STAR unique-mapping floor |
| `min_assigned_counts` | 1e6 | **below this the library is `DEGENERATE`** and biology is suppressed regardless of how large its CPMs look |
| `max_multi_loci_pct` | 40.0 | STAR multi-loci warning |
| `max_too_many_loci_pct` | 5.0 | STAR too-many-loci warning |
| `max_unmapped_too_short_pct` | 40.0 | |
| `max_mismatch_rate_pct` | 2.0 | |
| `bias_5p3p_low` / `_high` | 0.5 / 2.0 | qualimap 5′–3′ bias; healthy ≈ 1 |
| `min_exonic_pct` | 40.0 | |
| `max_intergenic_pct` | 30.0 | |
| `ratio_min_support_reads` | 10 | raw denominator reads a Gnai2:Gnao1 ratio needs before its **magnitude** is quotable |

`min_assigned_counts` earns its own note. CPM is a ratio: `pool100cells_S8`
reports 39,881 CPM of actin from only 60,744 assigned counts. A near-empty
library produces impressive-looking CPMs that mean nothing.

`ratio_min_support_reads` likewise. The recorded values "422.7:1" and "385:1"
for `target100cellsRep1_S6` are the *same measurement* with `Gnao1 = 7 reads`;
`ratio_low_support` marks it, the direction of the population call still holds,
and the magnitude must never be quoted as precise.

### `markers:` — panels

`pan_vno_neuron` (Trpc2), `v1r_population` (Gnai2), `v2r_population` (Gnao1),
`housekeeping_actin` (6 actins), `vr_families` (Vmn1r, Vmn2r),
`moe_receptor_families` (Olfr), `moe_transduction` (Adcy3, Cnga2, Gnal),
`tissue_shared` (**Omp — never a VNO marker**).

`Gnai2`/`Gnao1` are deliberately excluded from the tissue panels: they split
V1R vs V2R *within* the VNO and are broadly expressed elsewhere (`Gnai2` reads
165 CPM in the GFP− library that carries no VR signal at all). Tier-2 only.

---

## Cluster tiers, and the caveat

Two tiers, both in `ref/vr_clusters.tsv` and `ref/vr_gene_to_cluster.tsv`:
`cluster_id` at 200kb max gap (24 V1R, 37 V2R clusters) and `supercluster_id`
at 800kb.

**200kb is not a natural break.** The V1R inter-gene gap distribution has its
KDE minimum near 2Mb; 200kb was retained as the conservative choice and the
full threshold sweep is in `ref/vr_gtf_parse_report.txt`. The consequence is
concrete: a 217,366bp gap between `Vmn1r132` and `Vmn1r135` splits the chr7
V1R megacluster into `V1R_chr7_cl015` + `V1R_chr7_cl016`, which
`V1R_chr7_sc013` reunites. **Report both tiers for that region** — rows carry
`chr7_dual_tier_region = 1`. Historical labels `cl039` and `cl029` from earlier
notes are aliases: `cl029` = `V1R_chr6_cl008` (9 genes), and `cl039` =
`V1R_chr7_cl015` + `V1R_chr7_cl016`.

**V2R aggregation protects less than V1R**: 18 of 37 V2R clusters are
singletons and only 180 of 222 V2R genes sit in clusters of ≥5. Aggregation
cannot absorb ambiguity that has nowhere to go.

Two join rules, both non-negotiable and both enforced in code:

- **Join on `gene_id`, never `gene_name`.** `Vmn1r-ps5` and `Vmn2r118` each map
  to two distinct `gene_id`s on different chromosomes.
- **Filter `is_primary_assembly == 1`** (538 of 541 genes). This is what
  reconciles exactly to the 318 `Vmn1r` / 220 `Vmn2r` rows in the quant tables.

---

## Known limitations of individual-receptor calls

Read this before quoting any receptor name out of `vr_candidates.tsv`.

**Every candidate row carries `confirmation_status = tentative_unconfirmed`.**
That is by construction, not a backlog of work to close: 75bp reads inside an
85–95%-identical paralog cluster cannot be resolved to a single gene by
short-read quantification alone. Confirmation requires an orthogonal assay
(long reads, targeted amplicon, in-situ).

**Salmon's EM is not an independent witness.** STAR is not filtering the
ambiguity out for us either: retained multiple-loci reads run 1.66–36.35% of
input against 0.02–0.30% discarded as too-many-loci (median ratio 125×). The
ambiguity lands squarely in the EM, which is why cluster aggregation is
mandatory and why paralog-level claims need EM-independent evidence.

**Ranking rests on `bam_unique_mapq255`** (MAPQ 255 = STAR placed the read at
exactly one locus), counted mate-wise with and without duplicates. Salmon's
`aux_info/ambig_info.tsv` `UniqueCount` is a cross-check only: its "unique"
means unique to a *transcript*, so gene-level support is understated for the 28
of 538 VR genes with multiple transcripts. Full `eq_classes.txt` was not
written by this nf-core run.

**Several ranks rest on very few reads.** `top_over_background_floor` can be
`inf` against a zero background. Read the evidence columns before quoting a
rank.

**An even within-cluster split is necessary but NOT sufficient for an EM
artifact.** Within a detected even block, each member's unique MAPQ255
deduplicated reads are compared against `max(10 reads, 3 × median unique count
of cluster members outside the block)` — the cluster's own measured mismapping
floor. If exactly one co-dominant member clears it, that is redistribution. If
two or more clear it, that is **co-expression**, because monogenic choice is a
per-*cell* rule and a multi-cell pool legitimately captures two paralogs of one
cluster. Uniformity itself is tested by Monte Carlo from `Multinomial(N, 1/k)`
(4000 draws), not the asymptotic chi-square, which is invalid at N of tens.

`p_uniform` **saturates at 2.4994e-04 = 1/4001** — the permutation floor. Quote
it as `p < 2.5e-4`, a bound, never as a point estimate.

`em_flag` vocabulary, exhaustive: `no_signal`, `insufficient_signal`,
`no_redistribution_signature`, `single_paralog_only`,
`suspected_em_redistribution`. **Read the column; never recompute the verdict.**
A figure that applied its own visual evenness criterion instead of reading
`em_flag` mislabelled `target100cellsRep2_S7` / `V1R_chr7_cl013` as an EM
artifact — it is real co-expression, with `Vmn1r89` and `Vmn1r87` each carrying
~18.5k independent unique reads. The figure now reads the column.

**The one flagged artifact.** `target2cellsRep1_S3` / `V1R_chr7_cl016` is
`suspected_em_redistribution` (strong): EM fractions 0.500/0.498 across
`Vmn1r166`/`Vmn1r138` over 8,962 reads in the block, but unique MAPQ255
deduplicated reads are 51 vs 0. Preserve the asymmetry when reporting it — the
*call* is solid (a paralog holding 4,473 EM counts with zero single-locus reads
is not independently observed), while *which* paralog is the source rests on 51
unique reads, ~1.1% of `Vmn1r166`'s apparent expression. Related: the 99.9%
figure from earlier notes reproduces at the 800kb tier (`V1R_chr7_sc013` =
99.91% of that sample's VR signal), within which the structure is 66.8%
`cl016` / 33.1% `cl015`.

**When the two channels disagree, there is no call.** Same sample, second
cluster: `V1R_chr7_cl015`'s EM-dominant paralog `Vmn1r131` holds 99.9% of the
cluster's EM signal (4,452 counts) with **zero** unique reads, while
`Vmn1r103`/`Vmn1r104` carry 92 unique reads each and zero EM counts. Neither
channel can name this receptor; `evidence_contradiction` is populated and the
honest output is no individual call for that cluster.

**Pseudogene bleed is an unresolved mechanism, not a quantified correction.**
Largest instance: `V1R_chr17_cl021` in `target100cellsRep1_S6`, where
`Vmn1r-ps150` holds 38.0% of cluster signal with 2,824 unique reads — and that
cluster is *clean* on the redistribution test, so an even-split artifact does
not explain it. The two candidate mechanisms (EM leakage from an expressed
functional paralog vs genuine transcription of the pseudogene locus via
cluster-shared regulatory elements) are not separable from this data;
`pseudogene_mechanism` says so in every row.

Do not cite Dietschi et al. 2022 as a quantitative mouse prior for this. The
figures below were read out of the paper's own text and figure legends —
Dietschi Q. *et al.*, *Clustering of vomeronasal receptor genes is required for
transcriptional stability but not for choice*, Sci Adv 8(46) eabn7450,
doi:10.1126/sciadv.abn7450, PMC9668312, PMID 36383665 (open access) —
and each is quoted with the test it came from, because the distinction between
them is what matters here:

- **Pseudogene transcription inside vs outside mixed clusters** is a Wilcoxon
  rank sum test, significant in **rat only** (W = 1113, P = 0.002987) and **not
  significant in mouse** (W = 1214, P = 0.5704), on n = 8 mouse VNOs (141
  pseudogenes in mixed clusters, 16 outside) and n = 7 rat VNOs. Their proposed
  mechanism is regulatory — cluster-associated transcription-stabilizing
  elements — not multi-mapping. So this paper gives us **no** quantitative mouse
  expectation for pseudogene bleed, and cannot adjudicate our
  `V1R_chr17_cl021` case.
- What their work **does** support is the design choice behind our cluster
  tier: the effect of cluster identity on expression level, for functional
  genes of the mixed-cluster group, is ANOVA **P < 2 × 10⁻¹⁶**; and the effect
  of cluster membership on VSN abundance (choice probability) is ANOVA
  **P = 0.000232** (7 mice, 98 different VSN types — the legend prints this as
  "P = 7 mice", evidently a typo for *n*).
- On assignability they estimated, by counting copy numbers of every 100-nt
  sequence across an 8-kb V1rD transcript set, that **only 57%** of RNA-seq
  reads from V1rD transcripts could be assigned without ambiguity. Their
  response was STAR 2.7.1a with `--outFilterMultimapNmax 4` (retaining an
  estimated 91.6% of read alignments) and fractional counting of multimapped
  reads in featureCounts 1.6.5 via `-M --fraction`. Note this is a *different*
  strategy from ours: they distribute multimapped reads fractionally at the
  counting step, we aggregate to the cluster and require MAPQ255 single-locus
  evidence for any paralog-level claim. Their 57% is a useful order-of-magnitude
  check on how much ambiguity exists in this gene family, not a benchmark our
  numbers should match.

---

## Reusable packaging

The procedure is packaged as the **`vno-receptor-rnaseq`** skill (input
contract, CPM convention, tier hierarchy, cluster tiers, the EM unique-read
gate, flag vocabularies, and the failure modes above), with a `kernel.py`
sidecar of thin helpers: `vno_paths`, `vno_validate_run_dir`, `vno_read_table`,
`vno_cpm`, `vno_load_results`, `vno_clearance`, `vno_em_verdict`,
`vno_flag_vocabulary`, `vno_check_yaml_duplicate_keys`, `vno_run_plan`. The
heavy modules are not duplicated there — this work folder stays the source of
truth.

A **VNO Receptor Specialist** agent profile encodes the interpretation rules
(tissue check before sort validation, Omp never a VNO marker, conservative
multi-mapping stance, the priority hierarchy, mandatory uncertainty reporting).

## Conventions for anyone extending this

- Code goes in `bin/` as versioned, re-runnable scripts that read
  `config/project.yaml`. No notebook-style one-off cells; no sample name or
  path hardcoded inside a function.
- A consumer never re-derives an upstream verdict. Read the flag column.
- Deliverable figures load the `figure-style` skill and call
  `apply_figure_style()` first.
- `config/project.yaml` sections are owned by tracks. Append inside the
  existing block, then re-read the file and verify no top-level key is
  duplicated.
