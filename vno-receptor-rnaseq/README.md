# vno-receptor-rnaseq

Identify which vomeronasal receptor each mouse VNO sensory neuron expresses,
from an nf-core/rnaseq `star_salmon` results directory.

Each VNO neuron expresses **exactly one** vomeronasal receptor. The hard part is
not quantification: ~250 V1R (`Vmn1r*`) and ~120 V2R (`Vmn2r*`) genes sit in
genomic clusters of local duplicates sharing 85–95% nucleotide identity, so
short reads cannot be uniquely assigned within a cluster and Salmon's EM spreads
ambiguous reads across paralogs. Per-gene VR counts inside a cluster are not
trustworthy. This pipeline reports **how far down the chain of claims the data
actually supports you**, tier by tier, and refuses to name a receptor when the
evidence does not reach that far.

This repository is both a **Claude Science skill** and a **standalone command-line
pipeline**. Either works alone.

---

## Use it from the terminal

No Claude Science account needed. Runs on a laptop, a login node, or in a SLURM job.

```bash
git clone <this repo> && cd vno-receptor-rnaseq

# dependencies (pick one)
conda create -n vr python=3.11 pandas numpy matplotlib pyyaml scipy -c conda-forge
conda activate vr && conda install -c bioconda samtools
#   or: module load python3 && module load samtools          # HPC
#   or: python3 -m venv ~/vr-env && . ~/vr-env/bin/activate
#       pip install pandas numpy matplotlib pyyaml scipy     # + samtools separately

# check before spending a run
python3 bin/vr_analyze.py --selftest      # expect: 13 passed, 0 failed

# run it -- one command
python3 bin/vr_analyze.py \
    --results /path/to/your/nfcore/rnaseq/results \
    --gtf     /path/to/the/genes.gtf \
    --trial   myrun
```

`--results` is the nf-core/rnaseq output directory (the one containing
`star_salmon/`). `--gtf` is the annotation **that run used** — read it out of
`<results>/pipeline_info/params_*.json`. Outputs land in `./vr_out_myrun/results/`.

Full documentation: [`PIPELINE.md`](PIPELINE.md).
How to read every output: [`docs/VR_OUTPUT_GUIDE.md`](docs/VR_OUTPUT_GUIDE.md).
What each step does and why: [`docs/VR_PIPELINE_FLOWCHART.md`](docs/VR_PIPELINE_FLOWCHART.md).

### Exit codes carry meaning

| code | meaning | what to do |
|---|---|---|
| `0` | ran; at least one library reached a cluster-level call | read the report |
| `3` | ran; **nothing survived QC** | a *result*, not a crash — read the stop reasons |
| `4` | sample roles could not be inferred | rerun with `--target` / `--nontarget` regexes |
| `5` | missing dependency | install what it names |
| `6` | bad input | not an nf-core `star_salmon` tree, or the GTF is unreadable |

Do not wrap this in a retry-on-nonzero loop. `3` will never succeed on a retry,
and treating it as failure hides a real finding.

---

## Using it from Claude Science or Claude Code

`deploy/AGENT_USAGE.md` covers driving the pipeline in plain language rather
than remembering flags — skill install for Claude Science, and the
`CLAUDE.md`-driven CLI path for Claude Code. Neither needs an administrator.

## Deploying for a whole lab (zero install for end users)

`deploy/DEPLOY.md` covers making this available without per-user installs:

- **terminal users** — install as an Environment Modules module, so a colleague
  runs `module load vno-receptor` then `vr_analyze.py ...` with nothing to set
  up. `deploy/install_shared_module.sh` and the modulefile are included; the
  script was executed against a scratch prefix and the resulting module
  verified to load and run the pipeline (details in `deploy/DEPLOY.md`).
`deploy/DEPLOY.md` covers the shared-install route and says who to ask when the
modulefile directory is owner-only.

## Use it as a Claude Science skill

In a Claude Science session, in the `repl` tool:

```python
exec(open("/abs/path/to/vno-receptor-rnaseq/install.py").read())
```

Then in any session:

```python
skill({skill: "vno-receptor-rnaseq"})
```

The agent gains the pipeline commands, the interpretation rules, and helper
functions (`vno_load_results`, `vno_clearance`, `vno_em_verdict`,
`vno_install_check`, `vno_skill_version`, …) for reading the outputs. Re-run
`install.py` after `git pull` to update; it is idempotent.

Repository layout mirrors the skill contract — `SKILL.md` and `kernel.py` at the
root, so the repo root **is** the skill:

```
SKILL.md         the skill definition (frontmatter name + description, then guidance)
kernel.py        helper functions, loaded into the kernel when the skill loads
install.py       registers the two files above into your catalog
CLAUDE.md        project context for Claude Code (picked up automatically)
deploy/          AGENT_USAGE.md (Claude Science / Claude Code) + shared-module install
bin/             16 modules; vr_analyze.py is the one-command entry point
config/          project.template.yaml — every threshold, one place
ref/             mouse GRCm38 VR gene + cluster tables (prebuilt)
docs/            output guide, flowchart, QC module notes
PIPELINE.md      full pipeline reference (thresholds, outputs, limitations)
```

---

## What you need to supply

| input | where it comes from | notes |
|---|---|---|
| nf-core/rnaseq results | your own run, ≥3.18, `--aligner star_salmon` | needs `star_salmon/salmon.merged.gene_counts.tsv` and the per-sample BAMs |
| the reference GTF | whatever that run used | in `pipeline_info/params_*.json`; a *different* GTF gives wrong coordinates |
| which samples are which | you | `target` = the sorted population you want calls from; `nontarget` = negative control |

Sample roles are inferred from names (`target*`, `nontarget*`, `pool*`,
`GFPpos`, `GFPneg`, `*ctrl*`, `sorted*`). If a name matches neither vocabulary
the run stops with exit `4` and asks you to disambiguate:

```bash
--target 'sampleA|sampleB' --nontarget 'input|ctrl'
```

It refuses rather than guessing: reversing those two inverts the entire
sort-validation tier, and every downstream verdict would be confidently wrong.

---

## Scope

**Species.** The shipped `ref/` tables are **mouse GRCm38**, Ensembl chromosome
naming (`1`, `7`, `X` — not `chr1`). For another build, add `--force-ref` and the
VR reference is rebuilt from your GTF in one streaming pass. For another species,
the marker panels in `config/project.template.yaml` (`Trpc2`, `Omp`, `Gnai2`,
`Gnao1`, `Adcy3`, `Cnga2`, `Gnal`, `Olfr*`) and the `Vmn1r`/`Vmn2r` family
patterns need editing first — the analysis logic generalizes, the gene names do
not.

**What it will not do.** It does not confirm an individual receptor identity. At
75–100bp read lengths the unique-read channel carries roughly 1% of the signal
inside a dense cluster, so every individual call comes back
`tentative_unconfirmed` by construction. Cluster-level calls are the reliable
tier. Confirmation needs evidence that does not pass through the EM step —
longer reads, single-cell libraries, or targeted amplicon sequencing.

**`--no-bam`** runs without `samtools` and **degrades the result**: no
individual-receptor calls are possible and EM-redistribution flags fall back to
fraction-only evidence (measured: 20 named candidates instead of 25 on the
reference dataset). Cluster-level results are unaffected.

---

## Method notes

Two design points that are easy to get wrong and that the pipeline enforces
structurally rather than by convention:

**Tissue identity is checked before sort validation.** `Omp` marks mature
sensory neurons in *both* the main olfactory epithelium and the VNO, so an MOE
library can pass as a successful VNO sort on marker evidence alone. The
tissue-identity tier compares an MOE panel against a VNO-specific panel with an
absolute CPM floor, and a wrong-tissue library is stopped before any receptor
statement can be produced.

**An even within-cluster split is not sufficient evidence of an EM artifact.**
Paralogs at near-equal fractions usually mean one gene was expressed and EM split
the reads — but in a multi-cell pool two neurons expressing two paralogs of the
same cluster is expected biology, since monogenic choice is a per-*cell* rule.
The discriminator is whether both paralogs carry independent unique-read support.
Flagging on fractions alone would report real co-expression as a sequencing
artifact.

Cluster aggregation follows Dietschi et al. 2022, *Sci Adv* 8(46) eabn7450
([doi:10.1126/sciadv.abn7450](https://doi.org/10.1126/sciadv.abn7450)).

---

## Versioning

`kernel.py` carries `VNO_SKILL_VERSION`, readable as `vno_skill_version()`.
Because installs are per-user copies, that string is the only way two people can
tell whether they are running the same conventions. Bump it on every change.
