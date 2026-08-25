# VNO receptor RNA-seq — project context

Identify which vomeronasal receptor each mouse VNO sensory neuron expresses, from
an nf-core/rnaseq `star_salmon` results directory.

## The one command

```bash
python3 bin/vr_analyze.py --results <nfcore results dir> --gtf <genes.gtf> --trial <name>
```

It scaffolds a work folder, discovers samples, infers `target`/`nontarget` from
their names, writes and validates the config, runs every stage, and prints a
human report on stderr plus JSON on stdout. Outputs land in `./vr_out_<trial>/`.

Do not hand-assemble the six stages unless debugging one; `bin/run_pipeline.sh`
and the individual `bin/vr_*.py` modules exist for that case and are documented
in `PIPELINE.md`.

Before a real run: `python3 bin/vr_analyze.py --selftest` (expect 13 passed, 0
failed). `--dry-run` scaffolds and validates without running the stages — the
cheap way to confirm the inputs were understood.

## Reading the result

`--json-only` gives structured JSON on stdout with progress on stderr. Keys:
`libraries[]` (each with `highest_tier`, `stopped_because`), `clusters_called`,
`clusters_suppressed`, `called_clusters[]`, `candidates`, `candidate_rows_total`,
`top_candidates[]`, `flags{}`, `unique_read_channel`.

`clusters_called` **already excludes QC-suppressed libraries**; the withheld
count is `clusters_suppressed`. Do not recompute a call count from the raw tables
to cross-check it — that bypasses the gate the pipeline exists to enforce.

Full column-by-column interpretation: `docs/VR_OUTPUT_GUIDE.md`. Why each gate
exists: `docs/VR_PIPELINE_FLOWCHART.md`.

## Exit codes — read them, never retry blindly

| code | meaning |
|---|---|
| `0` | ran; at least one library reached a cluster-level call |
| `1` | a pipeline stage failed — a genuine error; the log path is printed |
| `3` | ran; **nothing survived QC** — this is a result, not a crash |
| `4` | sample roles could not be inferred; rerun with `--target`/`--nontarget` regexes |
| `5` | missing dependency (it names what) |
| `6` | bad input — not an nf-core `star_salmon` tree, or the GTF is unreadable |

**Check order matters when diagnosing an exit code.** Dependencies are checked
*before* inputs, so on a machine missing `samtools` a bad `--results` path still
reports `5` (missing dependency), not `6`. Resolve `5` first, then re-run — the
input error surfaces after. Verified: the same bad path returns `6` once the
dependency gate passes.

**Never wrap this in a retry-on-nonzero loop.** `3` will never succeed on a
retry; the correct response is to read the stop reasons and report them.

Exit `4` is deliberate: `target` vs `nontarget` cannot be guessed from every
naming scheme, and reversing them inverts the entire sort-validation tier, so it
refuses rather than defaulting. Ask the user which is which; do not guess.

## Dependencies

Python 3.9+ with pandas, numpy, matplotlib, pyyaml, scipy — plus `samtools` on
PATH. A system default python usually has none of the former. `samtools` is a
separate binary, not pip-installable.

`--no-bam` runs without samtools but **degrades the result**: no
individual-receptor calls are possible and EM-redistribution flags fall back to
fraction-only evidence (20 named candidates instead of 25 on the reference
dataset). Cluster-level results are unaffected. Say so if you use it; do not
reach for it silently to work around a missing binary.

## Interpretation rules — do not restate results in a way that breaks these

1. **Tissue identity is checked before sort validation.** `Omp` marks mature
   sensory neurons in *both* main olfactory epithelium and VNO, so an MOE library
   looks exactly like a failed VNO sort on marker evidence alone. Never present
   `Omp` as a VNO marker. The remedy for wrong tissue is wet-lab, not
   re-quantification — do not recommend a different annotation for it.
2. **A QC-failed library yields zero receptor statements**, regardless of how
   much VR signal it carries. One library in the reference dataset has 501 CPM of
   receptor signal and is correctly reported as unusable.
3. **Never report an individual receptor without its cluster context.** ~250 V1R
   (`Vmn1r*`) and ~120 V2R (`Vmn2r*`) genes sit in clusters of paralogs at 85–95%
   identity; Salmon's EM spreads ambiguous reads across them, so per-gene counts
   inside a cluster are not trustworthy. Cluster-level aggregation is the
   reliable tier.
4. **An even within-cluster split is not sufficient evidence of an artifact.**
   Monogenic choice is a per-*cell* rule: in a multi-cell pool, two paralogs each
   carrying independent unique-read support is real co-expression. The
   discriminator is the unique-read channel, not the fraction.
5. **Every individual-receptor call is `tentative_unconfirmed` by construction.**
   No evidence channel here can confirm one at these read lengths. If you find
   yourself writing "confirmed", something is wrong.

Read the pipeline's own flag columns; never recompute a verdict a module already
wrote, and never invent a flag value outside the vocabulary in
`vno_flag_vocabulary()` / `docs/VR_OUTPUT_GUIDE.md`.

## Configuration

Every threshold lives in `config/project.yaml` (generated) and
`config/project.template.yaml` (the shipped contract). Nothing is hardcoded in
`bin/`.

When editing a config, insert keys **inside** the existing blocks. A naive append
that creates a second top-level `markers:` or `thresholds:` key silently shadows
the first — YAML resolves duplicates last-wins with no error. This has broken the
project once. `vr_analyze.py` and `vr_init.py` both check for it.

## Scope

Shipped `ref/` tables are mouse GRCm38 with Ensembl chromosome naming (`1`, `7`,
`X` — not `chr1`). Another genome build needs `--force-ref` to rebuild from the
user's GTF. Another species needs the marker panels and family patterns in
`config/project.template.yaml` edited first — the logic generalizes, the gene
names do not.
