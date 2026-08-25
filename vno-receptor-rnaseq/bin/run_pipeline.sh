#!/usr/bin/env bash
#
# run_pipeline.sh -- end-to-end driver for the VNO vomeronasal-receptor
# RNA-seq pipeline (Stowers Lab / Natalie Cole; CCBB).
#
# Runs the whole chain in dependency order for one trial, or for every trial
# in config/project.yaml. The ONLY inputs are config/project.yaml and the
# nf-core/rnaseq result trees it points at; every table, report and figure
# under results/ is regenerated. Nothing is read from a previous run, so a
# cold start (empty results/) and a re-run over existing outputs must produce
# byte-identical tables.
#
# Stage order and why:
#   ref       build_vr_reference.py   VR gene annotation + cluster tiers from the
#                                     GRCm38 GTF. Trial-independent; skipped when
#                                     ref/vr_gene_annotation.tsv already exists
#                                     (--force-ref rebuilds, --refcheck rebuilds
#                                     into a scratch dir and diffs, which is the
#                                     cold-start test that does not disturb ref/).
#   qc        vr_sample_qc.py         technical QC + marker CPM + tissue/sort/
#                                     population verdicts. Always run with
#                                     --all-trials: it recomputes every trial
#                                     from source data (STAR logs, qualimap,
#                                     merged counts) and writes the COMBINED
#                                     results/sample_qc_all.tsv that the
#                                     quantification and report stages read.
#                                     Running it per trial would leave that
#                                     combined table holding one trial only, so
#                                     the result would depend on invocation
#                                     order -- exactly the failure mode this
#                                     driver exists to prevent.
#   quant     vr_quantify.py          cluster-level VR expression, within-cluster
#                                     fractions, EM-redistribution flags and
#                                     candidate receptors for ONE trial. Needs
#                                     samtools (BAM unique-read support).
#   report    vr_report.py            tier ladder + per-trial report. Run per
#                                     trial here; the cross-trial report
#                                     (vr_report_all.md) is written by `finalize`.
#   finalize  vr_report.py (all)      cross-trial report + combined tier tables,
#             vr_report.py --selftest tier-gate assertions,
#             vr_figures.py           VR deliverable figures,
#             vr_qc_figures.py        QC deliverable figures.
#
# Usage:
#   bin/run_pipeline.sh --all                       # trial2, trial1, finalize
#   bin/run_pipeline.sh --trial trial2              # one trial, no finalize
#   bin/run_pipeline.sh --trial trial2 --finalize   # one trial then finalize
#   bin/run_pipeline.sh --all --refcheck            # + reference reproducibility
#   bin/run_pipeline.sh --all --stage quant --stage report
#
# Options:
#   --config PATH   config/project.yaml (default: <work>/config/project.yaml)
#   --trial NAME    run one trial (repeatable)
#   --all           run every trial in the config, trial2 first
#   --stage NAME    restrict to these stages (repeatable): ref qc quant report finalize
#   --force-ref     rebuild ref/ even if it exists
#   --refcheck      rebuild the reference into a scratch dir and diff vs ref/
#   --no-bam        skip BAM unique-read support in quant (fast, weaker evidence)
#   --finalize      run the finalize stage after the named trials
#   --threads N     threads for quant (default 8)
#   --dry-run       print the commands without running them
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="$(cd "$HERE/.." && pwd)"
CONFIG="$WORK/config/project.yaml"
TRIALS=()
STAGES=()
FORCE_REF=0
REFCHECK=0
NO_BAM=0
FINALIZE=0
ALL=0
THREADS=8
DRY=0

die() { echo "[run_pipeline] ERROR: $*" >&2; exit 2; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)   CONFIG="$2"; shift 2 ;;
    --trial)    TRIALS+=("$2"); shift 2 ;;
    --all)      ALL=1; shift ;;
    --stage)    STAGES+=("$2"); shift 2 ;;
    --force-ref) FORCE_REF=1; shift ;;
    --refcheck) REFCHECK=1; shift ;;
    --no-bam)   NO_BAM=1; shift ;;
    --finalize) FINALIZE=1; shift ;;
    --threads)  THREADS="$2"; shift 2 ;;
    --dry-run)  DRY=1; shift ;;
    -h|--help)  sed -n '2,60p' "${BASH_SOURCE[0]}"; exit 0 ;;
    *)          die "unknown option: $1" ;;
  esac
done

[[ -f "$CONFIG" ]] || die "config not found: $CONFIG"

# ---- environment -----------------------------------------------------------
# Portable: works on an HPC with Environment Modules, in a conda/venv, or on a
# plain machine with the packages installed. Nothing here is site-specific.
#
# Interpreter selection, first hit wins:
#   1. $VR_PYTHON        -- explicit override
#   2. the module-provided python, IF Environment Modules exists and provides one
#   3. whatever python3 is on PATH (conda, venv, system)
# We only attempt `module load` when `module` actually exists, and a failed load
# is not fatal -- the import check below is the real gate.
if [[ -z "${VR_PYTHON:-}" ]] && { command -v module >/dev/null 2>&1 || [[ -n "${MODULESHOME:-}" ]]; }; then
  set +u
  # shellcheck disable=SC1091
  [[ -n "${MODULESHOME:-}" && -f "$MODULESHOME/init/bash" ]] && . "$MODULESHOME/init/bash"
  # Unversioned names first so this works on clusters with different pinnings.
  module load python/3.11.4 >/dev/null 2>&1 || module load python3 >/dev/null 2>&1 || module load python >/dev/null 2>&1 || true
  module load samtools/1.19 >/dev/null 2>&1 || module load samtools >/dev/null 2>&1 || true
  set -u
fi
PY="${VR_PYTHON:-python3}"
command -v "$PY" >/dev/null 2>&1 || die "python interpreter '$PY' not found (set \$VR_PYTHON to override)"
export PYTHONPATH="$WORK/bin${PYTHONPATH:+:$PYTHONPATH}"
export MPLBACKEND=Agg
export PYTHONHASHSEED=0          # keep any set/dict-order-sensitive output stable

# Real gate: can the interpreter import what the pipeline needs? The remedy
# depends on the platform, so report what is actually missing and offer all
# three install routes rather than assuming a module system.
if ! $PY - <<'PYEOF'
import sys
missing = []
for m in ("pandas", "numpy", "matplotlib", "yaml"):
    try:
        __import__(m)
    except ImportError:
        missing.append("pyyaml" if m == "yaml" else m)
if missing:
    sys.stderr.write("MISSING:" + ",".join(missing) + "\n")
    raise SystemExit(1)
import pandas, numpy
sys.stderr.write(f"[run_pipeline] python {sys.version.split()[0]} pandas {pandas.__version__} numpy {numpy.__version__}\n")
PYEOF
then
  cat >&2 <<'MSG'

  This python cannot import the packages the pipeline needs (see MISSING above).
  Remedy depends on your environment -- pick one:

    conda create -n vr python=3.11 pandas numpy matplotlib pyyaml scipy \
                 -c conda-forge && conda activate vr
    conda install -c bioconda samtools          # for the unique-read channel

    python3 -m venv ~/vr-env && . ~/vr-env/bin/activate
    pip install pandas numpy matplotlib pyyaml scipy

    module load python3                         # HPC with Environment Modules

  Or point the pipeline at a specific interpreter:
    VR_PYTHON=/path/to/python bin/run_pipeline.sh ...

MSG
  die "python environment is missing required packages"
fi

LOGDIR="$WORK/logs"
mkdir -p "$LOGDIR" "$WORK/results"
STAMP="$(date +%Y%m%d_%H%M%S)"
RUNLOG="$LOGDIR/run_pipeline_${STAMP}.log"

log() { echo "[run_pipeline $(date +%H:%M:%S)] $*" | tee -a "$RUNLOG"; }

run() {
  log "+ $*"
  if [[ "$DRY" == "1" ]]; then return 0; fi
  # shellcheck disable=SC2068
  if ! "$@" >>"$RUNLOG" 2>&1; then
    log "FAILED: $*  (see $RUNLOG)"
    tail -n 40 "$RUNLOG" >&2
    exit 1
  fi
}

want_stage() {
  [[ ${#STAGES[@]} -eq 0 ]] && return 0
  local s
  for s in "${STAGES[@]}"; do [[ "$s" == "$1" ]] && return 0; done
  return 1
}

# ---- trial list ------------------------------------------------------------
if [[ "$ALL" == "1" ]]; then
  # trial2 first: it is the trial that carries reportable VNO biology, so a
  # partial run still yields the informative half.
  CFG_TRIALS="$($PY -c 'import sys,vr_config; c=vr_config.load_config(sys.argv[1]); print(" ".join(sorted(vr_config.trials_of(c), key=lambda t:(t!="trial2",t))))' "$CONFIG")"
  read -r -a TRIALS <<<"$CFG_TRIALS"
  FINALIZE=1
fi
[[ ${#TRIALS[@]} -gt 0 ]] || die "give --trial NAME or --all"

# --- preflight: every named trial must exist in the config -------------------
# Catch an unregistered trial HERE, before ref/ and qc run, so the operator
# gets an actionable message instead of a KeyError traceback out of
# vr_quantify several minutes later.
for _t in "${TRIALS[@]}"; do
  if ! "$PY" - "$CONFIG" "$_t" <<'PYEOF'
import sys
import yaml
cfg = yaml.safe_load(open(sys.argv[1]))
trial = sys.argv[2]
trials = (cfg.get("trials") or {})
samples = (cfg.get("samples") or {})
ok = trial in trials and bool(samples.get(trial))
if not ok:
    known = sorted(set(trials) | set(samples))
    sys.stderr.write(
        "\n  trial %r is not registered in %s\n"
        "  known trials: %s\n\n"
        "  To add it, edit config/project.yaml INSIDE the existing blocks:\n"
        "    trials:\n"
        "      %s:\n"
        "        results: /path/to/nfcore/results/<dir>\n"
        "        fastq:   /path/to/fastqs\n"
        "        platform: Illumina\n"
        "    samples:\n"
        "      %s:\n"
        "        <sample_name>: {cell_type: target, n_cells: 2, prep_status: ok}\n\n"
        "  cell_type is 'target' (GFP+) or 'nontarget' (GFP-).\n"
        "  A naive append that creates a SECOND top-level trials:/samples: key\n"
        "  silently shadows the first -- YAML resolves last-wins with no error.\n"
        "  After editing, verify: grep -c '^trials:' config/project.yaml  (must be 1)\n"
        % (trial, sys.argv[1], ", ".join(known) or "(none)", trial, trial)
    )
    raise SystemExit(1)
PYEOF
  then
    die "trial '$_t' is not registered in $CONFIG -- see the message above"
  fi
done

log "work=$WORK config=$CONFIG trials=${TRIALS[*]} threads=$THREADS no_bam=$NO_BAM"
log "log file: $RUNLOG"

# ---- ref -------------------------------------------------------------------
REF_TABLE="$WORK/ref/vr_gene_annotation.tsv"
if want_stage ref; then
  if [[ -f "$REF_TABLE" && "$FORCE_REF" == "0" ]]; then
    log "stage ref: ref/ tables present -> skip (use --force-ref to rebuild)"
  else
    log "stage ref: building VR reference from the GRCm38 GTF (852MB parse)"
    run $PY "$WORK/bin/build_vr_reference.py" --config "$CONFIG" --outdir "$WORK/ref"
  fi
  if [[ "$REFCHECK" == "1" ]]; then
    SCRATCH="$WORK/results/refcheck_${STAMP}"
    mkdir -p "$SCRATCH"
    log "stage ref: cold-start reproducibility check -> $SCRATCH"
    run $PY "$WORK/bin/build_vr_reference.py" --config "$CONFIG" --outdir "$SCRATCH"
    REFDIFF="$WORK/results/refcheck_diff.txt"
    : > "$REFDIFF"
    for f in vr_gene_annotation.tsv vr_clusters.tsv vr_gene_to_cluster.tsv; do
      if diff -q "$WORK/ref/$f" "$SCRATCH/$f" >/dev/null 2>&1; then
        echo "IDENTICAL  $f" >> "$REFDIFF"
      else
        echo "DIFFERS    $f" >> "$REFDIFF"
        diff "$WORK/ref/$f" "$SCRATCH/$f" | head -20 >> "$REFDIFF" || true
      fi
    done
    cat "$REFDIFF" | tee -a "$RUNLOG"
  fi
fi

# ---- qc (all trials, once) -------------------------------------------------
if want_stage qc; then
  log "stage qc: technical QC + markers + verdicts for ALL trials (writes the combined sample_qc_all.tsv)"
  run $PY "$WORK/bin/vr_qc.py"      --config "$CONFIG" --all-trials
  run $PY "$WORK/bin/vr_markers.py" --config "$CONFIG" --all-trials
  run $PY "$WORK/bin/vr_sample_qc.py" --config "$CONFIG" --all-trials
  [[ -f "$WORK/results/sample_qc_all.tsv" ]] || die "qc stage did not write results/sample_qc_all.tsv"
fi

# ---- per-trial quant + report ---------------------------------------------
for T in "${TRIALS[@]}"; do
  if want_stage quant; then
    log "stage quant [$T]: cluster expression, within-cluster fractions, EM flags, candidates"
    QARGS=(--config "$CONFIG" --trial "$T" --threads "$THREADS")
    [[ "$NO_BAM" == "1" ]] && QARGS+=(--no-bam)
    run $PY "$WORK/bin/vr_quantify.py" "${QARGS[@]}"
  fi
  if want_stage report; then
    log "stage report [$T]: tier ladder + per-trial report"
    run $PY "$WORK/bin/vr_report.py" --config "$CONFIG" --trial "$T"
  fi
done

# ---- finalize --------------------------------------------------------------
if [[ "$FINALIZE" == "1" ]] && want_stage finalize; then
  log "stage finalize: tier-gate selftest -> results/tier_gate_selftest.txt"
  if [[ "$DRY" != "1" ]]; then
    if ! $PY "$WORK/bin/vr_report.py" --config "$CONFIG" --selftest \
         > "$WORK/results/tier_gate_selftest.txt" 2>&1; then
      log "FAILED: tier-gate selftest -- a gate assertion did not hold"
      tail -n 30 "$WORK/results/tier_gate_selftest.txt" >&2
      exit 1
    fi
    tail -n 3 "$WORK/results/tier_gate_selftest.txt" | tee -a "$RUNLOG"
  fi
  log "stage finalize: cross-trial report + combined tier tables"
  run $PY "$WORK/bin/vr_report.py" --config "$CONFIG"
  log "stage finalize: VR figures"
  run $PY "$WORK/bin/vr_figures.py" --config "$CONFIG"
  log "stage finalize: QC figures"
  run $PY "$WORK/bin/vr_qc_figures.py" --config "$CONFIG" --all-trials \
        --fig-dir "$WORK/results/figures"
fi

log "done. outputs under $WORK/results"
