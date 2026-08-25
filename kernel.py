"""Thin helpers for the VNO receptor RNA-seq pipeline.

The heavy modules live in the pipeline's own bin/ directory
(bin/vr_*.py); nothing here duplicates them. These helpers locate a run,
load its tables, apply the reporting gate, and expose the settled
conventions so a session never re-derives them.
"""

import os
import csv
import glob
import math

# Bump on every change to this skill or the pipeline it documents. Distribution
# is per-user copies (see skill/install_skill.py), so this is the only way a
# session can tell whether the local copy is current -- ask for
# vno_skill_version() when comparing notes with a collaborator.
VNO_SKILL_VERSION = "1.0.0"

# The pipeline is site-independent: it runs from any work folder on any machine
# with pandas/numpy/matplotlib/pyyaml (+ samtools for the unique-read channel).
# Work-folder resolution, first hit wins:
#   1. an explicit work_dir= argument to any helper here
#   2. $VR_WORK
#   3. $PWD, if it contains config/project.yaml
# There is deliberately NO site-specific default: a hardcoded path silently
# pointed every session at one cluster.
VNO_WORK_ENV_VAR = "VR_WORK"

# Optional remote target, read at call time via vno_ssh_target(). Set
# $VR_SSH_TARGET only if the data lives on a cluster; absent that the pipeline
# runs locally and nothing here needs SSH.
VNO_SSH_TARGET_ENV_VAR = "VR_SSH_TARGET"

# Dependencies, by platform. `samtools` is separate from the python stack and is
# required for the unique-read evidence channel -- without it vr_quantify stops
# at its preflight guard, and --no-bam degrades the result (no individual-
# receptor calls possible).
VNO_DEPS_PYTHON = ("pandas", "numpy", "matplotlib", "pyyaml", "scipy")
VNO_DEPS_BINARY = ("samtools",)
VNO_INSTALL_HINTS = {
    "conda": ("conda create -n vr python=3.11 pandas numpy matplotlib pyyaml scipy "
              "-c conda-forge && conda activate vr && "
              "conda install -c bioconda samtools"),
    "venv": ("python3 -m venv ~/vr-env && . ~/vr-env/bin/activate && "
             "pip install pandas numpy matplotlib pyyaml scipy   "
             "# samtools separately via your package manager"),
    "modules": "module load python3 && module load samtools   # HPC with Environment Modules",
}
VNO_CLUSTER_GAP_BP = 200000
VNO_SUPERCLUSTER_GAP_BP = 800000
VNO_PERMUTATION_FLOOR = 0.00024993751562109475
VNO_TIER_NAMES = ("tissue_identity", "sort_validation", "population_id",
                  "cluster_vr", "individual_vr")


def vno_skill_version():
    """Version of this skill copy, for comparing against a collaborator's.

    Distribution is per-user copies, so two people can be on different
    versions with no warning. Compare this string before trusting that two
    runs used the same conventions.
    """
    return VNO_SKILL_VERSION


def vno_ssh_target():
    """Remote target if one is configured, else None (local run)."""
    return os.environ.get(VNO_SSH_TARGET_ENV_VAR) or None


def vno_default_work_dir():
    """Resolve the work folder without assuming any particular site.

    Order: $VR_WORK, then $PWD if it holds config/project.yaml, then $PWD.
    Raises nothing -- callers that need a valid folder should follow up with
    vno_validate_work_dir().
    """
    env = os.environ.get(VNO_WORK_ENV_VAR)
    if env:
        return os.path.abspath(os.path.expanduser(env))
    cwd = os.getcwd()
    if os.path.isfile(os.path.join(cwd, "config", "project.yaml")):
        return cwd
    return cwd


def vno_validate_work_dir(work_dir=None):
    """Check a work folder is usable; return a dict of findings.

    Reports what is present rather than raising, so a session can tell the user
    exactly what to fix (missing config, missing bin/, missing ref tables).
    """
    wd = work_dir or vno_default_work_dir()
    p = vno_paths(wd)
    out = {"work_dir": wd, "ok": True, "missing": [], "notes": []}
    if not os.path.isfile(p["config"]):
        out["missing"].append(p["config"])
        out["notes"].append(
            "no config/project.yaml -- create the work folder with "
            "`python3 bin/vr_init.py --dest <dir> --trial <name> "
            "--results <star_salmon dir> --gtf <file>`")
    for key in ("bin", "results"):
        if not os.path.isdir(p[key]):
            out["missing"].append(p[key])
    ref = p["ref_tables"]["gene_annotation"]
    if not os.path.isfile(ref):
        out["notes"].append(
            "ref/vr_gene_annotation.tsv absent -- the `ref` stage will build it "
            "from the GTF on first run (one pass, a few minutes)")
    out["ok"] = not out["missing"]
    return out


def vno_install_check():
    """Report which dependencies are importable/on PATH, with a remedy.

    Platform-agnostic: no module system is assumed.
    """
    import importlib
    import shutil
    missing_py = []
    for mod in VNO_DEPS_PYTHON:
        name = "yaml" if mod == "pyyaml" else mod
        try:
            importlib.import_module(name)
        except ImportError:
            missing_py.append(mod)
    missing_bin = [b for b in VNO_DEPS_BINARY if shutil.which(b) is None]
    out = {
        "python_ok": not missing_py,
        "missing_python": missing_py,
        "binaries_ok": not missing_bin,
        "missing_binaries": missing_bin,
        "install_hints": VNO_INSTALL_HINTS if (missing_py or missing_bin) else {},
    }
    if missing_bin and not missing_py:
        out["degraded_mode"] = (
            "runnable with --no-bam, but the unique-read channel is lost: no "
            "individual-receptor calls, and EM-redistribution flags fall back "
            "to fraction-only evidence")
    return out


def vno_paths(work_dir=None):
    """Resolved layout of the work folder plus the pipeline module list."""
    if work_dir is None:
        work_dir = vno_default_work_dir()
    j = os.path.join
    return {
        "work": work_dir,
        "bin": j(work_dir, "bin"),
        "ref": j(work_dir, "ref"),
        "config": j(work_dir, "config", "project.yaml"),
        "results": j(work_dir, "results"),
        "figures": j(work_dir, "results", "figures"),
        "docs": j(work_dir, "docs"),
        "logs": j(work_dir, "logs"),
        "ref_tables": {
            "gene_annotation": j(work_dir, "ref", "vr_gene_annotation.tsv"),
            "clusters": j(work_dir, "ref", "vr_clusters.tsv"),
            "gene_to_cluster": j(work_dir, "ref", "vr_gene_to_cluster.tsv"),
            "parse_report": j(work_dir, "ref", "vr_gtf_parse_report.txt"),
        },
        "modules": ["build_vr_reference.py", "vr_config.py", "vr_qc.py",
                    "vr_markers.py", "vr_sample_qc.py", "vr_quantify.py",
                    "vr_report.py", "vr_figures.py", "vr_qc_figures.py"],
        "ssh_target": VNO_SSH_TARGET,
        "install_hints": VNO_INSTALL_HINTS,
    }


def vno_validate_run_dir(results_root):
    """Check an nf-core/rnaseq results dir against the input contract.

    `results_root` is the trial results directory (the parent of star_salmon/).
    Returns {ok, star_salmon, present, missing, samples, notes}. Run this
    before spending a cluster job on a new trial.
    """
    ss = os.path.join(results_root, "star_salmon")
    required = ["salmon.merged.gene_counts.tsv", "salmon.merged.gene_tpm.tsv",
                "tx2gene.tsv"]
    optional = ["salmon.merged.gene_counts_scaled.tsv"]
    present, missing, notes = [], [], []
    if not os.path.isdir(ss):
        return {"ok": False, "star_salmon": ss, "present": [],
                "missing": ["star_salmon/ (directory absent)"], "samples": [],
                "notes": ["not an nf-core star_salmon results directory"]}
    for f in required:
        (present if os.path.exists(os.path.join(ss, f)) else missing).append(f)
    for f in optional:
        if os.path.exists(os.path.join(ss, f)):
            present.append(f)
    samples = sorted(
        os.path.basename(p)[: -len(".Log.final.out")]
        for p in glob.glob(os.path.join(ss, "log", "*.Log.final.out")))
    if not samples:
        missing.append("log/<sample>.Log.final.out")
    for s in samples:
        if not os.path.exists(os.path.join(ss, s, "quant.sf")):
            notes.append("no quant.sf for " + s)
        if not os.path.isdir(os.path.join(ss, "qualimap", s)):
            notes.append("no qualimap dir for " + s)
        if not glob.glob(os.path.join(ss, s + "*.bam")):
            notes.append("no BAM for " + s
                         + " -- run vr_quantify.py --no-bam; no individual"
                           " receptor call is possible without MAPQ255 evidence")
    return {"ok": not missing, "star_salmon": ss, "present": present,
            "missing": missing, "samples": samples, "notes": notes}


def vno_read_table(path):
    """Read a pipeline TSV. Returns (DataFrame, provenance_comment).

    Pipeline outputs carry a leading '#' provenance line recording the CPM
    convention and source table; pandas would otherwise take it as the header.
    """
    import pandas as pd
    prov = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                prov.append(line.lstrip("#").strip())
            else:
                break
    df = pd.read_csv(path, sep="\t", comment="#")
    return df, " | ".join(prov)


def vno_cpm(counts, gene_id_col=None, gene_name_col=None):
    """The settled CPM convention: count / all-gene column sum * 1e6.

    `counts` is the UNSCALED salmon.merged.gene_counts.tsv as a DataFrame.
    Do not substitute the scaled or TPM table: they run 2.3-2.5x lower and
    every threshold in project.yaml is calibrated against this convention.
    """
    if gene_id_col is None:
        gene_id_col = "gene_id"
    if gene_name_col is None:
        gene_name_col = "gene_name"
    drop = [c for c in (gene_id_col, gene_name_col) if c in counts.columns]
    numeric = counts.drop(columns=drop)
    totals = numeric.sum(axis=0)
    cpm = numeric.divide(totals, axis=1) * 1e6
    for c in reversed(drop):
        cpm.insert(0, c, counts[c].values)
    return cpm


def vno_load_results(work_dir=None, trial=None):
    """Load every pipeline output table for a trial into one dict.

    Keys: sample_qc, technical_qc, marker_cpm, cluster_expression,
    within_cluster_fractions, artifact_flags, candidates, tier_status,
    tier_outcomes, plus 'provenance' and 'missing'.
    """
    if work_dir is None:
        work_dir = vno_default_work_dir()
    if trial is None:
        trial = "trial2"
    p = vno_paths(work_dir)
    tdir = os.path.join(p["results"], trial)
    wanted = {
        "sample_qc": os.path.join(tdir, "sample_qc.tsv"),
        "technical_qc": os.path.join(tdir, "technical_qc.tsv"),
        "marker_cpm": os.path.join(tdir, "marker_cpm.tsv"),
        "cluster_expression": os.path.join(tdir, "vr_cluster_expression.tsv"),
        "within_cluster_fractions": os.path.join(
            tdir, "vr_within_cluster_fractions.tsv"),
        "artifact_flags": os.path.join(tdir, "vr_artifact_flags.tsv"),
        "candidates": os.path.join(tdir, "vr_candidates.tsv"),
        "tier_status": os.path.join(tdir, "tier_status.tsv"),
        "tier_outcomes": os.path.join(tdir, "tier_outcomes.tsv"),
    }
    out = {"trial": trial, "dir": tdir, "provenance": {}, "missing": []}
    for key, path in wanted.items():
        if os.path.exists(path):
            df, prov = vno_read_table(path)
            out[key] = df
            if prov:
                out["provenance"][key] = prov
        else:
            out[key] = None
            out["missing"].append(os.path.basename(path))
    return out


def vno_clearance(sample_qc):
    """Apply the reporting clearance gate to a sample_qc DataFrame.

    Biology may be reported only where qc_overall == 'USABLE' AND
    suppress_biology is false. suppress_biology alone is NOT sufficient:
    libraries exist with suppress_biology=False and sort_verdict=FAIL.
    Returns a copy with a boolean `clearance` column added.
    """
    def _explicitly_false(v):
        # fail closed: a missing/blank suppress_biology is NOT clearance
        return str(v).strip().lower() in ("false", "0", "0.0", "no")
    df = sample_qc.copy()
    usable = df["qc_overall"].astype(str).str.strip().str.upper() == "USABLE"
    not_supp = df["suppress_biology"].map(_explicitly_false)
    df["clearance"] = usable & not_supp
    return df


def vno_em_verdict(row):
    """Read an artifact-flags row into an EM verdict. Never re-judges evenness.

    Consumes em_flag and the unique-read gate columns as the pipeline wrote
    them. A consumer that recomputes evenness from the fractions will
    eventually disagree with the pipeline -- that bug already happened once
    in this project's figures.
    """
    def _g(key, default=""):
        try:
            v = row[key]
        except (KeyError, IndexError, TypeError):
            return default
        if v is None:
            return default
        if isinstance(v, float) and math.isnan(v):
            return default
        return v

    flag = str(_g("em_flag")).strip()
    level = str(_g("em_flag_level")).strip()
    block = str(_g("even_block_size")).strip()
    supported = str(_g("n_block_members_unique_supported")).strip()
    try:
        n_supported = int(float(supported)) if supported else None
    except ValueError:
        n_supported = None
    try:
        block_n = int(float(block)) if block else 0
    except ValueError:
        block_n = 0

    if not flag:
        verdict = "not_evaluated"
        reading = ("em_flag is empty: this row was not tested (typically a"
                   " QC-suppressed sample or a non-cluster scope row). No EM"
                   " statement may be made from it.")
    elif flag == "suspected_em_redistribution":
        verdict = "em_redistribution"
        reading = ("one paralog cleared the unique-read gate inside the even"
                   " block: EM split the reads of a single expressed paralog")
    elif block_n >= 2 and n_supported is not None and n_supported >= 2:
        verdict = "co_expression"
        reading = ("two or more block members carry independent unique-read"
                   " support: real co-expression across the pool, not an"
                   " artifact (monogenic choice is a per-CELL rule)")
    elif flag == "no_redistribution_signature":
        verdict = "no_artifact"
        reading = "no even block detected; per-paralog structure stands as measured"
    elif flag == "single_paralog_only":
        verdict = "single_paralog"
        reading = "only one paralog detected in this cluster"
    elif flag in ("no_signal", "insufficient_signal"):
        verdict = "no_evidence"
        reading = "not enough signal in this cluster to test anything"
    else:
        verdict = "unknown_flag"
        reading = ("em_flag value %r is outside the pipeline vocabulary --"
                   " do not interpret it" % flag)

    p_unif = str(_g("even_block_p_uniform") or _g("p_uniform")).strip()
    p_note = ""
    if p_unif:
        try:
            if abs(float(p_unif) - VNO_PERMUTATION_FLOOR) < 1e-12:
                p_note = ("p is at the Monte Carlo floor (1/4001): quote as"
                          " 'p < 2.5e-4', a bound, not a point estimate")
        except ValueError:
            pass

    return {
        "verdict": verdict,
        "reading": reading,
        "em_flag": flag,
        "em_flag_level": level,
        "even_block_size": block_n,
        "n_block_members_unique_supported": n_supported,
        "block_unique_reads": str(_g("block_unique_reads")),
        "unique_support_threshold": _g("unique_support_threshold"),
        "unique_background_floor": _g("unique_background_floor"),
        "unique_channel_used": str(_g("unique_channel_used")),
        "p_uniform_note": p_note,
        "pseudogene_bleed_flag": str(_g("pseudogene_bleed_flag")),
        "confirmation_status": "tentative_unconfirmed",
    }


def vno_flag_vocabulary():
    """Exact allowed values per status/flag column. Do not invent others."""
    return {
        "tier_name": list(VNO_TIER_NAMES),
        "tier_status": ["PASS", "PASS_WITH_CAVEAT", "FAIL", "NO_DATA",
                        "SUPPRESSED"],
        "tissue_verdict": ["VNO", "VNO_dominant_mixed", "MOE",
                           "MOE_dominant_mixed", "ambiguous_mixed",
                           "no_tissue_signal"],
        "sort_verdict": ["PASS", "FAIL", "FAIL_WRONG_TISSUE"],
        "library_status": ["OK", "FAILED", "DEGENERATE"],
        "tech_verdict": ["PASS", "WARN", "FAIL"],
        "qc_overall": ["USABLE", "UNUSABLE"],
        "population_call": ["V1R_dominant", "V2R_dominant", "undetermined"],
        "em_flag": ["no_signal", "insufficient_signal",
                    "no_redistribution_signature", "single_paralog_only",
                    "suspected_em_redistribution"],
        "em_flag_level": ["none", "strong"],
        "pseudogene_bleed_flag": ["no_pseudogene_members",
                                  "no_apparent_pseudogene_expression",
                                  "apparent_pseudogene_expression",
                                  "pseudogene_only_cluster"],
        "evidence_type": ["bam_unique_mapq255",
                          "no_cluster_above_signal_threshold"],
        "confidence": ["moderate", "alternative_candidate", "unresolvable",
                       "no_call"],
        "confirmation_status": ["tentative_unconfirmed"],
        "interpretation_context": ["monogenic_expectation", "pooled_ambiguous"],
    }


def vno_check_yaml_duplicate_keys(path):
    """Detect duplicate top-level YAML keys (silent last-wins shadowing).

    project.yaml was broken once by an append that created a second
    top-level `markers:`/`thresholds:` block. YAML resolves last-wins with
    no error, so the older values silently vanished. Run this after every
    config edit. Returns {ok, duplicates, keys}.
    """
    seen = {}
    with open(path) as fh:
        for i, line in enumerate(fh, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            if line[0] in " \t-":
                continue
            if ":" not in line:
                continue
            key = line.split(":", 1)[0].strip()
            seen.setdefault(key, []).append(i)
    dups = {k: v for k, v in seen.items() if len(v) > 1}
    return {"ok": not dups, "duplicates": dups, "keys": sorted(seen)}


def vno_run_plan(work_dir=None, trial=None, rebuild_reference=False):
    """Ordered shell commands to run the pipeline on a trial."""
    if work_dir is None:
        work_dir = vno_default_work_dir()
    if trial is None:
        trial = "trial2"
    # No module load: bin/run_pipeline.sh detects conda / venv / Environment
    # Modules itself and fails with a platform-appropriate remedy.
    cmds = ["cd " + work_dir]
    if rebuild_reference:
        cmds.append("python bin/build_vr_reference.py --config "
                    "config/project.yaml --outdir ref")
    cmds += [
        "python bin/vr_sample_qc.py --all-trials",
        "python bin/vr_quantify.py --trial %s --threads 8" % trial,
        "python bin/vr_report.py --selftest",
        "python bin/vr_report.py",
        "python bin/vr_figures.py",
        "python bin/vr_qc_figures.py --all-trials --fig-dir results/figures",
    ]
    return {
        "trial": trial,
        "commands": cmds,
        "shell": "\n".join(cmds),
        "check_first": [
            "results/tier_status_all.tsv -> tier0_status is the verdict on"
            " tissue, NOT tissue_verdict alone (a GFP- nontarget legitimately"
            " reads no_tissue_signal and still passes tier 0)",
            "results/%s/marker_cpm.tsv -> tissue_verdict; MOE or"
            " MOE_dominant_mixed = wrong tissue, remedy is wet-lab" % trial,
            "results/%s/marker_cpm.tsv -> library_status not FAILED/DEGENERATE"
            % trial,
            "results/tier_status_all.tsv -> highest_tier_reported per sample",
            "results/%s/vr_cluster_expression.tsv where is_called==1" % trial,
            "results/%s/vr_artifact_flags.tsv -> em_flag BEFORE any gene name"
            % trial,
            "results/%s/vr_candidates.tsv -> all tentative_unconfirmed" % trial,
        ],
        "reminder": ("config edits go INSIDE existing blocks; run"
                     " vno_check_yaml_duplicate_keys() afterwards"),
    }
