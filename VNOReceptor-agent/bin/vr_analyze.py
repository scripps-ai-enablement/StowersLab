#!/usr/bin/env python3
"""vr_analyze.py -- one command, VNO receptor RNA-seq end to end.

    python3 vr_analyze.py --results /path/to/nfcore/results --gtf /path/to/genes.gtf

That is the whole interface. It scaffolds a work folder, discovers samples,
writes a validated config, runs every stage, and prints the answer plus a
machine-readable JSON summary. No config editing, no stage sequencing, no
environment setup beyond having the packages importable.

Designed to be driven by an agent or a human from a terminal, on a laptop or an
HPC login/compute node. Everything it decides is printed; everything it cannot
decide safely, it refuses to guess.

Exit codes
  0  ran, and at least one library reached a receptor-cluster call
  3  ran, but no library survived QC (a real result, not an error)
  4  sample roles could not be inferred -- rerun with --target/--nontarget
  5  missing dependency
  6  bad input (not an nf-core star_salmon tree, GTF unreadable)

Add --json-only for pure JSON on stdout (progress still goes to stderr).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
PKG = HERE.parent if HERE.name == "bin" else HERE

PY_DEPS = ("pandas", "numpy", "matplotlib", "yaml")
PY_PIP_NAMES = {"yaml": "pyyaml"}

# Sample-role heuristics. Deliberately conservative: a name that matches neither
# list is UNKNOWN, and UNKNOWN stops the run rather than defaulting. Getting
# target/nontarget backwards inverts the entire sort-validation tier.
# No \b anchors around these stems: real sample names run the label straight
# into a digit ("target100cells", "pool2cellsRep3"), and \btarget\b does not
# match there because a letter-digit junction is not a word boundary. Verified
# against the naming conventions in this project's two trials.
NONTARGET_PAT = re.compile(r"non[-_]?target|gfp[-_]?neg|negative|\bctrl\b|control", re.I)
TARGET_PAT = re.compile(r"target|gfp[-_]?pos|positive|pool|sorted", re.I)


def log(msg, *, quiet=False):
    if not quiet:
        sys.stderr.write(f"[vr_analyze] {msg}\n")
        sys.stderr.flush()


def die(code, msg, hint=None):
    sys.stderr.write(f"\n[vr_analyze] ERROR: {msg}\n")
    if hint:
        sys.stderr.write(hint.rstrip() + "\n")
    sys.stderr.write("\n")
    raise SystemExit(code)


# --------------------------------------------------------------- dependencies
def check_deps(allow_no_bam):
    import importlib
    missing = []
    for mod in PY_DEPS:
        try:
            importlib.import_module(mod)
        except ImportError:
            missing.append(PY_PIP_NAMES.get(mod, mod))
    if missing:
        die(5, f"missing python packages: {', '.join(missing)}",
            "  conda install -c conda-forge " + " ".join(missing) + "\n"
            "  or: pip install " + " ".join(missing) + "\n"
            "  or on an HPC: module load python3\n"
            "  or point at another interpreter: VR_PYTHON=/path/to/python")
    have_samtools = shutil.which("samtools") is not None
    if not have_samtools and not allow_no_bam:
        die(5, "samtools is not on PATH",
            "  samtools powers the unique-read evidence channel, which is what\n"
            "  separates an EM artifact from real co-expression.\n\n"
            "    conda install -c bioconda samtools\n"
            "    module load samtools        # HPC\n"
            "    apt install samtools / brew install samtools\n\n"
            "  Or rerun with --no-bam to proceed WITHOUT it. That degrades the\n"
            "  result: no individual-receptor calls are possible and\n"
            "  EM-redistribution flags fall back to fraction-only evidence.\n"
            "  Cluster-level results are unaffected.")
    return {"samtools": have_samtools}


# --------------------------------------------------------------------- inputs
def resolve_run_dir(p):
    """Accept the star_salmon dir, its parent, or a tree containing one."""
    p = pathlib.Path(os.path.abspath(os.path.expanduser(p)))
    if not p.is_dir():
        die(6, f"not a directory: {p}")
    if p.name == "star_salmon":
        return p, p.parent
    if (p / "star_salmon").is_dir():
        return p / "star_salmon", p
    hits = sorted(glob.glob(str(p / "**" / "star_salmon"), recursive=True))
    if len(hits) == 1:
        ss = pathlib.Path(hits[0])
        return ss, ss.parent
    if len(hits) > 1:
        die(6, f"found {len(hits)} star_salmon directories under {p}",
            "  Point --results at exactly one:\n    " + "\n    ".join(hits[:8]))
    die(6, f"no star_salmon directory at or under {p}",
        "  --results should be an nf-core/rnaseq output directory (the one\n"
        "  containing star_salmon/salmon.merged.gene_counts.tsv).")


def discover_samples(star_salmon):
    for name in ("salmon.merged.gene_counts.tsv", "salmon.merged.gene_tpm.tsv"):
        t = star_salmon / name
        if t.is_file():
            with open(t) as fh:
                header = fh.readline().rstrip("\n").split("\t")
            names = header[2:]
            if not names:
                die(6, f"{t} has no sample columns")
            return names, t
    die(6, f"no salmon.merged.gene_counts.tsv under {star_salmon}",
        "  Is this an nf-core/rnaseq >=3.18 star_salmon output directory?")


def assign_roles(samples, target_pats, nontarget_pats):
    """Return {sample: role} plus the list still unresolved."""
    roles, unknown = {}, []
    for s in samples:
        forced = None
        for pat in nontarget_pats:
            if re.search(pat, s, re.I):
                forced = "nontarget"
        for pat in target_pats:
            if re.search(pat, s, re.I):
                forced = "target"
        if forced:
            roles[s] = forced
            continue
        nt, tg = NONTARGET_PAT.search(s), TARGET_PAT.search(s)
        if nt and not tg:
            roles[s] = "nontarget"
        elif tg and not nt:
            roles[s] = "target"
        elif nt and tg:
            # "nontarget" contains "target" -- the negative reading wins.
            roles[s] = "nontarget"
        else:
            roles[s] = "UNKNOWN"
            unknown.append(s)
    return roles, unknown


def guess_n_cells(sample):
    m = re.search(r"(\d+)\s*cells?", sample, re.I)
    return m.group(1) if m else "null"


# --------------------------------------------------------------------- config
def carried_blocks():
    """The analysis contract, copied verbatim from the shipped template."""
    for cand in (PKG / "config" / "project.template.yaml",
                 HERE / "project.template.yaml",
                 PKG / "project.template.yaml"):
        if cand.is_file():
            text = cand.read_text()
            break
    else:
        die(6, "config/project.template.yaml not found beside the pipeline",
            "  It carries the threshold, marker and path blocks. Without it the\n"
            "  config would need hand-typed thresholds, which silently changes\n"
            "  verdicts. Re-unpack the pipeline package.")
    keep = ("paths", "markers", "thresholds", "qc_thresholds")
    out = []
    for part in re.split(r"(?m)^(?=[a-z_]+:)", text):
        m = re.match(r"([a-z_]+):", part)
        if m and m.group(1) in keep:
            out.append(part.rstrip() + "\n")
    got = {re.match(r"([a-z_]+):", b).group(1) for b in out}
    if set(keep) - got:
        die(6, f"template is missing blocks: {sorted(set(keep)-got)}")
    return "\n".join(out)


HEADER = """# Generated by vr_analyze.py -- do not hand-edit thresholds here.
# The analysis blocks below are copied verbatim from
# config/project.template.yaml. Edit the template to change defaults.

project: {project}
work: {work}

env:
  activate: ""

reference:
  genome_key: {genome_key}
  gtf: {gtf}
  fasta: {fasta}

trials:
  {trial}:
    results: {results_root}
    fastq: ""
    platform: {platform}

samples:
  {trial}:
{sample_lines}

"""


def write_config(work, trial, results_root, gtf, fasta, genome_key, platform,
                 roles, project):
    lines = []
    for s, role in roles.items():
        lines.append(f"    {s}: {{cell_type: {role}, n_cells: {guess_n_cells(s)}, "
                     f"prep_status: ok}}")
    cfg_dir = work / "config"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    path = cfg_dir / "project.yaml"
    path.write_text(HEADER.format(
        project=project, work=work, genome_key=genome_key,
        gtf=gtf, fasta=fasta or '""', trial=trial,
        results_root=results_root, platform=platform,
        sample_lines="\n".join(lines),
    ) + carried_blocks())
    return path


def validate_config(work, path, trial):
    """Load through the pipeline's own config module -- the real contract."""
    import yaml
    text = path.read_text()
    keys = re.findall(r"(?m)^([a-z_]+):", text)
    dupes = sorted({k for k in keys if keys.count(k) > 1})
    if dupes:
        die(6, f"generated config has duplicate top-level keys: {dupes}",
            "  YAML resolves duplicates last-wins with no error, silently\n"
            "  shadowing the earlier block. This is a bug in vr_analyze.py.")
    sys.path.insert(0, str(work / "bin"))
    import vr_config
    cfg = vr_config.load_config(str(path))
    P = vr_config.trial_paths(cfg, trial)
    vr_config.samples_of(cfg, trial)
    for key in ("gene_counts", "star_logs"):
        if not os.path.exists(P[key]):
            die(6, f"config resolves {key} to a path that does not exist:\n    {P[key]}",
                "  --results should point at the nf-core run ROOT (the directory\n"
                "  containing star_salmon/ and multiqc/), not deeper.")
    return cfg, P


# ------------------------------------------------------------------ execution
def stage_env(work):
    env = dict(os.environ)
    env["PYTHONPATH"] = str(work / "bin") + os.pathsep + env.get("PYTHONPATH", "")
    env["MPLBACKEND"] = "Agg"
    env["PYTHONHASHSEED"] = "0"
    env["VR_WORK"] = str(work)
    return env


def run_stage(name, argv, work, logf, quiet):
    t0 = time.time()
    log(f"{name} ...", quiet=quiet)
    with open(logf, "a") as fh:
        fh.write(f"\n=== {name}: {' '.join(argv)}\n")
        fh.flush()
        proc = subprocess.run(argv, cwd=str(work), env=stage_env(work),
                              stdout=fh, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    if proc.returncode != 0:
        tail = "".join(open(logf).readlines()[-25:])
        die(1, f"stage {name!r} failed (exit {proc.returncode})",
            f"  full log: {logf}\n\n{tail}")
    log(f"{name} done ({dt:.0f}s)", quiet=quiet)
    return dt


# -------------------------------------------------------------------- summary
def summarize(work, trial):
    import pandas as pd

    def rd(rel):
        p = work / "results" / rel
        if not p.is_file():
            return None
        with open(p) as fh:
            first = fh.readline()
        return pd.read_csv(p, sep="\t", skiprows=1 if first.startswith("#") else 0)

    ts = rd(f"{trial}/tier_status.tsv")
    ce = rd(f"{trial}/vr_cluster_expression.tsv")
    cd = rd(f"{trial}/vr_candidates.tsv")
    fl = rd(f"{trial}/vr_artifact_flags.tsv")

    out = {"trial": trial, "work_dir": str(work), "libraries": [],
           "clusters_called": 0, "candidates": 0, "flags": {}}
    if ts is not None:
        for _, r in ts.iterrows():
            out["libraries"].append({
                "sample": r.get("sample"),
                "cell_type": r.get("cell_type"),
                "tissue": r.get("tier0_status"),
                "qc": r.get("qc_overall"),
                "highest_tier": int(r["highest_tier_reported"])
                if pd.notna(r.get("highest_tier_reported")) else None,
                "highest_tier_name": r.get("highest_tier_name"),
                "stopped_because": (str(r.get("stop_reason"))[:200]
                                    if pd.notna(r.get("stop_reason")) else None),
            })
    # A summary MUST apply the same gate the report applies. `is_called` is set
    # by quantification without regard to QC, so a suppressed library can carry
    # called clusters -- surfacing those as results is exactly the failure mode
    # the tier ladder exists to prevent. Gate on suppress_biology, and report
    # the suppressed count separately so nothing is silently dropped.
    if ce is not None and "is_called" in ce:
        at_cluster = ce[(ce.is_called == 1) & (ce.tier == "cluster")]
        if "suppress_biology" in at_cluster:
            supp = at_cluster.suppress_biology.astype(str).str.lower().isin(
                ("true", "1", "yes"))
            called, blocked = at_cluster[~supp], at_cluster[supp]
        else:
            called, blocked = at_cluster, at_cluster.iloc[0:0]
        out["clusters_called"] = int(len(called))
        out["clusters_suppressed"] = int(len(blocked))
        if len(blocked):
            out["suppressed_note"] = (
                f"{len(blocked)} cluster row(s) in QC-failed libraries are NOT "
                "reported: they carry VR signal but no biology may be read from "
                "them. See tier_status.tsv stop_reason.")
        out["called_clusters"] = [
            {"sample": r["sample"], "cluster_id": r.cluster_id,
             "family": r.family, "cpm": round(float(r.cpm_sum), 2),
             "share_of_vr": round(float(r.share_of_sample_vr), 3)}
            for _, r in called.sort_values("cpm_sum", ascending=False).iterrows()
        ]
    if cd is not None:
        named = cd[cd.gene_name.notna()] if "gene_name" in cd else cd.iloc[0:0]
        # `candidates` counts NAMED candidate rows; `candidate_rows_total`
        # includes the no_call rows the pipeline emits for clusters it scored
        # but could not name. Reporting only the first would understate the
        # table; reporting only the second would overstate the answer.
        out["candidates"] = int(len(named))
        out["candidate_rows_total"] = int(len(cd))
        out["all_tentative"] = bool(
            (cd.confirmation_status == "tentative_unconfirmed").all()
        ) if "confirmation_status" in cd else None
        out["top_candidates"] = [
            {"sample": r["sample"], "cluster_id": r.cluster_id,
             "gene_name": r.gene_name, "confidence": r.confidence,
             "unique_reads": (int(r.unique_reads_bam_nodup)
                              if pd.notna(r.get("unique_reads_bam_nodup")) else None),
             "confirmation": r.confirmation_status}
            for _, r in named[named.get("rank", 1) == 1].iterrows()
        ] if "rank" in named else []
    if fl is not None and "em_flag" in fl:
        out["flags"] = {k: int(v) for k, v in fl.em_flag.value_counts().items()}
        if "pseudogene_bleed_flag" in fl:
            out["pseudogene_flags"] = {
                k: int(v) for k, v in fl.pseudogene_bleed_flag.value_counts().items()
            }
    return out


def print_report(s, quiet=False):
    if quiet:
        return
    w = sys.stderr.write
    w("\n" + "=" * 72 + "\n")
    w(f"  RESULT -- {s['trial']}\n")
    w("=" * 72 + "\n\n")
    w(f"  {len(s['libraries'])} libraries:\n")
    for L in s["libraries"]:
        tier = L["highest_tier"]
        ok = tier is not None and tier >= 3
        mark = "OK " if ok else "-- "
        # `highest_tier_reported` is the last tier that produced content; the
        # tier that FAILED is one further along. Printing only the former reads
        # as a contradiction next to a stop reason naming a later tier, so show
        # what was reached and let the stop reason name the gate that closed.
        reached = (f"reached tier {tier} ({L['highest_tier_name']})"
                   if tier is not None and tier >= 0
                   else "nothing reportable")
        w(f"    {mark}{L['sample']:<26} {str(L['cell_type']):<10} {reached}\n")
        if L["stopped_because"] and not ok:
            w(f"        {L['stopped_because'][:150]}\n")
    w(f"\n  Receptor clusters called: {s['clusters_called']}\n")
    if s.get("clusters_suppressed"):
        w(f"    ({s['clusters_suppressed']} more in QC-failed libraries, "
          "NOT reported -- no biology may be read from them)\n")
    for c in (s.get("called_clusters") or [])[:12]:
        w(f"    {c['sample']:<26} {c['cluster_id']:<18} "
          f"{c['cpm']:>10.1f} CPM  ({c['share_of_vr']:.0%} of VR)\n")
    w(f"\n  Individual candidates: {s['candidates']}")
    if s.get("all_tentative"):
        w("  -- ALL tentative_unconfirmed\n")
    else:
        w("\n")
    for c in (s.get("top_candidates") or [])[:12]:
        w(f"    {c['sample']:<26} {c['cluster_id']:<18} {str(c['gene_name']):<12} "
          f"{c['confidence']:<22} {c['unique_reads']} unique reads\n")
    if s.get("flags"):
        w("\n  EM-artifact flags: "
          + ", ".join(f"{k}={v}" for k, v in sorted(s["flags"].items())) + "\n")
    w(f"\n  Full outputs: {s['work_dir']}/results/\n")
    w(f"  Read first:   results/{s['trial']}/tier_status.tsv, then vr_report.md\n")
    w("  How to interpret every file: docs/VR_OUTPUT_GUIDE.md\n\n")


def selftest():
    """Assert the summary layer cannot report a QC-suppressed library.

    The pipeline enforces the tier ladder internally, but this script builds its
    own summary from the tables and could bypass it -- which it did, on the first
    version, surfacing called clusters from a failed library. These synthetic
    rows lock that shut.
    """
    import pandas as pd
    fails = []

    ce = pd.DataFrame([
        # a clean library with a real call
        {"sample": "good", "tier": "cluster", "cluster_id": "V1R_x_cl001",
         "family": "V1R", "is_called": 1, "cpm_sum": 900.0,
         "share_of_sample_vr": 0.9, "suppress_biology": False},
        # a suppressed library carrying VR signal -- must NOT be reported
        {"sample": "failed_lib", "tier": "cluster", "cluster_id": "V1R_x_cl010",
         "family": "V1R", "is_called": 1, "cpm_sum": 501.0,
         "share_of_sample_vr": 1.0, "suppress_biology": True},
        # supercluster duplicate of the clean row -- must not double-count
        {"sample": "good", "tier": "supercluster", "cluster_id": "V1R_x_sc001",
         "family": "V1R", "is_called": 1, "cpm_sum": 900.0,
         "share_of_sample_vr": 0.9, "suppress_biology": False},
    ])
    at_cluster = ce[(ce.is_called == 1) & (ce.tier == "cluster")]
    supp = at_cluster.suppress_biology.astype(str).str.lower().isin(
        ("true", "1", "yes"))
    called, blocked = at_cluster[~supp], at_cluster[supp]

    if len(called) != 1:
        fails.append(f"expected 1 reportable cluster, got {len(called)}")
    if "failed_lib" in set(called["sample"]):
        fails.append("FAIL: a suppressed library reached the reported set")
    if len(blocked) != 1:
        fails.append(f"expected 1 suppressed cluster, got {len(blocked)}")
    if "V1R_x_sc001" in set(called.cluster_id):
        fails.append("supercluster row leaked into the cluster-tier count")

    # role inference: the naming conventions that must resolve, and the ones
    # that must NOT be guessed.
    cases = [
        ("nontarget100cells_S8", "nontarget"), ("target100cellsRep1_S6", "target"),
        ("pool2cellsRep3_S7", "target"), ("GFPpos_rep1", "target"),
        ("GFPneg_ctrl", "nontarget"), ("input_control", "nontarget"),
        ("sorted_pos_1", "target"), ("sample_A", "UNKNOWN"), ("lib7", "UNKNOWN"),
    ]
    got, _ = assign_roles([c[0] for c in cases], [], [])
    for name, want in cases:
        if got[name] != want:
            fails.append(f"role({name}) = {got[name]!r}, expected {want!r}")

    for f in fails:
        sys.stderr.write(f"  [FAIL] {f}\n")
    n = len(fails)
    sys.stderr.write(
        f"\n  selftest: {4 + len(cases) - n} passed, {n} failed\n")
    return 1 if n else 0


# ----------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="vr_analyze.py",
        description="VNO receptor RNA-seq, end to end, in one command.")
    ap.add_argument("--results", required=True,
                    help="nf-core/rnaseq output directory (contains star_salmon/)")
    ap.add_argument("--gtf", required=True, help="reference GTF used for that run")
    ap.add_argument("--fasta", default="", help="reference FASTA (optional)")
    ap.add_argument("--out", default=None,
                    help="work folder to create (default: ./vr_out_<trial>)")
    ap.add_argument("--trial", default="run1", help="label for this dataset")
    ap.add_argument("--genome-key", default="unspecified")
    ap.add_argument("--platform", default="Illumina")
    ap.add_argument("--project", default="VNO receptor RNA-seq")
    ap.add_argument("--target", action="append", default=[], metavar="REGEX",
                    help="sample-name pattern to force role=target (repeatable)")
    ap.add_argument("--nontarget", action="append", default=[], metavar="REGEX",
                    help="sample-name pattern to force role=nontarget (repeatable)")
    ap.add_argument("--threads", type=int, default=0,
                    help="0 = auto (min(8, cpu_count))")
    ap.add_argument("--no-bam", action="store_true",
                    help="skip the unique-read channel; DEGRADES the result")
    ap.add_argument("--force-ref", action="store_true",
                    help="rebuild the VR reference from the GTF even if present")
    ap.add_argument("--json-only", action="store_true",
                    help="only JSON on stdout; progress to stderr")
    ap.add_argument("--dry-run", action="store_true",
                    help="scaffold and validate, then stop before the stages")
    ap.add_argument("--selftest", action="store_true",
                    help="assert the summary gate and role inference, then exit")
    if argv is None:
        argv = sys.argv[1:]
    if "--selftest" in argv:
        raise SystemExit(selftest())
    args = ap.parse_args(argv)

    quiet = args.json_only
    t_start = time.time()

    caps = check_deps(args.no_bam)
    use_bam = caps["samtools"] and not args.no_bam
    if not use_bam:
        log("running WITHOUT the unique-read channel: no individual-receptor "
            "calls will be possible", quiet=quiet)

    star_salmon, results_root = resolve_run_dir(args.results)
    gtf = os.path.abspath(os.path.expanduser(args.gtf))
    if not os.path.isfile(gtf):
        die(6, f"GTF not readable: {gtf}")
    log(f"input: {results_root}", quiet=quiet)

    samples, quant = discover_samples(star_salmon)
    roles, unknown = assign_roles(samples, args.target, args.nontarget)
    log(f"discovered {len(samples)} samples in {quant.name}", quiet=quiet)
    for s in samples:
        log(f"    {s:<30} -> {roles[s]}", quiet=quiet)
    if unknown:
        die(4, f"cannot infer target/nontarget for: {', '.join(unknown)}",
            "  target vs nontarget cannot be guessed from these names, and\n"
            "  getting it backwards inverts the entire sort-validation tier, so\n"
            "  this stops rather than defaulting.\n\n"
            "  Say which is which with regexes, e.g.:\n"
            f"    --target '{unknown[0]}' \\\n"
            "    --nontarget 'ctrl|input'\n\n"
            "  target    = the sorted population you want receptor calls from\n"
            "  nontarget = the negative-control population")

    work = pathlib.Path(os.path.abspath(os.path.expanduser(
        args.out or f"./vr_out_{args.trial}")))
    for sub in ("bin", "config", "ref", "results", "docs", "logs"):
        (work / sub).mkdir(parents=True, exist_ok=True)
    # Self-contained: copy the pipeline in so the work folder stands alone.
    for src in sorted((PKG / "bin").glob("*.py")) + sorted((PKG / "bin").glob("*.sh")):
        if ".bak" not in src.name:
            shutil.copy2(src, work / "bin" / src.name)
    for src in sorted((PKG / "config").glob("*.template.yaml")):
        shutil.copy2(src, work / "config" / src.name)
    if (PKG / "ref").is_dir() and not args.force_ref:
        for src in (PKG / "ref").iterdir():
            if src.is_file():
                shutil.copy2(src, work / "ref" / src.name)
    for src in (PKG / "docs").glob("*.md") if (PKG / "docs").is_dir() else []:
        shutil.copy2(src, work / "docs" / src.name)
    log(f"work folder: {work}", quiet=quiet)

    cfg_path = write_config(work, args.trial, results_root, gtf, args.fasta,
                            args.genome_key, args.platform, roles, args.project)
    validate_config(work, cfg_path, args.trial)
    log(f"config written and validated: {cfg_path}", quiet=quiet)

    if args.dry_run:
        log("--dry-run: stopping before the stages", quiet=quiet)
        print(json.dumps({"work_dir": str(work), "config": str(cfg_path),
                          "samples": roles, "dry_run": True}, indent=1))
        return 0

    threads = args.threads or min(8, os.cpu_count() or 4)
    py = os.environ.get("VR_PYTHON") or sys.executable
    B = str(work / "bin")
    C = str(cfg_path)
    logf = work / "logs" / f"vr_analyze_{time.strftime('%Y%m%d_%H%M%S')}.log"

    need_ref = args.force_ref or not (work / "ref" / "vr_gene_annotation.tsv").is_file()
    if need_ref:
        run_stage("reference (streaming GTF parse, a few minutes)",
                  [py, f"{B}/build_vr_reference.py", "--config", C,
                   "--outdir", str(work / "ref")], work, logf, quiet)
    else:
        log("reference tables present, skipping the GTF parse "
            "(--force-ref to rebuild)", quiet=quiet)

    run_stage("qc and gating", [py, f"{B}/vr_sample_qc.py", "--config", C,
                               "--all-trials"], work, logf, quiet)
    qargs = [py, f"{B}/vr_quantify.py", "--config", C, "--trial", args.trial,
             "--threads", str(threads)]
    if not use_bam:
        qargs.append("--no-bam")
    run_stage("quantification", qargs, work, logf, quiet)
    run_stage("gate self-test", [py, f"{B}/vr_report.py", "--config", C,
                                 "--selftest"], work, logf, quiet)
    run_stage("report", [py, f"{B}/vr_report.py", "--config", C,
                         "--trial", args.trial], work, logf, quiet)
    for fig, extra in ((f"{B}/vr_figures.py", []),
                       (f"{B}/vr_qc_figures.py",
                        ["--all-trials", "--fig-dir", str(work / "results" / "figures")])):
        run_stage(f"figures ({pathlib.Path(fig).stem})",
                  [py, fig, "--config", C] + extra, work, logf, quiet)

    s = summarize(work, args.trial)
    s["elapsed_s"] = round(time.time() - t_start, 1)
    s["unique_read_channel"] = use_bam
    s["log"] = str(logf)
    print_report(s, quiet=quiet)
    print(json.dumps(s, indent=1, default=str))

    usable = [L for L in s["libraries"]
              if L["highest_tier"] is not None and L["highest_tier"] >= 3]
    if not usable:
        log("no library reached a cluster-level call -- see stop reasons above. "
            "This is a result, not a crash.", quiet=quiet)
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
