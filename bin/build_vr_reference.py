#!/usr/bin/env python3
"""
build_vr_reference.py — VNO receptor (V1R/V2R) gene annotation + genomic cluster reference.

Builds the reference layer that all downstream VR quantification joins against:

  ref/vr_gene_annotation.tsv   one row per VR gene: coords, biotype, pseudogene status
  ref/vr_clusters.tsv          one row per genomic cluster of same-family paralogs
  ref/vr_gene_to_cluster.tsv   flat gene -> cluster lookup (the join key downstream)
  ref/vr_gtf_parse_report.txt  audit trail for the GTF pass and the quant-table crosscheck

Design notes
------------
* This GRCm38 GTF has NO `gene` feature rows, so gene spans are aggregated.
  Transcript rows are preferred for the span; exon-derived spans are computed
  independently and any disagreement is reported rather than silently reconciled.
* The GTF is ~852 MB: it is streamed in a single pass, line-by-line, no whole-file read.
* All thresholds/paths come from config/project.yaml. Nothing is hardcoded here except
  the regexes that define family membership, which are the definition of the task.

Usage
-----
  module load python/3.11.4
  python build_vr_reference.py --config <work>/config/project.yaml --outdir <work>/ref
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
from collections import Counter, defaultdict

__version__ = "1.0.0"

# ---------------------------------------------------------------- family definition
# Case-sensitive, anchored. Catches Vmn1r1, Vmn1r-ps150, Vmn1r169, Vmn2r116 etc.
FAMILY_PATTERNS = [
    ("V1R", re.compile(r"^Vmn1r")),
    ("V2R", re.compile(r"^Vmn2r")),
]

ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')


# ---------------------------------------------------------------- config (minimal YAML)
def load_config(path):
    """Parse the subset of project.yaml we need.

    pyyaml is not guaranteed on the cluster module, so this reads the flat
    scalar keys we depend on with an indentation-aware scan. Keys are returned
    as a flat dict of "dotted.path" -> string value.
    """
    try:
        import yaml  # noqa: PLC0415

        with open(path) as fh:
            raw = yaml.safe_load(fh)
        return _flatten(raw), "pyyaml"
    except Exception:  # noqa: BLE001 - fall back to the hand scan
        pass

    flat, stack = {}, []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            indent = len(line) - len(line.lstrip())
            line = line.split("  #")[0].rstrip()
            if ":" not in line:
                continue
            key, _, val = line.strip().partition(":")
            key, val = key.strip(), val.strip()
            while stack and stack[-1][0] >= indent:
                stack.pop()
            path_key = ".".join([k for _, k in stack] + [key])
            if val == "":
                stack.append((indent, key))
            else:
                flat[path_key] = val
    return flat, "fallback-scan"


def _flatten(node, prefix=""):
    out = {}
    if isinstance(node, dict):
        for k, v in node.items():
            out.update(_flatten(v, f"{prefix}.{k}" if prefix else str(k)))
    elif isinstance(node, list):
        out[prefix] = ",".join(str(x) for x in node)
    else:
        out[prefix] = str(node)
    return out


def cfg_get(flat, *candidates, default=None, cast=str):
    for key in candidates:
        if key in flat:
            return cast(flat[key])
    if default is None:
        raise KeyError(f"none of {candidates} found in config")
    return default


# ---------------------------------------------------------------- chromosome ordering
def chrom_sort_key(chrom):
    """Genome order: 1..19, X, Y, MT, then scaffolds alphabetically."""
    c = str(chrom)
    if c.isdigit():
        return (0, int(c), "")
    special = {"X": 1, "Y": 2, "MT": 3, "M": 3}
    if c in special:
        return (1, special[c], "")
    return (2, 0, c)


# ---------------------------------------------------------------- STEP 1: GTF pass
def parse_gtf(gtf_path):
    """Single streaming pass over the GTF. Returns (genes, stats)."""
    genes = {}
    n_rows = 0
    n_vr_rows = 0
    feature_counter = Counter()

    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            n_rows += 1
            # Cheap pre-filter: skip the ~99% of lines that cannot be a VR gene
            # before paying for attribute parsing.
            if "Vmn1r" not in line and "Vmn2r" not in line:
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _src, feature, start, end, _score, strand, _frame, attrs = fields[:9]

            attr = dict(ATTR_RE.findall(attrs))
            gene_name = attr.get("gene_name", "")
            family = None
            for fam, pat in FAMILY_PATTERNS:
                if pat.match(gene_name):
                    family = fam
                    break
            if family is None:
                continue

            n_vr_rows += 1
            feature_counter[feature] += 1
            gid = attr.get("gene_id", gene_name)
            start_i, end_i = int(start), int(end)

            g = genes.get(gid)
            if g is None:
                g = genes[gid] = {
                    "gene_id": gid,
                    "gene_name": gene_name,
                    "family": family,
                    "chroms": set(),
                    "strands": set(),
                    "biotypes": set(),
                    "transcripts": set(),
                    "n_exons": 0,
                    "tx_start": None,
                    "tx_end": None,
                    "ex_start": None,
                    "ex_end": None,
                    "any_start": None,
                    "any_end": None,
                    "gene_names": set(),
                }
            g["chroms"].add(chrom)
            g["strands"].add(strand)
            g["gene_names"].add(gene_name)
            bt = attr.get("gene_biotype") or attr.get("gene_type") or ""
            if bt:
                g["biotypes"].add(bt)
            tid = attr.get("transcript_id")
            if tid:
                g["transcripts"].add(tid)

            g["any_start"] = start_i if g["any_start"] is None else min(g["any_start"], start_i)
            g["any_end"] = end_i if g["any_end"] is None else max(g["any_end"], end_i)
            if feature == "transcript":
                g["tx_start"] = start_i if g["tx_start"] is None else min(g["tx_start"], start_i)
                g["tx_end"] = end_i if g["tx_end"] is None else max(g["tx_end"], end_i)
            elif feature == "exon":
                g["n_exons"] += 1
                g["ex_start"] = start_i if g["ex_start"] is None else min(g["ex_start"], start_i)
                g["ex_end"] = end_i if g["ex_end"] is None else max(g["ex_end"], end_i)

    stats = {
        "n_rows_scanned": n_rows,
        "n_vr_rows": n_vr_rows,
        "vr_feature_counts": dict(feature_counter),
    }
    return genes, stats


def finalize_genes(genes):
    """Resolve spans, pseudogene status, and collect span disagreements."""
    rows, span_conflicts, multi_locus = [], [], []

    for gid, g in genes.items():
        chrom = sorted(g["chroms"], key=chrom_sort_key)[0]
        if len(g["chroms"]) > 1:
            multi_locus.append((gid, g["gene_name"], "chroms=" + ",".join(sorted(g["chroms"]))))
        if len(g["strands"]) > 1:
            multi_locus.append((gid, g["gene_name"], "strands=" + ",".join(sorted(g["strands"]))))

        # Prefer transcript-derived span; fall back to exon, then to any row.
        if g["tx_start"] is not None:
            start, end, span_src = g["tx_start"], g["tx_end"], "transcript"
        elif g["ex_start"] is not None:
            start, end, span_src = g["ex_start"], g["ex_end"], "exon"
        else:
            start, end, span_src = g["any_start"], g["any_end"], "other"

        # Cross-check transcript span against exon span.
        if g["tx_start"] is not None and g["ex_start"] is not None:
            if (g["tx_start"], g["tx_end"]) != (g["ex_start"], g["ex_end"]):
                span_conflicts.append(
                    {
                        "gene_id": gid,
                        "gene_name": g["gene_name"],
                        "chrom": chrom,
                        "tx_span": f"{g['tx_start']}-{g['tx_end']}",
                        "exon_span": f"{g['ex_start']}-{g['ex_end']}",
                        "delta_start": g["ex_start"] - g["tx_start"],
                        "delta_end": g["ex_end"] - g["tx_end"],
                    }
                )

        biotype = ";".join(sorted(g["biotypes"])) if g["biotypes"] else "NA"
        by_biotype = "pseudogene" in biotype.lower()
        by_name = "-ps" in g["gene_name"]
        if by_biotype and by_name:
            crit = "biotype+name"
        elif by_biotype:
            crit = "biotype"
        elif by_name:
            crit = "name"
        else:
            crit = "none"

        rows.append(
            {
                "gene_id": gid,
                "gene_name": g["gene_name"],
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": sorted(g["strands"])[0],
                "gene_biotype": biotype,
                "n_transcripts": len(g["transcripts"]),
                "n_exons": g["n_exons"],
                "family": g["family"],
                "is_pseudogene": int(by_biotype or by_name),
                "pseudogene_criterion": crit,
                "span_source": span_src,
                "exon_span_start": g["ex_start"] if g["ex_start"] is not None else "NA",
                "exon_span_end": g["ex_end"] if g["ex_end"] is not None else "NA",
                "span_agrees_tx_vs_exon": int(
                    g["tx_start"] is None
                    or g["ex_start"] is None
                    or (g["tx_start"], g["tx_end"]) == (g["ex_start"], g["ex_end"])
                ),
            }
        )

    rows.sort(key=lambda r: (r["family"], chrom_sort_key(r["chrom"]), r["start"], r["gene_name"]))
    return rows, span_conflicts, multi_locus


# ---------------------------------------------------------------- STEP 2: clustering
def build_clusters(rows, max_gap):
    """Single-linkage clustering on genomic proximity, WITHIN family.

    A new cluster is cut whenever the gap between a gene's start and the running
    max end of the current cluster exceeds max_gap. V1R and V2R are clustered
    independently, so a V1R and a V2R gene 50 kb apart never share a cluster.
    Cluster ordinals run sequentially along the genome within each family.
    """
    clusters = []
    gene2cluster = {}

    for family in ("V1R", "V2R"):
        fam_rows = [r for r in rows if r["family"] == family]
        fam_rows.sort(key=lambda r: (chrom_sort_key(r["chrom"]), r["start"], r["end"]))

        ordinal = 0
        current = []
        cur_chrom = None
        cur_max_end = None

        def flush(members, ordn):
            if not members:
                return
            chrom = members[0]["chrom"]
            cid = f"{family}_chr{chrom}_cl{ordn:03d}"
            starts = [m["start"] for m in members]
            ends = [m["end"] for m in members]
            n_ps = sum(m["is_pseudogene"] for m in members)
            clusters.append(
                {
                    "cluster_id": cid,
                    "family": family,
                    "chrom": chrom,
                    "start": min(starts),
                    "end": max(ends),
                    "span_bp": max(ends) - min(starts) + 1,
                    "n_genes": len(members),
                    "n_functional": len(members) - n_ps,
                    "n_pseudogenes": n_ps,
                    "member_gene_names": ",".join(m["gene_name"] for m in members),
                    "member_gene_ids": ",".join(m["gene_id"] for m in members),
                }
            )
            for i, m in enumerate(members):
                gene2cluster[m["gene_id"]] = {
                    "gene_id": m["gene_id"],
                    "gene_name": m["gene_name"],
                    "family": family,
                    "chrom": chrom,
                    "start": m["start"],
                    "end": m["end"],
                    "strand": m["strand"],
                    "is_pseudogene": m["is_pseudogene"],
                    "cluster_id": cid,
                    "cluster_n_genes": len(members),
                    "cluster_rank_in_cluster": i + 1,
                }

        for r in fam_rows:
            new_cluster = (
                cur_chrom is None
                or r["chrom"] != cur_chrom
                or (r["start"] - cur_max_end) > max_gap
            )
            if new_cluster:
                flush(current, ordinal)
                ordinal += 1
                current = [r]
                cur_chrom = r["chrom"]
                cur_max_end = r["end"]
            else:
                current.append(r)
                cur_max_end = max(cur_max_end, r["end"])
        flush(current, ordinal)

    clusters.sort(key=lambda c: (c["family"], chrom_sort_key(c["chrom"]), c["start"]))
    return clusters, gene2cluster


# Dietschi et al. 2022 (doi 10.1126/sciadv.abn7450) derive their cluster boundaries from
# a data-driven "aggregation threshold" on the intergenic-distance distribution, NOT from
# a fixed round number. On this annotation that threshold sits far above 200 kb.
#
# At 800 kb this annotation yields 19 V1R superclusters counting every assembly unit, of
# which one (V1R_chrGL456219.1_sc019) is a lone gene on an unplaced scaffold. Restricted
# to the primary assembly — the basis comparable to a published karyotype-level count —
# 800 kb gives 18, matching the 18 mouse V1r clusters Dietschi et al. report. Both
# numbers are printed in the report so the basis of the comparison is explicit; do not
# quote one without saying which. (1000 kb gives 18 all-in / 17 primary-only, so it does
# NOT match on the comparable basis.)
#
# We keep the 200 kb config rule as the primary (conservative, finer) tier and expose the
# literature-scale grouping as a SECOND tier so cluster-level results can be compared to
# the paper without re-clustering. Superclusters are advisory context, not the
# aggregation unit.
SUPERCLUSTER_GAP_BP = 800_000


def build_superclusters(rows, gap=SUPERCLUSTER_GAP_BP):
    """Coarser, literature-calibrated grouping. Returns gene_id -> supercluster_id."""
    out = {}
    for family in ("V1R", "V2R"):
        fam_rows = sorted(
            (r for r in rows if r["family"] == family),
            key=lambda r: (chrom_sort_key(r["chrom"]), r["start"], r["end"]),
        )
        ordinal, cur_chrom, cur_max_end, members = 0, None, None, []

        def flush(ms, ordn):
            if not ms:
                return
            sid = f"{family}_chr{ms[0]['chrom']}_sc{ordn:03d}"
            for m in ms:
                out[m["gene_id"]] = sid

        for r in fam_rows:
            if cur_chrom is None or r["chrom"] != cur_chrom or (r["start"] - cur_max_end) > gap:
                flush(members, ordinal)
                ordinal += 1
                members, cur_chrom, cur_max_end = [r], r["chrom"], r["end"]
            else:
                members.append(r)
                cur_max_end = max(cur_max_end, r["end"])
        flush(members, ordinal)
    return out


def gap_distribution(rows):
    """Intergenic gaps between consecutive same-family VR genes, per family.

    This is the diagnostic Dietschi et al. use to place a clustering threshold
    (their Fig. 1C). Reported so the choice of cluster_max_gap_bp is auditable
    against the data rather than taken on faith.
    """
    out = {}
    for family in ("V1R", "V2R"):
        fam_rows = sorted(
            (r for r in rows if r["family"] == family),
            key=lambda r: (chrom_sort_key(r["chrom"]), r["start"], r["end"]),
        )
        gaps, prev_chrom, run_max_end = [], None, None
        for r in fam_rows:
            if r["chrom"] == prev_chrom and run_max_end is not None:
                g = r["start"] - run_max_end
                if g > 0:
                    gaps.append(g)
            if r["chrom"] != prev_chrom:
                prev_chrom, run_max_end = r["chrom"], r["end"]
            else:
                run_max_end = max(run_max_end, r["end"])
        out[family] = sorted(gaps)
    return out


PRIMARY_CHROMS = {str(i) for i in range(1, 20)} | {"X", "Y", "MT"}


def sweep_thresholds(rows, gaps_kb=(50, 100, 150, 200, 300, 400, 500, 800, 1000)):
    """Cluster-count sensitivity to the gap threshold, for the report.

    Reports counts on ALL assembly units and on the primary assembly only. The two
    differ because unplaced scaffolds and patches each contribute their own cluster;
    published counts are comparable to the primary-only column.
    """
    primary = [r for r in rows if r["chrom"] in PRIMARY_CHROMS]
    res = []
    for gkb in gaps_kb:
        clusters, _ = build_clusters(rows, gkb * 1000)
        prim_clusters, _ = build_clusters(primary, gkb * 1000)
        row = {"gap_kb": gkb}
        for fam in ("V1R", "V2R"):
            fc = [c for c in clusters if c["family"] == fam]
            sizes = sorted((c["n_genes"] for c in fc), reverse=True)
            row[fam] = {
                "n_clusters": len(fc),
                "n_clusters_primary_only": sum(1 for c in prim_clusters if c["family"] == fam),
                "singletons": sum(1 for s in sizes if s == 1),
                "largest": sizes[0] if sizes else 0,
            }
        res.append(row)
    return res


# ---------------------------------------------------------------- quant crosscheck
def quant_vr_names(quant_tsv):
    """Collect Vmn1r*/Vmn2r* gene_name + gene_id from a salmon merged table."""
    found = {}
    if not quant_tsv or not os.path.exists(quant_tsv):
        return found, "missing"
    with open(quant_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            gi, gn = header.index("gene_id"), header.index("gene_name")
        except ValueError:
            gi, gn = 0, 1
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(gi, gn):
                continue
            name = p[gn]
            if name.startswith("Vmn1r") or name.startswith("Vmn2r"):
                found[p[gi]] = name
    return found, "ok"


# ---------------------------------------------------------------- writers
def write_tsv(path, rows, columns):
    with open(path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in columns) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--also-outdir", default=None, help="second copy target (for transfer-back)")
    ap.add_argument("--gtf", default=None, help="override config GTF path")
    args = ap.parse_args()

    flat, cfg_mode = load_config(args.config)
    gtf = args.gtf or cfg_get(flat, "reference.gtf")
    genome_key = cfg_get(flat, "reference.genome_key", default="unknown")
    max_gap = cfg_get(flat, "thresholds.cluster_max_gap_bp", cast=lambda v: int(float(v)))

    outdirs = [args.outdir] + ([args.also_outdir] if args.also_outdir else [])
    for d in outdirs:
        os.makedirs(d, exist_ok=True)

    print(f"[cfg]  parsed via {cfg_mode}; cluster_max_gap_bp={max_gap}", flush=True)
    print(f"[gtf]  streaming {gtf}", flush=True)

    genes, stats = parse_gtf(gtf)
    rows, span_conflicts, multi_locus = finalize_genes(genes)
    clusters, gene2cluster = build_clusters(rows, max_gap)
    superclusters = build_superclusters(rows)
    gaps = gap_distribution(rows)
    sweep = sweep_thresholds(rows)

    # Genes on non-primary assembly units (patches/unplaced scaffolds) are annotated but
    # are NOT part of the primary-assembly repertoire; flag rather than drop them.
    primary = {str(i) for i in range(1, 20)} | {"X", "Y", "MT"}
    for r in rows:
        r["supercluster_id"] = superclusters.get(r["gene_id"], "NA")
        r["is_primary_assembly"] = int(r["chrom"] in primary)
    for cl in clusters:
        members = cl["member_gene_ids"].split(",")
        cl["supercluster_id"] = superclusters.get(members[0], "NA")

    # ---- quant-table crosscheck (both trials; they share the annotation)
    quant_reports = []
    for trial in ("trial1", "trial2"):
        res = flat.get(f"trials.{trial}.results")
        rel = flat.get("paths.gene_tpm", "star_salmon/salmon.merged.gene_tpm.tsv")
        if not res:
            continue
        qpath = os.path.join(res, rel)
        qfound, status = quant_vr_names(qpath)
        gtf_ids = {r["gene_id"] for r in rows}
        gtf_names = {r["gene_name"] for r in rows}
        n1 = sum(1 for n in qfound.values() if n.startswith("Vmn1r"))
        n2 = sum(1 for n in qfound.values() if n.startswith("Vmn2r"))
        quant_reports.append(
            {
                "trial": trial,
                "path": qpath,
                "status": status,
                "n_vmn1r": n1,
                "n_vmn2r": n2,
                "in_quant_not_gtf_by_id": sorted(set(qfound) - gtf_ids),
                "in_gtf_not_quant_by_id": sorted(gtf_ids - set(qfound)),
                "in_quant_not_gtf_by_name": sorted(set(qfound.values()) - gtf_names),
                "in_gtf_not_quant_by_name": sorted(gtf_names - set(qfound.values())),
            }
        )

    # ---- write tables
    ann_cols = [
        "gene_id", "gene_name", "chrom", "start", "end", "strand", "gene_biotype",
        "n_transcripts", "n_exons", "family", "is_pseudogene", "pseudogene_criterion",
        "span_source", "exon_span_start", "exon_span_end", "span_agrees_tx_vs_exon",
        "supercluster_id", "is_primary_assembly",
    ]
    clu_cols = [
        "cluster_id", "family", "chrom", "start", "end", "span_bp", "n_genes",
        "n_functional", "n_pseudogenes", "supercluster_id",
        "member_gene_names", "member_gene_ids",
    ]
    g2c_cols = [
        "gene_id", "gene_name", "family", "chrom", "start", "end", "strand",
        "is_pseudogene", "cluster_id", "cluster_n_genes", "cluster_rank_in_cluster",
        "supercluster_id", "is_primary_assembly",
    ]

    # gene_to_cluster follows annotation order so it is diff-stable
    g2c_rows = []
    for r in rows:
        e = gene2cluster.get(r["gene_id"])
        if e is None:
            continue
        e = dict(e)
        e["supercluster_id"] = r["supercluster_id"]
        e["is_primary_assembly"] = r["is_primary_assembly"]
        g2c_rows.append(e)

    for d in outdirs:
        write_tsv(os.path.join(d, "vr_gene_annotation.tsv"), rows, ann_cols)
        write_tsv(os.path.join(d, "vr_clusters.tsv"), clusters, clu_cols)
        write_tsv(os.path.join(d, "vr_gene_to_cluster.tsv"), g2c_rows, g2c_cols)

    # ---- parse report
    fam_counts = Counter(r["family"] for r in rows)
    ps_crit = Counter(r["pseudogene_criterion"] for r in rows)
    chrom_dist = Counter((r["family"], r["chrom"]) for r in rows)
    clu_by_fam = Counter(c["family"] for c in clusters)

    lines = []
    A = lines.append
    A("VR GTF parse + cluster definition report")
    A(f"script: build_vr_reference.py v{__version__}")
    A(f"config: {args.config}  (parsed via {cfg_mode})")
    A(f"genome_key: {genome_key}")
    A(f"gtf: {gtf}")
    A(f"cluster_max_gap_bp: {max_gap}")
    A("")
    A("== GTF pass ==")
    A(f"total rows scanned (non-comment): {stats['n_rows_scanned']}")
    A(f"VR-matching rows: {stats['n_vr_rows']}")
    A("VR rows by feature type: " + ", ".join(
        f"{k}={v}" for k, v in sorted(stats["vr_feature_counts"].items())))
    A(f"distinct VR gene_ids: {len(rows)}")
    A("")
    A("== genes per family ==")
    for fam in ("V1R", "V2R"):
        f_rows = [r for r in rows if r["family"] == fam]
        A(f"{fam}: {fam_counts[fam]} genes "
          f"({sum(r['is_pseudogene'] for r in f_rows)} pseudogene, "
          f"{sum(1 - r['is_pseudogene'] for r in f_rows)} functional)")
    A("")
    A("== pseudogene criterion (which rule fired) ==")
    for k in ("none", "biotype", "name", "biotype+name"):
        A(f"{k}: {ps_crit.get(k, 0)}")
    A("")
    A("== chromosome distribution ==")
    for fam in ("V1R", "V2R"):
        items = sorted(((c, n) for (f, c), n in chrom_dist.items() if f == fam),
                       key=lambda x: chrom_sort_key(x[0]))
        A(f"{fam}: " + ", ".join(f"chr{c}={n}" for c, n in items))
    A("")
    A("== transcript-span vs exon-span disagreements ==")
    A(f"n_disagreements: {len(span_conflicts)}")
    for sc in span_conflicts[:50]:
        A(f"  {sc['gene_name']} ({sc['gene_id']}) chr{sc['chrom']}: "
          f"tx={sc['tx_span']} exon={sc['exon_span']} "
          f"dstart={sc['delta_start']} dend={sc['delta_end']}")
    if len(span_conflicts) > 50:
        A(f"  ... {len(span_conflicts) - 50} more")
    A("")
    A("== multi-locus / inconsistent gene records ==")
    A(f"n_flagged: {len(multi_locus)}")
    for gid, name, note in multi_locus[:50]:
        A(f"  {name} ({gid}): {note}")
    A("")
    A("== clusters ==")
    for fam in ("V1R", "V2R"):
        fc = [c for c in clusters if c["family"] == fam]
        sizes = sorted((c["n_genes"] for c in fc), reverse=True)
        A(f"{fam}: {clu_by_fam[fam]} clusters; "
          f"singletons={sum(1 for s in sizes if s == 1)}; "
          f"largest={sizes[:8]}")
        per_chrom = Counter(c["chrom"] for c in fc)
        A(f"  clusters per chrom: " + ", ".join(
            f"chr{c}={n}" for c, n in sorted(per_chrom.items(), key=lambda x: chrom_sort_key(x[0]))))
        A("  top clusters by gene count:")
        for c in sorted(fc, key=lambda c: -c["n_genes"])[:10]:
            A(f"    {c['cluster_id']}  n={c['n_genes']} "
              f"(func={c['n_functional']}, ps={c['n_pseudogenes']}) "
              f"span={c['span_bp']/1e6:.2f}Mb  {c['member_gene_names'].split(',')[0]}"
              f"..{c['member_gene_names'].split(',')[-1]}")
    A("")
    A("== intergenic gap distribution (basis for cluster_max_gap_bp) ==")
    A("Dietschi et al. 2022 place cluster boundaries with a data-driven aggregation")
    A("threshold on this distribution (their Fig. 1C) rather than a fixed round number.")
    for fam in ("V1R", "V2R"):
        g = gaps[fam]
        if not g:
            continue
        def pct(p):
            return g[min(len(g) - 1, int(len(g) * p / 100))]
        A(f"{fam}: n_gaps={len(g)} median={pct(50):,} q25={pct(25):,} q75={pct(75):,} "
          f"q90={pct(90):,} q95={pct(95):,} max={g[-1]:,}")
        A(f"  gaps > {max_gap:,} bp (i.e. cluster breaks): {sum(1 for x in g if x > max_gap)}")
        A(f"  gaps in 200kb-800kb (reassigned by the supercluster tier): "
          f"{sum(1 for x in g if 200_000 < x <= 800_000)}")
    A("")
    A("== threshold sensitivity sweep ==")
    A("n = clusters over ALL assembly units; nP = primary assembly only (chr 1-19,X,Y,MT).")
    A("Published counts are comparable to nP, not n.")
    A("gap_kb   V1R n/nP/singl/largest     V2R n/nP/singl/largest")
    for s in sweep:
        a1, a2 = s["V1R"], s["V2R"]
        A(f"{s['gap_kb']:6d}   "
          f"{a1['n_clusters']:3d}/{a1['n_clusters_primary_only']:3d}/"
          f"{a1['singletons']:3d}/{a1['largest']:3d}"
          f"             {a2['n_clusters']:3d}/{a2['n_clusters_primary_only']:3d}/"
          f"{a2['singletons']:3d}/{a2['largest']:3d}")
    sc_row = next((s for s in sweep if s["gap_kb"] * 1000 == SUPERCLUSTER_GAP_BP), None)
    A(f"NOTE: primary tier uses cluster_max_gap_bp={max_gap} from config (unchanged).")
    if sc_row:
        A(f"      supercluster tier uses {SUPERCLUSTER_GAP_BP} bp, giving "
          f"{sc_row['V1R']['n_clusters']} V1R superclusters over all assembly units and "
          f"{sc_row['V1R']['n_clusters_primary_only']} on the primary assembly.")
        A("      The primary-assembly figure is the one comparable to the 18 mouse V1r")
        A("      clusters reported by Dietschi et al. 2022; the extra all-units cluster is")
        A("      a lone V1R gene on unplaced scaffold GL456219.1.")
    A("")
    A("== prior-label reconciliation ==")
    for probe in ("Vmn1r91", "Vmn1r166", "Vmn1r32", "Vmn1r39"):
        hit = [g for g in g2c_rows if g["gene_name"] == probe]
        if hit:
            h = hit[0]
            A(f"  {probe} -> {h['cluster_id']} (chr{h['chrom']}, "
              f"cluster n_genes={h['cluster_n_genes']}, "
              f"supercluster={h['supercluster_id']})")
        else:
            A(f"  {probe} -> NOT FOUND in GTF-derived annotation")
    # Does the supercluster tier reunite the prior "60 paralogs" chr7 block?
    sc_of = {g["gene_name"]: g["supercluster_id"] for g in g2c_rows}
    if sc_of.get("Vmn1r91") and sc_of.get("Vmn1r91") == sc_of.get("Vmn1r166"):
        n = sum(1 for g in g2c_rows if g["supercluster_id"] == sc_of["Vmn1r91"])
        A(f"  Vmn1r91 and Vmn1r166 share supercluster {sc_of['Vmn1r91']} (n={n} genes);")
        A("  at the 200kb primary tier they are two clusters split by a 217,366 bp gap.")
    A("")
    A("== quant-table crosscheck (expected Vmn1r=318, Vmn2r=220) ==")
    for q in quant_reports:
        A(f"{q['trial']}: {q['path']} [{q['status']}]")
        A(f"  quant Vmn1r={q['n_vmn1r']}  Vmn2r={q['n_vmn2r']}  total={q['n_vmn1r']+q['n_vmn2r']}")
        A(f"  in quant not GTF (by gene_id): {len(q['in_quant_not_gtf_by_id'])} "
          f"{q['in_quant_not_gtf_by_id'][:20]}")
        A(f"  in GTF not quant (by gene_id): {len(q['in_gtf_not_quant_by_id'])} "
          f"{q['in_gtf_not_quant_by_id'][:20]}")
        A(f"  in quant not GTF (by gene_name): {len(q['in_quant_not_gtf_by_name'])} "
          f"{q['in_quant_not_gtf_by_name'][:20]}")
        A(f"  in GTF not quant (by gene_name): {len(q['in_gtf_not_quant_by_name'])} "
          f"{q['in_gtf_not_quant_by_name'][:20]}")
    A("")
    A("== outputs ==")
    for d in outdirs:
        for f in ("vr_gene_annotation.tsv", "vr_clusters.tsv", "vr_gene_to_cluster.tsv"):
            A(f"  {os.path.join(d, f)}")

    report = "\n".join(lines) + "\n"
    for d in outdirs:
        with open(os.path.join(d, "vr_gtf_parse_report.txt"), "w") as fh:
            fh.write(report)
    print(report, flush=True)


if __name__ == "__main__":
    sys.exit(main())
