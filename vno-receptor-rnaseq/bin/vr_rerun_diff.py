#!/usr/bin/env python3
"""
vr_rerun_diff.py -- compare a freshly regenerated results/ tree against a
backed-up one and report every difference.

A re-run that does not reproduce the committed tables is a bug in the
pipeline, not a nuisance, so this compares deliberately strictly and then
classifies what it finds:

  IDENTICAL      byte-identical (after dropping '#' provenance/comment lines,
                 which carry no computed value)
  NUMERIC_ONLY   same shape, same keys, same strings; only float cells differ,
                 all within tolerance -> reproducible up to float formatting
  DIFFERS        anything else: changed cell values, changed row/column sets,
                 changed row counts
  MISSING_NEW    present in the backup, absent from the re-run
  EXTRA_NEW      present in the re-run, absent from the backup

Markdown reports and figures are compared separately: reports by content hash
with a diff of the non-timestamp lines (report headers embed a generation
timestamp, so a hash difference there is expected and is reported as
TIMESTAMP_ONLY when the only differing lines match a date/time pattern);
figures by byte size only, since matplotlib PNG bytes are not required to be
reproducible across runs and pixel equality is not the property under test.

Usage:
  vr_rerun_diff.py --old results.prerun_backup --new results --out results/rerun_diff.tsv
"""
from __future__ import annotations

import argparse
import hashlib
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np

RTOL = 1e-9
ATOL = 1e-12
TIMESTAMP_RE = re.compile(
    r"\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}"          # 2026-08-20 14:03
    r"|\d{4}-\d{2}-\d{2}"                          # bare date
    r"|generated|Generated|timestamp|run at|UTC"
)


def sha(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def content_lines(path: str) -> List[str]:
    """File lines with '#' comment lines dropped (provenance, not data)."""
    with open(path, "r", errors="replace") as fh:
        return [ln.rstrip("\n") for ln in fh if not ln.startswith("#")]


def walk(root: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames
                       if not d.startswith("refcheck_") and d != "__pycache__"]
        for fn in filenames:
            if fn.startswith("."):
                continue
            full = os.path.join(dirpath, fn)
            out[os.path.relpath(full, root)] = full
    return out


def read_table(path: str) -> Optional[pd.DataFrame]:
    try:
        df = pd.read_csv(path, sep="\t", comment="#", dtype=str,
                         keep_default_na=False)
    except Exception:
        return None
    return df


def key_columns(df: pd.DataFrame) -> List[str]:
    prefer = ["trial", "sample", "cluster_id", "supercluster_id", "gene_id",
              "gene_name", "family", "rank", "tier", "cluster"]
    return [c for c in prefer if c in df.columns]


def compare_tables(old: str, new: str) -> Tuple[str, str]:
    a, b = read_table(old), read_table(new)
    if a is None or b is None:
        return ("DIFFERS", "unparseable as TSV by one side")
    if list(a.columns) != list(b.columns):
        only_a = [c for c in a.columns if c not in b.columns]
        only_b = [c for c in b.columns if c not in a.columns]
        if not only_a and not only_b:
            return ("DIFFERS", f"column ORDER changed: {list(a.columns)[:6]} vs {list(b.columns)[:6]}")
        return ("DIFFERS",
                f"columns changed; only_backup={only_a[:8]} only_rerun={only_b[:8]}")
    if len(a) != len(b):
        return ("DIFFERS", f"row count {len(a)} -> {len(b)}")

    keys = key_columns(a)
    if keys:
        a = a.sort_values(keys, kind="mergesort").reset_index(drop=True)
        b = b.sort_values(keys, kind="mergesort").reset_index(drop=True)
        for k in keys:
            if not a[k].equals(b[k]):
                n = int((a[k] != b[k]).sum())
                return ("DIFFERS", f"key column '{k}' differs in {n} rows")

    numeric_cols, str_diffs, num_diffs = [], [], []
    for col in a.columns:
        s_a, s_b = a[col], b[col]
        if s_a.equals(s_b):
            continue
        fa = pd.to_numeric(s_a, errors="coerce")
        fb = pd.to_numeric(s_b, errors="coerce")
        both_num = fa.notna() & fb.notna()
        # a column is "numeric" if every disagreeing cell parses as a float
        disagree = s_a != s_b
        if bool((disagree & ~both_num).any()):
            idx = list(a.index[disagree & ~both_num])[:3]
            ex = "; ".join(
                f"row{i} {keys and a.loc[i, keys[-1]] or i}: "
                f"'{s_a.loc[i]}' -> '{s_b.loc[i]}'" for i in idx)
            str_diffs.append(f"{col}: {int(disagree.sum())} cells [{ex}]")
            continue
        numeric_cols.append(col)
        close = np.isclose(fa[disagree], fb[disagree], rtol=RTOL, atol=ATOL,
                           equal_nan=True)
        if not bool(close.all()):
            bad = fa[disagree][~close]
            i = bad.index[0]
            num_diffs.append(
                f"{col}: {int((~close).sum())}/{int(disagree.sum())} cells beyond "
                f"tol (e.g. row {i}: {fa.loc[i]!r} -> {fb.loc[i]!r})")

    if str_diffs or num_diffs:
        return ("DIFFERS", " | ".join((str_diffs + num_diffs)[:4]))
    if numeric_cols:
        return ("NUMERIC_ONLY",
                f"float formatting only in {len(numeric_cols)} col(s): "
                f"{numeric_cols[:6]}")
    return ("IDENTICAL", "")


def compare_text(old: str, new: str) -> Tuple[str, str]:
    la, lb = content_lines(old), content_lines(new)
    if la == lb:
        return ("IDENTICAL", "")
    sa, sb = set(la), set(lb)
    only_a = [x for x in la if x not in sb]
    only_b = [x for x in lb if x not in sa]
    changed = only_a + only_b
    if changed and all(TIMESTAMP_RE.search(x) for x in changed):
        return ("TIMESTAMP_ONLY", f"{len(changed)} timestamp/provenance lines")
    ex = []
    for x in only_a[:2]:
        ex.append(f"- {x[:120]}")
    for x in only_b[:2]:
        ex.append(f"+ {x[:120]}")
    return ("DIFFERS",
            f"{len(only_a)} lines only in backup, {len(only_b)} only in rerun :: "
            + " ".join(ex))


def compare_binary(old: str, new: str) -> Tuple[str, str]:
    za, zb = os.path.getsize(old), os.path.getsize(new)
    if sha(old) == sha(new):
        return ("IDENTICAL", "")
    rel = abs(za - zb) / max(za, 1)
    if rel < 0.02:
        return ("BINARY_EQUIV",
                f"bytes differ, size {za} -> {zb} ({rel*100:.2f}%); PNG bytes "
                f"are not required to be reproducible")
    return ("DIFFERS", f"size {za} -> {zb} ({rel*100:.1f}% change)")


def classify(rel: str, old: str, new: str) -> Tuple[str, str]:
    ext = os.path.splitext(rel)[1].lower()
    if ext in (".tsv", ".csv"):
        return compare_tables(old, new)
    if ext in (".png", ".pdf", ".svg", ".jpg"):
        return compare_binary(old, new)
    return compare_text(old, new)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old", required=True, help="backed-up results tree")
    ap.add_argument("--new", required=True, help="regenerated results tree")
    ap.add_argument("--out", required=True, help="TSV report path")
    a = ap.parse_args(argv)

    A, B = walk(a.old), walk(a.new)
    rows = []
    for rel in sorted(set(A) | set(B)):
        if rel in A and rel not in B:
            rows.append((rel, "MISSING_NEW", "in backup, not regenerated"))
        elif rel in B and rel not in A:
            rows.append((rel, "EXTRA_NEW", "produced by re-run, absent from backup"))
        else:
            verdict, detail = classify(rel, A[rel], B[rel])
            rows.append((rel, verdict, detail))

    df = pd.DataFrame(rows, columns=["path", "verdict", "detail"])
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    df.to_csv(a.out, sep="\t", index=False)

    counts = df["verdict"].value_counts().to_dict()
    print(f"[rerun_diff] {len(df)} paths compared -> {a.out}")
    for k in sorted(counts):
        print(f"[rerun_diff]   {k}: {counts[k]}")
    bad = df[df["verdict"].isin(["DIFFERS", "MISSING_NEW"])]
    if len(bad):
        print("\n[rerun_diff] NON-REPRODUCIBLE PATHS:")
        for _, r in bad.iterrows():
            print(f"  {r['verdict']:<12} {r['path']}\n      {r['detail']}")
    else:
        print("[rerun_diff] every backed-up output was reproduced")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
