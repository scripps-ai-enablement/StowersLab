#!/usr/bin/env python3
"""
vr_report.py -- tiered reporting engine for the VNO receptor pipeline.

The point of this module is that the priority hierarchy is enforced
STRUCTURALLY, not by convention. A tier's content is produced by a callable
that `TierGate.emit()` refuses to invoke unless every upstream tier passed for
that sample. There is no code path in the report writer that can reach a
receptor call for a library that failed tissue identity, sort validation or
library viability -- not because the writer is careful, but because the
producer function is never called and the gate substitutes an explicit
suppression line naming the failing tier and its reason.

Tier ladder
-----------
    tier 0  tissue_identity   is this VNO at all?           (Phase-0 finding)
    tier 1  sort_validation   is this a viable GFP+/GFP- library? (Rule 1 + 4)
    tier 2  population_id     V1R vs V2R via Gnai2:Gnao1    (Rule 3)
    tier 3  cluster_vr        cluster-level VR calls         (reliable tier)
    tier 4  individual_vr     individual receptor candidates (always tentative)

Verdicts are READ from the QC layer's flag columns (tissue_verdict,
sort_verdict, population_call, library_status, suppress_biology,
blocking_flags) in results/sample_qc_all.tsv. This module never re-derives a
verdict from CPM values -- one implementation of each rule, in vr_markers.py.

Tiers 3 and 4 consume the quantification track's tables:
    results/<trial>/vr_cluster_expression.tsv
    results/<trial>/vr_within_cluster_fractions.tsv
    results/<trial>/vr_artifact_flags.tsv
    results/<trial>/vr_candidates.tsv
If a table is absent the tier reports status NO_DATA with the path it looked
for, and the rest of the report still builds. That is deliberate: the report
runs on QC alone and fills in the VR tiers when they land.

Usage
-----
    python vr_report.py                      # all trials + combined report
    python vr_report.py --trial trial2
    python vr_report.py --selftest           # synthetic-record gate tests
"""

from __future__ import annotations

import argparse
import datetime as _dt
import os
import sys
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Sequence

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from vr_config import (  # noqa: E402
    load_config,
    marker_genes,
    samples_of,
    threshold,
    trial_paths,
    trials_of,
)

__all__ = [
    "TIERS",
    "TIER_NAMES",
    "TierOutcome",
    "TierSuppressed",
    "TierGate",
    "VR_TABLES",
    "load_qc_records",
    "load_vr_tables",
    "build_sample_report",
    "tier_status_table",
    "write_trial_report",
    "write_combined_report",
    "selftest",
]

# --------------------------------------------------------------------------
# Tier ladder
# --------------------------------------------------------------------------

TIERS: List[int] = [0, 1, 2, 3, 4]
TIER_NAMES: Dict[int, str] = {
    0: "tissue_identity",
    1: "sort_validation",
    2: "population_id",
    3: "cluster_vr",
    4: "individual_vr",
}
TIER_TITLES: Dict[int, str] = {
    0: "Tier 0 -- tissue identity (VNO vs main olfactory epithelium)",
    1: "Tier 1 -- sort validation and library viability",
    2: "Tier 2 -- population identification (V1R vs V2R)",
    3: "Tier 3 -- cluster-level VR calls (reliable tier)",
    4: "Tier 4 -- individual receptor candidates (tentative)",
}

# Statuses. PASS and PASS_CAVEAT are the only ones that let a downstream tier run.
PASS = "PASS"
PASS_CAVEAT = "PASS_WITH_CAVEAT"
FAIL = "FAIL"
SUPPRESSED = "SUPPRESSED"
NO_DATA = "NO_DATA"
_OPENING = (PASS, PASS_CAVEAT)

# Documented filenames from the quantification track, relative to
# results/<trial>/. Tier 3/4 degrade to NO_DATA when one is absent.
VR_TABLES: Dict[str, str] = {
    "cluster_expression": "vr_cluster_expression.tsv",
    "within_cluster_fractions": "vr_within_cluster_fractions.tsv",
    "artifact_flags": "vr_artifact_flags.tsv",
    "candidates": "vr_candidates.tsv",
}

# QC columns this module reads. Missing ones are treated as unknown, which
# fails closed (a tier cannot pass on an absent verdict).
QC_VERDICT_COLUMNS = [
    "tissue_verdict",
    "sort_verdict",
    "population_call",
    "library_status",
    "suppress_biology",
    "blocking_flags",
    "qc_overall",
]

# Tissue verdicts that establish VNO identity.
_TISSUE_VNO = {"VNO", "VNO_dominant_mixed"}
_TISSUE_WRONG = {"MOE"}
# Library statuses that invalidate every CPM ratio in the sample.
_LIBRARY_DEAD = {"FAILED", "DEGENERATE"}
_POPULATION_CALLS = {"V1R_dominant", "V2R_dominant", "mixed", "ambiguous"}


def _truthy(v: Any) -> bool:
    """Tolerant bool for TSV round-trips ('True'/'true'/1/True)."""
    if isinstance(v, bool):
        return v
    if v is None:
        return False
    if isinstance(v, (int, float)):
        try:
            if pd.isna(v):
                return False
        except (TypeError, ValueError):
            pass
        return bool(v)
    return str(v).strip().lower() in {"true", "t", "yes", "y", "1"}


def _s(rec: Dict[str, Any], key: str, default: str = "") -> str:
    v = rec.get(key, default)
    if v is None:
        return default
    try:
        if isinstance(v, float) and pd.isna(v):
            return default
    except (TypeError, ValueError):
        pass
    return str(v)


def _f(rec: Dict[str, Any], key: str) -> Optional[float]:
    v = rec.get(key)
    if v is None or v == "":
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if pd.isna(f) else f


# --------------------------------------------------------------------------
# Tier outcome + gate
# --------------------------------------------------------------------------


@dataclass
class TierOutcome:
    """One tier's result for one sample."""

    tier: int
    name: str
    status: str
    reason: str
    blocked_by: Optional[int] = None
    lines: List[str] = field(default_factory=list)
    data: Dict[str, Any] = field(default_factory=dict)

    @property
    def opens_downstream(self) -> bool:
        return self.status in _OPENING

    def as_row(self, trial: str, sample: str) -> Dict[str, Any]:
        return {
            "trial": trial,
            "sample": sample,
            "tier": self.tier,
            "tier_name": self.name,
            "status": self.status,
            "blocked_by_tier": self.blocked_by if self.blocked_by is not None else "",
            "blocked_by_name": TIER_NAMES.get(self.blocked_by, "")
            if self.blocked_by is not None
            else "",
            "reason": self.reason,
        }


class TierSuppressed(Exception):
    """Raised by TierGate.require() when an upstream tier did not pass."""

    def __init__(self, tier: int, blocked_by: int, reason: str):
        self.tier = tier
        self.blocked_by = blocked_by
        self.reason = reason
        super().__init__(
            f"tier {tier} ({TIER_NAMES.get(tier, '?')}) suppressed: upstream tier "
            f"{blocked_by} ({TIER_NAMES.get(blocked_by, '?')}) did not pass -- {reason}"
        )


class TierGate:
    """
    Structural enforcement of the priority hierarchy for ONE sample.

    Tiers 0-2 are decided from the QC layer's verdict columns at construction
    time. Tiers 3-4 are decided by `emit()`: the producer callable supplied for
    a tier is invoked ONLY when every lower tier opened. When it is not
    invoked, the gate records an explicit SUPPRESSED outcome naming the lowest
    failing tier and its reason.

    A caller cannot bypass the gate by "just calling the producer": the report
    writer only ever reaches tier content through `emit()`, and `require()`
    raises for any imperative code that tries to compute a tier directly.
    """

    def __init__(self, record: Dict[str, Any], cfg: Optional[Dict[str, Any]] = None):
        self.record = dict(record)
        self.cfg = cfg
        self.trial = _s(self.record, "trial", "NA")
        self.sample = _s(self.record, "sample", "NA")
        self.cell_type = _s(self.record, "cell_type", "unknown")
        self.outcomes: Dict[int, TierOutcome] = {}
        self._seal_qc_tiers()

    # -- QC-derived tiers -------------------------------------------------

    def _seal_qc_tiers(self) -> None:
        for tier, decider in ((0, self._decide_tier0), (1, self._decide_tier1),
                              (2, self._decide_tier2)):
            blocked = self.first_failure(below=tier)
            if blocked is not None:
                self.outcomes[tier] = TierOutcome(
                    tier=tier,
                    name=TIER_NAMES[tier],
                    status=SUPPRESSED,
                    reason=(
                        f"not evaluated -- tier {blocked} ({TIER_NAMES[blocked]}) "
                        f"{self._blocked_reason(blocked)}"
                    ),
                    blocked_by=blocked,
                )
            else:
                self.outcomes[tier] = decider()

    def _decide_tier0(self) -> TierOutcome:
        v = _s(self.record, "tissue_verdict", "unknown")
        mk = lambda st, rs: TierOutcome(0, TIER_NAMES[0], st, rs,  # noqa: E731
                                       data={"tissue_verdict": v})
        if v in _TISSUE_VNO:
            caveat = v == "VNO_dominant_mixed"
            return mk(
                PASS_CAVEAT if caveat else PASS,
                f"tissue_verdict={v}"
                + (
                    " -- VNO dominant with a real minor MOE component; treat as a "
                    "sort-purity note, not a tissue failure"
                    if caveat
                    else " -- VNO-specific panel establishes vomeronasal identity"
                ),
            )
        if v in _TISSUE_WRONG:
            return mk(
                FAIL,
                "tissue_verdict=MOE -- the main-olfactory marker panel is above the "
                "tissue floor while the VNO-specific panel (Trpc2, Vmn1r*, Vmn2r*) is "
                "not. This is the WRONG TISSUE, not a failed VNO sort; Trpc2=0 is the "
                "expected value for main olfactory epithelium. No VR biology is "
                "reportable and the remedy is wet-lab (dissection / sort gate), not "
                "re-quantification.",
            )
        if v == "no_tissue_signal":
            if self.cell_type == "nontarget":
                return mk(
                    PASS_CAVEAT,
                    "tissue_verdict=no_tissue_signal -- expected for a GFP- nontarget "
                    "library: absence of the VNO-specific panel is the designed "
                    "outcome, so tissue identity is not a gate for this sample. It "
                    "serves as the sort-specificity control, not as VNO tissue "
                    "evidence.",
                )
            return mk(
                FAIL,
                "tissue_verdict=no_tissue_signal -- both marker panels sit below the "
                "absolute CPM floor, so this target library carries no tissue "
                "information either way. VNO identity is UNCONFIRMED (this is not "
                "positive evidence of the wrong tissue), and an unconfirmed tissue "
                "cannot support a receptor call.",
            )
        return mk(FAIL, f"tissue_verdict={v!r} is not a recognised verdict -- failing "
                        "closed rather than assuming VNO identity")

    def _decide_tier1(self) -> TierOutcome:
        sort_v = _s(self.record, "sort_verdict", "unknown")
        lib = _s(self.record, "library_status", "unknown")
        suppress = _truthy(self.record.get("suppress_biology"))
        data = {"sort_verdict": sort_v, "library_status": lib,
                "suppress_biology": suppress}
        mk = lambda st, rs: TierOutcome(1, TIER_NAMES[1], st, rs, data=data)  # noqa: E731

        # Library viability first: CPM is a ratio, so a dead library's sort
        # measurement is not merely failing, it is meaningless.
        if lib in _LIBRARY_DEAD:
            why = _s(self.record, "library_reason") or f"library_status={lib}"
            return mk(
                FAIL,
                f"library_status={lib} -- {why} Every CPM in this sample is a ratio "
                "taken against a near-empty or failed library, so no downstream "
                "number is interpretable however large it looks.",
            )
        if suppress:
            flags = _s(self.record, "blocking_flags")
            return mk(FAIL, "QC layer set suppress_biology=True"
                            + (f" (blocking_flags: {flags})" if flags else ""))
        if sort_v == "PASS":
            note = _s(self.record, "sort_reason") or "sort_verdict=PASS"
            return mk(PASS, note)
        if sort_v.startswith("FAIL"):
            note = _s(self.record, "sort_reason") or f"sort_verdict={sort_v}"
            return mk(FAIL, f"sort_verdict={sort_v} -- {note}")
        return mk(FAIL, f"sort_verdict={sort_v!r} is not a recognised verdict -- "
                        "failing closed")

    def _decide_tier2(self) -> TierOutcome:
        call = _s(self.record, "population_call", "undetermined")
        note = _s(self.record, "population_note")
        low = _truthy(self.record.get("ratio_low_support"))
        ratio_str = _s(self.record, "gnai2_gnao1_ratio_str")
        data = {"population_call": call, "ratio_low_support": low,
                "gnai2_gnao1_ratio_str": ratio_str}
        mk = lambda st, rs: TierOutcome(2, TIER_NAMES[2], st, rs, data=data)  # noqa: E731
        if call in _POPULATION_CALLS:
            if low:
                return mk(
                    PASS_CAVEAT,
                    f"population_call={call}; Gnai2:Gnao1 = {ratio_str or 'NA'} but the "
                    "denominator has too few raw reads for the ratio MAGNITUDE to be "
                    "quotable (Poisson-unstable). The direction of the call stands; "
                    f"the number does not. {note}".strip(),
                )
            return mk(PASS, f"population_call={call}; Gnai2:Gnao1 = "
                            f"{ratio_str or 'NA'}. {note}".strip())
        return mk(FAIL, f"population_call={call} -- {note or 'no population signal'}")

    # -- gate mechanics ---------------------------------------------------

    def first_failure(self, below: int) -> Optional[int]:
        """
        Lowest tier < `below` that did not open the next tier. A tier that has
        not been evaluated yet counts as a failure: permission for tier N can
        never be granted before tier N-1 has actually been decided.
        """
        for t in TIERS:
            if t >= below:
                break
            oc = self.outcomes.get(t)
            if oc is None or not oc.opens_downstream:
                return t
        return None

    def _blocked_reason(self, blocked: int) -> str:
        up = self.outcomes.get(blocked)
        if up is None:
            return (f"NOT_EVALUATED: tier {blocked} ({TIER_NAMES[blocked]}) has not "
                    "been decided yet, so no downstream tier can be permitted")
        return f"{up.status}: {up.reason}"

    def allowed(self, tier: int) -> bool:
        return self.first_failure(below=tier) is None

    def require(self, tier: int) -> None:
        """Imperative guard: raise TierSuppressed unless every upstream tier opened."""
        blocked = self.first_failure(below=tier)
        if blocked is not None:
            raise TierSuppressed(tier, blocked, self._blocked_reason(blocked))

    def emit(self, tier: int, producer: Callable[[Dict[str, Any]], TierOutcome]) -> TierOutcome:
        """
        Record tier `tier`. `producer` is called ONLY when every upstream tier
        opened; otherwise it is not called at all and a SUPPRESSED outcome is
        recorded naming the lowest failing tier. A producer that raises
        TierSuppressed is converted to the same suppression record rather than
        propagating.
        """
        if tier in self.outcomes and tier <= 2:
            return self.outcomes[tier]  # QC tiers are sealed at construction
        blocked = self.first_failure(below=tier)
        if blocked is not None:
            oc = TierOutcome(
                tier=tier,
                name=TIER_NAMES[tier],
                status=SUPPRESSED,
                reason=(
                    f"suppressed by tier {blocked} ({TIER_NAMES[blocked]}) "
                    f"{self._blocked_reason(blocked)}"
                ),
                blocked_by=blocked,
            )
        else:
            try:
                oc = producer(self.record)
            except TierSuppressed as exc:  # defence in depth
                oc = TierOutcome(tier, TIER_NAMES[tier], SUPPRESSED, str(exc),
                                 blocked_by=exc.blocked_by)
        self.outcomes[tier] = oc
        return oc

    # -- summary ----------------------------------------------------------

    def highest_tier_reported(self) -> int:
        """Highest tier whose status is PASS/PASS_WITH_CAVEAT; -1 if none."""
        best = -1
        for t in TIERS:
            oc = self.outcomes.get(t)
            if oc is not None and oc.opens_downstream:
                best = t
            else:
                break
        return best

    def stop_reason(self) -> str:
        top = self.highest_tier_reported()
        nxt = top + 1
        oc = self.outcomes.get(nxt)
        if oc is None:
            return "all tiers reported"
        return f"stopped at tier {nxt} ({oc.name}) [{oc.status}]: {oc.reason}"

    def summary_row(self) -> Dict[str, Any]:
        top = self.highest_tier_reported()
        return {
            "trial": self.trial,
            "sample": self.sample,
            "cell_type": self.cell_type,
            "qc_overall": _s(self.record, "qc_overall"),
            "highest_tier_reported": top,
            "highest_tier_name": TIER_NAMES.get(top, "none"),
            "stopped_at_tier": top + 1 if (top + 1) in self.outcomes else "",
            "stopped_at_name": TIER_NAMES.get(top + 1, ""),
            "stop_status": self.outcomes[top + 1].status if (top + 1) in self.outcomes else "",
            "stop_reason": self.stop_reason(),
            **{f"tier{t}_status": self.outcomes[t].status if t in self.outcomes else ""
               for t in TIERS},
        }


# --------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------


def load_qc_records(cfg, trial: Optional[str] = None,
                    qc_path: Optional[str] = None) -> List[Dict[str, Any]]:
    """
    Read results/sample_qc_all.tsv (or a per-trial sample_qc.tsv) into records.
    Comment lines beginning '#' carry the CPM-convention provenance note and
    are skipped by the reader but recoverable via `read_provenance_note`.
    """
    work = cfg.get("work")
    path = qc_path or os.path.join(work, "results", "sample_qc_all.tsv")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"QC table not found: {path} -- run bin/vr_sample_qc.py first")
    df = pd.read_csv(path, sep="\t", comment="#", dtype=str, keep_default_na=False)
    missing = [c for c in QC_VERDICT_COLUMNS if c not in df.columns]
    if missing:
        raise KeyError(f"{path} is missing required verdict columns: {missing}")
    if trial:
        df = df[df["trial"] == trial]
    recs = df.to_dict("records")
    # stable order: config sample order within trial
    order: Dict[str, int] = {}
    for t in trials_of(cfg):
        for i, s in enumerate(samples_of(cfg, t)):
            order[f"{t}/{s}"] = i
    recs.sort(key=lambda r: (r.get("trial", ""),
                             order.get(f"{r.get('trial')}/{r.get('sample')}", 999)))
    return recs


def read_provenance_note(path: str) -> str:
    """First '#' comment line of a table (the CPM-convention note)."""
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    return line.lstrip("#").strip()
                break
    except OSError:
        pass
    return ""


def load_vr_tables(cfg, trial: str) -> Dict[str, Any]:
    """
    Load the quantification track's tables for one trial. Returns
    {key: DataFrame | None, "_missing": {key: path}, "_present": {key: path}}.
    Absent tables are None -- tiers 3/4 then report NO_DATA with the path.
    """
    out_dir = trial_paths(cfg, trial)["out_dir"]
    res: Dict[str, Any] = {"_missing": {}, "_present": {}}
    for key, fname in VR_TABLES.items():
        p = os.path.join(out_dir, fname)
        if os.path.exists(p):
            try:
                res[key] = pd.read_csv(p, sep="\t", comment="#")
                res["_present"][key] = p
            except Exception as exc:  # unreadable table is not a silent pass
                res[key] = None
                res["_missing"][key] = f"{p} (unreadable: {exc})"
        else:
            res[key] = None
            res["_missing"][key] = p
    return res


def _rows_for_sample(df: Optional[pd.DataFrame], sample: str) -> Optional[pd.DataFrame]:
    if df is None or len(df) == 0:
        return None
    for col in ("sample", "sample_id", "Sample"):
        if col in df.columns:
            sub = df[df[col].astype(str) == sample]
            return sub
    # wide matrix: clusters x samples
    if sample in df.columns:
        idcol = df.columns[0]
        sub = df[[idcol, sample]].copy()
        sub.columns = ["cluster_id", "value"]
        sub["sample"] = sample
        return sub
    return None


_SC_CACHE: Dict[str, Dict[str, str]] = {}


def _supercluster_map(cfg) -> Dict[str, str]:
    """
    cluster_id -> supercluster_id from ref/vr_clusters.tsv. Empty dict if the
    reference table is absent, so the dual-tier note degrades to silence rather
    than raising.
    """
    path = os.path.join(cfg.get("work", ""), "ref", "vr_clusters.tsv")
    if path in _SC_CACHE:
        return _SC_CACHE[path]
    m: Dict[str, str] = {}
    if os.path.exists(path):
        try:
            df = pd.read_csv(path, sep="\t", usecols=["cluster_id", "supercluster_id"])
            m = dict(zip(df["cluster_id"].astype(str),
                         df["supercluster_id"].astype(str)))
        except Exception:
            m = {}
    _SC_CACHE[path] = m
    return m


def _pick(df: pd.DataFrame, *names: str) -> Optional[str]:
    for n in names:
        if n in df.columns:
            return n
    return None


# --------------------------------------------------------------------------
# Tier 3 / tier 4 producers
# --------------------------------------------------------------------------


def make_cluster_producer(cfg, vr: Dict[str, Any], trial: str) -> Callable:
    """Producer for tier 3: cluster-level VR calls."""

    def produce(rec: Dict[str, Any]) -> TierOutcome:
        sample = _s(rec, "sample")
        df = vr.get("cluster_expression")
        if df is None:
            return TierOutcome(
                3, TIER_NAMES[3], NO_DATA,
                "quantification track has not written "
                f"{vr['_missing'].get('cluster_expression')} yet -- tier 3 will "
                "populate on the next run of this report; QC tiers above are final.",
            )
        # The table carries BOTH cluster (200kb) and supercluster (800kb) rows in
        # a `tier` column. Mixing them would double-count every region, so the
        # 200kb rows drive the call and the 800kb rows supply the dual-tier note.
        d = df.copy()
        if "tier" in d.columns:
            cl_rows = d[d["tier"].astype(str) == "cluster"]
            sc_rows = d[d["tier"].astype(str) == "supercluster"]
        else:
            cl_rows, sc_rows = d, d.iloc[0:0]
        sub = _rows_for_sample(cl_rows, sample)
        if sub is None or len(sub) == 0:
            return TierOutcome(
                3, TIER_NAMES[3], NO_DATA,
                f"no rows for sample {sample} in "
                f"{vr['_present'].get('cluster_expression')}",
            )
        valcol = _pick(sub, "cpm_sum", "cluster_cpm", "cpm", "cluster_sum_cpm", "value")
        idcol = _pick(sub, "cluster_id", "cluster")
        if valcol is None or idcol is None:
            return TierOutcome(
                3, TIER_NAMES[3], NO_DATA,
                f"cluster table present but lacks an identifiable cluster/value "
                f"column pair (saw: {list(sub.columns)})",
            )
        sub = sub.copy()
        sub[valcol] = pd.to_numeric(sub[valcol], errors="coerce")
        sub = sub.sort_values(valcol, ascending=False)

        # `is_called` is the quantification track's own call flag. Use it when
        # present rather than re-deriving a threshold here: one implementation
        # of the calling rule, in the track that owns it.
        callcol = _pick(sub, "is_called")
        if callcol is not None:
            called = sub[pd.to_numeric(sub[callcol], errors="coerce").fillna(0) > 0]
            basis = "is_called flag from the quantification track"
        else:
            called = sub[sub[valcol] > 0]
            basis = "non-zero cluster CPM (no is_called column in the table)"
        n_nz = int(len(called))

        n_cells = _s(rec, "n_cells")
        ct = _s(rec, "cell_type")
        key = ("nontarget_100cell" if ct == "nontarget"
               else "target_2cell" if n_cells == "2" else "target_100cell")
        exp = (threshold(cfg, "expected_clusters", {}) or {}).get(key)

        lines = [
            f"Clusters called: **{n_nz}** ({basis})"
            + (f"; expected {exp[0]}-{exp[1]} for a {key.replace('_', ' ')}"
               if isinstance(exp, (list, tuple)) and len(exp) == 2 else ""),
        ]
        if isinstance(exp, (list, tuple)) and len(exp) == 2:
            lo, hi = exp
            if n_nz < lo:
                lines.append(
                    f"* **Below the expected range** ({n_nz} < {lo}). Read depth and "
                    "capture efficiency bound what can be detected, so absence here "
                    "is not evidence that a receptor is unexpressed.")
            elif n_nz > hi:
                if ct == "nontarget":
                    lines.append(
                        f"* **Above the expected range** ({n_nz} > {hi}). A GFP- "
                        "control should carry no VR cluster signal, so any cluster "
                        "here is contamination or carry-over -- see the sort-purity "
                        "line below.")
                else:
                    lines.append(
                        f"* **Above the expected range** ({n_nz} > {hi}). Each VNO "
                        "neuron expresses exactly one receptor, so excess clusters "
                        "indicate more cells than intended, ambient RNA, or EM spread "
                        "ACROSS clusters. Treat the set as a candidate list, not a "
                        "co-expression result.")
        if ct == "nontarget":
            cap = threshold(cfg, "nontarget_total_vr_cpm_max", default=100)
            tot = float(called[valcol].sum()) if n_nz else 0.0
            lines.append(
                f"* **Sort-purity check:** GFP- control, {n_nz} cluster(s) called, "
                f"total {tot:.2f} CPM against a {cap} CPM ceiling. "
                + ("Within the ceiling." if tot <= float(cap) else
                   "**Above the ceiling** -- VR-expressing cells entered the GFP- "
                   "fraction, weakening every target/nontarget contrast.")
                + " No receptor is called from a nontarget library either way.")

        # dual-tier: 200kb clusters sharing an 800kb supercluster
        sc = _supercluster_map(cfg)
        if sc and n_nz:
            from collections import defaultdict
            by_sc = defaultdict(list)
            for cid in called[idcol].astype(str):
                s = sc.get(cid)
                if s:
                    by_sc[s].append(cid)
            for s, members in ((k, v) for k, v in by_sc.items() if len(v) > 1):
                lines.append(
                    f"* **Dual-tier note:** {', '.join(f'`{m}`' for m in members)} are "
                    f"all called and all belong to supercluster `{s}` (800 kb tier). "
                    f"Under the 800 kb definition this is ONE genomic region, not "
                    f"{len(members)} independent calls. The 200 kb split is the "
                    "conservative choice, not a natural break -- both tiers stand.")
        # the track's own chr7 dual-tier marker, if it set one
        d7 = _pick(sub, "chr7_dual_tier_region")
        if d7 is not None:
            hits = called[pd.to_numeric(called[d7], errors="coerce").fillna(0) > 0]
            if len(hits):
                lines.append(
                    "* The quantification track marked "
                    + ", ".join(f"`{c}`" for c in hits[idcol].astype(str))
                    + " as lying in the chr7 dual-tier region (historical \"cluster "
                    "039\", 60 paralogs, split by the 200 kb rule at the 217,366 bp "
                    "Vmn1r132-Vmn1r135 gap and reunited by supercluster "
                    "`V1R_chr7_sc013`).")

        if n_nz:
            shcol = _pick(called, "share_of_sample_vr")
            dcol = _pick(called, "n_member_genes_detected")
            mcol = _pick(called, "n_member_genes")
            lines += ["", "| cluster | family | CPM | share of sample VR | paralogs "
                      "detected / in cluster |", "|---|---|---:|---:|---:|"]
            for _, rr in called.head(10).iterrows():
                fam = str(rr.get("family", ""))
                sh = pd.to_numeric(pd.Series([rr[shcol]]), errors="coerce").iloc[0] \
                    if shcol else float("nan")
                lines.append(
                    f"| `{rr[idcol]}` | {fam} | {rr[valcol]:.1f} | "
                    + (f"{sh:.3f}" if sh == sh else "NA") + " | "
                    + (f"{rr[dcol]:.0f} / {rr[mcol]:.0f}" if dcol and mcol else "NA")
                    + " |")

        # Per-cluster QC of the within-cluster structure. TWO INDEPENDENT
        # questions are asked of every called cluster and they have different
        # answers, so each cluster gets exactly ONE entry carrying both:
        #
        #   (1) redistribution: did the EM split one transcript's reads across
        #       paralogs? Answered by `em_flag` (Monte-Carlo test against
        #       uniformity + unique-read gate).
        #   (2) pseudogene bleed: does a pseudogene inside a functional cluster
        #       carry signal? Answered by `pseudogene_bleed_flag`, and its
        #       mechanism is UNRESOLVED -- it is an open question, not a verdict.
        #
        # Listing a cluster under both a "judged not real" heading and a "found
        # clean" heading, as an earlier version did, is a contradiction: clean on
        # (1) and flagged on (2) is a coherent and common state, and the reader
        # needs to see it as one row, not two conflicting ones.
        flags = _rows_for_sample(vr.get("artifact_flags"), sample)
        redist_lines: List[str] = []      # (1) fired: the EM is misleading here
        pseudo_lines: List[str] = []      # (2) fired: open mechanism question
        clean_lines: List[str] = []       # neither fired
        if flags is not None and len(flags):
            f = flags.copy()
            if "scope" in f.columns:
                f = f[f["scope"].astype(str) == "cluster"]
            emc = _pick(f, "em_flag")
            lvl = _pick(f, "em_flag_level")
            note = _pick(f, "em_flag_note")
            ccol2 = _pick(f, "cluster_id")
            rs = _pick(f, "read_support")
            ebs = _pick(f, "even_block_size")
            psb = _pick(f, "pseudogene_bleed_flag")
            psm = _pick(f, "pseudogene_mechanism")
            psg = _pick(f, "pseudogene_genes_detected")
            pss = _pick(f, "pseudogene_share_of_cluster")
            quiet = {"no_signal", "insufficient_signal", "none", "nan", ""}
            alarm_flags = {"suspected_em_redistribution"}

            def _emverdict(fr) -> str:
                em = str(fr[emc]) if emc else ""
                interp = str(fr.get("interpretation", "") or "")
                out = f"redistribution test: {em}" if em else "redistribution test: n/a"
                if interp and interp != "nan":
                    out += f" ({interp})"
                return out

            for _, fr in f.iterrows():
                em = str(fr[emc]) if emc else ""
                cid = f"`{fr[ccol2]}`" if ccol2 else "cluster"
                pb = str(fr[psb]) if psb else ""
                pseudo_fired = (pb == "apparent_pseudogene_expression")
                if em.lower() in quiet and not pseudo_fired:
                    continue

                if em in alarm_flags:
                    bits = [cid, f"**{em}**"]
                    if lvl and str(fr[lvl]).lower() not in quiet:
                        bits.append(f"level {fr[lvl]}")
                    if ebs and pd.notna(fr[ebs]):
                        try:
                            if float(fr[ebs]) > 1:
                                bits.append(
                                    f"even block of {float(fr[ebs]):.0f} paralogs")
                        except (TypeError, ValueError):
                            pass
                    if rs and pd.notna(fr[rs]):
                        try:
                            bits.append(f"{float(fr[rs]):.0f} reads of support")
                        except (TypeError, ValueError):
                            pass
                    if note and str(fr[note]) not in ("", "nan"):
                        bits.append(str(fr[note]))
                    if pseudo_fired:
                        bits.append("ALSO carries pseudogene signal (see the "
                                    "pseudogene section)")
                    redist_lines.append(" — ".join(b for b in bits if b))
                elif pseudo_fired:
                    # Clean on redistribution, but a pseudogene carries signal.
                    # State BOTH so the entry cannot be read as a redistribution
                    # call, and say plainly that the mechanism is unresolved.
                    genes = str(fr[psg]) if psg else ""
                    share = ""
                    if pss and pd.notna(fr[pss]):
                        try:
                            share = f", {float(fr[pss]) * 100:.1f}% of cluster signal"
                        except (TypeError, ValueError):
                            share = ""
                    mech = str(fr[psm]) if psm else ""
                    pseudo_lines.append(
                        f"{cid} — pseudogene signal present"
                        + (f" ({genes}{share})" if genes and genes != "nan"
                           else share)
                        + f". {_emverdict(fr)} — so this is NOT a redistribution "
                          "call. **Mechanism unresolved**: EM leakage from an "
                          "expressed paralog in the same cluster, or genuine "
                          "transcription of the pseudogene locus (cluster-shared "
                          "regulatory elements). This module does not adjudicate. "
                          "Dietschi et al. 2022 is NOT a quantitative prior for "
                          "mouse: their pseudogene result was significant in rat "
                          "only (P=0.003), not mouse (W=1214, P=0.5704), and their "
                          "proposed mechanism is regulatory, not multi-mapping."
                        + (f" Track note: {mech}"
                           if mech and mech != "nan"
                           and "Dietschi" not in mech else ""))
                else:
                    clean_lines.append(f"{cid} — {_emverdict(fr)}")

        if redist_lines:
            lines += ["", "**EM-redistribution findings** — the within-cluster split "
                      "here is an artifact of read assignment, not per-paralog "
                      "biology:"]
            lines += [f"- {fl}" for fl in redist_lines]
        elif clean_lines or pseudo_lines:
            lines += ["", "**No EM-redistribution artifact was found in any cluster "
                      "called for this sample.**"]
        if pseudo_lines:
            lines += ["", "Pseudogene signal inside functional clusters — an OPEN "
                      "mechanism question, not a redistribution verdict (each of "
                      "these passed the redistribution test):"]
            lines += [f"- {pl}" for pl in pseudo_lines]
        if clean_lines:
            lines += ["", "Clusters checked and clean on both tests (reportable "
                      "results, not warnings):"]
            lines += [f"- {cl}" for cl in clean_lines]

        return TierOutcome(
            3, TIER_NAMES[3], PASS,
            f"{n_nz} cluster(s) called; cluster-level aggregation is the reliable "
            "readout at 75bp",
            lines=lines,
            data={"n_clusters_called": n_nz,
                  "top_clusters": called.head(10)[[idcol, valcol]].to_dict("records"),
                  "expected_range": exp, "n_redistribution_findings": len(redist_lines),
                  "n_pseudogene_open": len(pseudo_lines),
                  "n_clusters_clean": len(clean_lines)},
        )

    return produce


def make_individual_producer(cfg, vr: Dict[str, Any], trial: str) -> Callable:
    """
    Producer for tier 4: individual receptor candidates. Every row carries its
    cluster context, its read support, and a confirmation status. The marker
    cross-check (Rule 3) is READ from the quantification track's
    `marker_consistency` column when present and only derived here as a
    fallback, so the rule has one implementation. Contradictions are surfaced,
    never smoothed over.
    """

    def produce(rec: Dict[str, Any]) -> TierOutcome:
        sample = _s(rec, "sample")
        cand = vr.get("candidates")
        frac = vr.get("within_cluster_fractions")
        if cand is None and frac is None:
            return TierOutcome(
                4, TIER_NAMES[4], NO_DATA,
                "quantification track has not written "
                f"{vr['_missing'].get('candidates')} / "
                f"{vr['_missing'].get('within_cluster_fractions')} yet -- no "
                "individual receptor call is made. Cluster-level results above stand "
                "on their own.",
            )
        sub = _rows_for_sample(cand, sample) if cand is not None else None
        if sub is None or len(sub) == 0:
            return TierOutcome(
                4, TIER_NAMES[4], NO_DATA,
                f"the candidate table carries no rows for {sample}: the "
                "quantification track found no cluster above its signal threshold "
                "from which to nominate a paralog. Cluster-level results stand.",
            )
        gene = _pick(sub, "gene_name", "gene", "receptor")
        gid = _pick(sub, "gene_id")
        clus = _pick(sub, "cluster_id", "cluster")
        scl = _pick(sub, "supercluster_id")
        fr = _pick(sub, "em_frac_of_cluster", "within_cluster_fraction", "fraction")
        rd = _pick(sub, "unique_reads_bam", "em_counts", "reads", "read_support")
        uq = _pick(sub, "unique_share_of_cluster")
        conf = _pick(sub, "confidence")
        cstat = _pick(sub, "confirmation_status")
        mcons = _pick(sub, "marker_consistency")
        emf = _pick(sub, "em_flag")
        fam = _pick(sub, "family")
        ncand = _pick(sub, "n_candidates_reported")
        evid = _pick(sub, "evidence_type")
        notes = _pick(sub, "notes")
        dom = _pick(sub, "is_dominant_cluster")

        # rows the track itself marks as "no call" are not candidates
        if conf is not None:
            real = sub[~sub[conf].astype(str).isin(["no_call", "nan", ""])]
        else:
            real = sub
        if len(real) == 0:
            reason = (str(sub.iloc[0][evid]) if evid else "no candidate nominated")
            return TierOutcome(
                4, TIER_NAMES[4], NO_DATA,
                f"the quantification track nominated no candidate for {sample} "
                f"(evidence_type = {reason}). No individual receptor is named; the "
                "cluster-level result is the statement that stands.",
            )

        pop = _s(rec, "population_call")
        lines = [
            "Every row is a TENTATIVE within-cluster assignment. At 75bp, paralogs "
            "inside one cluster share 85-95% nucleotide identity and reads are not "
            "uniquely assignable; Salmon's EM distributes them. **The cluster column, "
            "not the gene column, is the defensible unit.** The confirmation status "
            "and read support below are what a reader should weigh, not the gene name "
            "alone.",
            "",
            "| candidate | family | cluster (200 kb) | supercluster (800 kb) | EM "
            "fraction | unique reads | unique share | EM flag | confidence | "
            "confirmation | marker cross-check |",
            "|---|---|---|---|---:|---:|---:|---|---|---|---|",
        ]
        contradictions = 0

        def _num(v, spec="{:.3f}"):
            try:
                f = float(v)
                return "NA" if pd.isna(f) else spec.format(f)
            except (TypeError, ValueError):
                return "NA"

        for _, rr in real.iterrows():
            fam_v = str(rr[fam]) if fam else ""
            if mcons is not None and str(rr[mcons]) not in ("", "nan"):
                mc = str(rr[mcons])
                xcheck = {"consistent": "consistent (track check)",
                          "inconsistent": "**CONTRADICTION** (track check)"}.get(mc, mc)
                if mc.lower().startswith(("inconsist", "contradict", "conflict")):
                    contradictions += 1
            else:
                if fam_v.upper().startswith("V1R"):
                    ok = pop == "V1R_dominant"
                    xcheck = ("consistent (Gnai2>Gnao1)" if ok else
                              f"**CONTRADICTION**: V1R candidate in a {pop} sample")
                elif fam_v.upper().startswith("V2R"):
                    ok = pop == "V2R_dominant"
                    xcheck = ("consistent (Gnao1>Gnai2)" if ok else
                              f"**CONTRADICTION**: V2R candidate in a {pop} sample")
                else:
                    ok, xcheck = True, "n/a"
                if not ok:
                    contradictions += 1
            lines.append(
                "| {g}{gi} | {f} | {c}{d} | {sc} | {fr} | {rd} | {uq} | {ef} | {cf} "
                "| {cs} | {xc} |".format(
                    g=f"*{rr[gene]}*" if gene else "?",
                    gi=f"<br>`{rr[gid]}`" if gid else "",
                    f=fam_v or "?",
                    c=f"`{rr[clus]}`" if clus else "**unknown cluster — NOT "
                                                  "reportable**",
                    d=" (dominant)" if dom and str(rr[dom]) in ("1", "True", "true")
                      else "",
                    sc=f"`{rr[scl]}`" if scl and str(rr[scl]) not in ("", "nan")
                       else "—",
                    fr=_num(rr[fr]) if fr else "NA",
                    rd=_num(rr[rd], "{:.0f}") if rd else "NA",
                    uq=_num(rr[uq]) if uq else "NA",
                    ef=str(rr[emf]) if emf else "—",
                    cf=str(rr[conf]) if conf else "—",
                    cs=str(rr[cstat]) if cstat else "unconfirmed",
                    xc=xcheck,
                )
            )

        # candidate-set framing: the project accepts 2-3 candidates per cell
        if ncand is not None:
            per_cluster = (real.groupby(real[clus].astype(str))[ncand].first()
                           if clus else None)
            if per_cluster is not None and len(per_cluster):
                bits = ", ".join(f"`{k}`: {int(float(v))}"
                                 for k, v in per_cluster.items()
                                 if str(v) not in ("", "nan"))
                if bits:
                    lines += ["", "Candidates nominated per called cluster — "
                              + bits + ". The project's stated tolerance is 2-3 "
                              "candidates per cell, so a cluster reporting more than "
                              "one paralog is the expected outcome of 75bp reads, not "
                              "a failure of the assay."]
        ev = sorted({str(rr) for rr in real[evid]}) if evid else []
        if ev:
            lines += ["", "Evidence basis: " + ", ".join(f"`{e}`" for e in ev)
                      + ". Unique-read evidence (MAPQ-filtered BAM) is stronger than "
                      "an EM fraction, because it does not depend on the "
                      "redistribution step that makes within-cluster calls unsafe."]
        nn = sorted({str(rr) for rr in real[notes]} - {"", "nan"}) if notes else []
        if nn:
            lines += ["", "Track notes: " + " / ".join(nn)]

        if contradictions:
            lines += [
                "",
                f"**{contradictions} candidate(s) contradict the tier-2 population "
                f"call ({pop}).** Rule 3 requires Gnai2 > Gnao1 for a V1R call. The "
                "contradiction is reported, not resolved: either the population call "
                "or the within-cluster assignment is wrong, and the cluster-level "
                "result (tier 3) is the statement that survives.",
            ]
        # confirmation-status honesty: nothing here is confirmed
        if cstat is not None:
            statuses = sorted({str(v) for v in real[cstat]} - {"", "nan"})
            if statuses and all(s.startswith("tentative") or "unconfirmed" in s
                                for s in statuses):
                lines += [
                    "",
                    "Every candidate in this table carries confirmation status "
                    + ", ".join(f"`{s}`" for s in statuses)
                    + ". **No individual receptor identity is confirmed by this "
                      "pipeline.** Confirmation requires evidence that does not go "
                      "through the EM step — longer reads, targeted amplicon "
                      "sequencing of the cluster, or in-situ/immunostaining.",
                ]
        return TierOutcome(
            4, TIER_NAMES[4], PASS_CAVEAT,
            f"{len(real)} tentative candidate(s) across "
            f"{real[clus].nunique() if clus else 0} cluster(s); "
            f"{contradictions} marker contradiction(s); no identity confirmed",
            lines=lines,
            data={"n_candidates": int(len(real)),
                  "n_contradictions": contradictions},
        )

    return produce


# --------------------------------------------------------------------------
# Report assembly
# --------------------------------------------------------------------------


def build_sample_report(cfg, rec: Dict[str, Any], vr: Dict[str, Any]) -> TierGate:
    """Run the full tier ladder for one sample record. Returns the sealed gate."""
    gate = TierGate(rec, cfg)
    trial = _s(rec, "trial")
    gate.emit(3, make_cluster_producer(cfg, vr, trial))
    gate.emit(4, make_individual_producer(cfg, vr, trial))
    return gate


def tier_status_table(gates: Sequence[TierGate]) -> pd.DataFrame:
    return pd.DataFrame([g.summary_row() for g in gates])


def tier_outcome_table(gates: Sequence[TierGate]) -> pd.DataFrame:
    rows = []
    for g in gates:
        for t in TIERS:
            oc = g.outcomes.get(t)
            if oc is not None:
                rows.append(oc.as_row(g.trial, g.sample))
    return pd.DataFrame(rows)


def _threshold_provenance(cfg) -> List[List[str]]:
    """(threshold, value, config section, what it gates) rows for the header."""
    T = lambda k, sec="thresholds": threshold(cfg, k, default=None, section=sec)  # noqa: E731
    rows = [
        ["target_trpc2_min", T("target_trpc2_min"), "thresholds",
         "Rule 1 sort validation: target Trpc2 CPM must exceed this"],
        ["target_trpc2_concern", T("target_trpc2_concern"), "thresholds",
         "Rule 1: target Trpc2 below this is a failed sort"],
        ["nontarget_trpc2_max", T("nontarget_trpc2_max"), "thresholds",
         "Rule 1: GFP- nontarget Trpc2 CPM must stay below this"],
        ["failed_lib_actin_cpm_max", T("failed_lib_actin_cpm_max"), "thresholds",
         "Rule 4 failed-library gate (with Trpc2 below its own max)"],
        ["failed_lib_trpc2_cpm_max", T("failed_lib_trpc2_cpm_max"), "thresholds",
         "Rule 4 failed-library gate"],
        ["v1r_dominant_ratio_min", T("v1r_dominant_ratio_min"), "thresholds",
         "Rule 3 population call: Gnai2:Gnao1 above this = V1R-dominant"],
        ["v2r_dominant_ratio_max", T("v2r_dominant_ratio_max"), "thresholds",
         "Rule 3 population call: Gnai2:Gnao1 below this = V2R-dominant"],
        ["nontarget_total_vr_cpm_max", T("nontarget_total_vr_cpm_max"), "thresholds",
         "sort-purity ceiling on VR signal in the GFP- library"],
        ["cluster_max_gap_bp", T("cluster_max_gap_bp"), "thresholds",
         "genomic cluster definition (see cluster caveat below)"],
        ["tissue_panel_floor_cpm", T("tissue_panel_floor_cpm"), "thresholds",
         "tier 0: absolute CPM floor below which a marker panel is noise"],
        ["tissue_dominance_ratio", T("tissue_dominance_ratio"), "thresholds",
         "tier 0: factor by which one tissue panel must exceed the other"],
        ["min_assigned_counts", T("min_assigned_counts", "qc_thresholds"),
         "qc_thresholds",
         "library_status=DEGENERATE below this all-gene assigned-count total"],
        ["ratio_min_support_reads", T("ratio_min_support_reads", "qc_thresholds"),
         "qc_thresholds",
         "raw Gnao1 reads required before the ratio MAGNITUDE is quotable"],
        ["min_input_reads", T("min_input_reads", "qc_thresholds"), "qc_thresholds",
         "library viability floor"],
        ["min_uniquely_mapped_pct", T("min_uniquely_mapped_pct", "qc_thresholds"),
         "qc_thresholds", "technical QC: unique mapping rate floor"],
    ]
    return [[str(a), str(b), str(c), str(d)] for a, b, c, d in rows]


CLUSTER_CAVEAT = """\
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
expectation.\
"""


def report_header(cfg, trial: Optional[str], recs: Sequence[Dict[str, Any]],
                  vr_by_trial: Dict[str, Dict[str, Any]]) -> List[str]:
    """Provenance block: inputs, reference, CPM convention, thresholds, caveats."""
    work = cfg["work"]
    ref = cfg["reference"]
    qc_path = os.path.join(work, "results", "sample_qc_all.tsv")
    L: List[str] = []
    scope = trial if trial else "all trials"
    L += [
        f"# VNO receptor tiered analysis report -- {scope}",
        "",
        f"*Generated {_dt.datetime.now().strftime('%Y-%m-%d %H:%M')} by "
        f"`bin/vr_report.py`. Project: {cfg.get('project', 'NA')}.*",
        "",
        "## Provenance",
        "",
        "### Reference",
        "",
        f"* Genome key: **{ref['genome_key']}**",
        f"* GTF: `{ref['gtf']}`",
        f"* FASTA: `{ref['fasta']}`",
        "* Chromosome naming is Ensembl (`1`, `7`, `X`), not `chr1`. This GTF has no "
        "`gene` feature rows; gene spans in the VR reference were derived by "
        "aggregating exon/transcript rows per `gene_id`.",
        "* VR reference build: `bin/build_vr_reference.py` -> `ref/vr_gene_annotation.tsv`, "
        "`ref/vr_clusters.tsv`, `ref/vr_gene_to_cluster.tsv`; parse report in "
        "`ref/vr_gtf_parse_report.txt`.",
        "* 541 VR genes were annotated (319 V1R, 222 V2R); 538 are on the primary "
        "assembly. Filtering `is_primary_assembly == 1` reconciles exactly to the 318 "
        "Vmn1r / 220 Vmn2r rows present in the merged quantification table.",
        "* **Joins are on `gene_id`, never `gene_name`**: *Vmn1r-ps5* and *Vmn2r118* "
        "each map to two distinct gene_ids on different chromosomes.",
        "",
        "### Input tables",
        "",
        f"* QC verdicts: `{qc_path}`",
    ]
    for t in (trials_of(cfg) if not trial else [trial]):
        p = trial_paths(cfg, t)
        L.append(f"* {t} ({p['platform']}) quantification: `{p['gene_counts']}`")
        vr = vr_by_trial.get(t, {})
        for k, path in (vr.get("_present") or {}).items():
            L.append(f"  * VR table `{k}`: `{path}`")
        for k, path in (vr.get("_missing") or {}).items():
            L.append(f"  * VR table `{k}`: **NOT PRESENT** (looked for `{path}`) -- "
                     "tiers 3/4 report NO_DATA")
    note = read_provenance_note(qc_path)
    L += [
        "",
        "### CPM convention",
        "",
        "* CPM = gene count / (column sum over **all** genes in the **unscaled** "
        "`salmon.merged.gene_counts.tsv`) x 1e6.",
        "* The scaled and TPM tables give values 2.3-2.5x lower and must NOT be used "
        "for CPM in this project.",
    ]
    if note:
        L.append(f"* As recorded in the QC table header: {note}")
    L += [
        "",
        "### Thresholds in force",
        "",
        "| threshold | value | config section | gates |",
        "|---|---:|---|---|",
    ]
    for name, val, sec, desc in _threshold_provenance(cfg):
        L.append(f"| `{name}` | {val} | `{sec}` | {desc} |")
    L += [
        "",
        f"All values are read from `{cfg.get('_config_path')}`; none are hardcoded in "
        "the report module.",
        "",
        "### Cluster definition and its limits",
        "",
        CLUSTER_CAVEAT,
        "",
        "### Multi-mapping context",
        "",
        "* STAR retained 1.66-36.35% of input reads as multi-locus alignments versus "
        "0.02-0.30% discarded as too-many-loci (median ratio 125x). STAR is **not** "
        "filtering VR paralog ambiguity; it lands in Salmon's EM. Cluster-level "
        "aggregation is therefore mandatory, not optional.",
        "* **Even within-cluster splits in a MULTI-cell pool are expected biology, "
        "not automatically an artifact.** Monogenic receptor choice is a per-CELL "
        "rule. In a 100-cell pool, two neurons each expressing a different paralog "
        "of the same cluster produce a ~50/50 within-cluster split that is REAL "
        "co-expression across the pool. Only the unique-read channel separates that "
        "from EM redistribution: redistribution divides ONE transcript's reads, so "
        "unique (MAPQ-filtered, deduplicated) support collapses onto a single member, "
        "whereas genuine co-expression reproduces the split in unique reads too. The "
        "quantification track therefore gates its redistribution call on both a "
        "Monte-Carlo test against uniformity AND the requirement that only one "
        "co-dominant member clears the unique-read threshold. Neither this report nor "
        "the figures re-derive that judgement; both read `em_flag` from "
        "`vr_artifact_flags.tsv`. This matters for the planned stimulus-response "
        "experiment, which pools responsive cells and will make co-expression the "
        "common case.\n"
        "* Equal-fraction signature: paralogs within one cluster at near-equal "
        "fractions (e.g. 33/33/33) are a CANDIDATE redistribution artifact -- one "
        "expressed paralog whose reads the EM split. It is confirmed as an artifact "
        "only when the unique-read channel fails to reproduce the split (see the "
        "preceding caveat); in a single-cell or 2-cell library, where monogenic choice "
        "makes co-expression unlikely a priori, the artifact reading is the more "
        "probable one.",
        "* For calibration, Dietschi et al. 2022 report that only 57% of V1rD reads "
        "were unambiguously assignable under STAR `--outFilterMultimapNmax 4` with "
        "featureCounts `-M --fraction`.",
        "",
        "### Tier ladder",
        "",
        "A tier is reported only when **every** upstream tier passed for that sample. "
        "This is enforced by `TierGate.emit()`, which does not invoke a tier's content "
        "producer when an upstream tier failed and substitutes an explicit suppression "
        "line instead.",
        "",
    ]
    for t in TIERS:
        L.append(f"{t}. **{TIER_NAMES[t]}** -- {TIER_TITLES[t].split('-- ', 1)[-1]}")
    L.append("")
    return L


def _verdict_table(recs: Sequence[Dict[str, Any]]) -> List[str]:
    L = [
        "## Per-sample QC verdicts (read from the QC layer, not re-derived)",
        "",
        "| trial | sample | type | cells | QC | tissue | sort | population | library |"
        " suppress_biology | Trpc2 CPM | actin CPM |",
        "|---|---|---|---:|---|---|---|---|---|---|---:|---:|",
    ]
    for r in recs:
        tr = _f(r, "Trpc2_cpm")
        ac = _f(r, "actin_sum_cpm")
        L.append(
            "| {t} | `{s}` | {ct} | {n} | {q} | {ti} | {so} | {po} | {li} | {sb} | "
            "{tp} | {aa} |".format(
                t=_s(r, "trial"), s=_s(r, "sample"), ct=_s(r, "cell_type"),
                n=_s(r, "n_cells"), q=_s(r, "qc_overall"),
                ti=_s(r, "tissue_verdict"), so=_s(r, "sort_verdict"),
                po=_s(r, "population_call"), li=_s(r, "library_status"),
                sb=_s(r, "suppress_biology"),
                tp=f"{tr:.2f}" if tr is not None else "NA",
                aa=f"{ac:.1f}" if ac is not None else "NA",
            )
        )
    L.append("")
    return L


def _sample_section(gate: TierGate) -> List[str]:
    rec = gate.record
    top = gate.highest_tier_reported()
    L = [
        f"### `{gate.sample}` ({_s(rec, 'cell_type')}, {_s(rec, 'n_cells')} cells, "
        f"prep {_s(rec, 'prep_status')})",
        "",
        f"**Highest tier reported: {top if top >= 0 else 'none'}"
        + (f" ({TIER_NAMES[top]})" if top >= 0 else "") + ".** "
        + gate.stop_reason(),
        "",
    ]
    for t in TIERS:
        oc = gate.outcomes.get(t)
        if oc is None:
            continue
        L.append(f"**{TIER_TITLES[t]}** — `{oc.status}`")
        L.append("")
        if oc.status == SUPPRESSED:
            L.append(f"> SUPPRESSED. {oc.reason}")
        elif oc.status == NO_DATA:
            L.append(f"> NO DATA. {oc.reason}")
        elif oc.status == FAIL:
            L.append(f"> FAIL. {oc.reason}")
        else:
            L.append(f"{oc.reason}")
        if oc.lines:
            L.append("")
            L += oc.lines
        L.append("")
    return L


def write_trial_report(cfg, trial: str, gates: Sequence[TierGate],
                       recs: Sequence[Dict[str, Any]],
                       vr: Dict[str, Any], out_dir: Optional[str] = None) -> str:
    out_dir = out_dir or trial_paths(cfg, trial)["out_dir"]
    os.makedirs(out_dir, exist_ok=True)
    L = report_header(cfg, trial, recs, {trial: vr})
    L += _verdict_table(recs)
    L += ["## Tier outcomes by sample", ""]
    for g in gates:
        L += _sample_section(g)
    path = os.path.join(out_dir, "vr_report.md")
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")
    return path


def write_combined_report(cfg, all_gates: Dict[str, List[TierGate]],
                          all_recs: Dict[str, List[Dict[str, Any]]],
                          vr_by_trial: Dict[str, Dict[str, Any]],
                          out_dir: Optional[str] = None) -> str:
    out_dir = out_dir or os.path.join(cfg["work"], "results")
    os.makedirs(out_dir, exist_ok=True)
    flat_recs = [r for t in sorted(all_recs) for r in all_recs[t]]
    flat_gates = [g for t in sorted(all_gates) for g in all_gates[t]]
    L = report_header(cfg, None, flat_recs, vr_by_trial)

    # headline
    status = tier_status_table(flat_gates)
    reportable = status[status["highest_tier_reported"] >= 1]["sample"].tolist()
    L += [
        "## Headline",
        "",
        f"* {len(flat_gates)} libraries assessed; **{len(reportable)}** passed tier 1 "
        "(tissue identity + sort validation + library viability) and are the only "
        "libraries from which any biology is reported.",
        "",
        "**Trial 1 is the wrong tissue.** All four trial-1 libraries are main "
        "olfactory epithelium, not VNO: *Trpc2* is 0.00 CPM in every one while the "
        "main-olfactory panel is high (e.g. `pool2cellsRep3_S7`: *Adcy3* 2374, "
        "*Cnga2* 585, *Gnal* 4277, Olfr sum 1113 CPM). Both trials were quantified "
        "against the same annotation, which yields 17,918-37,230 *Trpc2* counts in "
        "trial-2 VNO libraries, so this is not a quantification artifact. *Omp* is "
        "tissue-shared -- OMP-Cre labels mature main-olfactory AND VNO neurons -- so "
        "it cannot be used as a VNO-specific marker, and high *Omp* is exactly what "
        "makes a main-olfactory library resemble a successful GFP+ VNO sort. No "
        "trial-1 VR biology is reportable; the remedy is wet-lab, not "
        "re-quantification.",
        "",
        "**A VR family signal can survive a library that fails every other check.** "
        "`target2cellsRep2_S4` carries a Vmn1r family sum of 501.5 CPM in a library "
        "whose actin sum is 0.88 CPM and whose *Trpc2* is 0.16 CPM. That is a failed "
        "library, and the VR signal in it must never be read as a receptor call. It is "
        "the single clearest reason this report gates structurally rather than by "
        "convention.",
        "",
    ]
    L += _verdict_table(flat_recs)
    L += ["## Tier reached, all samples", "",
          "| trial | sample | highest tier | stopped at | status | reason |",
          "|---|---|---|---|---|---|"]
    for g in flat_gates:
        row = g.summary_row()
        L.append("| {t} | `{s}` | {h} ({hn}) | {st} {sn} | {ss} | {rr} |".format(
            t=row["trial"], s=row["sample"], h=row["highest_tier_reported"],
            hn=row["highest_tier_name"], st=row["stopped_at_tier"],
            sn=row["stopped_at_name"], ss=row["stop_status"],
            rr=row["stop_reason"].replace("|", "/")[:240]))
    L.append("")
    for trial in sorted(all_gates):
        L += [f"## {trial} -- tier outcomes by sample", ""]
        for g in all_gates[trial]:
            L += _sample_section(g)
    path = os.path.join(out_dir, "vr_report_all.md")
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")
    return path


# --------------------------------------------------------------------------
# Self-test: synthetic records against the gate
# --------------------------------------------------------------------------


def _synthetic_records() -> Dict[str, Dict[str, Any]]:
    """Synthetic sample records exercising each suppression path."""
    return {
        "clean_target": dict(
            trial="tX", sample="clean_target", cell_type="target", n_cells="100",
            tissue_verdict="VNO", sort_verdict="PASS", population_call="V1R_dominant",
            library_status="OK", suppress_biology="False", blocking_flags="",
            qc_overall="USABLE", ratio_low_support="False",
            gnai2_gnao1_ratio_str="120.0:1"),
        "wrong_tissue": dict(
            trial="tX", sample="wrong_tissue", cell_type="target", n_cells="2",
            tissue_verdict="MOE", sort_verdict="FAIL_WRONG_TISSUE",
            population_call="V2R_dominant", library_status="OK",
            suppress_biology="True", blocking_flags="WRONG_TISSUE_MOE",
            qc_overall="UNUSABLE"),
        "failed_library": dict(
            trial="tX", sample="failed_library", cell_type="target", n_cells="2",
            tissue_verdict="VNO", sort_verdict="FAIL", population_call="V1R_dominant",
            library_status="FAILED", suppress_biology="True",
            library_reason="actin sum 0.9 CPM < 100 AND Trpc2 0.16 CPM < 10.",
            blocking_flags="library_FAILED", qc_overall="UNUSABLE"),
        "degenerate_library": dict(
            trial="tX", sample="degenerate_library", cell_type="target", n_cells="100",
            tissue_verdict="MOE", sort_verdict="FAIL_WRONG_TISSUE",
            population_call="undetermined", library_status="DEGENERATE",
            suppress_biology="True", qc_overall="UNUSABLE"),
        "nontarget_control": dict(
            trial="tX", sample="nontarget_control", cell_type="nontarget",
            n_cells="100", tissue_verdict="no_tissue_signal", sort_verdict="PASS",
            population_call="V1R_dominant", library_status="OK",
            suppress_biology="False", ratio_low_support="True",
            gnai2_gnao1_ratio_str="infinity", qc_overall="USABLE"),
        "target_no_tissue_signal": dict(
            trial="tX", sample="target_no_tissue_signal", cell_type="target",
            n_cells="2", tissue_verdict="no_tissue_signal", sort_verdict="FAIL",
            population_call="V1R_dominant", library_status="OK",
            suppress_biology="False", qc_overall="UNUSABLE"),
        "sort_fail_only": dict(
            trial="tX", sample="sort_fail_only", cell_type="target", n_cells="2",
            tissue_verdict="VNO", sort_verdict="FAIL",
            sort_reason="target Trpc2 0.37 CPM <= 100; no VNO-neuron signal",
            population_call="V1R_dominant", library_status="OK",
            suppress_biology="False", qc_overall="UNUSABLE"),
        "population_undetermined": dict(
            trial="tX", sample="population_undetermined", cell_type="target",
            n_cells="100", tissue_verdict="VNO", sort_verdict="PASS",
            population_call="undetermined", library_status="OK",
            suppress_biology="False", qc_overall="USABLE",
            population_note="both Gnai2 and Gnao1 are 0 CPM"),
        "unknown_verdicts": dict(
            trial="tX", sample="unknown_verdicts", cell_type="target", n_cells="2",
            tissue_verdict="", sort_verdict="", population_call="",
            library_status="", suppress_biology="", qc_overall=""),
    }


def selftest(verbose: bool = True) -> Dict[str, Any]:
    """
    Feed synthetic records to the gate and assert the suppression behaviour.
    Returns {'passed': n, 'failed': n, 'checks': [...]}. Exit code of
    `--selftest` is non-zero if any assertion fails.
    """
    recs = _synthetic_records()
    checks: List[Dict[str, str]] = []
    failures: List[str] = []

    def check(name: str, cond: bool, detail: str) -> None:
        checks.append({"check": name, "result": "PASS" if cond else "FAIL",
                       "detail": detail})
        if not cond:
            failures.append(f"{name}: {detail}")

    # A producer that records whether it was ever invoked. The central claim of
    # the design is that this is NOT called for a suppressed sample.
    class Spy:
        def __init__(self, tier: int):
            self.tier = tier
            self.calls = 0

        def __call__(self, rec: Dict[str, Any]) -> TierOutcome:
            self.calls += 1
            return TierOutcome(self.tier, TIER_NAMES[self.tier], PASS,
                               "synthetic VR payload", lines=["RECEPTOR CALL EMITTED"])

    # 1. clean target: every tier opens, producers run, tier 4 content present
    g = TierGate(recs["clean_target"])
    s3, s4 = Spy(3), Spy(4)
    g.emit(3, s3)
    g.emit(4, s4)
    check("clean_target reaches tier 4", g.highest_tier_reported() == 4,
          f"highest={g.highest_tier_reported()}")
    check("clean_target invokes both VR producers", s3.calls == 1 and s4.calls == 1,
          f"tier3 calls={s3.calls}, tier4 calls={s4.calls}")

    # 2. wrong tissue: tier 0 FAIL, tiers 1-4 SUPPRESSED, producers NEVER called
    g = TierGate(recs["wrong_tissue"])
    s3, s4 = Spy(3), Spy(4)
    g.emit(3, s3)
    g.emit(4, s4)
    check("wrong_tissue fails tier 0", g.outcomes[0].status == FAIL,
          g.outcomes[0].status)
    check("wrong_tissue suppresses tiers 1-4",
          all(g.outcomes[t].status == SUPPRESSED for t in (1, 2, 3, 4)),
          str({t: g.outcomes[t].status for t in TIERS}))
    check("wrong_tissue never invokes a VR producer", s3.calls == 0 and s4.calls == 0,
          f"tier3 calls={s3.calls}, tier4 calls={s4.calls}")
    check("wrong_tissue suppression names tier 0",
          all(g.outcomes[t].blocked_by == 0 for t in (1, 2, 3, 4)),
          str({t: g.outcomes[t].blocked_by for t in (1, 2, 3, 4)}))
    check("wrong_tissue reports no tier", g.highest_tier_reported() == -1,
          str(g.highest_tier_reported()))
    check("wrong_tissue reason states wrong tissue, not failed sort",
          "WRONG TISSUE" in g.outcomes[0].reason, g.outcomes[0].reason[:80])

    # 3. failed library: tier 0 passes (tissue is VNO) but tier 1 FAILs on
    #    library_status, so no receptor call can escape. This is the
    #    target2cellsRep2_S4 shape: real VR family signal in a dead library.
    g = TierGate(recs["failed_library"])
    s3, s4 = Spy(3), Spy(4)
    g.emit(3, s3)
    g.emit(4, s4)
    check("failed_library passes tier 0", g.outcomes[0].opens_downstream,
          g.outcomes[0].status)
    check("failed_library fails tier 1", g.outcomes[1].status == FAIL,
          g.outcomes[1].status)
    check("failed_library suppresses tiers 2-4",
          all(g.outcomes[t].status == SUPPRESSED for t in (2, 3, 4)),
          str({t: g.outcomes[t].status for t in (2, 3, 4)}))
    check("failed_library never invokes a VR producer",
          s3.calls == 0 and s4.calls == 0, f"{s3.calls},{s4.calls}")
    check("failed_library tier-1 reason cites library_status",
          "FAILED" in g.outcomes[1].reason, g.outcomes[1].reason[:80])

    # 4. degenerate library
    g = TierGate(recs["degenerate_library"])
    s4 = Spy(4)
    g.emit(3, Spy(3))
    g.emit(4, s4)
    check("degenerate_library emits no receptor call", s4.calls == 0, str(s4.calls))

    # 5. nontarget control: no_tissue_signal is EXPECTED, so tier 0 opens with a
    #    caveat and the sort-specificity check is reportable.
    g = TierGate(recs["nontarget_control"])
    check("nontarget no_tissue_signal opens tier 0 with caveat",
          g.outcomes[0].status == PASS_CAVEAT, g.outcomes[0].status)
    check("nontarget passes tier 1", g.outcomes[1].status == PASS,
          g.outcomes[1].status)
    check("nontarget low ratio support downgrades tier 2 to caveat",
          g.outcomes[2].status == PASS_CAVEAT, g.outcomes[2].status)
    check("nontarget tier 2 caveat says the magnitude is not quotable",
          "Poisson-unstable" in g.outcomes[2].reason, g.outcomes[2].reason[:80])

    # 6. target with no_tissue_signal: tissue UNCONFIRMED -> tier 0 fails closed
    g = TierGate(recs["target_no_tissue_signal"])
    check("target no_tissue_signal fails tier 0", g.outcomes[0].status == FAIL,
          g.outcomes[0].status)
    check("target no_tissue_signal is not called wrong-tissue",
          "not positive evidence" in g.outcomes[0].reason,
          g.outcomes[0].reason[:80])

    # 7. sort failure alone (tissue fine, library fine)
    g = TierGate(recs["sort_fail_only"])
    s3 = Spy(3)
    g.emit(3, s3)
    check("sort_fail_only stops at tier 1", g.highest_tier_reported() == 0,
          str(g.highest_tier_reported()))
    check("sort_fail_only suppresses tier 3 producer", s3.calls == 0, str(s3.calls))

    # 8. population undetermined: tiers 0-1 report, tier 2 fails, VR tiers suppressed
    g = TierGate(recs["population_undetermined"])
    s3 = Spy(3)
    g.emit(3, s3)
    check("population_undetermined reports through tier 1",
          g.highest_tier_reported() == 1, str(g.highest_tier_reported()))
    check("population_undetermined suppresses cluster tier", s3.calls == 0,
          str(s3.calls))
    check("population_undetermined suppression names tier 2",
          g.outcomes[3].blocked_by == 2, str(g.outcomes[3].blocked_by))

    # 9. unknown/empty verdicts fail closed rather than defaulting to pass
    g = TierGate(recs["unknown_verdicts"])
    s4 = Spy(4)
    g.emit(3, Spy(3))
    g.emit(4, s4)
    check("empty verdicts fail closed at tier 0", g.outcomes[0].status == FAIL,
          g.outcomes[0].status)
    check("empty verdicts emit no receptor call", s4.calls == 0, str(s4.calls))

    # 10. require() raises for imperative bypass attempts
    g = TierGate(recs["wrong_tissue"])
    raised = False
    try:
        g.require(4)
    except TierSuppressed as exc:
        raised = exc.blocked_by == 0
    check("require(4) raises TierSuppressed on a wrong-tissue sample", raised,
          "TierSuppressed raised with blocked_by=0" if raised else "no raise")
    ok = False
    gc = TierGate(recs["clean_target"])
    gc.emit(3, Spy(3))          # tier 3 must actually be decided first
    try:
        gc.require(4)
        ok = True
    except TierSuppressed:
        ok = False
    check("require(4) permits a clean sample once tier 3 has passed", ok, "no raise")

    # 11. a producer that itself raises TierSuppressed is converted, not propagated
    g = TierGate(recs["clean_target"])

    def naughty(rec):
        raise TierSuppressed(4, 1, "producer-side guard")

    oc = g.emit(4, naughty)
    check("producer-raised TierSuppressed becomes a SUPPRESSED outcome",
          oc.status == SUPPRESSED, oc.status)

    # 12. out-of-order emission: asking for tier 4 before tier 3 has been
    #     decided must NOT be permitted. Tier 3 is genuinely unevaluated, so
    #     "not yet decided" has to fail closed the same way a FAIL does --
    #     otherwise a caller could skip the cluster tier and go straight to a
    #     receptor call.
    g = TierGate(recs["clean_target"])
    s4 = Spy(4)
    oc = g.emit(4, s4)
    check("tier 4 emitted before tier 3 is suppressed", oc.status == SUPPRESSED,
          oc.status)
    check("out-of-order tier 4 does not invoke its producer", s4.calls == 0,
          str(s4.calls))
    check("out-of-order suppression names tier 3 as unevaluated",
          oc.blocked_by == 3 and "NOT_EVALUATED" in oc.reason, oc.reason[:90])
    raised_oo = False
    try:
        TierGate(recs["clean_target"]).require(4)
    except TierSuppressed as exc:
        raised_oo = exc.blocked_by == 3
    check("require(4) before tier 3 raises rather than KeyError", raised_oo,
          "TierSuppressed(blocked_by=3)" if raised_oo else "did not raise correctly")

    # 13. emit() cannot overwrite a sealed QC tier
    g = TierGate(recs["wrong_tissue"])
    oc = g.emit(0, lambda rec: TierOutcome(0, TIER_NAMES[0], PASS, "forged pass"))
    check("emit() cannot forge a pass on a sealed QC tier", oc.status == FAIL,
          oc.status)

    res = {"passed": sum(1 for c in checks if c["result"] == "PASS"),
           "failed": len(failures), "checks": checks, "failures": failures}
    if verbose:
        for c in checks:
            print(f"[{c['result']}] {c['check']} -- {c['detail']}")
        print(f"\n{res['passed']} passed, {res['failed']} failed")
    return res


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def run(cfg, trials: Optional[Sequence[str]] = None) -> Dict[str, Any]:
    trials = list(trials) if trials else trials_of(cfg)
    all_recs: Dict[str, List[Dict[str, Any]]] = {}
    all_gates: Dict[str, List[TierGate]] = {}
    vr_by_trial: Dict[str, Dict[str, Any]] = {}
    written: List[str] = []
    for t in trials:
        recs = load_qc_records(cfg, trial=t)
        vr = load_vr_tables(cfg, t)
        gates = [build_sample_report(cfg, r, vr) for r in recs]
        all_recs[t], all_gates[t], vr_by_trial[t] = recs, gates, vr
        written.append(write_trial_report(cfg, t, gates, recs, vr))
        out_dir = trial_paths(cfg, t)["out_dir"]
        st = tier_status_table(gates)
        st.to_csv(os.path.join(out_dir, "tier_status.tsv"), sep="\t", index=False)
        written.append(os.path.join(out_dir, "tier_status.tsv"))
        oc = tier_outcome_table(gates)
        oc.to_csv(os.path.join(out_dir, "tier_outcomes.tsv"), sep="\t", index=False)
        written.append(os.path.join(out_dir, "tier_outcomes.tsv"))
    written.append(write_combined_report(cfg, all_gates, all_recs, vr_by_trial))
    flat = [g for t in sorted(all_gates) for g in all_gates[t]]
    res_dir = os.path.join(cfg["work"], "results")
    tier_status_table(flat).to_csv(os.path.join(res_dir, "tier_status_all.tsv"),
                                   sep="\t", index=False)
    tier_outcome_table(flat).to_csv(os.path.join(res_dir, "tier_outcomes_all.tsv"),
                                    sep="\t", index=False)
    written += [os.path.join(res_dir, "tier_status_all.tsv"),
                os.path.join(res_dir, "tier_outcomes_all.tsv")]
    return {"written": written, "gates": all_gates, "records": all_recs,
            "vr": vr_by_trial}


def _main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--trial", action="append", default=None)
    ap.add_argument("--selftest", action="store_true",
                    help="run the synthetic-record gate tests and exit")
    args = ap.parse_args(argv)
    if args.selftest:
        return 0 if selftest()["failed"] == 0 else 1
    cfg = load_config(args.config)
    res = run(cfg, args.trial)
    for p in res["written"]:
        print(p)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
