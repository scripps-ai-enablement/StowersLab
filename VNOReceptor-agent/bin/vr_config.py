#!/usr/bin/env python3
"""
vr_config.py -- shared config/path contract loader for the vr_analysis pipeline.

Every module in bin/ reads paths, sample metadata and numeric thresholds through
this loader. Nothing downstream hardcodes a sample name, a directory or a cutoff.

Usage
-----
    from vr_config import load_config, trial_paths, samples_of, threshold
    cfg = load_config()                      # finds ../config/project.yaml
    p   = trial_paths(cfg, "trial2")         # dict of resolved absolute paths
    for s, meta in samples_of(cfg, "trial2").items(): ...

Notes
-----
* The on-disk results directories are spelled "trail1"/"trail2" (sic). That
  typo lives only in config/project.yaml (trials.<t>.results); no module should
  ever construct it.
* Thresholds are looked up by dotted key with an explicit default so that a
  config missing a newer key degrades loudly (the default is reported) rather
  than silently.
"""

from __future__ import annotations

import os
from typing import Any, Dict, Iterable, Mapping, Optional

import yaml

__all__ = [
    "DEFAULT_CONFIG_PATH",
    "load_config",
    "trial_paths",
    "samples_of",
    "trials_of",
    "threshold",
    "marker_genes",
    "TECH_THRESHOLD_DEFAULTS",
]

# Canonical location of the work folder's config. Overridable via --config or
# the VR_ANALYSIS_CONFIG environment variable.
def _resolve_work_dir():
    """Locate the work folder without assuming any particular site.

    Resolution order, first hit wins:
      1. $VR_WORK              -- explicit override, works anywhere
      2. <this file>/..        -- the normal case: bin/ inside the work folder
      3. $PWD                  -- running from the work folder root
    No site-specific path is ever used as a fallback; if none of the above
    contains config/project.yaml the caller must pass --config explicitly.
    """
    env = os.environ.get("VR_WORK")
    if env:
        return os.path.abspath(os.path.expanduser(env))
    here = os.path.dirname(os.path.abspath(__file__))
    parent = os.path.dirname(here)
    if os.path.isfile(os.path.join(parent, "config", "project.yaml")):
        return parent
    if os.path.isfile(os.path.join(os.getcwd(), "config", "project.yaml")):
        return os.getcwd()
    return parent


_WORK_DEFAULT = _resolve_work_dir()
DEFAULT_CONFIG_PATH = os.path.join(_WORK_DEFAULT, "config", "project.yaml")

# Technical-QC thresholds owned by the QC module. These are *defaults*; if the
# config carries a `qc_thresholds:` section its values win. Kept here so the
# module is runnable against a config that predates the section.
TECH_THRESHOLD_DEFAULTS: Dict[str, float] = {
    # library viability
    "min_input_reads": 1_000_000,
    "min_uniquely_mapped_pct": 50.0,
    # alignment-quality warnings
    "max_multi_loci_pct": 40.0,
    "max_too_many_loci_pct": 5.0,
    "max_unmapped_too_short_pct": 40.0,
    "max_mismatch_rate_pct": 2.0,
    # coverage evenness: qualimap 5'-3' bias, healthy ~= 1
    "bias_5p3p_low": 0.5,
    "bias_5p3p_high": 2.0,
    # genomic origin
    "min_exonic_pct": 40.0,
    "max_intergenic_pct": 30.0,
}


def _resolve_config_path(path: Optional[str] = None) -> str:
    if path:
        return os.path.abspath(path)
    env = os.environ.get("VR_ANALYSIS_CONFIG")
    if env:
        return os.path.abspath(env)
    # ../config/project.yaml relative to this file (works inside the work folder)
    here = os.path.dirname(os.path.abspath(__file__))
    local = os.path.join(os.path.dirname(here), "config", "project.yaml")
    if os.path.exists(local):
        return local
    return DEFAULT_CONFIG_PATH


def load_config(path: Optional[str] = None) -> Dict[str, Any]:
    """Load project.yaml and stash the resolved path under `_config_path`."""
    cfgpath = _resolve_config_path(path)
    with open(cfgpath) as fh:
        cfg = yaml.safe_load(fh)
    if not isinstance(cfg, dict):
        raise ValueError(f"{cfgpath}: expected a YAML mapping at top level")
    cfg["_config_path"] = cfgpath
    for req in ("trials", "paths", "samples", "markers", "thresholds"):
        if req not in cfg:
            raise KeyError(f"{cfgpath}: missing required top-level key '{req}'")
    return cfg


def trials_of(cfg: Mapping[str, Any]) -> list:
    """Trial names in a stable order (trial1 before trial2)."""
    return sorted(cfg["trials"].keys())


def trial_paths(cfg: Mapping[str, Any], trial: str) -> Dict[str, str]:
    """Resolve every relative entry in `paths:` against this trial's results dir."""
    if trial not in cfg["trials"]:
        raise KeyError(f"unknown trial {trial!r}; known: {sorted(cfg['trials'])}")
    tinfo = cfg["trials"][trial]
    root = tinfo["results"]
    out = {"trial": trial, "results_root": root, "platform": tinfo.get("platform", "NA")}
    for key, rel in cfg["paths"].items():
        out[key] = os.path.join(root, rel)
    out["work"] = cfg.get("work", _WORK_DEFAULT)
    out["out_dir"] = os.path.join(out["work"], "results", trial)
    return out


def samples_of(cfg: Mapping[str, Any], trial: str) -> Dict[str, Dict[str, Any]]:
    """Sample -> metadata dict (cell_type, n_cells, prep_status) for one trial."""
    smap = cfg["samples"].get(trial)
    if not smap:
        raise KeyError(f"no samples configured for trial {trial!r}")
    out: Dict[str, Dict[str, Any]] = {}
    for name, meta in smap.items():
        m = dict(meta or {})
        m.setdefault("cell_type", "unknown")
        m.setdefault("n_cells", None)
        m.setdefault("prep_status", "unknown")
        m["trial"] = trial
        m["sample"] = name
        out[name] = m
    return out


def threshold(cfg: Mapping[str, Any], key: str, default: Any = None,
              section: str = "thresholds") -> Any:
    """
    Fetch a numeric threshold. `section` selects the config block; for keys in
    the QC module's own `qc_thresholds:` block the TECH_THRESHOLD_DEFAULTS table
    supplies the fallback so the module runs against an older config.
    """
    block = cfg.get(section) or {}
    if key in block:
        return block[key]
    if section == "qc_thresholds" and key in TECH_THRESHOLD_DEFAULTS:
        return TECH_THRESHOLD_DEFAULTS[key]
    if default is not None:
        return default
    raise KeyError(f"threshold {section}.{key} not in config and no default given")


def marker_genes(cfg: Mapping[str, Any], group: Optional[str] = None) -> Any:
    """`markers` block, or one named list from it."""
    m = cfg["markers"]
    if group is None:
        return m
    if group not in m:
        raise KeyError(f"marker group {group!r} not in config; known: {sorted(m)}")
    val = m[group]
    return list(val) if isinstance(val, Iterable) and not isinstance(val, str) else [val]
