# Provenance notes for the cold-start validation run

## results/vr_cluster_map.png has no producer in the pipeline (ORPHAN)

The cold-start re-run reproduced 31 of 36 backed-up outputs byte-identically,
3 as timestamp-only differences, and flagged exactly one MISSING_NEW path:
`vr_cluster_map.png`.

Cause: no module under `bin/` writes that filename. `grep -rn cluster_map bin/`
matches only `_supercluster_map()` in `vr_report.py`, an unrelated dict helper.
`build_vr_reference.py` contains no plotting code at all. The file was therefore
produced out-of-band during phase 0 -- an exploratory figure that was never
wired into a module -- and it is not regenerable from `config/project.yaml`
plus the data.

This is a real finding of the validation exercise, not a run failure: it means
the committed `results/` tree contained one artifact that the pipeline could not
reproduce. Resolution options, in order of preference:

1. If the figure is wanted as a deliverable, add a producer to `bin/vr_figures.py`
   so it is regenerated like the other five figures.
2. If it was exploratory, leave it out of `results/` (the backup at
   `results.prerun_backup/vr_cluster_map.png` retains it).

Nothing else in the tree is unaccounted for. The three TIMESTAMP_ONLY paths
(`vr_report.md` x2, `vr_report_all.md`) differ only in their generation-timestamp
lines; their computed content is identical. `refcheck_diff.txt` is EXTRA_NEW
because the reference-reproducibility check is new in this run.
