# QC / tissue-identity / sort-validation / population-ID modules

Track: library viability, tissue identity, sort validation, population
identification. These modules GATE everything downstream: no VR cluster or
receptor call may be reported for a sample whose `suppress_biology` is True.

## Modules (bin/)
| file | role |
|---|---|
| `vr_config.py`      | loads `config/project.yaml`; all paths, sample metadata and thresholds come from here |
| `vr_qc.py`          | STAR / qualimap / MultiQC parsing -> `technical_qc.tsv` |
| `vr_markers.py`     | CPM + tissue identity, Rule 1 sort validation, Rule 3 actin sanity, Rule 4 population call -> `marker_cpm.tsv` |
| `vr_sample_qc.py`   | joins the two, derives `blocking_flags` -> `sample_qc.tsv`, `sample_qc_summary.txt` |
| `vr_qc_figures.py`  | the three deliverable figures |

## Run
    module load python/3.11.4
    python bin/vr_sample_qc.py --all-trials
    python bin/vr_qc_figures.py --all-trials --fig-dir results/figures

## Priority hierarchy
    (0) TISSUE IDENTITY  <- gates everything below it
    (1) sort validation (Rule 1, Trpc2)
    (2) population ID via Gnai2:Gnao1 (Rule 4)
    (3) cluster-level VR calls        <- the reliable tier
    (4) individual receptor ID        <- always flagged tentative

## (0) Tissue identity — why it comes first

**OMP-Cre drives GFP in mature MAIN OLFACTORY (MOE) neurons as well as VNO
neurons.** `Omp` is therefore a TISSUE-SHARED marker: it marks mature olfactory
neurons of *both* epithelia and cannot discriminate them. A GFP+ sort taken from
the wrong tissue — or from VNO contaminated by MOE — yields Omp-high,
Trpc2-zero cells that superficially look like a successful sort.

Reading a high `Omp` as evidence of a VNO population is a real error mode, and
it is the reason this check runs before Rule 1: asking "did the sort enrich VNO
neurons?" is only meaningful once the tissue is VNO.

    MOE panel           Olfr* summed, Adcy3, Cnga2, Gnal
    VNO-specific panel  Trpc2, Vmn1r* summed, Vmn2r* summed
    SHARED (excluded)   Omp
    NOT tissue markers  Gnai2 / Gnao1 -- these split V1R vs V2R WITHIN the VNO
                        and are broadly expressed elsewhere (Gnai2 reads 165 CPM
                        in a GFP- library with no VR signal at all)

Panels are compared on their **maximum** member, not their sum, so one strongly
expressed marker establishes the tissue without dilution by silent family
members. An **absolute CPM floor** (`tissue_panel_floor_cpm`, 100) is required
rather than a bare ratio: below the floor a panel is noise, and a ratio of two
noise values carries no information.

| verdict | rule |
|---|---|
| `VNO` | VNO panel >= floor, MOE panel < floor |
| `MOE` | MOE panel >= floor, VNO panel < floor |
| `VNO_dominant_mixed` | both >= floor, VNO > 3x MOE |
| `MOE_dominant_mixed` | both >= floor, MOE > 3x VNO |
| `ambiguous_mixed` | both >= floor, within 3x |
| `no_tissue_signal` | both < floor — no tissue information either way |

## Downstream contract
Read `results/sample_qc_all.tsv` and honour these columns:

* `suppress_biology` (bool) — **hard gate**. True => report NO biology for that
  sample, whatever VR reads appear in it. Set by a FAILED/DEGENERATE library
  *or* by a wrong-tissue verdict.
* `tissue_verdict` — see the table above. Check this BEFORE `sort_verdict`.
* `qc_overall` — `USABLE` | `USE_WITH_CAUTION` | `UNUSABLE`.
* `blocking_flags` — `;`-separated; the literal string `none` when clean.
  Entries prefixed `warn_` are advisory; all others are gates.
* `sort_verdict` — `PASS` | `CONCERN` | `FAIL` | `FAIL_WRONG_TISSUE` | `NA`.
* `library_status` — `OK` | `SUSPECT` | `FAILED` | `DEGENERATE`.
* `population_call` — `V1R_dominant` | `V2R_dominant` | `mixed` | `undetermined`.
  Check `ratio_low_support` before quoting a ratio magnitude.

## Four findings that shape how these verdicts must be read

### 1. Trial 1 is MAIN OLFACTORY EPITHELIUM, not VNO

Three of four trial-1 libraries are MOE; the fourth carries no tissue signal at
all. `pool2cellsRep3_S7` — which the project brief calls "the cleanest 2-cell
sample" — is Adcy3 2374 / Cnga2 585 / Gnal 4277 / Olfr 1113 CPM with Trpc2 = 0:
a textbook mature MOE neuron. Its Omp of 4028 CPM is a real mature-neuron
signal, but it is *MOE* Omp. `pool100cells_S8` puts 48,729 CPM — about 4.9% of
its entire library — into olfactory receptors with zero VR signal.

Trial 2's targets are the mirror image: Trpc2 1029-1476 CPM, Vmn1r sum
772-6036 CPM, MOE panel low.

**Trial-1 VR biology is not reportable, and the remedy is wet-lab (dissection /
sort gate), NOT bioinformatic.** Re-quantifying against a different annotation
will not change it:

* Both trials used the SAME genome key (`mouse-ensembl-grcm38-r91`), the same
  aligner (star_salmon) and the same trimmer.
* Both `tx2gene.tsv` files contain the same single Trpc2 transcript
  (`ENSMUST00000124189` -> `ENSMUSG00000100254`).
* On that identical annotation, trial 2 assigns Trpc2 **17,918 / 19,846 /
  37,230** raw counts in `target2cellsRep1_S3` / `target100cellsRep1_S6` /
  `target100cellsRep2_S7`.

An annotation yielding ~37k counts in one run cannot "drop out" in another run
using the same annotation. Trial 1 vs trial 2 is a TISSUE difference.

### 2. Multi-mapping loss is dominated by the RETAINED channel

Reads mapped to multiple loci run 1.66-36.35% of input; reads discarded as
too-many-loci run 0.02-0.30% — a median 125x difference. VR paralog ambiguity is
therefore not filtered out at alignment; it lands in Salmon's EM. Cluster-level
aggregation is mandatory, exactly as the project brief specifies.

### 3. A large CPM can come from a near-empty library

CPM is a ratio. `pool100cells_S8` has only 60,744 counts assigned across all
genes, so its 39,881 CPM actin sum is "all of a tiny pie" and it slipped past
the actin arm of the failed-library gate. The `min_assigned_counts` gate catches
this as `library_status = DEGENERATE` and suppresses biology. Always read
`library_total_counts` alongside any CPM.

### 4. Ratio magnitudes need read support

The Gnai2:Gnao1 ratio is denominator-free, so normalisation cannot shift it —
but at small read counts Poisson noise can. `target100cellsRep1_S6` rests on 7
raw Gnao1 reads, where a one-read shift moves the ratio ~15% (422.7:1 here vs
385:1 in an earlier hand check — the same measurement). The `ratio_low_support`
flag marks these; the direction of the call is robust, the magnitude is not.

## Sort purity note
`target100cellsRep1_S6` and `target100cellsRep2_S7` are `VNO_dominant_mixed`:
they carry a real minor MOE component (Cnga2 180/284, Olfr 177/69 CPM) on top of
a dominant VNO signal. This is reported as sort-purity information, not
suppressed — both remain usable.

## Thresholds
Rule 1/3/4 cutoffs and the tissue keys live in `markers:` / `thresholds:`. The
technical-QC cutoffs this track added live in `qc_thresholds:` (append-only).
