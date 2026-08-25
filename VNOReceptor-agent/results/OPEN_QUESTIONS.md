# Open bioinformatics questions — status and recommendations

Generated 2026-08-20 14:50 from the tables under `results/`. Every recommendation below is grounded in a number from those tables; where a number does not exist yet, the entry says what would have to be computed rather than appealing to general principle.

---

## Q1. STAR multi-mapping diagnostic — **ANSWERED**

**Question.** Is VR paralog ambiguity being discarded by STAR's multi-mapping filters, or is it reaching the quantifier?

**Status: answered, and the answer determines the whole method.** It reaches the quantifier.

| library | retained multi-loci | discarded too-many-loci | ratio |
|---|---|---|---|
| pool100cells_S8 | 1.66% | 0.30% | 6× |
| pool2cellsRep1_S5 | 36.35% | 0.24% | 151× |
| pool2cellsRep2_S6 | 13.90% | 0.20% | 70× |
| pool2cellsRep3_S7 | 4.24% | 0.03% | 141× |
| nontarget100cells_S8 | 11.45% | 0.07% | 164× |
| target100cellsRep1_S6 | 10.92% | 0.12% | 91× |
| target100cellsRep2_S7 | 8.54% | 0.14% | 61× |
| target2cellsRep1_S3 | 7.57% | 0.07% | 108× |
| target2cellsRep2_S4 | 22.71% | 0.12% | 189× |
| target2cellsRep3_S5 | 18.28% | 0.02% | 914× |

Retained multiple-loci reads are **1.66–36.35%** of input; reads discarded as too-many-loci are **0.02–0.30%**. The retained channel dominates in 10/10 libraries, with a median ratio of ~125×.

**Conclusion.** STAR is *not* filtering VR paralog ambiguity out of the data. With `--outFilterMultimapNmax` at the nf-core default, a read matching a handful of 85–95%-identical paralogs is retained as a multi-mapper and handed to Salmon, where the EM distributes it. The ambiguity is therefore an *estimation* problem in the quantifier, not a *loss* problem in the aligner. Two consequences: (a) cluster-level aggregation is mandatory rather than merely prudent, and (b) any fix must act on the estimation step or on the underlying read information — tightening alignment filters would only convert a biased estimate into a missing one.

**Recommendation.** No further work. This question is closed; the diagnostic columns (`pct_multi_loci`, `pct_too_many_loci`, `dominant_multimap_channel`) are in `results/sample_qc_all.tsv` for every future run.

---

## Q2. Salmon selective alignment with a genome decoy

**Question.** Would re-quantifying with selective alignment against a decoy-aware index improve VR assignment?

**Status: scoped, expected payoff low for this specific problem.**

Decoys solve a particular failure: a read originating from an unannotated genomic or intronic region that finds a spurious best hit in the transcriptome. The decoy gives that read somewhere honest to go, suppressing false transcript assignment.

**That is not the failure mode here.** The ambiguity in this dataset is genuine sequence identity *between annotated paralogous transcripts*. A 75bp read from *Vmn1r166* matching *Vmn1r138* at 90% identity is not mis-assigned because it lacks a genomic home; it is ambiguous because two real transcripts explain it about equally well. A decoy cannot break that tie — both candidates remain in the index, and the EM still has to apportion the read.

**What a decoy *would* help with, quantified.** The relevant signal is the intergenic fraction, which is where genomic-origin reads show up. In the four cleared libraries it runs 2.7–5.0%. The one library with a genuinely high intergenic fraction is `pool100cells_S8` at 44.3% — and that is a trial-1 MOE library already excluded for tissue. So the population where a decoy would act is either small or already gated out.

**Expected payoff.** Low for paralog resolution; modest and worth having for general quantification hygiene if the data are ever reprocessed for another purpose. It would *not* change any cluster-level call in this report, and it would not resolve the `V1R_chr7_cl015` contradiction, because both competing genes there are annotated VR transcripts.

**Recommendation.** Do not re-run on this account alone. If the lab reprocesses for an unrelated reason, enable it then. Cost is one index build plus one quantification pass per sample; the honest expectation is that cluster-level shares move by less than the Monte-Carlo noise already reported.

---

## Q3. Sub-cluster aggregation by sequence identity

**Question.** Should clusters be defined by sequence identity rather than genomic proximity?

**Status: the proximity choice is already known to be arbitrary; the identity alternative is untested.**

Three facts from the reference layer bear on this. (1) **200kb is not a natural break**: the V1R inter-gene gap distribution has its KDE minimum near 2Mb, an order of magnitude away, so the threshold is a conservative convention rather than a discovered boundary. (2) The convention has a visible consequence — a 217,366bp gap, 17kb over the rule, splits the chr7 region so that the same signal reads as 99.91% of one supercluster or 66.8%/33.1% of two clusters. (3) V2R clustering is fragmented: 18 of 37 clusters are singletons and only 180/222 genes sit in clusters of ≥5, so proximity aggregation protects V2R much less than V1R.

**Why sequence identity is the better-motivated grouping.** The thing being defended against is reads that cannot be assigned between two sequences. That is a property of the *sequences*, and genomic proximity is only a proxy for it — a good proxy, since local duplication produces both adjacency and similarity, but a proxy that can fail in both directions: adjacent-but-divergent paralogs get pooled unnecessarily, and distant-but-similar ones stay split.

**Concrete evidence that the proxy is failing here.** In `V1R_chr7_cl015`/`V1R_chr7_cl016` the two evidence channels disagree across the cluster boundary: *Vmn1r131* (cl015) holds 4,452 EM counts with zero unique reads, while *Vmn1r166* (cl016) holds 4,489 EM counts with only 51 unique reads. A grouping that put the actual sequence-similar set in one bin might make this region interpretable as a single unit instead of two partly-contradictory ones.

**What would have to be computed.** Three steps, all with tools already on the cluster: (1) extract VR CDS/transcript sequences for the 538 primary-assembly genes from the GRCm38 FASTA (`seqkit`, minutes); (2) all-pairs identity within family — `blastn` 538×538 is trivial at this scale — and, better matched to the actual failure, a *k*-mer-sharing matrix at k = 31 and at the read length, since mappability is about shared *k*-mers, not global identity; (3) single-linkage or spectral clustering on that matrix, then re-derive cluster-level tables and compare the resulting calls against the current ones. Step (3) is the actual test: does the grouping change any call, or only the labels?

**Expected payoff.** Moderate and mostly *diagnostic*. With three usable GFP+ libraries and two clusters called in the 2-cell pool, a regrouping will not manufacture new biology. Its value is that it would tell you whether the current cluster-level calls are robust to the boundary convention — which is currently unknown and is the single largest unexamined assumption in the pipeline.

**Recommendation.** Worth doing, and cheap (a day of compute at most). Do it as a *sensitivity analysis* on the existing calls rather than as a replacement clustering: keep the genomic tiers as the reported ones, add an identity-based tier, and report where the three disagree. Prioritise it above Q2 and Q4.

---

## Q4. De novo assembly on dominant-cluster reads

**Question.** Could assembling the reads from a dominant cluster reconstruct the expressed transcript directly?

**Status: feasible only for the deepest cases, and it inherits the same ambiguity.**

The read support that would be available, per candidate:

| library | cluster | gene | unique MAPQ255 (all / dedup) | EM counts |
|---|---|---|---|---|
| target100cellsRep1_S6 | V1R_chr6_cl008 | *Vmn1r35* | 33,670 / 15,329 | 17,097 |
| target100cellsRep1_S6 | V1R_chr7_cl013 | *Vmn1r89* | 33,139 / 18,750 | 14,767 |
| target100cellsRep1_S6 | V1R_chr7_cl019 | *Vmn1r185* | 29,712 / 17,083 | 13,050 |
| target100cellsRep1_S6 | V1R_chr7_cl018 | *Vmn1r184* | 24,549 / 13,560 | 10,847 |
| target100cellsRep1_S6 | V1R_chr7_cl017 | *Vmn1r178* | 25,992 / 15,908 | 7,322 |
| target100cellsRep2_S7 | V1R_chr6_cl008 | *Vmn1r37* | 53,841 / 17,159 | 41,188 |
| target100cellsRep2_S7 | V1R_chr7_cl013 | *Vmn1r89* | 35,153 / 18,468 | 15,971 |
| target100cellsRep2_S7 | V1R_chr17_cl021 | *Vmn1r236* | 16,723 / 10,646 | 5,717 |
| target2cellsRep1_S3 | V1R_chr7_cl016 | *Vmn1r166* | 70 / 51 | 4,489 |
| target2cellsRep1_S3 | V1R_chr7_cl015 | *Vmn1r103* | 186 / 92 | 0 |

**Depth is adequate in the 100-cell pools.** The top candidates there carry 16,723–53,841 unique reads, and a V1R coding sequence is roughly ~900bp (e.g. *Vmn1r91* spans 923bp of genome). Even at 30% duplication that is coverage in the hundreds — assembly is not depth-limited there.

**Depth is marginal in the 2-cell pool, which is the case that matters most.** `target2cellsRep1_S3`'s flagged cluster has 51 deduplicated unique read-mates for *Vmn1r166* — at ~76bp per mate that is roughly 4× nominal coverage of a ~900bp CDS if they were evenly spread, and they are not. An assembly from ambiguous reads plus 51 anchors is not going to produce a confident contig.

**The deeper problem is that assembly does not escape the ambiguity.** If you assemble the *cluster's* reads, the ambiguous majority is exactly the input, and an assembler faced with two 90%-identical templates will either collapse them into one consensus (losing the distinction you wanted) or fragment at every divergent site. If you assemble only the unambiguous reads, you have thrown away 98.9% of the signal in the case you care about. The information content of the read set is what it is; assembly re-arranges it rather than adding to it.

**Expected payoff.** Low for the 2-cell pools, where the question is live. Possibly useful in the 100-cell pools as a *confirmation* device: an assembled contig that matches one paralog's divergent sites and not the other's would independently corroborate a call — but in those pools the unique-read evidence is already strong (16k–54k reads), so it would confirm what is not in doubt.

**Recommendation.** Deprioritise. If tried, spend the effort on `target100cellsRep2_S7`'s `V1R_chr7_cl013` — where *Vmn1r89* and *Vmn1r87* both carry ~32–35k unique reads and the interesting question is whether two transcripts are genuinely present — rather than on the 2-cell sample. Longer reads (Q6) are the better route to the 2-cell case.

---

## Q5. Variant calling within the dominant cluster

**Question.** Could SNVs separate paralogs that reads cannot be assigned between?

**Status: the right idea, and there is a specific test case for it.**

This is the most promising of the four open computational routes, because it attacks the actual problem: paralogs differ at *specific positions*, and a read overlapping such a position is diagnostic even when the rest of it is ambiguous. Unlike assembly, this extracts information that the EM discards.

**The constraint is geometric.** Per-mate reads are ~75–76bp (STAR's 143–149 "average input read length" is the sum of both mates; see Q6), and the observed mismatch rate in the cleared libraries is 0.21–0.23%, i.e. sequencing error is low enough that a genuine paralog-diagnostic mismatch is distinguishable from noise. At 85–95% identity, paralogs differ every ~10–20bp on average, so even a 75bp read should span several diagnostic sites — the information is present in the reads today. What is missing is a pipeline step that uses it. (At 2×150bp it would span proportionally more, and crucially would link them, which is why Q5 and Q6 compound.)

**Concrete test case: `V1R_chr7_cl015` in `target2cellsRep1_S3`.** This is the cluster where the two evidence channels contradict each other:

- *Vmn1r131*: 4,452 EM counts (99.9% of the cluster), **0** unique MAPQ255 reads.
- *Vmn1r103* / *Vmn1r104*: 92 deduplicated unique reads each, **0** EM counts.

A variant-based approach makes a falsifiable prediction here. Pile up *all* cluster-assigned reads against each candidate paralog's sequence and genotype the diagnostic sites: if the expressed transcript is *Vmn1r131*, the pileup should match *Vmn1r131*'s alleles at sites where it differs from *Vmn1r103*/*104*, and the ~4,452 EM-assigned reads become real evidence. If instead the pileup matches *Vmn1r103*/*104*, then the EM assignment is wrong and the 92 unique reads were the honest signal. Either outcome resolves a contradiction the current pipeline can only report.

**What it would take.** The BAMs already exist and `bcftools` 1.20 and `samtools` 1.19 are on the cluster. Steps: (1) build a paralog-diagnostic site catalogue per cluster by aligning member CDS sequences to each other (this is the same sequence work as Q3, so the two share infrastructure); (2) pile up cluster reads including multi-mappers — deliberately, since multi-mappers are the reads carrying the information — at those sites; (3) compute a per-paralog likelihood from the allele counts and report it as a third evidence channel alongside EM and unique reads. Roughly a week's work. The main methodological care needed is that MAPQ filtering must be *disabled* for this pileup, which is the opposite of normal practice and would need a guard so it never leaks into the standard path.

**Expected payoff.** High relative to the others, and it is the only route that could upgrade a call from `tentative_unconfirmed` to confirmed using existing data. Honest caveat: in the 2-cell case it operates on ~4,500 cluster reads, so it will produce a likelihood ratio rather than a certainty, and a diagnostic site falling in a region where both paralogs are identical contributes nothing.

**Recommendation.** Highest-priority computational work of the four. Do Q3's sequence comparison first (shared infrastructure), then Q5. Treat `V1R_chr7_cl015` as the acceptance test: a method that cannot adjudicate that cluster is not yet working.

---

## Q6. The planned 150bp PE upgrade — what it does and does not buy

**Read the length field carefully.** STAR reports "average input read length" as the sum of BOTH mates for paired-end data, so the 143–149 in the table below is the fragment's sequenced total, not the per-read length. Measured directly from the aligned BAMs, the per-mate length is **75–76bp** (modal 76bp, `samtools view -f 64/-f 128`), confirming the briefing's 2×75bp. The planned upgrade to 2×150bp is therefore a genuine doubling of per-read length, not a marginal change.

| library | STAR avg input read length (both mates) | avg mapped length | per-mate length |
|---|---|---|---|
| pool100cells_S8 | 143bp | 123.8bp | ~72bp |
| pool2cellsRep1_S5 | 144bp | 141.4bp | ~72bp |
| pool2cellsRep2_S6 | 145bp | 140.5bp | ~72bp |
| pool2cellsRep3_S7 | 149bp | 148.0bp | ~74bp |
| nontarget100cells_S8 | 149bp | 148.6bp | ~74bp |
| target100cellsRep1_S6 | 149bp | 149.2bp | ~74bp |
| target100cellsRep2_S7 | 149bp | 149.4bp | ~74bp |
| target2cellsRep1_S3 | 149bp | 148.8bp | ~74bp |
| target2cellsRep2_S4 | 146bp | 145.6bp | ~73bp |
| target2cellsRep3_S5 | 146bp | 144.0bp | ~73bp |

The per-mate column is simply the STAR field halved, and it lands slightly below the 76bp nominal because STAR's figure is a mean over reads that have already been adapter- and quality-trimmed — the trimmed fraction runs 1.4–5.8% across these libraries. The direct BAM measurement resolves the distinction: the modal SEQ length is 76bp with a 74–76bp spread, i.e. untrimmed reads sit at the nominal length and the sub-76 average reflects trimming, not a shorter protocol.

**What longer reads improve.**

1. *Unique-read support per paralog* — the quantity every individual call in this report rests on. A longer read spans more divergent sites, so a larger fraction of reads become uniquely placeable. This is the direct fix for the case that is currently weakest: *Vmn1r166*'s attribution rests on 51 unique reads out of 4,489 apparent counts (~1.1%). Raising that fraction is exactly what changes `tentative_unconfirmed` to confirmed.
2. *Ability to span two or more diagnostic sites in one read*, which is what makes Q5's variant approach strong rather than marginal: a read carrying two linked diagnostic alleles is far more informative than two reads each carrying one.
3. *Mate-pair span*, which helps when one mate lands in a conserved stretch and the other in a divergent one.

**What longer reads do not fix.**

1. **Regions where paralogs are identical.** Where two VR genes share an exact stretch longer than the read, no read length below that stretch helps. This is a property of the duplication history, not the assay.
2. **Sample count.** Only 3 usable GFP+ libraries exist. Longer reads on the same 3 libraries give better-resolved calls for the same three cells' worth of choices; they do not add biological replication. Given that `Rep1`/`Rep2`/`Rep3` are independent libraries from different cells, more libraries — not longer reads — is what buys statistical statements about receptor-choice frequency.
3. **The trial-1 tissue problem.** Nothing about read length touches it.
4. **Library failure.** `target2cellsRep2_S4` failed at 0.88 CPM actin; sequencing it longer would produce longer reads from the same failed prep.

**Recommendation.** Proceed with 2×150bp — it is a real doubling of per-mate length from the current 75–76bp, and it acts directly on the quantity every individual call in this report is limited by. Expect it to raise unique-read support substantially and to make Q5's variant approach materially stronger (a 150bp read spanning two or three diagnostic sites carries linkage information that two 75bp reads do not).

Two caveats on allocation, though. First, longer reads will not rescue a region where paralogs are identical over a stretch longer than the read, and they do nothing for the trial-1 tissue problem or for a failed prep. Second, with only 3 usable GFP+ libraries, read length and library count buy different things: longer reads give better-resolved calls for the same three cells' worth of choices, while more libraries are what would support any statement about receptor-choice frequency across the population. If the budget allows only one, the ordering depends on the question — **resolve which receptor these cells chose** favours 150bp; **characterise how choice is distributed** favours more libraries. Both are downstream of fixing the tissue sampling.

---

## Priority summary

| # | question | status | payoff | recommendation |
|---|---|---|---|---|
| Q1 | STAR multi-mapping diagnostic | **answered** | — | closed; ambiguity reaches the EM, cluster aggregation mandatory |
| Q5 | variant calling in cluster | scoped | **high** | **do first** (after Q3's sequence work); acceptance test = `V1R_chr7_cl015` |
| Q3 | sequence-identity sub-clustering | scoped | moderate, diagnostic | do as a sensitivity analysis; shares infrastructure with Q5 |
| Q2 | Salmon decoy selective alignment | scoped | **low** for paralogs | skip unless reprocessing anyway |
| Q4 | de novo assembly | scoped | low | deprioritise; longer reads are the better route |
| Q6 | 150bp PE upgrade | planned | small (already ~149bp) | proceed, but more libraries beats longer reads |

**Above all of these:** the trial-1 tissue finding and the 3-usable-GFP+-library count are wet-lab constraints, and no computational work in this table substitutes for fixing them.

