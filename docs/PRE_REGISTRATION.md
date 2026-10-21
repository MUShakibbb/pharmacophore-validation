# Pre-registration

Committed to git **before any fit score existed**. The commit order in this repository is the evidence:
the pre-registration commit precedes the results commit. That is checkable, which is the entire point.

Target: **p38α MAPK** (CHEMBL260), ATP site, type-I (hinge-binding) binding mode.

---

## Primary metric

**Güner–Henry (GH).** Chosen over MCC because GH penalises low specificity more heavily, and a
high-recall / low-specificity filter saturates the screening library and is operationally useless.

**Secondary:** MCC, ROC-AUC, EF1%, BEDROC.

**GH is reported at 1:50 only.** At 1:50, GH ≥ 0.5 would require roughly recall 0.6 at FPR 1% —
an enrichment factor around 27, beyond what any pharmacophore achieves. The pass bar is therefore
expressed as AUC and EF1%, not GH.

> **Post-hoc note (recorded, not changed).** The decoy pool was sized at 1:50 against an earlier split
> with 1,479 test actives. After the split was corrected to 1,231 actives, the **realised ratio was
> 1:60**. The registration said 1:50; the run was at 1:60. Both are recorded here rather than silently
> reconciled, because a ratio that drifts after registration is exactly the kind of thing that makes a
> score incomparable.

## The bar — all must hold on the scaffold-test split at 1:50, with bootstrap lower bounds

| Metric | Point estimate | 95% CI lower bound |
|---|---|---|
| ROC-AUC (fit-score ranking) | ≥ 0.75 | ≥ 0.70 |
| EF1% | ≥ 5 | ≥ 3 |
| Operating point | recall ≥ 0.40, decoy pass ≤ 10% | — |

## Gates

| Gate | When | Continue only if |
|---|---|---|
| **G1** | End of Day 1, **dev set only** | TPR/FPR ≥ 2.0 at recall ≥ 0.40 **and** dev AUC ≥ 0.65 |
| **G2** | End of Day 2 | dev AUC ≥ 0.70 **and** TPR/FPR ≥ 3 at recall ≥ 0.40 |

**Failing either gate means stop.** No extra tiers, no tolerance searching, no "one more configuration".
**At most three configurations total.**

## The test set

Touched **once**, on Day 3. Not used for tuning, not used for model selection, not peeked at.

The test set is stratified by binding mode: ~22% of test actives carry an aryl-urea / type-II-like motif
— the chemotype the model was deliberately built to exclude. The **type-I subset is primary**; the full
set is reported as secondary. Neither is dropped.

## Configurations, in order

| # | Configuration |
|---|---|
| C0 | Baseline — 2 isotropic HBA features (Met109, Lys53), both required |
| C1 | + directional H-bond vector at Met109 (45° cone); Lys53 made optional |
| C2 | + exclusion volumes from the 8 type-I ensemble structures (exclude a point if ≥6 of 8 structures have a protein atom within 2.2 Å) |

## Evaluation protocol

- **Decoys:** stratified 10,000-decoy subsample from the full 73,950, stratified on the six DUD-E
  matching properties. FPR standard error ≈ 0.4%.
- **k selection:** tuned on **training actives plus training-matched decoys only**.
- **No cross-ratio comparison.** All GH values at 1:50.

## What was registered as the expected failure mode

Recall collapsing faster than decoy pass, as exclusion volumes suppress flexible molecules that lose
conformer coverage at 8 conformers.

**This did not happen.** See `docs/FINDINGS.md` — the registered signal was mis-specified.

## Deviations, logged at the time

- Dev decoys: only 3,600 of 73,950 matched the current training actives, so 6,400 new train-matched
  decoys were generated.
- One post-registration code change: `n_jobs` 13 → 8 (memory pressure).
- The clash tolerance was calibrated on the co-crystal ligands **before** registration.
