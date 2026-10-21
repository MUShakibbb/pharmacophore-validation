# pharmacophore-validation

**A pre-registered negative result, and the four checks that produced it.**

I set out to build a structure-based pharmacophore for a kinase target. I pre-registered a
performance gate, ran it, and **failed it**. This repository is the record: the gate, the numbers,
the stop, and the reusable checks that made the failure legible instead of mysterious.

The model is not the contribution. **The measurement discipline is.**

---

## What happened, in one paragraph

I derived a structure-based pharmacophore from an 11-structure crystal ensemble and a ligand-based
model from 4,929 ChEMBL actives, then pre-registered a gate before computing any fit score: **TPR/FPR
≥ 2.0 and ROC-AUC ≥ 0.65** on a scaffold-split held-out set at the decoy ratio. Three
configurations — baseline, hinge-vector, exclusion-volumes — reached **TPR/FPR 1.17–1.28 and AUC
0.54–0.56**. I stopped. No further tiers, no tolerance search, no test-set peek.

**The interesting part is *why* it failed**, and that is what this repo is about.

---

## The four findings

### 1. A 2D feature-presence prefilter predicts failure before you dock anything

A pharmacophore can only match a molecule that *possesses* its required feature types. That is
conformation-independent, so counting feature types per molecule is an **exact** prefilter — no false
negatives.

| Model | Required features | Chemically | Decoys rejected |
|---|---|---|---|
| Target A | 2 × HBA | 2 acceptors | **6.1%** |
| Target B | 3 × HBA + 1 × HBD | 3 acceptors + 1 donor | **3.6%** |

Against an expected 60–85% rejection. The features were HBA/HBD — **the most common feature types in
drug-like space.** Almost every decoy possesses them.

> **If your feature-presence prefilter cannot reject anything, your pharmacophore is chemically
> unspecific and cannot discriminate.** This costs seconds and runs before any conformer generation.

### 2. Güner–Henry is ratio-dependent, and comparing across ratios is invalid

The same recall and FPR give **wildly different GH** depending on the active:decoy ratio:

| Ratio | Actives | Decoys | GH |
|---|---|---|---|
| ~1:3 | 1,562 | 500 | **0.369** |
| ~1:60 | 1,231 | 73,950 | **0.086** |

I initially reported the 0.369. It was an artefact of a tuning run on a 500-decoy subsample. **At the
realised full-set ratio it was 0.086.**

Note the ratio itself drifted: the decoy pool was sized at 1:50 against an earlier split with 1,479
test actives. Correcting the split to 1,231 actives changed the realised ratio to **1:60** — a fact
that only surfaces if the ratio is carried alongside the score.

`gh_score.py` in this repo **refuses to compare scores computed at different ratios.**

### 3. Property-matched decoys can be property-confounded

DUD-E-style matching on six properties does not guarantee the sets are jointly matched. A random
forest on **physicochemical properties alone** — no 3D information, no receptor — scored:

| Set | Property-only AUC |
|---|---|
| Development | **0.829** |
| (reported elsewhere) | 0.80–0.83 |

**If a property-only model separates your actives from your decoys, your AUC partly measures property
separation, not model quality.** `decoy_sanity.py` runs this check in a few minutes and should be a
mandatory first step.

### 4. Public benchmarks leak into public training data

A "frozen external benchmark" I intended to use had **31 of its 40 actives present in ChEMBL** — the
same database my actives came from. It was not external for anything trained on ChEMBL.

`benchmark_overlap.py` checks InChIKey overlap between a benchmark and a training source before you
trust the word "external".

---

## What the gate actually showed

| Config | k | TPR/FPR | Recall | AUC | Verdict |
|---|---|---|---|---|---|
| C0 baseline (2 isotropic HBA) | 2 | 1.211 | 0.622 | **0.543** | fail |
| C1 + hinge H-bond vector, gatekeeper optional | 3 | 1.169 | 0.515 | **0.549** | fail |
| C2 + exclusion volumes | 1 | 1.279 | 0.529 | **0.563** | fail |
| *required* | — | *≥ 2.0* | *≥ 0.40* | *≥ 0.65* | — |

Two honest observations:

- **The baseline was a coin flip** (AUC 0.543). It matched plenty (recall 0.62) and discriminated
  nothing.
- **Geometry helped, marginally** — TPR/FPR 1.211 → 1.279, AUC 0.543 → 0.563. The right direction,
  nowhere near enough.
- **The registered failure signal did not fire.** I had predicted conformer collapse; the actual
  failure was different. A mis-specified diagnostic is a finding too.

---

## Repository layout

```
src/
  gh_score.py           corrected Güner–Henry, with a ratio guard
  decoy_sanity.py       property-only baseline — the check that should run first
  benchmark_overlap.py  InChIKey overlap between a benchmark and a training source
results/
  dev_results.json      the actual gate failure, verbatim
  lbp_hypothesis_scores.csv
docs/
  PRE_REGISTRATION.md   the gate, as committed before any score existed
  FINDINGS.md           longer write-up of the four findings
tests/
  test_gh_score.py      reproduces the published numbers
```

## Quick start

```bash
pip install -r requirements.txt
pytest tests/ -q                       # reproduces the GH numbers in this README
python src/decoy_sanity.py --actives actives.csv --decoys decoys.csv
python src/benchmark_overlap.py --benchmark bench.csv --training chembl_actives.csv
```

## On the negative result

I am aware that a failed model is a weak thing to publish. Three things make it worth writing down:

1. **The gate was registered before the score existed.** The git history shows the pre-registration
   commit preceding the results commit. That is checkable, which is the point.
2. **Each failure has a named cause** with a number attached — not "it didn't work".
3. **Three of the four checks would have saved the project**, and none of them is standard practice.

The reusable lessons are the ratio guard, the property-only baseline, and the benchmark overlap
check. They apply to any virtual screen, not just this one.

## Limitations

- **This does not show that no pharmacophore could work here.** It shows that *this model form*, at 8
  conformers, against a property-confounded decoy set, did not reach the gate.
- **The conformer count was low.** 8 per molecule. More conformers would change the result in an
  unknown direction — more conformers also raise decoy pass rates. That was never tested.
- **One target, one split, one decoy protocol.** No claim beyond that.
- **Retrospective ChEMBL enrichment says little about natural-product transfer**, which was the
  original application. Nothing here addresses that.

## Licence

MIT — see `LICENSE`.
