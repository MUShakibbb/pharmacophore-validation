# Findings

Four checks, each of which caught something real in this project. Written up as a methods note: what
the check does, what it found here, and what it would catch elsewhere.

---

## 1. Feature-presence prefilter

### What it does

A pharmacophore match requires the molecule to *possess* the feature types the model specifies. Feature
presence is a property of the molecular graph, not of any conformation. So counting feature types per
molecule is an **exact** necessary-condition test: a molecule that lacks the required types can never
match, in any conformation. No false negatives — unlike a conformer-based prefilter, which can reject a
flexible molecule that would have matched given better sampling.

Cost: one pass over the library with a feature factory. Seconds to minutes.

### What it found

| Model | Required | Chemically | Decoy rejection |
|---|---|---|---|
| A | 2 × HBA | 2 acceptors | **6.1%** |
| B | 3 × HBA + 1 × HBD | 3 acceptors + 1 donor | **3.6%** |

An 60–85% rejection was expected. Both models required only hydrogen-bond acceptors and donors — the
most abundant feature types in drug-like chemistry. Nearly every decoy possesses them, so the
prefilter could not reject.

### The rule

> **If the feature-presence prefilter cannot reject anything, the pharmacophore is chemically
> unspecific and will not discriminate.**

This is a cheap falsification test that runs before any conformer generation or docking. It would have
predicted the failure of model B from its feature list alone.

### What it would catch elsewhere

Pharmacophores built from a single co-crystal, where the features are whatever that one ligand happened
to touch. Those models routinely require only acceptors and aromatics, and routinely fail to
discriminate.

---

## 2. Ratio-dependent metrics

### What it does

Güner–Henry is the conventional score for a *binary* screening filter. It is a function of the
active:decoy ratio as well as the model's performance:

```
GH = (Ha / (4 · Ht · A)) · (3A + Ht) · (1 − (Ht − Ha) / (D − A))
```

with `Ha` hits among actives, `Ht` total hits, `A` actives, `D` database size. Change the ratio at fixed
recall and FPR and the score moves.

### What it found

The same model, same recall (0.653), same false-positive rate (0.516):

| Ratio | Actives | Decoys | GH |
|---|---|---|---|
| ~1:3 | 1,562 | 500 | **0.369** |
| ~1:60 | 1,231 | 73,950 | **0.086** |

I reported the 0.369 first. It came from a tuning run on a 500-decoy subsample, not the full-set
evaluation. **A four-fold inflation from the ratio alone.**

The ratio also drifted underneath me: the pool was sized at 1:50 against an earlier split with 1,479
test actives. Correcting the split to 1,231 changed the realised ratio to **1:60**. Nobody would have
noticed, because the score was being reported without its ratio.

### The rule

> **Never compare GH across active:decoy ratios.** Recompute both at the same ratio, or don't compare.

`gh_score.py` makes the ratio part of the score's identity and raises `RatioMismatchError` on
cross-ratio comparison.

### What it would catch elsewhere

Any paper reporting GH or enrichment factors without stating the ratio. This is common, and it makes
reported scores incomparable even when the models are directly comparable.

---

## 3. Property-only baseline

### What it does

Train a classifier on six physicochemical descriptors — MW, cLogP, HBD, HBA, rotatable bonds, formal
charge — and no structural or 3D information. Cross-validate. Report AUC.

If that classifier separates actives from decoys, the two sets differ in properties. Any 3D model's AUC
is then partly measuring that property separation rather than binding.

### What it found

| Set | Property-only AUC |
|---|---|
| Development split | **0.829** |

The decoys were matched on those same six properties by a DUD-E-style nearest-neighbour procedure, and
they were still separable. Matching marginally on each property does not produce a jointly matched
set — the joint distribution can differ even when every marginal agrees.

### The rule

> **Run the property-only baseline before any 3D model.** If AUC > ~0.60, the decoy set is confounded
> and no 3D AUC from it can be interpreted.

0.60 is a judgement call, not a law. Pre-register whatever threshold you use.

### What it would catch elsewhere

Most published decoy sets. DUD-E itself has been criticised on exactly this point, and the criticism is
easy to test — which is why not testing it is hard to excuse.

---

## 4. Benchmark overlap

### What it does

InChIKey overlap between a benchmark's actives and a training source. Uses the 14-character
connectivity block by default, so stereoisomers and salt forms of a training compound still count as
overlap — the conservative choice.

### What it found

A benchmark I intended to use as an "independent external" set had **31 of its 40 actives present in
ChEMBL** — the same database my actives came from.

It was not external for anything trained on ChEMBL. It would have been external only for a model
trained on something else.

### The rule

> **A benchmark is external only with respect to a named training source.** Check the overlap, and state
> the fraction.

### What it would catch elsewhere

Almost every "external validation" claim where both benchmark and training data trace back to ChEMBL —
which is most of them.

---

## A fifth observation, from the failure itself

The registered failure signal did not fire.

I predicted that exclusion volumes would suppress recall faster than decoy pass — flexible molecules
losing conformer coverage. The measured signal was `failure_path: false`. The actual failure was
different: at the most constrained configuration, only **0.59% of actives** could match at all, because
eight rigid conformers cannot satisfy a four-point alignment plus exclusion volumes.

**A mis-specified diagnostic is itself a finding.** It says the mental model of the failure was wrong,
and it tells the next person what to check first — conformer coverage, not recall collapse.

---

## What this does and does not show

**Shows:** this model form, at eight conformers per molecule, against a property-confounded decoy set,
did not reach a pre-registered TPR/FPR ≥ 2.0 and AUC ≥ 0.65 gate.

**Does not show:** that no pharmacophore could work on this target. Three limits were identified and
only one was tested:

| Limit | Status |
|---|---|
| Conformer sampling (8 per molecule) | **untested** — more conformers would change the result in an unknown direction, since more conformers also raise decoy pass |
| Decoy confounding (property-only AUC 0.829) | **demonstrated**, unfixed |
| Feature specificity (prefilter rejects 3.6%) | **demonstrated**, inherent to the feature set |

Two of the three were never addressed. The negative result is about the configuration that was run, not
about the method in general.
