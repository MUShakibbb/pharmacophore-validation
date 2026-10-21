"""Property-only baseline: is your decoy set separable without any 3D information?

If a random forest trained on six physicochemical descriptors — molecular weight,
cLogP, H-bond donors, H-bond acceptors, rotatable bonds, formal charge — can tell
your actives from your decoys, then the two sets differ in ways that have nothing
to do with binding. Any AUC your 3D model achieves is then partly measuring
property separation, and the benchmark cannot resolve model quality.

This check costs minutes and runs before any docking or conformer generation.

Threshold
---------
AUC <= 0.60  the decoy set is adequately matched; proceed
AUC  > 0.60  confounded; fix the decoy set before trusting any AUC

0.60 is a judgement call, not a law. It is the level at which the property
signal stops dominating a typical 3D model's 0.7–0.8. State whatever threshold
you use before running, and treat it as pre-registered.

Usage
-----
    python src/decoy_sanity.py --actives actives.csv --decoys decoys.csv
    python src/decoy_sanity.py --actives a.csv --decoys d.csv --smiles-col smiles
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, Lipinski
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

RDLogger.DisableLog("rdApp.*")

DESCRIPTORS = {
    "mw": Descriptors.MolWt,
    "clogp": Crippen.MolLogP,
    "hbd": Lipinski.NumHDonors,
    "hba": Lipinski.NumHAcceptors,
    "rotb": Lipinski.NumRotatableBonds,
    "charge": lambda m: Chem.GetFormalCharge(m),
}


def featurise(smiles: list[str]) -> tuple[np.ndarray, list[int]]:
    """Six DUD-E properties per molecule. Unparseable SMILES are dropped."""
    rows, kept = [], []
    for i, smi in enumerate(smiles):
        mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        if mol is None:
            continue
        rows.append([fn(mol) for fn in DESCRIPTORS.values()])
        kept.append(i)
    if not rows:
        raise ValueError("no valid SMILES could be parsed from the input")
    return np.asarray(rows, dtype=float), kept


def property_only_auc(
    actives: list[str],
    decoys: list[str],
    n_splits: int = 5,
    seed: int = 0,
) -> dict:
    """Cross-validated AUC of a property-only classifier.

    Cross-validation is stratified and repeated over folds; the reported AUC is
    the mean of held-out folds, so it is not an in-sample fit.
    """
    Xa, _ = featurise(actives)
    Xd, _ = featurise(decoys)
    X = np.vstack([Xa, Xd])
    y = np.concatenate([np.ones(len(Xa)), np.zeros(len(Xd))])

    folds = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    aucs = []
    for train_idx, test_idx in folds.split(X, y):
        clf = RandomForestClassifier(
            n_estimators=300, min_samples_leaf=3, random_state=seed, n_jobs=-1
        )
        clf.fit(X[train_idx], y[train_idx])
        proba = clf.predict_proba(X[test_idx])[:, 1]
        aucs.append(roc_auc_score(y[test_idx], proba))

    return {
        "auc_mean": float(np.mean(aucs)),
        "auc_std": float(np.std(aucs)),
        "auc_folds": [float(a) for a in aucs],
        "n_actives": int(len(Xa)),
        "n_decoys": int(len(Xd)),
        "n_dropped": int(len(actives) + len(decoys) - len(X)),
    }


def read_smiles(path: str, column: str | None = None) -> list[str]:
    """Read SMILES from a .smi/.txt (one per line) or a CSV with a header."""
    with open(path, encoding="utf-8") as fh:
        head = fh.readline()
        fh.seek(0)
        if "," in head:
            import csv

            reader = csv.DictReader(fh)
            col = column or _guess_column(reader.fieldnames or [])
            return [row[col] for row in reader if row.get(col)]
        return [line.split()[0] for line in fh if line.strip()]


def _guess_column(fieldnames: list[str]) -> str:
    for candidate in ("smiles", "SMILES", "canonical_smiles", "Smiles"):
        if candidate in fieldnames:
            return candidate
    raise ValueError(
        f"no SMILES column found in {fieldnames}; pass --smiles-col explicitly"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--actives", required=True, help="actives: .smi or .csv")
    parser.add_argument("--decoys", required=True, help="decoys: .smi or .csv")
    parser.add_argument("--smiles-col", default=None, help="CSV column holding SMILES")
    parser.add_argument("--threshold", type=float, default=0.60,
                        help="AUC above which the decoy set is called confounded")
    parser.add_argument("--splits", type=int, default=5)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)

    actives = read_smiles(args.actives, args.smiles_col)
    decoys = read_smiles(args.decoys, args.smiles_col)
    result = property_only_auc(actives, decoys, args.splits, args.seed)

    print(f"actives: {result['n_actives']}   decoys: {result['n_decoys']}"
          f"   dropped (unparseable): {result['n_dropped']}")
    print(f"property-only AUC: {result['auc_mean']:.3f} "
          f"+/- {result['auc_std']:.3f}   (folds: "
          f"{', '.join(f'{a:.3f}' for a in result['auc_folds'])})")

    if result["auc_mean"] > args.threshold:
        print(f"\nFAIL  AUC {result['auc_mean']:.3f} > {args.threshold:.2f}: "
              "the decoy set is separable on properties alone.")
        print("      Any 3D model's AUC on this set is partly property "
              "separation. Fix the decoys before trusting it.")
        return 1

    print(f"\nPASS  AUC {result['auc_mean']:.3f} <= {args.threshold:.2f}: "
          "properties alone do not separate the sets.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
