"""Is your "external" benchmark actually external?

A frozen benchmark is only external for models that never saw it. If its actives
overlap the training source — same ChEMBL release, same vendor catalogue, same
publication — it is an internal test wearing an external label, and the
enrichment it reports is optimistic.

This checks InChIKey overlap between a benchmark's actives and a training set.

The 14-character prefix (the connectivity block) is used by default, so
stereoisomers and salt forms of a training compound still count as overlap. That
is the conservative choice: if the same scaffold connectivity was in training,
the benchmark is not measuring generalisation.

Usage
-----
    python src/benchmark_overlap.py --benchmark bench.smi --training train.smi
    python src/benchmark_overlap.py --benchmark b.csv --training t.csv --key-col inchikey
"""

from __future__ import annotations

import argparse
import csv
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def to_key(value: str, full: bool = False) -> str | None:
    """InChIKey for a SMILES or an existing InChIKey. None if unusable."""
    value = (value or "").strip()
    if not value:
        return None
    # Already a key? 14 chars = connectivity block, 27 = full.
    if len(value) in (14, 27) and value.isupper() and "-" in value:
        return value if full else value.split("-")[0]
    mol = Chem.MolFromSmiles(value)
    if mol is None:
        return None
    key = Chem.MolToInchiKey(mol)
    return key if full else key.split("-")[0]


def read_keys(path: str, column: str | None, full: bool = False) -> set[str]:
    keys: set[str] = set()
    with open(path, encoding="utf-8") as fh:
        head = fh.readline()
        fh.seek(0)
        if "," in head:
            reader = csv.DictReader(fh)
            col = column or _guess_column(reader.fieldnames or [])
            values = [row.get(col, "") for row in reader]
        else:
            values = [line.split()[0] for line in fh if line.strip()]
    for value in values:
        key = to_key(value, full=full)
        if key:
            keys.add(key)
    return keys


def _guess_column(fieldnames: list[str]) -> str:
    for candidate in ("inchikey", "InChIKey", "smiles", "SMILES", "canonical_smiles"):
        if candidate in fieldnames:
            return candidate
    raise ValueError(
        f"no InChIKey or SMILES column found in {fieldnames}; pass --key-col"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--training", required=True)
    parser.add_argument("--key-col", default=None)
    parser.add_argument("--full-key", action="store_true",
                        help="use the full 27-char key (stereo-sensitive)")
    parser.add_argument("--threshold", type=float, default=0.05,
                        help="overlap fraction above which the benchmark is not external")
    args = parser.parse_args(argv)

    bench = read_keys(args.benchmark, args.key_col, args.full_key)
    train = read_keys(args.training, args.key_col, args.full_key)

    if not bench:
        print("ERROR: no usable keys in the benchmark file")
        return 2

    overlap = bench & train
    fraction = len(overlap) / len(bench)

    print(f"benchmark actives : {len(bench)}")
    print(f"training keys     : {len(train)}")
    print(f"overlap           : {len(overlap)} ({fraction:.1%})")

    if overlap:
        print(f"shared keys (first 10): {sorted(overlap)[:10]}")

    if fraction > args.threshold:
        print(f"\nFAIL  {fraction:.1%} overlap > {args.threshold:.0%}: this benchmark "
              "is not external for this training source.")
        print("      Enrichment measured on it is an internal test, not "
              "generalisation. Say so, or choose another benchmark.")
        return 1

    print(f"\nPASS  {fraction:.1%} overlap <= {args.threshold:.0%}: "
          "the benchmark is acceptably external.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
