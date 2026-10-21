"""Güner–Henry score for binary pharmacophore filters, with a ratio guard.

The Güner–Henry (GH) score is the conventional metric for a *binary* screening
filter — one that returns hit / no-hit rather than a ranking. For a ranked
screen, EF1% / BEDROC / ROC-AUC are the appropriate metrics instead; GH is
reported at a fixed active:decoy ratio and only means something at that ratio.

Reference
---------
Güner, O. F. & Henry, D. R. (2000). Metric for analyzing hit lists and
pharmacophore searches. *J. Chem. Inf. Comput. Sci.* — the score is

    GH = (Ha / (4 · Ht · A)) · (3A + Ht) · (1 − (Ht − Ha) / (D − A))

    Ha = hits among actives        Ht = total hits
    A  = actives in the database   D  = database size (actives + decoys)

Why this module exists
----------------------
GH depends on the active:decoy ratio. Computing it at 3:1 and at 1:50 gives
very different numbers for the same model, so scores from different ratios are
not comparable — yet they get compared all the time. This module makes the
ratio part of the score's identity and raises if you try to compare across it.

Convention: every comparison here is actives-vs-decoys, where "hit" means
"passed the filter". A decoy hit is a false positive.

Typical reading (Güner & Henry 2000, and common practice since):
    GH >= 0.5   useful
    GH >= 0.7   good
    GH >= 0.9   excellent
"""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt


class RatioMismatchError(ValueError):
    """Raised when scores computed at different active:decoy ratios are compared."""


@dataclass(frozen=True)
class Contingency:
    """A binary screen's outcome against a known-truth actives/decoy set.

    Attributes
    ----------
    tp : actives that passed the filter
    fn : actives that did not
    fp : decoys that passed the filter
    tn : decoys that did not
    """

    tp: int
    fn: int
    fp: int
    tn: int

    def __post_init__(self) -> None:
        for name in ("tp", "fn", "fp", "tn"):
            value = getattr(self, name)
            if not isinstance(value, int):
                raise TypeError(f"{name} must be int, got {type(value).__name__}")
            if value < 0:
                raise ValueError(f"{name} must be >= 0, got {value}")
        if self.n_actives == 0:
            raise ValueError("no actives: the screen set contains zero positives")
        if self.n_decoys == 0:
            raise ValueError("no decoys: a ratio is undefined without negatives")

    @property
    def n_actives(self) -> int:
        return self.tp + self.fn

    @property
    def n_decoys(self) -> int:
        return self.fp + self.tn

    @property
    def n_database(self) -> int:
        return self.n_actives + self.n_decoys

    @property
    def ratio(self) -> float:
        """Active:decoy ratio as actives per decoy (1:50 -> 0.02)."""
        return self.n_actives / self.n_decoys

    @property
    def recall(self) -> float:
        return self.tp / self.n_actives

    @property
    def false_positive_rate(self) -> float:
        return self.fp / self.n_decoys

    @property
    def specificity(self) -> float:
        return 1.0 - self.false_positive_rate

    @property
    def precision(self) -> float:
        hits = self.tp + self.fp
        return self.tp / hits if hits else 0.0

    @property
    def tpr_fpr(self) -> float:
        """Enrichment as a ratio of rates. Undefined (inf) at zero decoy pass."""
        fpr = self.false_positive_rate
        return self.recall / fpr if fpr else float("inf")

    @property
    def mcc(self) -> float:
        """Matthews correlation coefficient. 0.0 when a margin is empty."""
        tp, fn, fp, tn = self.tp, self.fn, self.fp, self.tn
        denom = sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        return (tp * tn - fp * fn) / denom if denom else 0.0


@dataclass(frozen=True)
class GHScore:
    """A GH score bound to the ratio it was computed at."""

    value: float
    ratio: float
    contingency: Contingency

    @property
    def ratio_label(self) -> str:
        """Human-readable ratio, e.g. '1:50' or '1:3.1'."""
        if self.ratio >= 1.0:
            return f"{self.ratio:.0f}:1"
        inverse = 1.0 / self.ratio
        return f"1:{inverse:.0f}" if inverse >= 10 else f"1:{inverse:.1f}"

    def compare(self, other: "GHScore") -> int:
        """Compare two scores, refusing to compare across ratios.

        Returns -1, 0 or 1. Raises RatioMismatchError if the ratios differ by
        more than a rounding tolerance, because a GH comparison across ratios
        is not a comparison of models.
        """
        if abs(self.ratio - other.ratio) > 1e-9 * max(self.ratio, other.ratio):
            raise RatioMismatchError(
                f"cannot compare GH across ratios: {self.ratio_label} vs "
                f"{other.ratio_label}. GH is ratio-dependent; recompute both at "
                "the same active:decoy ratio before comparing."
            )
        if self.value < other.value:
            return -1
        return 0 if self.value == other.value else 1

    def __str__(self) -> str:
        return f"GH={self.value:.3f} (at {self.ratio_label})"


def gh_score(contingency: Contingency) -> GHScore:
    """Compute the Güner–Henry score for a binary screen.

    Raises ValueError when the score is undefined: no actives in the database,
    no hits at all, or every compound a hit (the (D − A) denominator collapses).
    """
    tp, fn, fp, tn = contingency.tp, contingency.fn, contingency.fp, contingency.tn
    ha = tp                      # hits among actives
    ht = tp + fp                 # total hits
    a = tp + fn                  # actives in database
    d = tp + fn + fp + tn        # database size

    if ht == 0:
        raise ValueError("no hits: GH is undefined when the filter passes nothing")
    if d == a:
        raise ValueError("no decoys: GH is undefined without a negative class")

    value = (ha / (4.0 * ht * a)) * (3.0 * a + ht) * (1.0 - (ht - ha) / (d - a))
    return GHScore(value=value, ratio=contingency.ratio, contingency=contingency)


def gh_at_fixed_ratio(
    recall: float,
    false_positive_rate: float,
    n_actives: int,
    n_decoys: int,
) -> GHScore:
    """GH from rates rather than counts — the form used to show ratio dependence.

    Useful for answering "what would this model score at the pre-registered
    ratio?" without re-running the screen.
    """
    if not 0.0 <= recall <= 1.0:
        raise ValueError(f"recall must be in [0, 1], got {recall}")
    if not 0.0 <= false_positive_rate <= 1.0:
        raise ValueError(f"false_positive_rate must be in [0, 1], got {false_positive_rate}")
    tp = round(recall * n_actives)
    fp = round(false_positive_rate * n_decoys)
    return gh_score(Contingency(tp=tp, fn=n_actives - tp, fp=fp, tn=n_decoys - fp))
