"""Tests reproduce the numbers reported in the README.

If these pass, the GH implementation agrees with the published values, and the
ratio-dependence claim is demonstrated rather than asserted.
"""

import pytest

from src.gh_score import (
    Contingency,
    GHScore,
    RatioMismatchError,
    gh_at_fixed_ratio,
    gh_score,
)


def test_gh_reproduces_reported_full_run():
    """The real 1:50 gate run: TP=30 FP=1845 FN=93 TN=9355 -> GH 0.061."""
    c = Contingency(tp=30, fn=93, fp=1845, tn=9355)
    score = gh_score(c)

    assert round(score.value, 3) == 0.061
    assert round(c.recall, 3) == 0.244
    assert round(c.false_positive_rate, 3) == 0.165
    assert round(c.tpr_fpr, 2) == 1.48
    assert round(c.mcc, 3) == 0.022
    assert round(c.precision, 4) == 0.0160


def test_gh_is_ratio_dependent():
    """Same recall and FPR, two ratios: 0.369 at ~3:1 versus 0.087 at 1:50."""
    at_three_to_one = gh_at_fixed_ratio(recall=0.653, false_positive_rate=0.516,
                                        n_actives=1562, n_decoys=500)
    at_one_to_fifty = gh_at_fixed_ratio(recall=0.653, false_positive_rate=0.516,
                                        n_actives=1231, n_decoys=73950)

    assert round(at_three_to_one.value, 3) == 0.369
    assert round(at_one_to_fifty.value, 3) == 0.087

    # A four-fold difference from the ratio alone. This is why cross-ratio
    # comparison is refused below.
    assert at_three_to_one.value / at_one_to_fifty.value > 4.0


def test_compare_refuses_cross_ratio():
    """The ratio guard is the point of the module, so it gets its own test."""
    a = gh_at_fixed_ratio(0.65, 0.52, 1562, 500)     # ~3:1
    b = gh_at_fixed_ratio(0.65, 0.52, 1231, 73950)   # 1:50

    with pytest.raises(RatioMismatchError) as excinfo:
        a.compare(b)

    # The message must say what to do, not just that it failed.
    assert "same active:decoy ratio" in str(excinfo.value)


def test_compare_allows_same_ratio():
    a = gh_at_fixed_ratio(0.65, 0.10, 1000, 50000)
    b = gh_at_fixed_ratio(0.65, 0.20, 1000, 50000)

    assert a.compare(b) == 1      # lower decoy pass -> higher GH
    assert b.compare(a) == -1
    assert a.compare(a) == 0


def test_ratio_label():
    assert gh_at_fixed_ratio(0.5, 0.1, 100, 5000).ratio_label == "1:50"
    assert gh_at_fixed_ratio(0.5, 0.1, 3000, 1000).ratio_label == "3:1"


def test_gh_rejects_degenerate_input():
    with pytest.raises(ValueError, match="no hits"):
        gh_score(Contingency(tp=0, fn=100, fp=0, tn=10000))

    with pytest.raises(ValueError, match="no actives"):
        gh_score(Contingency(tp=0, fn=0, fp=10, tn=100))

    with pytest.raises(ValueError, match="no decoys"):
        gh_score(Contingency(tp=5, fn=5, fp=0, tn=0))


def test_gh_at_fixed_ratio_validates_rates():
    with pytest.raises(ValueError, match="recall"):
        gh_at_fixed_ratio(1.5, 0.1, 100, 1000)
    with pytest.raises(ValueError, match="false_positive_rate"):
        gh_at_fixed_ratio(0.5, -0.1, 100, 1000)


def test_contingency_rejects_negative_counts():
    with pytest.raises(ValueError, match="tp"):
        Contingency(tp=-1, fn=10, fp=1, tn=100)


def test_score_str_is_unambiguous():
    """1,231 actives against 73,950 decoys is 1:60, not the nominal 1:50.

    The decoy pool was sized at 1:50 against an earlier split that had 1,479
    test actives. Correcting the split changed the realised ratio, which is
    exactly the kind of drift a ratio label has to make visible.
    """
    s = gh_at_fixed_ratio(0.65, 0.516, 1231, 73950)
    assert s.ratio_label == "1:60"
    assert f"GH={s.value:.3f}" in str(s)
    assert "1:60" in str(s)
