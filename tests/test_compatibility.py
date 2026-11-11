import numpy as np

from spur.compatibility import PureGainCompatibilityChecker
from spur.enumeration import RouteEnumerator
from spur.route import RouteType
from spur.simulation import sample_route_jmd


def test_pure_gain_route_jmd_is_pure_gain_compatible():
    rng = np.random.default_rng(1)
    sample1_cn = (3, 1)
    sample2_cn = (3, 1)
    route = RouteEnumerator().enumerate(
        sample1_cn,
        sample2_cn,
        max_truncal_cn=4,
        max_branch_events=3,
        route_type=RouteType.PURE_GAIN,
    )[0]

    _, q = sample_route_jmd(route, sample1_cn, sample2_cn, rng)

    checker = PureGainCompatibilityChecker(
        max_truncal_cn=4,
        route_filter="max_branch_events",
        route_filter_threshold=3,
    )
    assert checker.is_compatible(q, sample1_cn, sample2_cn)


def test_some_mixed_gain_loss_jmd_is_conclusive():
    rng = np.random.default_rng(2)
    sample1_cn = (4, 2)
    sample2_cn = (3, 1)
    routes = RouteEnumerator().enumerate(
        sample1_cn,
        sample2_cn,
        max_truncal_cn=6,
        max_branch_events=4,
        route_type=RouteType.MIXED_GAIN_LOSS,
    )
    checker = PureGainCompatibilityChecker(
        max_truncal_cn=6,
        route_filter="max_branch_events",
        route_filter_threshold=4,
    )

    found_conclusive = False
    for route in routes[:200]:
        _, q = sample_route_jmd(route, sample1_cn, sample2_cn, rng)
        if not checker.is_compatible(q, sample1_cn, sample2_cn):
            found_conclusive = True
            break

    assert found_conclusive


def test_compatibility_checker_applies_route_filter_to_candidates():
    sample1_cn = (4, 2)
    sample2_cn = (3, 1)
    checker = PureGainCompatibilityChecker(
        max_truncal_cn=6,
        route_filter="total_events",
        route_filter_threshold=8,
    )

    candidates = checker._candidates(sample1_cn, sample2_cn)

    assert candidates
    assert all(candidate.route.total_events <= 8 for candidate in candidates)


def test_compatibility_checker_requires_filter_threshold_when_filter_is_enabled():
    checker = PureGainCompatibilityChecker(
        max_truncal_cn=6,
        route_filter="total_events",
    )

    try:
        checker._candidates((4, 2), (3, 1))
    except ValueError as exc:
        assert "--route-filter-threshold" in str(exc)
    else:
        raise AssertionError("Expected missing route filter threshold to raise ValueError.")
