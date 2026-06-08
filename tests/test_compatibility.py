import numpy as np

from spur_codex.compatibility import PureGainCompatibilityChecker
from spur_codex.enumeration import RouteEnumerator
from spur_codex.route import RouteType
from spur_codex.simulation import sample_route_jmd


def test_pure_gain_route_jmd_is_pure_gain_compatible():
    rng = np.random.default_rng(1)
    sample1_cn = (3, 1)
    sample2_cn = (3, 1)
    route = RouteEnumerator().enumerate(
        sample1_cn,
        sample2_cn,
        max_mrca_cn=4,
        max_branch_events=3,
        route_type=RouteType.PURE_GAIN,
    )[0]

    _, q = sample_route_jmd(route, sample1_cn, sample2_cn, rng)

    checker = PureGainCompatibilityChecker(max_mrca_cn=4, max_branch_events=3)
    assert checker.is_compatible(q, sample1_cn, sample2_cn)


def test_some_mixed_gain_loss_jmd_is_conclusive():
    rng = np.random.default_rng(2)
    sample1_cn = (4, 2)
    sample2_cn = (3, 1)
    routes = RouteEnumerator().enumerate(
        sample1_cn,
        sample2_cn,
        max_mrca_cn=6,
        max_branch_events=4,
        route_type=RouteType.MIXED_GAIN_LOSS,
    )
    checker = PureGainCompatibilityChecker(max_mrca_cn=6, max_branch_events=4)

    found_conclusive = False
    for route in routes[:200]:
        _, q = sample_route_jmd(route, sample1_cn, sample2_cn, rng)
        if not checker.is_compatible(q, sample1_cn, sample2_cn):
            found_conclusive = True
            break

    assert found_conclusive
