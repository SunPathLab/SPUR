from spur_codex.enumeration import BranchEnumerator, RouteEnumerator, TruncalEnumerator
from spur_codex.route import RouteType


def test_branch_enumerator_uses_copy_fate_event_filter():
    truncal = TruncalEnumerator().enumerate_for_state(3, 1)[0]
    branches = BranchEnumerator().enumerate(
        truncal=truncal,
        target_total_cn=3,
        target_minor_cn=1,
        max_events=1,
    )

    assert branches
    assert all(branch.n_events <= 1 for branch in branches)


def test_route_enumerator_splits_all_routes_into_pure_gain_and_mixed_gain_loss():
    enumerator = RouteEnumerator()
    all_routes = enumerator.enumerate((4, 2), (3, 1), max_truncal_cn=6, max_branch_events=3)
    pure_gain_routes = enumerator.enumerate(
        (4, 2),
        (3, 1),
        max_truncal_cn=6,
        max_branch_events=3,
        route_type=RouteType.PURE_GAIN,
    )
    mixed_gain_loss_routes = enumerator.enumerate(
        (4, 2),
        (3, 1),
        max_truncal_cn=6,
        max_branch_events=3,
        route_type=RouteType.MIXED_GAIN_LOSS,
    )

    assert len(all_routes) == len(pure_gain_routes) + len(mixed_gain_loss_routes)
    assert all(route.is_pure_gain for route in pure_gain_routes)
    assert all(route.is_mixed_gain_loss for route in mixed_gain_loss_routes)


def test_updated_mixed_gain_loss_counts_for_4_2_vs_3_1_thresholds():
    enumerator = RouteEnumerator()
    expected = {
        3: (2, 425, 427),
        4: (2, 1433, 1435),
        5: (2, 2330, 2332),
    }

    for threshold, (n_pure_gain, n_mixed_gain_loss, n_all) in expected.items():
        pure_gain_routes = enumerator.enumerate(
            (4, 2),
            (3, 1),
            max_truncal_cn=6,
            max_branch_events=threshold,
            route_type=RouteType.PURE_GAIN,
        )
        mixed_gain_loss_routes = enumerator.enumerate(
            (4, 2),
            (3, 1),
            max_truncal_cn=6,
            max_branch_events=threshold,
            route_type=RouteType.MIXED_GAIN_LOSS,
        )

        assert len(pure_gain_routes) == n_pure_gain
        assert len(mixed_gain_loss_routes) == n_mixed_gain_loss
        assert len(pure_gain_routes) + len(mixed_gain_loss_routes) == n_all
