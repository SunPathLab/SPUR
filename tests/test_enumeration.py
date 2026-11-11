from spur.enumeration import BranchEnumerator, RouteEnumerator, TruncalEnumerator
from spur.route import RouteType


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


def test_truncal_enumerator_keeps_loh_and_non_loh_modes_separate():
    enumerator = TruncalEnumerator()

    loh_truncals = enumerator.enumerate_loh_up_to(4)
    non_loh_truncals = enumerator.enumerate_up_to(4)

    assert loh_truncals
    assert {truncal.truncal_total_cn for truncal in loh_truncals} == {1, 2, 3, 4}
    assert all(truncal.truncal_minor_cn == 0 for truncal in loh_truncals)
    assert all(truncal.truncal_minor_cn > 0 for truncal in non_loh_truncals)


def test_empty_truncal_minor_allele_has_one_empty_fate_vector():
    assert BranchEnumerator._allele_vectors(0, 0) == [()]
    assert BranchEnumerator._allele_vectors(0, 1) == []


def test_route_enumerator_uses_loh_truncals_only_for_loh_loh_pair():
    routes = RouteEnumerator().enumerate(
        (3, 0),
        (2, 0),
        max_truncal_cn=3,
        max_branch_events=5,
    )

    assert routes
    assert all(route.truncal_minor_cn == 0 for route in routes)


def test_route_enumerator_uses_non_loh_truncals_for_loh_non_loh_pair():
    routes = RouteEnumerator().enumerate(
        (2, 0),
        (2, 1),
        max_truncal_cn=3,
        max_branch_events=5,
    )

    assert routes
    assert all(route.truncal_minor_cn > 0 for route in routes)


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
