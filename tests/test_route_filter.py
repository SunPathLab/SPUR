import numpy as np

from spur.enumeration import RouteEnumerator
from spur.route import RouteType
from spur.route_filter import collect_h_by_truncal, filter_routes, group_by_truncal, split_routes


def test_split_routes_matches_route_type_filters():
    routes = RouteEnumerator().enumerate((4, 2), (3, 1), max_truncal_cn=6, max_branch_events=3)
    split = split_routes(routes)

    assert split.pure_gain == filter_routes(routes, RouteType.PURE_GAIN)
    assert split.mixed_gain_loss == filter_routes(routes, RouteType.MIXED_GAIN_LOSS)
    assert len(routes) == len(split.pure_gain) + len(split.mixed_gain_loss)


def test_group_by_truncal_uses_canonical_key():
    routes = RouteEnumerator().enumerate((3, 1), (3, 1), max_truncal_cn=4, max_branch_events=2)
    groups = group_by_truncal(routes)

    assert len(groups) > 1
    assert sum(len(group.routes) for group in groups) == len(routes)
    for group in groups:
        assert group.truncal_key == group.truncal.canonical_key
        assert "major=" in group.truncal_key
        assert "minor=" in group.truncal_key
        assert "L" in group.truncal_key
        assert all(route.truncal_key == group.truncal_key for route in group.routes)


def test_route_h_matches_group_h_tensor_slice():
    routes = RouteEnumerator().enumerate((3, 1), (3, 1), max_truncal_cn=4, max_branch_events=2)
    group = group_by_truncal(routes)[0]

    H = group.h_tensor()

    assert H.shape[0] == len(group.routes)
    assert np.array_equal(group.routes[0].h(), H[0])


def test_collect_h_by_truncal_returns_one_tensor_per_group():
    routes = RouteEnumerator().enumerate((3, 1), (3, 1), max_truncal_cn=4, max_branch_events=2)
    groups = group_by_truncal(routes)
    by_key = collect_h_by_truncal(routes)

    assert set(by_key) == {group.truncal_key for group in groups}
    for group in groups:
        assert by_key[group.truncal_key].shape[0] == len(group.routes)
