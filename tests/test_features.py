import numpy as np

from spur_codex.enumeration import RouteEnumerator
from spur_codex.features import RouteScoreResult
from spur_codex.route import RouteType


def test_route_score_result_score_features_uses_mixed_gain_loss_names():
    pure_gain_routes = RouteEnumerator().enumerate(
        (3, 1),
        (3, 1),
        max_truncal_cn=3,
        max_branch_events=2,
        route_type=RouteType.PURE_GAIN,
    )
    mixed_gain_loss_routes = RouteEnumerator().enumerate(
        (3, 1),
        (3, 1),
        max_truncal_cn=3,
        max_branch_events=2,
        route_type=RouteType.MIXED_GAIN_LOSS,
    )[:3]
    result = RouteScoreResult(
        routes=tuple(pure_gain_routes + mixed_gain_loss_routes),
        scores=np.array([-10.0, -12.0, -9.0, -11.0, -13.0]),
    )

    features = result.score_features()

    assert features["best_pg"] == -10.0
    assert features["best_mixed_gain_loss"] == -9.0
    assert features["delta_best_mixed_gain_loss_minus_pg"] == 1.0
    assert features["n_pg_routes"] == 2
    assert features["n_mixed_gain_loss_routes"] == 3


def test_route_score_result_score_features_handles_pure_gain_only_routes():
    routes = RouteEnumerator().enumerate(
        (3, 1),
        (3, 1),
        max_truncal_cn=3,
        max_branch_events=2,
        route_type=RouteType.PURE_GAIN,
    )
    result = RouteScoreResult(
        routes=tuple(routes),
        scores=np.array([-10.0, -12.0]),
    )

    features = result.score_features()

    assert features["best_pg"] == -10.0
    assert features["n_pg_routes"] == 2
    assert features["n_mixed_gain_loss_routes"] == 0
    assert np.isnan(features["best_mixed_gain_loss"])
    assert np.isnan(features["delta_best_mixed_gain_loss_minus_pg"])
    assert np.isnan(features["mixed_gain_loss_evidence"])
    assert np.isnan(features["frac_mixed_gain_loss_ge_best_pg"])


def test_route_score_result_score_features_handles_mixed_gain_loss_only_routes():
    routes = RouteEnumerator().enumerate(
        (3, 1),
        (3, 1),
        max_truncal_cn=3,
        max_branch_events=2,
        route_type=RouteType.MIXED_GAIN_LOSS,
    )[:3]
    result = RouteScoreResult(
        routes=tuple(routes),
        scores=np.array([-9.0, -11.0, -13.0]),
    )

    features = result.score_features()

    assert np.isnan(features["best_pg"])
    assert features["best_mixed_gain_loss"] == -9.0
    assert features["n_pg_routes"] == 0
    assert features["n_mixed_gain_loss_routes"] == 3
    assert np.isnan(features["percentile_rank"])
    assert np.isnan(features["score_ratio"])
    assert np.isnan(features["pg_evidence"])


def test_route_score_result_route_score_table_uses_route_keys():
    routes = RouteEnumerator().enumerate(
        (3, 1),
        (3, 1),
        max_truncal_cn=3,
        max_branch_events=2,
        route_type=RouteType.PURE_GAIN,
    )
    scores = np.array([-10.0, -12.0])
    result = RouteScoreResult(routes=tuple(routes), scores=scores)

    table = result.route_score_table(seed=7)

    assert table["truncal"] == [route.truncal_key for route in routes]
    assert table["branch1"] == [route.branch1.canonical_key for route in routes]
    assert table["branch2"] == [route.branch2.canonical_key for route in routes]
    assert np.array_equal(table["score_seed_7"], scores)
