import numpy as np

from spur_codex.enumeration import RouteEnumerator
from spur_codex.route import RouteType
from spur_codex.scoring import Scorer, binom_pmf_tensor, score_thb_shared_t_shared_b


def test_score_thb_shared_t_shared_b_returns_one_score_per_route():
    T = np.array(
        [
            [0.7, 0.3],
            [0.2, 0.8],
        ],
        dtype=float,
    )
    H = np.array(
        [
            [[1.0, 0.0], [0.0, 1.0]],
            [[0.0, 1.0], [1.0, 0.0]],
        ],
        dtype=float,
    )
    B = np.full((2, 3, 2), 0.5, dtype=float)

    scores = score_thb_shared_t_shared_b(T, H, B)

    assert scores.shape == (2,)
    assert np.all(np.isfinite(scores))
    assert np.allclose(scores, np.log(0.5) * 3)


def test_binom_pmf_tensor_shape():
    vaf_samples = np.array(
        [
            [0.2, 0.3, 0.25, 0.35],
            [0.4, 0.5, 0.45, 0.55],
        ]
    )

    tensor = binom_pmf_tensor(
        alt1=[2, 3],
        depth1=[10, 10],
        alt2=[4, 5],
        depth2=[10, 10],
        vaf_samples=vaf_samples,
    )

    assert tensor.shape == (2, 2, 2)
    assert np.all(tensor >= 0)


def test_score_thb_shared_t_shared_b_prefers_likely_route():
    T = np.array([[1.0]], dtype=float)
    H = np.array(
        [
            [[1.0, 0.0]],
            [[0.0, 1.0]],
        ],
        dtype=float,
    )
    B = np.array(
        [
            [
                [0.9, 0.1],
                [0.9, 0.1],
            ]
        ],
        dtype=float,
    )

    log_scores = score_thb_shared_t_shared_b(T, H, B)

    expected_log_scores = np.array([np.log(0.9) * 2, np.log(0.1) * 2])

    assert np.allclose(log_scores, expected_log_scores)
    assert log_scores[0] > log_scores[1]


def test_score_thb_shared_t_shared_b_uses_requested_float_dtype():
    T = np.array([[0.7, 0.3]], dtype=np.float64)
    H = np.array([[[1.0, 0.0], [0.0, 1.0]]], dtype=np.float64)
    B = np.full((1, 2, 2), 0.5, dtype=np.float64)

    scores = score_thb_shared_t_shared_b(T, H, B, dtype=np.float32)

    assert scores.dtype == np.dtype(np.float32)


def test_scorer_score_returns_route_score_result_for_alt_depth_arrays():
    routes = RouteEnumerator().enumerate(
        (3, 1),
        (3, 1),
        max_truncal_cn=3,
        max_branch_events=2,
        route_type=RouteType.PURE_GAIN,
    )

    result = Scorer().score(
        alt1=[2, 3],
        dep1=[10, 10],
        alt2=[2, 3],
        dep2=[10, 10],
        routes=routes,
        sample1_pu=0.8,
        sample2_pu=0.8,
        n_samples=4,
        burn_in=2,
        skips=1,
        rng=np.random.default_rng(1),
        dtype=np.float32,
    )

    assert result.routes == tuple(routes)
    assert result.scores.shape == (len(routes),)
    assert result.scores.dtype == np.dtype(np.float32)
    assert np.all(np.isfinite(result.scores))
