import numpy as np
from scipy.special import logsumexp
from scipy.stats import binom

from spur.enumeration import RouteEnumerator
from spur.route import RouteType
from spur.route_filter import group_by_truncal
from spur.sampling import JointVAFSampler, TimingSampler
from spur.scoring import Scorer, binom_pmf_tensor, score_thb_shared_t_shared_b
from spur.tree import leaf_to_root_paths


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


def _plain_route_jmd(route, timing):
    """Compute a route's JMD from copy fates and tree paths, without H."""
    major_paths = leaf_to_root_paths(route.truncal.tree.major_tree)
    paths = list(major_paths)
    if route.truncal.tree.minor_tree is not None:
        minor_event_offset = 2 * route.truncal.truncal_major_cn - 1
        paths.extend(
            [event + minor_event_offset for event in path]
            for path in leaf_to_root_paths(route.truncal.tree.minor_tree)
        )

    multiplicity1 = np.zeros(len(timing), dtype=int)
    multiplicity2 = np.zeros(len(timing), dtype=int)
    for copy_index, path in enumerate(paths):
        for event_index in path:
            multiplicity1[event_index] += route.branch1.copy_fate_vector[copy_index]
            multiplicity2[event_index] += route.branch2.copy_fate_vector[copy_index]

    major1 = route.sample1_total_cn - route.sample1_minor_cn
    major2 = route.sample2_total_cn - route.sample2_minor_cn
    q = np.zeros(major1 * major2, dtype=float)
    for event_index, event_weight in enumerate(timing):
        m1 = multiplicity1[event_index]
        m2 = multiplicity2[event_index]
        if 1 <= m1 <= major1 and 1 <= m2 <= major2:
            joint_index = (m1 - 1) * major2 + (m2 - 1)
            q[joint_index] += event_weight

    total = q.sum()
    return q / total if total > 0 else q


def _plain_route_score(route, timing_samples, vaf_samples, alt1, dep1, alt2, dep2):
    """Direct two-sample likelihood calculation, with explicit loops over JMs."""
    n_joint = vaf_samples.shape[1] // 2
    n_variants = len(alt1)
    per_timing_log_likelihoods = []
    for timing_index, timing in enumerate(timing_samples):
        q = _plain_route_jmd(route, timing)
        variant_log_likelihoods = np.empty(n_variants, dtype=float)
        for variant in range(n_variants):
            jm_log_likelihoods = np.empty(n_joint, dtype=float)
            for jm in range(n_joint):
                p1 = vaf_samples[timing_index, jm]
                p2 = vaf_samples[timing_index, n_joint + jm]
                jm_log_likelihoods[jm] = (
                    binom.logpmf(alt1[variant], dep1[variant], p1)
                    + binom.logpmf(alt2[variant], dep2[variant], p2)
                )
            supported = q > 0
            variant_log_likelihoods[variant] = logsumexp(
                np.log(q[supported]) + jm_log_likelihoods[supported]
            )
        per_timing_log_likelihoods.append(variant_log_likelihoods.sum())
    return logsumexp(per_timing_log_likelihoods) - np.log(len(timing_samples))


def test_matrix_scorer_matches_independent_two_sample_route_likelihood(monkeypatch):
    # Use a seeded generator only to create common random draws. Both scoring
    # calculations below receive these exact VAF and timing samples.
    routes = RouteEnumerator().enumerate(
        (4, 2), (3, 1), max_truncal_cn=4, max_branch_events=3,
        route_type=RouteType.MIXED_GAIN_LOSS,
    )
    route_rng = np.random.default_rng(20261103)

    alt1 = np.array([4, 2, 7, 3, 5, 1], dtype=int)
    dep1 = np.array([30, 24, 32, 27, 29, 21], dtype=int)
    alt2 = np.array([6, 3, 2, 8, 4, 5], dtype=int)
    dep2 = np.array([31, 25, 22, 34, 28, 30], dtype=int)

    n_samples = 9
    sample_rng = np.random.default_rng(66254)
    vaf_samples = JointVAFSampler().sample_joint_vafs(
        sample1_total_cn=4, sample1_minor_cn=2,
        sample2_total_cn=3, sample2_minor_cn=1,
        sample1_pu=0.8, sample2_pu=0.8,
        n_samples=n_samples, scale=0.01, rng=sample_rng,
    )
    timing_by_key = {}
    for group in group_by_truncal(routes):
        timing_by_key[group.truncal_key] = TimingSampler().sample_timing_vectors(
            group.truncal, rng=sample_rng, n_samples=n_samples,
            burn_in=10, skips=3,
        )

    # Some enumerated routes can have zero shared-mutation mass for a timing
    # draw; the matrix scorer's normalization is undefined for those rows.
    # Keep randomly selected routes whose shared mass is positive for every
    # common timing draw, which is the domain of the score being compared.
    valid_routes = [
        route for route in routes
        if all(
            _plain_route_jmd(route, timing).sum() > 0
            for timing in timing_by_key[route.truncal_key]
        )
    ]
    assert valid_routes
    selected_routes = tuple(
        valid_routes[i]
        for i in route_rng.choice(len(valid_routes), size=min(7, len(valid_routes)), replace=False)
    )
    assert selected_routes

    # Have Scorer run its ordinary matrix pipeline, but inject the shared draws
    # so the comparison is not dependent on matching RNG call order.
    monkeypatch.setattr(
        "spur.scoring.JointVAFSampler.sample_joint_vafs",
        lambda self, **kwargs: vaf_samples.copy(),
    )
    monkeypatch.setattr(
        "spur.scoring.TimingSampler.sample_timing_vectors",
        lambda self, truncal, **kwargs: timing_by_key[truncal.canonical_key].copy(),
    )
    matrix_result = Scorer().score(
        alt1, dep1, alt2, dep2, selected_routes,
        sample1_pu=0.8, sample2_pu=0.8,
        n_samples=n_samples, scale=0.01, burn_in=10, skips=3,
        rng=np.random.default_rng(1), dtype=np.float64,
    )

    direct_scores = np.array([
        _plain_route_score(
            route,
            timing_by_key[route.truncal_key],
            vaf_samples,
            alt1, dep1, alt2, dep2,
        )
        for route in matrix_result.routes
    ])
    assert np.allclose(matrix_result.scores, direct_scores, rtol=1e-10, atol=1e-10)
