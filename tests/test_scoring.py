import numpy as np

from spur_codex.scoring import THBScorer, binom_pmf_tensor, posterior_probabilities


def test_thb_scorer_returns_one_score_per_route():
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

    scores = THBScorer().score(T, H, B)

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
        n_joint_groups=2,
    )

    assert tensor.shape == (2, 2, 2)
    assert np.all(tensor >= 0)


def test_posterior_calculation_prefers_likely_route():
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

    log_scores = THBScorer().score(T, H, B)
    posterior = posterior_probabilities(log_scores)

    expected_log_scores = np.array([np.log(0.9) * 2, np.log(0.1) * 2])
    expected_posterior = np.array([0.9**2, 0.1**2]) / (0.9**2 + 0.1**2)

    assert np.allclose(log_scores, expected_log_scores)
    assert np.allclose(posterior, expected_posterior)
    assert posterior[0] > posterior[1]
