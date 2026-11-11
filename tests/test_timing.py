import numpy as np

from spur.enumeration import TruncalEnumerator
from spur.sampling import TimingSampler


def test_hit_and_run_repeats_fixed_timing_when_null_space_is_empty():
    truncal = TruncalEnumerator().enumerate_for_state(2, 1)[0]
    rng = np.random.default_rng(123)

    samples = TimingSampler().sample_timing_vectors(truncal, rng, n_samples=4)

    assert samples.shape == (4, 2)
    assert np.allclose(samples, np.ones((4, 2)))


def test_hit_and_run_timing_samples_satisfy_constraints():
    truncal = TruncalEnumerator().enumerate_for_state(3, 1)[0]
    matrix = truncal.constraint_matrix
    rng = np.random.default_rng(123)

    samples = TimingSampler().sample_timing_vectors(truncal, rng, n_samples=8, burn_in=5, skips=2, dtype=None)

    assert samples.shape == (8, matrix.shape[1])
    assert np.all(samples >= 0)
    assert np.all(samples <= 1)
    assert np.allclose(samples @ matrix.T, np.ones((8, matrix.shape[0])), atol=1e-8)


def test_sample_timing_vectors_can_return_one_hit_and_run_sample():
    truncal = TruncalEnumerator().enumerate_for_state(3, 1)[0]
    matrix = truncal.constraint_matrix
    rng = np.random.default_rng(123)

    timing = TimingSampler().sample_timing_vectors(truncal, rng, n_samples=1, burn_in=5, skips=2)[0]

    assert timing.shape == (matrix.shape[1],)
    assert np.all(timing >= 0)
    assert np.all(timing <= 1)
    assert np.allclose(matrix @ timing, np.ones(matrix.shape[0]), atol=1e-6)


def test_sample_timing_vectors_supports_wgd_root_equality():
    truncal = TruncalEnumerator().enumerate_for_state(4, 2)[0]
    matrix = np.vstack(
        (
            truncal.constraint_matrix,
            np.array([1, 0, 0, -1, 0, 0]),
        )
    )
    values = np.array([1, 1, 1, 1, 0])

    samples = TimingSampler().sample_timing_vectors_from_constraints(
        matrix,
        values,
        np.random.default_rng(123),
        n_samples=8,
        burn_in=5,
        skips=2,
        dtype=None,
    )

    assert np.allclose(samples @ matrix.T, values, atol=1e-8)
