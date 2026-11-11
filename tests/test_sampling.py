import numpy as np

from spur.sampling import JointVAFSampler, vafEst


def test_vafEst_matches_old_spur_formula_for_all_multiplicities():
    vafs = vafEst(nt=4, nb=2, pu=0.8)

    expected = [
        0.8 * 1 / (0.8 * 4 + 2 * (1 - 0.8)),
        0.8 * 2 / (0.8 * 4 + 2 * (1 - 0.8)),
    ]
    assert np.allclose(vafs, expected)


def test_vafEst_supports_single_multiplicity():
    assert np.isclose(vafEst(nt=4, nb=2, pu=0.8, ncmut=2), 0.8 * 2 / (0.8 * 4 + 2 * (1 - 0.8)))


def test_joint_vaf_sampler_uses_total_minor_inputs():
    samples = JointVAFSampler().sample_joint_vafs(
        3,
        1,
        3,
        1,
        sample1_pu=0.8,
        sample2_pu=0.8,
        n_samples=2,
        rng=np.random.default_rng(1),
    )
    assert samples.shape == (2, 8)
