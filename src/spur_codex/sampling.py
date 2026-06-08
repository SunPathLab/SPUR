from __future__ import annotations

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import nnls
from scipy.stats import truncnorm

from .hit_and_run import hit_and_run
from .route import Truncal



def vafEst(nt, nb, pu, pa=1, ncmut="all"):
    vafs = []
    epu = pa * pu
    if ncmut == "all":
        ncm = list(range(1, nt - nb + 1))
        for i in ncm:
            vafs.append(epu * ncm[i - 1] / (epu * nt + 2 * (1 - epu)))
    else:
        vafs = epu * ncmut / (epu * nt + 2 * (1 - epu))

    return vafs

class TimingSampler:
    def sample_timing_vectors(self,
        truncal: Truncal,
        rng: np.random.Generator,
        n_samples: int = 1000,
        burn_in: int = 200,
        skips: int = 50,
        dtype=np.float32,
    ) -> np.ndarray:
        matrix = truncal.constraint_matrix
        init = self._feasible_timing(matrix)
        basis = null_space(matrix)
        samples = hit_and_run(basis, init, rng, n_samples=n_samples, burn_in=burn_in, skips=skips)
        samples = np.clip(samples, 0.0, 1.0)
        if dtype is not None:
            samples = samples.astype(dtype)
        return samples

    def _feasible_timing(self, matrix: np.ndarray) -> np.ndarray:
        init, residual = nnls(matrix, np.ones(matrix.shape[0]))
        if not np.allclose(matrix @ init, np.ones(matrix.shape[0]), atol=1e-8):
            raise ValueError(f"Could not find feasible timing; NNLS residual={residual}.")
        return init



class JointVAFSampler:
    def sample_joint_vafs(self, sample1_total_cn: int, sample1_minor_cn: int,
                    sample2_total_cn: int, sample2_minor_cn: int, 
                    sample1_pu: float, sample2_pu: float, sample1_prevalence: float = 1, sample2_prevalence: float = 1, multiplicity="all", n_samples: int = 1000, scale: float = 0.1, rng: np.random.Generator = None) -> np.ndarray:
        
        expected_sample1_vafs = vafEst(sample1_total_cn, sample1_minor_cn, sample1_pu, sample1_prevalence, multiplicity)
        expected_sample2_vafs = vafEst(sample2_total_cn, sample2_minor_cn, sample2_pu, sample2_prevalence, multiplicity)
        
        vafs_samples = np.zeros((n_samples, len(expected_sample1_vafs)*len(expected_sample2_vafs)*2))
        n_joint_mult_groups = len(expected_sample1_vafs)*len(expected_sample2_vafs)
        for i in range(len(expected_sample1_vafs)):
            for j in range(len(expected_sample2_vafs)):
                a_trunc = 0
                b_trunc = 1
                a_1, b_1 = (a_trunc - expected_sample1_vafs[i])/scale, (b_trunc-expected_sample1_vafs[i])/scale
                a_2, b_2 = (a_trunc - expected_sample2_vafs[j])/scale, (b_trunc-expected_sample2_vafs[j])/scale
                joint_mult_idx = i*len(expected_sample2_vafs) + j
                vafs_samples[:, joint_mult_idx] = \
                                truncnorm.rvs(a=a_1,b=b_1, loc=expected_sample1_vafs[i], scale=scale, size=n_samples, random_state=rng)
                vafs_samples[:, n_joint_mult_groups + joint_mult_idx] = \
                                truncnorm.rvs(a=a_2,b=b_2, loc=expected_sample2_vafs[j], scale=scale, size=n_samples, random_state=rng)
        
        return vafs_samples
