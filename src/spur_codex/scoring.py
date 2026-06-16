from __future__ import annotations
from typing import Sequence

import numpy as np
from scipy.special import gammaln, logsumexp, xlog1py, xlogy

from spur_codex.route import Route
from spur_codex.route_filter import group_by_truncal
from spur_codex.sampling import TimingSampler, JointVAFSampler
from spur_codex.features import RouteScoreResult


class Scorer:
    def score(self, alt1, dep1, alt2, dep2,
              routes: Sequence[Route],
              sample1_pu: float, sample2_pu: float, 
              sample1_prevalence: float = 1, sample2_prevalence: float = 1, n_samples: int = 1000, scale: float = 0.1, 
              burn_in: int = 200, skips: int = 50, 
              rng: np.random.Generator = None,
              l_chunk: int = 256, r_chunk: int = 1024, eps: float = 1e-30,
              dtype=np.float64) -> RouteScoreResult:
        dtype = np.dtype(dtype)
        if not routes:
            return RouteScoreResult(routes=tuple(), scores=np.asarray([], dtype=dtype))
        
        
        
        # First get B, which is the product of binomial PMFs for the observed alt/depth in the two samples, across a grid of joint multiplicity groups
        joint_vafs_samples = JointVAFSampler().sample_joint_vafs(
            sample1_total_cn=routes[0].sample1_total_cn,
            sample1_minor_cn=routes[0].sample1_minor_cn,
            sample2_total_cn=routes[0].sample2_total_cn,
            sample2_minor_cn=routes[0].sample2_minor_cn,
            sample1_pu=sample1_pu,
            sample2_pu=sample2_pu,
            sample1_prevalence=sample1_prevalence,
            sample2_prevalence=sample2_prevalence,
            multiplicity="all",
            n_samples=n_samples,
            scale=scale,
            rng=rng
        )
        
        B = binom_pmf_tensor(alt1, dep1, alt2, dep2, joint_vafs_samples)
        
        # first split routes into groups according to their truncal tree
        groups = group_by_truncal(routes)
        
        # For each group of routes, get the shared T, H, and score the routes
        all_routes = []
        all_scores = []
        for group in groups:
            # first get T, which is the timing samples for the truncal tree shared by the group
            T = TimingSampler().sample_timing_vectors(group.truncal, 
                                                      rng=rng, 
                                                      n_samples=n_samples, 
                                                      burn_in=burn_in, 
                                                      skips=skips)
            
            # then get H, which is the timing-joint multiplicity contribution tensor
            H = group.h_tensor()
            
            scores = score_thb_shared_t_shared_b(T, H, B, l_chunk=l_chunk, r_chunk=r_chunk, eps=eps, dtype=dtype)
            all_routes.extend(group.routes)
            all_scores.extend(scores)
        return RouteScoreResult(routes=tuple(all_routes), scores=np.asarray(all_scores, dtype=dtype))


def score_thb_shared_t_shared_b(T: np.ndarray, H: np.ndarray, B: np.ndarray, l_chunk: int = 256, r_chunk: int = 1024, eps: float = 1e-30, dtype=np.float64) -> np.ndarray:
    dtype = np.dtype(dtype)
    T = np.asarray(T, dtype=dtype)
    H = np.asarray(H, dtype=dtype)
    B = np.asarray(B, dtype=dtype)
    eps = max(eps, np.finfo(dtype).tiny)

    n_timing, n_events = T.shape
    n_routes, n_events_h, n_joint_groups = H.shape
    n_timing_b, _, n_joint_groups_b = B.shape
    if n_events != n_events_h:
        raise ValueError("T and H disagree on the number of timing/event columns.")
    if n_timing != n_timing_b:
        raise ValueError("T and B disagree on the number of timing samples.")
    if n_joint_groups != n_joint_groups_b:
        raise ValueError("H and B disagree on the number of joint multiplicity groups.")

    out = np.full(n_routes, -np.inf, dtype=dtype)
    for route_start in range(0, n_routes, r_chunk):
        route_end = min(route_start + r_chunk, n_routes)
        for timing_start in range(0, n_timing, l_chunk):
            timing_end = min(timing_start + l_chunk, n_timing)
            Q = np.einsum("le,reh->lhr", T[timing_start:timing_end], H[route_start:route_end])
            Q = Q / Q.sum(axis=1, keepdims=True)
            route_likelihood = np.matmul(B[timing_start:timing_end], Q)
            route_likelihood = np.log(np.maximum(route_likelihood, eps))
            sample_scores = route_likelihood.sum(axis=1)
            out[route_start:route_end] = np.logaddexp(
                out[route_start:route_end],
                logsumexp(sample_scores, axis=0),
            )
    return (out - np.asarray(np.log(n_timing), dtype=dtype)).astype(dtype, copy=False)


def binom_pmf_tensor(alt1, depth1, alt2, depth2, vaf_samples: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    alt1 = np.asarray(alt1, dtype=np.float64)[None, :, None]
    depth1 = np.asarray(depth1, dtype=np.float64)[None, :, None]
    alt2 = np.asarray(alt2, dtype=np.float64)[None, :, None]
    depth2 = np.asarray(depth2, dtype=np.float64)[None, :, None]

    n_joint_groups = vaf_samples.shape[1]//2 # lengt of jm indices = number of jm
    
    p1 = np.clip(vaf_samples[:, None, :n_joint_groups], eps, 1.0 - eps)
    p2 = np.clip(vaf_samples[:, None, n_joint_groups : 2 * n_joint_groups], eps, 1.0 - eps)

    logp1 = _binom_logpmf(alt1, depth1, p1)
    logp2 = _binom_logpmf(alt2, depth2, p2)
    return np.exp(logp1 + logp2)


def _binom_logpmf(k: np.ndarray, n: np.ndarray, p: np.ndarray) -> np.ndarray:
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1) + xlogy(k, p) + xlog1py(n - k, -p)
