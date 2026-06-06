from __future__ import annotations

import numpy as np
from scipy.special import gammaln, logsumexp, xlog1py, xlogy


class THBScorer:
    def score(self, T: np.ndarray, H: np.ndarray, B: np.ndarray, l_chunk: int = 256, r_chunk: int = 1024, eps: float = 1e-30) -> np.ndarray:
        return score_thb(T, H, B, l_chunk=l_chunk, r_chunk=r_chunk, eps=eps)

    def posterior_probabilities(self, log_scores: np.ndarray) -> np.ndarray:
        return posterior_probabilities(log_scores)


def score_thb(T: np.ndarray, H: np.ndarray, B: np.ndarray, l_chunk: int = 256, r_chunk: int = 1024, eps: float = 1e-30) -> np.ndarray:
    n_timing, n_events = T.shape
    n_routes, n_events_h, n_joint_groups = H.shape
    n_timing_b, _, n_joint_groups_b = B.shape
    if n_events != n_events_h:
        raise ValueError("T and H disagree on the number of timing/event columns.")
    if n_timing != n_timing_b:
        raise ValueError("T and B disagree on the number of timing samples.")
    if n_joint_groups != n_joint_groups_b:
        raise ValueError("H and B disagree on the number of joint multiplicity groups.")

    out = np.full(n_routes, -np.inf, dtype=float)
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
    return out - np.log(n_timing)


def posterior_probabilities(log_scores: np.ndarray) -> np.ndarray:
    log_scores = np.asarray(log_scores, dtype=float)
    return np.exp(log_scores - logsumexp(log_scores))


def binom_pmf_tensor(alt1, depth1, alt2, depth2, vaf_samples: np.ndarray, n_joint_groups: int, eps: float = 1e-12) -> np.ndarray:
    alt1 = np.asarray(alt1, dtype=np.float64)[None, :, None]
    depth1 = np.asarray(depth1, dtype=np.float64)[None, :, None]
    alt2 = np.asarray(alt2, dtype=np.float64)[None, :, None]
    depth2 = np.asarray(depth2, dtype=np.float64)[None, :, None]

    p1 = np.clip(vaf_samples[:, None, :n_joint_groups], eps, 1.0 - eps)
    p2 = np.clip(vaf_samples[:, None, n_joint_groups : 2 * n_joint_groups], eps, 1.0 - eps)

    logp1 = _binom_logpmf(alt1, depth1, p1)
    logp2 = _binom_logpmf(alt2, depth2, p2)
    return np.exp(logp1 + logp2)


def _binom_logpmf(k: np.ndarray, n: np.ndarray, p: np.ndarray) -> np.ndarray:
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1) + xlogy(k, p) + xlog1py(n - k, -p)
