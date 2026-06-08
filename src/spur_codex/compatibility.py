from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.optimize import nnls

from .encoding import RouteEncoder
from .enumeration import RouteEnumerator
from .route import Route
from .route_filter import split_routes


@dataclass(frozen=True)
class PureGainFeasibilityCandidate:
    constraint_matrix: np.ndarray
    operator_matrix: np.ndarray
    support: np.ndarray
    route: Route


class PureGainCompatibilityChecker:
    def __init__(self, max_mrca_cn: int, max_branch_events: int, atol: float = 1e-8):
        self.max_mrca_cn = max_mrca_cn
        self.max_branch_events = max_branch_events
        self.atol = atol
        self.encoder = RouteEncoder()

    def is_compatible(self, q: np.ndarray, sample1_cn: tuple[int, int], sample2_cn: tuple[int, int]) -> bool:
        q = np.asarray(q, dtype=float).reshape(-1)
        q_sum = q.sum()
        if q_sum <= 0:
            return False
        q = q / q_sum

        for candidate in self._candidates(sample1_cn, sample2_cn):
            if np.any((q > self.atol) & ~candidate.support):
                continue
            if self._candidate_is_compatible(q, candidate):
                return True
        return False

    def _candidate_is_compatible(self, q: np.ndarray, candidate: PureGainFeasibilityCandidate) -> bool:
        n_timings = candidate.constraint_matrix.shape[1]
        test_matrix = np.vstack(
            (
                candidate.constraint_matrix,
                candidate.operator_matrix - np.outer(q, np.ones(n_timings)),
            )
        )
        test_b = np.hstack(
            (
                np.ones(candidate.constraint_matrix.shape[0]),
                np.zeros_like(q),
            )
        )
        _, residual = nnls(test_matrix, test_b)
        return bool(residual <= self.atol)

    @lru_cache(maxsize=None)
    def _candidates(self, sample1_cn: tuple[int, int], sample2_cn: tuple[int, int]) -> tuple[PureGainFeasibilityCandidate, ...]:
        all_routes = RouteEnumerator().enumerate(
            sample1_cn=sample1_cn,
            sample2_cn=sample2_cn,
            max_mrca_cn=self.max_mrca_cn,
            max_branch_events=self.max_branch_events,
        )
        routes = split_routes(all_routes).pure_gain

        candidates = []
        n_joint_groups = (sample1_cn[0] - sample1_cn[1]) * (sample2_cn[0] - sample2_cn[1])
        for route in routes:
            route_constraint = route.truncal.constraint_matrix
            H = self.encoder.h_tensor([route])
            if H is None:
                continue
            operator = H[0].T
            if operator.shape[0] != n_joint_groups:
                raise ValueError("Unexpected pure-gain operator shape.")
            candidates.append(
                PureGainFeasibilityCandidate(
                    constraint_matrix=route_constraint,
                    operator_matrix=operator,
                    support=operator.sum(axis=1) > 0,
                    route=route,
                )
            )
        return tuple(candidates)
