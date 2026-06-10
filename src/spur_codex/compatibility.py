from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.optimize import nnls

from .encoding import RouteEncoder
from .enumeration import RouteEnumerator
from .route import Route
from .route_filter import RouteFilterName, apply_route_filter, split_routes


@dataclass(frozen=True)
class PureGainFeasibilityCandidate:
    constraint_matrix: np.ndarray
    h_matrix: np.ndarray
    support: np.ndarray
    route: Route


class PureGainCompatibilityChecker:
    def __init__(
        self,
        max_truncal_cn: int | None = None,
        route_filter: RouteFilterName = "none",
        route_filter_threshold: int | None = None,
        atol: float = 1e-8,
    ):
        if max_truncal_cn is None:
            raise TypeError("max_truncal_cn is required.")
        self.max_truncal_cn = max_truncal_cn
        self.route_filter = route_filter
        self.route_filter_threshold = route_filter_threshold
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
        H_transposed = candidate.h_matrix.T
        test_matrix = np.vstack(
            (
                candidate.constraint_matrix,
                H_transposed - np.outer(q, np.ones(n_timings)),
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
            max_truncal_cn=self.max_truncal_cn,
            max_branch_events=self._enumeration_branch_event_cap(sample1_cn, sample2_cn),
        )
        filtered_routes = apply_route_filter(
            all_routes,
            route_filter=self.route_filter,
            threshold=self.route_filter_threshold,
        )
        routes = split_routes(filtered_routes).pure_gain

        candidates = []
        n_joint_groups = (sample1_cn[0] - sample1_cn[1]) * (sample2_cn[0] - sample2_cn[1])
        for route in routes:
            route_constraint = route.truncal.constraint_matrix
            H = self.encoder.h_tensor([route])
            if H is None:
                continue
            h_matrix = H[0]
            if h_matrix.shape[1] != n_joint_groups:
                raise ValueError("Unexpected pure-gain H shape.")
            candidates.append(
                PureGainFeasibilityCandidate(
                    constraint_matrix=route_constraint,
                    h_matrix=h_matrix,
                    support=h_matrix.sum(axis=0) > 0,
                    route=route,
                )
            )
        return tuple(candidates)

    def _enumeration_branch_event_cap(self, sample1_cn: tuple[int, int], sample2_cn: tuple[int, int]) -> int:
        if self.route_filter != "none":
            if self.route_filter_threshold is None:
                raise ValueError("--route-filter-threshold is required when --route-filter is not 'none'.")
            return self.route_filter_threshold
        return self.max_truncal_cn + max(sample1_cn[0], sample2_cn[0])
