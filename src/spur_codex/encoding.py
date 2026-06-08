from __future__ import annotations

from typing import Sequence

import numpy as np

from .route import Route

class RouteEncoder:
    def _copy_fate_matrices(self, routes: Sequence[Route]) -> tuple[np.ndarray, np.ndarray]:
        if not routes:
            return np.array([]), np.array([])
        branch1 = np.asarray([route.branch1.copy_fate_vector for route in routes], dtype=float)
        branch2 = np.asarray([route.branch2.copy_fate_vector for route in routes], dtype=float)
        return branch1, branch2

    def _multiplicity_matrices(self, routes: Sequence[Route]) -> tuple[np.ndarray, np.ndarray]:
        if not routes:
            return np.array([]), np.array([])
        truncal = routes[0].truncal
        for route in routes:
            if route.truncal_key != truncal.canonical_key:
                raise ValueError("Multiplicity matrices can only be built for routes sharing one truncal tree.")

        matrix = truncal.constraint_matrix
        branch1, branch2 = self._copy_fate_matrices(routes)
        return (matrix.T @ branch1.T).T, (matrix.T @ branch2.T).T

    def h_tensor(self, routes: Sequence[Route], dtype=np.float32) -> np.ndarray | None:
        if not routes:
            return None
        first = routes[0]
        for route in routes:
            if route.sample1_total_cn != first.sample1_total_cn or route.sample1_minor_cn != first.sample1_minor_cn:
                raise ValueError("H tensor routes must share sample 1 total/minor CN.")
            if route.sample2_total_cn != first.sample2_total_cn or route.sample2_minor_cn != first.sample2_minor_cn:
                raise ValueError("H tensor routes must share sample 2 total/minor CN.")
        matrix1, matrix2 = self._multiplicity_matrices(routes)
        return build_paired_h(
            matrix1,
            matrix2,
            sample1_total_cn=first.sample1_total_cn,
            sample1_minor_cn=first.sample1_minor_cn,
            sample2_total_cn=first.sample2_total_cn,
            sample2_minor_cn=first.sample2_minor_cn,
            dtype=dtype,
        )


def build_paired_h(
    matrix1: np.ndarray,
    matrix2: np.ndarray,
    sample1_total_cn: int,
    sample1_minor_cn: int,
    sample2_total_cn: int,
    sample2_minor_cn: int,
    dtype=np.float32,
) -> np.ndarray | None:
    if matrix1.shape[0] == 0 or matrix2.shape[0] == 0:
        return None
    if matrix1.shape != matrix2.shape:
        raise ValueError("Paired H construction requires matrix1 and matrix2 to have the same shape.")

    sample1_major_cn = sample1_total_cn - sample1_minor_cn
    sample2_major_cn = sample2_total_cn - sample2_minor_cn
    n_routes, n_events = matrix1.shape
    n_joint_groups = sample1_major_cn * sample2_major_cn
    H = np.zeros((n_routes, n_events, n_joint_groups), dtype=dtype)

    m1 = matrix1.astype(np.int32, copy=False)
    m2 = matrix2.astype(np.int32, copy=False)
    valid = (m1 > 0) & (m2 > 0)
    joint_idx = (m1 - 1) * sample2_major_cn + (m2 - 1)

    route_idx, event_idx = np.nonzero(valid)
    H[route_idx, event_idx, joint_idx[route_idx, event_idx]] = 1
    return H
