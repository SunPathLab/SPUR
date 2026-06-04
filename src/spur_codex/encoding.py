from __future__ import annotations

from typing import Sequence

import numpy as np

from .objects import Route
from .tree import constraint_matrix


class RouteEncoder:
    def copy_fate_matrices(self, routes: Sequence[Route]) -> tuple[np.ndarray, np.ndarray]:
        if not routes:
            return np.array([]), np.array([])
        branch1 = np.asarray([route.branch1.copy_fate_vector for route in routes], dtype=float)
        branch2 = np.asarray([route.branch2.copy_fate_vector for route in routes], dtype=float)
        return branch1, branch2

    def multiplicity_matrices(self, routes: Sequence[Route]) -> tuple[np.ndarray, np.ndarray]:
        if not routes:
            return np.array([]), np.array([])
        truncal = routes[0].truncal
        for route in routes:
            if route.truncal != truncal:
                raise ValueError("Multiplicity matrices can only be built for routes sharing one truncal tree.")

        matrix = constraint_matrix(
            truncal.tree.major_tree,
            truncal.tree.minor_tree,
            truncal.major_cn,
            truncal.minor_cn,
        )
        branch1, branch2 = self.copy_fate_matrices(routes)
        return (matrix.T @ branch1.T).T, (matrix.T @ branch2.T).T

    def paired_multiplicity_matrices(self, routes: Sequence[Route]) -> tuple[np.ndarray, np.ndarray]:
        return self.multiplicity_matrices(routes)

    def timing_contribution_tensor(self, routes: Sequence[Route], sample1_major_cn: int, sample2_major_cn: int, dtype=np.float32) -> np.ndarray | None:
        matrix1, matrix2 = self.paired_multiplicity_matrices(routes)
        return build_paired_j(matrix1, matrix2, sample1_major_cn, sample2_major_cn, dtype=dtype)


def build_paired_j(matrix1: np.ndarray, matrix2: np.ndarray, sample1_major_cn: int, sample2_major_cn: int, dtype=np.float32) -> np.ndarray | None:
    if matrix1.shape[0] == 0 or matrix2.shape[0] == 0:
        return None
    if matrix1.shape != matrix2.shape:
        raise ValueError("Paired J construction requires matrix1 and matrix2 to have the same shape.")

    n_routes, n_events = matrix1.shape
    n_joint_groups = sample1_major_cn * sample2_major_cn
    J = np.zeros((n_routes, n_events, n_joint_groups), dtype=dtype)

    m1 = matrix1.astype(np.int32, copy=False)
    m2 = matrix2.astype(np.int32, copy=False)
    valid = (m1 > 0) & (m2 > 0)
    joint_idx = (m1 - 1) * sample2_major_cn + (m2 - 1)

    route_idx, event_idx = np.nonzero(valid)
    J[route_idx, event_idx, joint_idx[route_idx, event_idx]] = 1
    return J
