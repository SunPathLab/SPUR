from __future__ import annotations

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import linprog

from .objects import Truncal
from .tree import constraint_matrix


def sample_timing(truncal: Truncal, rng: np.random.Generator, n_steps: int = 200, step_scale: float = 0.25) -> np.ndarray:
    matrix = constraint_matrix(
        truncal.tree.major_tree,
        truncal.tree.minor_tree,
        truncal.truncal_major_cn,
        truncal.truncal_minor_cn,
    )
    init = _feasible_timing(matrix)
    basis = null_space(matrix)
    if basis.shape[1] == 0:
        return init

    current = init
    for _ in range(n_steps):
        direction = basis @ rng.normal(size=basis.shape[1])
        norm = np.linalg.norm(direction)
        if norm == 0:
            continue
        direction = direction / norm
        lo, hi = _nonnegative_step_bounds(current, direction)
        if not np.isfinite(lo) or not np.isfinite(hi) or lo >= hi:
            continue
        width = hi - lo
        lower = max(lo, -step_scale * width)
        upper = min(hi, step_scale * width)
        current = current + rng.uniform(lower, upper) * direction
        current = np.maximum(current, 0.0)
    return current


def _feasible_timing(matrix: np.ndarray) -> np.ndarray:
    result = linprog(
        c=np.zeros(matrix.shape[1]),
        A_eq=matrix,
        b_eq=np.ones(matrix.shape[0]),
        bounds=[(0, None)] * matrix.shape[1],
        method="highs",
    )
    if not result.success:
        raise ValueError(f"Could not find feasible timing: {result.message}")
    return result.x


def _nonnegative_step_bounds(point: np.ndarray, direction: np.ndarray) -> tuple[float, float]:
    lower = -np.inf
    upper = np.inf
    for value, delta in zip(point, direction):
        if delta > 0:
            lower = max(lower, -value / delta)
        elif delta < 0:
            upper = min(upper, -value / delta)
    return lower, upper
