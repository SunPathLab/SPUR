from __future__ import annotations

import numpy as np

from .encoding import RouteEncoder
from .objects import Route
from .timing import sample_timing


def route_jmd(route: Route, sample1_cn: tuple[int, int], sample2_cn: tuple[int, int], timing: np.ndarray) -> np.ndarray:
    H = RouteEncoder().h_tensor([route], sample1_cn[0], sample2_cn[0])
    if H is None:
        raise ValueError("Route produced no timing contribution tensor.")
    q = timing @ H[0]
    q_sum = q.sum()
    if q_sum <= 0:
        raise ValueError("Route/timing produced zero shared-SSNV JMD mass.")
    return q / q_sum


def sample_route_jmd(route: Route, sample1_cn: tuple[int, int], sample2_cn: tuple[int, int], rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    timing = sample_timing(route.truncal, rng)
    q = route_jmd(route, sample1_cn, sample2_cn, timing)
    return timing, q
