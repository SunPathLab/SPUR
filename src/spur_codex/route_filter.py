from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Callable, Iterable, Sequence

import numpy as np

from .encoding import RouteEncoder
from .objects import Route, RouteType


RoutePredicate = Callable[[Route], bool]


@dataclass(frozen=True)
class RouteSplit:
    pure_gain: tuple[Route, ...]
    mixed_gain_loss: tuple[Route, ...]


@dataclass(frozen=True)
class RouteGroup:
    truncal_key: str
    routes: tuple[Route, ...]

    @property
    def truncal(self):
        return self.routes[0].truncal

    def h_tensor(self, sample1_major_cn: int, sample2_major_cn: int, dtype=np.float32) -> np.ndarray:
        tensor = RouteEncoder().h_tensor(
            self.routes,
            sample1_major_cn=sample1_major_cn,
            sample2_major_cn=sample2_major_cn,
            dtype=dtype,
        )
        if tensor is None:
            raise ValueError(f"Route group {self.truncal_key} produced no H tensor.")
        return tensor


def filter_routes(
    routes: Iterable[Route],
    route_type: RouteType = RouteType.ALL,
    predicate: RoutePredicate | None = None,
) -> tuple[Route, ...]:
    filtered = []
    for route in routes:
        if not route.matches_type(route_type):
            continue
        if predicate is not None and not predicate(route):
            continue
        filtered.append(route)
    return tuple(filtered)


def split_routes(routes: Iterable[Route]) -> RouteSplit:
    pure_gain = []
    mixed_gain_loss = []
    for route in routes:
        if route.is_pure_gain:
            pure_gain.append(route)
        else:
            mixed_gain_loss.append(route)
    return RouteSplit(pure_gain=tuple(pure_gain), mixed_gain_loss=tuple(mixed_gain_loss))


def group_by_truncal(routes: Iterable[Route]) -> tuple[RouteGroup, ...]:
    grouped: dict[str, list[Route]] = defaultdict(list)
    for route in routes:
        grouped[route.truncal_key].append(route)
    return tuple(RouteGroup(truncal_key=key, routes=tuple(value)) for key, value in grouped.items())


def collect_h_by_truncal(
    routes: Sequence[Route],
    sample1_major_cn: int,
    sample2_major_cn: int,
    dtype=np.float32,
) -> dict[str, np.ndarray]:
    return {
        group.truncal_key: group.h_tensor(
            sample1_major_cn=sample1_major_cn,
            sample2_major_cn=sample2_major_cn,
            dtype=dtype,
        )
        for group in group_by_truncal(routes)
    }
