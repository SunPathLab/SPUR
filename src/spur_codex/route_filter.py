from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Callable, Iterable, Literal, Sequence

import numpy as np

from .encoding import RouteEncoder
from .route import Route, RouteType


RoutePredicate = Callable[[Route], bool]
RouteFilterName = Literal[
    "none",
    "total_events",
    "max_branch_events",
    "total_branch_events",
    "single_sample_total_events",
]


class RouteFilter:
    @staticmethod
    def filter_by_total_events(routes: tuple[Route, ...], threshold: int) -> tuple[Route, ...]:
        return tuple(route for route in routes if route.total_events <= threshold)

    @staticmethod
    def filter_by_max_branch_events(routes: tuple[Route, ...], threshold: int) -> tuple[Route, ...]:
        return tuple(route for route in routes if np.max((route.branch1.n_events, route.branch2.n_events)) <= threshold)

    @staticmethod
    def filter_by_total_branch_events(routes: tuple[Route, ...], threshold: int) -> tuple[Route, ...]:
        return tuple(route for route in routes if route.branch1.n_events + route.branch2.n_events <= threshold)
    
    @staticmethod
    def filter_by_single_sample_total_events(routes: tuple[Route, ...], threshold: int) -> tuple[Route, ...]:
        return tuple(route for route in routes if np.max((route.truncal.n_events + route.branch1.n_events,
                                                          route.truncal.n_events + route.branch2.n_events)) <= threshold)




def apply_route_filter(
    routes: Iterable[Route],
    route_filter: RouteFilterName = "none",
    threshold: int | None = None,
) -> tuple[Route, ...]:
    routes = tuple(routes)
    if route_filter == "none":
        return routes
    if threshold is None:
        raise ValueError("--route-filter-threshold is required when --route-filter is not 'none'.")

    route_filter_map = {
        "total_events": RouteFilter.filter_by_total_events,
        "max_branch_events": RouteFilter.filter_by_max_branch_events,
        "total_branch_events": RouteFilter.filter_by_total_branch_events,
        "single_sample_total_events": RouteFilter.filter_by_single_sample_total_events,
    }

    try:
        filter_method = route_filter_map[route_filter]
    except KeyError as exc:
        raise ValueError(f"Unsupported route filter: {route_filter}") from exc
    return filter_method(routes, threshold)


@dataclass(frozen=True)
class RouteSplit:
    pure_gain: tuple[Route, ...]
    mixed_gain_loss: tuple[Route, ...]


@dataclass(frozen=True)
class RouteGroupByTruncal:
    truncal_key: str
    routes: tuple[Route, ...]

    @property
    def truncal(self):
        return self.routes[0].truncal

    def h_tensor(self, dtype=np.float32) -> np.ndarray:
        tensor = RouteEncoder().h_tensor(
            self.routes,
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


def group_by_truncal(routes: Iterable[Route]) -> tuple[RouteGroupByTruncal, ...]:
    grouped: dict[str, list[Route]] = defaultdict(list)
    for route in routes:
        grouped[route.truncal_key].append(route)
    return tuple(RouteGroupByTruncal(truncal_key=key, routes=tuple(value)) for key, value in grouped.items())


def collect_h_by_truncal(
    routes: Sequence[Route],
    dtype=np.float32,
) -> dict[str, np.ndarray]:
    return {
        group.truncal_key: group.h_tensor(dtype=dtype)
        for group in group_by_truncal(routes)
    }
