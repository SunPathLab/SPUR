from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Tuple

from .tree import TreeNode


class RouteType(Enum):
    ALL = "all"
    PURE_GAIN = "pure_gain"
    MIXED_GAIN_LOSS = "mixed_gain_loss"


@dataclass(frozen=True)
class GainHistoryTree:
    major_tree: TreeNode
    minor_tree: TreeNode


@dataclass(frozen=True)
class Truncal:
    major_cn: int
    minor_cn: int
    tree: GainHistoryTree
    tree_hash: str

    @property
    def n_events(self) -> int:
        return (self.major_cn - 1) + (self.minor_cn - 1)


@dataclass(frozen=True)
class Branch:
    copy_fate_vector: Tuple[int, ...]
    branch_index: int

    @property
    def has_deletion(self) -> bool:
        return any(copy_count == 0 for copy_count in self.copy_fate_vector)

    @property
    def is_pure_gain(self) -> bool:
        return not self.has_deletion

    @property
    def n_events(self) -> int:
        return sum(abs(copy_count - 1) for copy_count in self.copy_fate_vector)


@dataclass(frozen=True)
class Route:
    truncal: Truncal
    branch1: Branch
    branch2: Branch

    @property
    def is_mixed_gain_loss(self) -> bool:
        return self.branch1.has_deletion or self.branch2.has_deletion

    @property
    def is_pure_gain(self) -> bool:
        return not self.is_mixed_gain_loss

    def matches_type(self, route_type: RouteType) -> bool:
        if route_type == RouteType.ALL:
            return True
        if route_type == RouteType.PURE_GAIN:
            return self.is_pure_gain
        if route_type == RouteType.MIXED_GAIN_LOSS:
            return self.is_mixed_gain_loss
        raise ValueError(f"Unsupported route type: {route_type}")
