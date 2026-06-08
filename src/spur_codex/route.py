from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Tuple

import numpy as np

from .tree import TreeNode, leaf_to_root_paths, union_tree_shape


# route type
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
    truncal_total_cn: int
    truncal_minor_cn: int
    tree: GainHistoryTree

    @property
    def truncal_major_cn(self) -> int:
        return self.truncal_total_cn - self.truncal_minor_cn

    @property
    def n_events(self) -> int:
        return (self.truncal_major_cn - 1) + (self.truncal_minor_cn - 1)

    @property
    def canonical_key(self) -> str:
        shape = union_tree_shape(self.tree.major_tree, self.tree.minor_tree)
        return f"truncal_{self.truncal_total_cn}:{self.truncal_minor_cn}_{shape}"
    
    @property
    def constraint_matrix(self) -> np.ndarray:
        major_paths = leaf_to_root_paths(self.tree.major_tree)
        minor_paths = leaf_to_root_paths(self.tree.minor_tree)

        major_block = np.zeros((self.truncal_major_cn, 2 * self.truncal_major_cn - 1), dtype=int)
        for row, path in enumerate(major_paths):
            for col in path:
                major_block[row, col] = 1

        minor_block = np.zeros((self.truncal_minor_cn, 2 * self.truncal_minor_cn - 1), dtype=int)
        for row, path in enumerate(minor_paths):
            for col in path:
                minor_block[row, col] = 1

        top_right = np.zeros((major_block.shape[0], minor_block.shape[1]), dtype=int)
        bottom_left = np.zeros((minor_block.shape[0], major_block.shape[1]), dtype=int)
        return np.block([[major_block, top_right], [bottom_left, minor_block]])


@dataclass(frozen=True)
class Branch:
    copy_fate_vector: Tuple[int, ...]
    branch_index: int
    sample_total_cn: int
    sample_minor_cn: int

    @property
    def has_deletion(self) -> bool:
        return any(copy_count == 0 for copy_count in self.copy_fate_vector)

    @property
    def is_pure_gain(self) -> bool:
        return not self.has_deletion

    @property
    def n_events(self) -> int:
        return sum(abs(copy_count - 1) for copy_count in self.copy_fate_vector)

    @property
    def canonical_key(self) -> str:
        return ",".join(map(str, self.copy_fate_vector))
    

@dataclass(frozen=True)
class Route:
    truncal: Truncal
    branch1: Branch
    branch2: Branch
    
    @property
    def truncal_total_cn(self) -> int:
        return self.truncal.truncal_total_cn
    
    @property
    def truncal_minor_cn(self) -> int:
        return self.truncal.truncal_minor_cn
    
    @property
    def sample1_total_cn(self) -> int:
        return self.branch1.sample_total_cn
    
    @property
    def sample1_minor_cn(self) -> int:
        return self.branch1.sample_minor_cn
    
    @property
    def sample2_total_cn(self) -> int:
        return self.branch2.sample_total_cn
    
    @property
    def sample2_minor_cn(self) -> int:
        return self.branch2.sample_minor_cn

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

    @property
    def truncal_key(self) -> str:
        return self.truncal.canonical_key
    
    @property
    def total_events(self) -> int:
        return self.truncal.n_events + self.branch1.n_events + self.branch2.n_events
    

    def h(self, dtype=np.float32) -> np.ndarray:
        from .encoding import RouteEncoder

        tensor = RouteEncoder().h_tensor(
            [self],
            dtype=dtype,
        )
        if tensor is None:
            raise ValueError("Route produced no H tensor.")
        return tensor[0]
