from __future__ import annotations

from itertools import combinations_with_replacement
from typing import Iterable, List, Tuple

from .route import Branch, GainHistoryTree, Route, RouteType, Truncal
from .tree import all_full_binary_trees, assign_node_ids


class TruncalEnumerator:
    def enumerate_for_state(self, total_cn: int, minor_cn: int) -> List[Truncal]:
        if total_cn < 0 or minor_cn < 0:
            raise ValueError("Truncal copy numbers must be nonnegative.")
        major_cn = total_cn - minor_cn
        if major_cn < minor_cn:
            raise ValueError("Minor copy number cannot exceed major copy number.")

        major_trees = all_full_binary_trees(major_cn)
        
        # LOH or nonLOH
        minor_trees = all_full_binary_trees(minor_cn) if minor_cn > 0 else [None]
        
        # assign node IDs to all trees, so that we can extract paths from leaves to root later
        for tree in major_trees + minor_trees:
            if tree is not None:
                assign_node_ids(tree)

        # for each (major_tree,minor_tree),
        # construct a gain history tree,
        # (Computation note: since the object is immutable, no need to copy the tree for each truncal.)
        # Then for each gain history tree, construct a truncal object
        truncals: List[Truncal] = []
        for major_tree in major_trees:
            for minor_tree in minor_trees:
                gain_tree = GainHistoryTree(major_tree=major_tree, minor_tree=minor_tree)
                truncals.append(
                    Truncal(
                        truncal_total_cn=total_cn,
                        truncal_minor_cn=minor_cn,
                        tree=gain_tree,
                    )
                )
        return truncals

    def enumerate_up_to(self, max_truncal_cn: int) -> List[Truncal]:
        truncals: List[Truncal] = []
        for total_cn in range(2, max_truncal_cn + 1):
            for minor_cn in range(1, total_cn // 2 + 1):
                truncals.extend(self.enumerate_for_state(total_cn, minor_cn))
        return truncals

    def enumerate_loh_up_to(self, max_truncal_cn: int) -> List[Truncal]:
        truncals: List[Truncal] = []
        for total_cn in range(1, max_truncal_cn + 1):
            truncals.extend(self.enumerate_for_state(total_cn, minor_cn=0))
        return truncals


class BranchEnumerator:
    def enumerate(self, truncal: Truncal, target_total_cn: int, target_minor_cn: int, max_events: int) -> List[Branch]:
        if target_total_cn < 0 or target_minor_cn < 0:
            raise ValueError("Target copy numbers must be nonnegative.")
        target_major_cn = target_total_cn - target_minor_cn
        if target_major_cn < target_minor_cn:
            raise ValueError("Minor copy number cannot exceed major copy number.")

        major_vectors = self._allele_vectors(truncal.truncal_major_cn, target_major_cn)
        minor_vectors = self._allele_vectors(truncal.truncal_minor_cn, target_minor_cn)

        branches: List[Branch] = []
        next_index = 0
        for major_vector in major_vectors:
            for minor_vector in minor_vectors:
                # combine the major and minor vectors to form the copy fate vector
                copy_fate_vector = major_vector + minor_vector
                
                branch = Branch(
                    copy_fate_vector=copy_fate_vector,
                    branch_index=next_index,
                    sample_total_cn=target_total_cn,
                    sample_minor_cn=target_minor_cn,
                )
                
                next_index += 1
                
                if branch.n_events <= max_events:
                    branches.append(branch)
        return branches

    @staticmethod
    def _allele_vectors(n_truncal_copies: int, target_cn: int) -> List[Tuple[int, ...]]:
        if n_truncal_copies == 0:
            return [()] if target_cn == 0 else []
        if target_cn == 0:
            return [tuple(0 for _ in range(n_truncal_copies))]

        vectors: List[Tuple[int, ...]] = []
        for selected in combinations_with_replacement(range(n_truncal_copies), target_cn):
            counts = [0] * n_truncal_copies
            for idx in selected:
                counts[idx] += 1
            vectors.append(tuple(counts))
        return vectors


class RouteEnumerator:
    def __init__(self, truncal_enumerator: TruncalEnumerator | None = None, branch_enumerator: BranchEnumerator | None = None):
        self.truncal_enumerator = truncal_enumerator or TruncalEnumerator()
        self.branch_enumerator = branch_enumerator or BranchEnumerator()

    def enumerate_for_truncal(
        self,
        truncal: Truncal,
        sample1_cn: tuple[int, int],
        sample2_cn: tuple[int, int],
        max_branch_events: int,
        route_type: RouteType = RouteType.ALL,
    ) -> List[Route]:
        branch1_routes = self.branch_enumerator.enumerate(truncal, sample1_cn[0], sample1_cn[1], max_branch_events)
        branch2_routes = self.branch_enumerator.enumerate(truncal, sample2_cn[0], sample2_cn[1], max_branch_events)

        routes: List[Route] = []
        for branch1 in branch1_routes:
            for branch2 in branch2_routes:
                route = Route(truncal=truncal, branch1=branch1, branch2=branch2)
                if route.matches_type(route_type):
                    routes.append(route)
        return routes

    def enumerate(
        self,
        sample1_cn: tuple[int, int],
        sample2_cn: tuple[int, int],
        max_truncal_cn: int | None = None,
        max_branch_events: int = 5,
        route_type: RouteType = RouteType.ALL,
    ) -> List[Route]:
        if max_truncal_cn is None:
            raise TypeError("max_truncal_cn is required.")

        if sample1_cn[1] == 0 and sample2_cn[1] == 0:
            truncals = self.truncal_enumerator.enumerate_loh_up_to(max_truncal_cn)
        else:
            truncals = self.truncal_enumerator.enumerate_up_to(max_truncal_cn)

        routes: List[Route] = []
        for truncal in truncals:
            routes.extend(
                self.enumerate_for_truncal(
                    truncal=truncal,
                    sample1_cn=sample1_cn,
                    sample2_cn=sample2_cn,
                    max_branch_events=max_branch_events,
                    route_type=route_type,
                )
            )
        return routes
