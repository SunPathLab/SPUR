from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import networkx as nx
import numpy as np


@dataclass
class TreeNode:
    left: Optional["TreeNode"] = None
    right: Optional["TreeNode"] = None
    id: Optional[int] = None

    @property
    def is_leaf(self) -> bool:
        return self.left is None and self.right is None


def all_full_binary_trees(n_leaves: int) -> List[TreeNode]:
    if n_leaves <= 0:
        return []
    if n_leaves == 1:
        return [TreeNode()]

    trees: List[TreeNode] = []
    for left_leaves in range(1, n_leaves):
        right_leaves = n_leaves - left_leaves
        if right_leaves < left_leaves:
            continue
        left_trees = all_full_binary_trees(left_leaves)
        right_trees = all_full_binary_trees(right_leaves)
        for i, left in enumerate(left_trees):
            for j, right in enumerate(right_trees):
                if left_leaves == right_leaves and j > i:
                    break
                trees.append(TreeNode(left=left, right=right))
    return trees


def assign_node_ids(root: TreeNode) -> None:
    queue = [root]
    next_id = 0
    while queue:
        node = queue.pop(0)
        if node.left is not None:
            queue.append(node.left)
        if node.right is not None:
            queue.append(node.right)
        node.id = next_id
        next_id += 1


def leaf_to_root_paths(root: TreeNode) -> List[List[int]]:
    if root.id is None:
        raise ValueError("Tree node IDs must be assigned before path extraction.")

    paths: List[List[int]] = []

    def visit(node: TreeNode, path: List[int]) -> None:
        if node.id is None:
            raise ValueError("Tree node IDs must be assigned before path extraction.")
        next_path = path + [node.id]
        if node.is_leaf:
            paths.append(next_path)
            return
        if node.left is not None:
            visit(node.left, next_path)
        if node.right is not None:
            visit(node.right, next_path)

    visit(root, [])
    return paths


def tree_hash(root: TreeNode) -> str:
    if root.id is None:
        raise ValueError("Tree node IDs must be assigned before hashing.")
    graph = nx.Graph()
    stack = [root]
    while stack:
        node = stack.pop()
        graph.add_node(node.id)
        for child in (node.left, node.right):
            if child is not None:
                graph.add_edge(node.id, child.id)
                stack.append(child)
    return nx.weisfeiler_lehman_graph_hash(graph)


def union_tree_hash(major_tree: TreeNode, minor_tree: TreeNode) -> str:
    major_graph = _tree_graph(major_tree)
    minor_graph = _tree_graph(minor_tree)
    union = nx.union(major_graph, minor_graph, rename=("major", "minor"))
    return nx.weisfeiler_lehman_graph_hash(union)


def tree_shape(root: TreeNode) -> str:
    if root.left is None and root.right is None:
        return "L"
    children = []
    if root.left is not None:
        children.append(tree_shape(root.left))
    if root.right is not None:
        children.append(tree_shape(root.right))
    if len(children) != 2:
        raise ValueError("Expected a full binary tree node to have zero or two children.")
    left, right = sorted(children)
    return f"({left},{right})"


def union_tree_shape(major_tree: TreeNode, minor_tree: TreeNode) -> str:
    return f"major={tree_shape(major_tree)};minor={tree_shape(minor_tree)}"


def constraint_matrix(major_tree: TreeNode, minor_tree: TreeNode, major_cn: int, minor_cn: int) -> np.ndarray:
    major_paths = leaf_to_root_paths(major_tree)
    minor_paths = leaf_to_root_paths(minor_tree)

    major_block = np.zeros((major_cn, 2 * major_cn - 1), dtype=float)
    for row, path in enumerate(major_paths):
        for col in path:
            major_block[row, col] = 1.0

    minor_block = np.zeros((minor_cn, 2 * minor_cn - 1), dtype=float)
    for row, path in enumerate(minor_paths):
        for col in path:
            minor_block[row, col] = 1.0

    top_right = np.zeros((major_block.shape[0], minor_block.shape[1]), dtype=float)
    bottom_left = np.zeros((minor_block.shape[0], major_block.shape[1]), dtype=float)
    return np.block([[major_block, top_right], [bottom_left, minor_block]])


def _tree_graph(root: TreeNode) -> nx.Graph:
    if root.id is None:
        raise ValueError("Tree node IDs must be assigned before graph construction.")
    graph = nx.Graph()
    stack = [root]
    while stack:
        node = stack.pop()
        graph.add_node(node.id)
        for child in (node.left, node.right):
            if child is not None:
                graph.add_edge(node.id, child.id)
                stack.append(child)
    return graph
