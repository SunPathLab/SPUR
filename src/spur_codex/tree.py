from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

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
