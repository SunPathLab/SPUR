import numpy as np

from spur.route import GainHistoryTree, Truncal
from spur.tree import (
    TreeNode,
    all_full_binary_trees,
    assign_node_ids,
    leaf_to_root_paths,
    tree_shape,
    union_tree_shape,
)


# Verifies that the generator generates correct number of non-isomorphic full-binary-tree
def test_all_full_binary_trees_returns_unlabeled_shapes():
    counts = {n_leaves: len(all_full_binary_trees(n_leaves)) for n_leaves in range(1, 6)}
    assert counts == {
        1: 1,
        2: 1,
        3: 1,
        4: 2,
        5: 3,
    }


# Verifies that canonical shape strings ignore left/right child ordering.
def test_tree_shape_is_canonical_under_left_right_swap():
    left_heavy = TreeNode(left=TreeNode(left=TreeNode(), right=TreeNode()), right=TreeNode())
    right_heavy = TreeNode(left=TreeNode(), right=TreeNode(left=TreeNode(), right=TreeNode()))

    assert tree_shape(left_heavy) == tree_shape(right_heavy)
    assert tree_shape(TreeNode()) == "L"


# Verifies that shape serialization rejects malformed trees, since SPUR route
# trees should be full binary trees.
def test_tree_shape_rejects_non_full_binary_tree():
    root = TreeNode(left=TreeNode())

    try:
        tree_shape(root)
    except ValueError as exc:
        assert "full binary tree" in str(exc)
    else:
        raise AssertionError("Expected ValueError for a non-full binary tree")


# Verifies the node ID convention used by path extraction and constraint matrix
# construction: breadth-first IDs from the root.
def test_assign_node_ids_uses_breadth_first_order_and_paths():
    root = TreeNode(
        left=TreeNode(),
        right=TreeNode(left=TreeNode(), right=TreeNode()),
    )

    assign_node_ids(root)

    assert root.id == 0
    assert root.left.id == 1
    assert root.right.id == 2
    assert root.right.left.id == 3
    assert root.right.right.id == 4
    assert leaf_to_root_paths(root) == [[0, 1], [0, 2, 3], [0, 2, 4]]


# Verifies that path extraction fails loudly before node IDs are assigned.
def test_leaf_to_root_paths_requires_node_ids():
    try:
        leaf_to_root_paths(TreeNode())
    except ValueError as exc:
        assert "IDs must be assigned" in str(exc)
    else:
        raise AssertionError("Expected ValueError before assigning node IDs")


# Verifies that the truncal constraint matrix keeps major and minor allele tree
# blocks separate and places path indicators in the expected columns.
def test_truncal_constraint_matrix_has_major_and_minor_blocks():
    # Major tree: 2 leaves
    # Minor tree: 1 leaf
    major_tree = TreeNode(left=TreeNode(), right=TreeNode())
    minor_tree = TreeNode()
    assign_node_ids(major_tree)
    assign_node_ids(minor_tree)

    truncal = Truncal(
        truncal_total_cn=3,
        truncal_minor_cn=1,
        tree=GainHistoryTree(major_tree=major_tree, minor_tree=minor_tree),
    )
    matrix = truncal.constraint_matrix

    expected = np.array(
        [
            [1, 1, 0, 0],
            [1, 0, 1, 0],
            [0, 0, 0, 1],
        ]
    )
    assert np.array_equal(matrix, expected)


# Verifies that tree shape strings can identify trees without requiring
# hash-specific node IDs.
def test_union_shape_is_the_tree_identifier_without_assigned_ids():
    major_tree = TreeNode(left=TreeNode(), right=TreeNode())
    minor_tree = TreeNode()

    assert union_tree_shape(major_tree, minor_tree) == "major=(L,L);minor=L"


# Verifies that a truncal LOH history can represent the absent minor tree.
def test_union_shape_represents_empty_minor_tree():
    major_tree = TreeNode(left=TreeNode(), right=TreeNode())

    assert union_tree_shape(major_tree, None) == "major=(L,L);minor=empty"


# Verifies that truncal LOH counts one homolog deletion plus retained-allele gains.
def test_truncal_loh_event_count_includes_homolog_deletion():
    truncal = Truncal(
        truncal_total_cn=2,
        truncal_minor_cn=0,
        tree=GainHistoryTree(
            major_tree=TreeNode(left=TreeNode(), right=TreeNode()),
            minor_tree=None,
        ),
    )

    assert truncal.n_events == 2
