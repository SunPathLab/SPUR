import importlib.util
from pathlib import Path

import numpy as np
import pytest


SCRIPT_PATH = Path(__file__).parents[1] / "scripts" / "simulate_parsimony_route.py"
SPEC = importlib.util.spec_from_file_location("simulate_parsimony_route", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
get_non_wgd_parsimony_route_constraint_matrix = (
    MODULE.get_non_wgd_parsimony_route_constraint_matrix
)
get_non_wgd_parsimony_branches = MODULE.get_non_wgd_parsimony_branches
get_wgd_parsimony_branches = MODULE.get_wgd_parsimony_branches
get_wgd_parsimony_route_constraint_matrix = (
    MODULE.get_wgd_parsimony_route_constraint_matrix
)


# Verifies that each supported non-WGD route uses the truncal CN state
# transcribed from the parsimony-route diagrams.
@pytest.mark.parametrize(
    ("sample_cn", "expected"),
    [
        (((3, 1), (2, 1)), np.eye(2, dtype=int)),
        (
            ((4, 2), (3, 1)),
            np.array(
                [
                    [1, 1, 0, 0],
                    [1, 0, 1, 0],
                    [0, 0, 0, 1],
                ]
            ),
        ),
        (
            ((4, 1), (3, 1)),
            np.array(
                [
                    [1, 1, 0, 0],
                    [1, 0, 1, 0],
                    [0, 0, 0, 1],
                ]
            ),
        ),
        (((4, 1), (2, 1)), np.eye(2, dtype=int)),
        (
            ((3, 1), (3, 1)),
            np.array(
                [
                    [1, 1, 0, 0],
                    [1, 0, 1, 0],
                    [0, 0, 0, 1],
                ]
            ),
        ),
        (
            ((4, 2), (4, 2)),
            np.array(
                [
                    [1, 1, 0, 0, 0, 0],
                    [1, 0, 1, 0, 0, 0],
                    [0, 0, 0, 1, 1, 0],
                    [0, 0, 0, 1, 0, 1],
                ]
            ),
        ),
    ],
)
def test_non_wgd_parsimony_constraint_matrix(sample_cn, expected):
    sample1_cn, sample2_cn = sample_cn

    actual = get_non_wgd_parsimony_route_constraint_matrix(
        sample1_cn[0],
        sample1_cn[1],
        sample2_cn[0],
        sample2_cn[1],
    )

    assert np.array_equal(actual, expected)


# Verifies that minimum-event branch enumeration reproduces the fate vectors
# drawn for each supported non-WGD parsimony route.
@pytest.mark.parametrize(
    ("sample_cn", "truncal_cn", "expected_branch1", "expected_branch2"),
    [
        (((3, 1), (2, 1)), (2, 1), {(2, 1)}, {(1, 1)}),
        (((4, 2), (3, 1)), (3, 1), {(1, 1, 2)}, {(1, 1, 1)}),
        (
            ((4, 1), (3, 1)),
            (3, 1),
            {(2, 1, 1), (1, 2, 1)},
            {(1, 1, 1)},
        ),
        (((4, 1), (2, 1)), (2, 1), {(3, 1)}, {(1, 1)}),
        (((3, 1), (3, 1)), (3, 1), {(1, 1, 1)}, {(1, 1, 1)}),
        (((4, 2), (4, 2)), (4, 2), {(1, 1, 1, 1)}, {(1, 1, 1, 1)}),
    ],
)
def test_non_wgd_parsimony_branches(
    sample_cn,
    truncal_cn,
    expected_branch1,
    expected_branch2,
):
    sample1_cn, sample2_cn = sample_cn
    truncal, sample1_branches, sample2_branches = get_non_wgd_parsimony_branches(
        sample1_cn[0],
        sample1_cn[1],
        sample2_cn[0],
        sample2_cn[1],
    )

    assert (truncal.truncal_total_cn, truncal.truncal_minor_cn) == truncal_cn
    assert {branch.copy_fate_vector for branch in sample1_branches} == expected_branch1
    assert {branch.copy_fate_vector for branch in sample2_branches} == expected_branch2


# Verifies that WGD routes use the same minimum-event branch enumeration after
# selecting their WGD-derived truncal state.
@pytest.mark.parametrize(
    ("sample_cn", "truncal_cn", "expected_branch1", "expected_branch2"),
    [
        (
            ((3, 1), (2, 1)),
            (3, 1),
            {(1, 1, 1)},
            {(1, 0, 1), (0, 1, 1)},
        ),
        (
            ((4, 2), (3, 1)),
            (4, 2),
            {(1, 1, 1, 1)},
            {(1, 1, 1, 0), (1, 1, 0, 1)},
        ),
        (
            ((5, 2), (4, 2)),
            (4, 2),
            {(2, 1, 1, 1), (1, 2, 1, 1)},
            {(1, 1, 1, 1)},
        ),
        (
            ((5, 2), (3, 1)),
            (4, 2),
            {(2, 1, 1, 1), (1, 2, 1, 1)},
            {(1, 1, 1, 0), (1, 1, 0, 1)},
        ),
        (
            ((4, 1), (3, 1)),
            (3, 1),
            {(2, 1, 1), (1, 2, 1)},
            {(1, 1, 1)},
        ),
        (((3, 1), (3, 1)), (3, 1), {(1, 1, 1)}, {(1, 1, 1)}),
        (((4, 2), (4, 2)), (4, 2), {(1, 1, 1, 1)}, {(1, 1, 1, 1)}),
        (((4, 1), (4, 1)), (4, 1), {(1, 1, 1, 1)}, {(1, 1, 1, 1)}),
        (
            ((5, 2), (5, 2)),
            (5, 2),
            {(1, 1, 1, 1, 1)},
            {(1, 1, 1, 1, 1)},
        ),
    ],
)
def test_wgd_parsimony_branches(
    sample_cn,
    truncal_cn,
    expected_branch1,
    expected_branch2,
):
    sample1_cn, sample2_cn = sample_cn
    truncal, sample1_branches, sample2_branches = get_wgd_parsimony_branches(
        sample1_cn[0],
        sample1_cn[1],
        sample2_cn[0],
        sample2_cn[1],
    )

    assert (truncal.truncal_total_cn, truncal.truncal_minor_cn) == truncal_cn
    assert {branch.copy_fate_vector for branch in sample1_branches} == expected_branch1
    assert {branch.copy_fate_vector for branch in sample2_branches} == expected_branch2


# Verifies the additional root-timing equality for visible WGD truncal states.
@pytest.mark.parametrize(
    ("sample_cn", "expected_wgd_row"),
    [
        (((4, 2), (3, 1)), np.array([1, 0, 0, -1, 0, 0])),
        (((5, 2), (5, 2)), np.array([1, 0, 0, 0, 0, -1, 0, 0])),
    ],
)
def test_wgd_parsimony_constraint_matrix_adds_root_equality(
    sample_cn,
    expected_wgd_row,
):
    sample1_cn, sample2_cn = sample_cn
    matrix = get_wgd_parsimony_route_constraint_matrix(
        sample1_cn[0],
        sample1_cn[1],
        sample2_cn[0],
        sample2_cn[1],
    )

    assert np.array_equal(matrix[-1], expected_wgd_row)
