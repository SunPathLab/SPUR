import numpy as np
from spur_codex.encoding import RouteEncoder, build_paired_h
from spur_codex.enumeration import RouteEnumerator, TruncalEnumerator


def test_build_paired_h_uses_sample1_major_order():
    matrix1 = np.array([[1, 1, 2, 2]])
    matrix2 = np.array([[1, 2, 1, 2]])

    H = build_paired_h(matrix1, matrix2, sample1_major_cn=2, sample2_major_cn=2)

    assert H.shape == (1, 4, 4)
    assert np.array_equal(np.nonzero(H[0])[1], np.array([0, 1, 2, 3]))


def test_route_encoder_builds_paired_route_tensor_for_one_truncal():
    truncal = TruncalEnumerator().enumerate_for_state(2, 1)[0]
    routes = RouteEnumerator().enumerate_for_truncal(
        truncal=truncal,
        sample1_cn=(2, 1),
        sample2_cn=(2, 1),
        max_branch_events=3,
    )

    encoder = RouteEncoder()
    matrix1, matrix2 = encoder.paired_multiplicity_matrices(routes)
    H = encoder.h_tensor(routes, sample1_major_cn=2, sample2_major_cn=2)

    assert matrix1.shape == matrix2.shape
    assert matrix1.shape[0] == len(routes)
    assert H.shape == (len(routes), matrix1.shape[1], 4)


def test_route_encoder_rejects_mixed_truncals():
    routes = RouteEnumerator().enumerate((2, 1), (2, 1), max_mrca_cn=4, max_branch_events=2)
    truncal_keys = {route.truncal_key for route in routes}
    assert len(truncal_keys) > 1

    try:
        RouteEncoder().paired_multiplicity_matrices(routes)
    except ValueError as exc:
        assert "sharing one truncal" in str(exc)
    else:
        raise AssertionError("Expected ValueError for routes from mixed truncal trees")
