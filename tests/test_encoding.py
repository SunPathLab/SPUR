import numpy as np
from spur.encoding import RouteEncoder, build_paired_h
from spur.enumeration import RouteEnumerator, TruncalEnumerator


def test_build_paired_h_uses_sample1_total_minor_order():
    matrix1 = np.array([[1, 1, 2, 2]])
    matrix2 = np.array([[1, 2, 1, 2]])

    H = build_paired_h(
        matrix1,
        matrix2,
        sample1_total_cn=4,
        sample1_minor_cn=2,
        sample2_total_cn=4,
        sample2_minor_cn=2,
    )

    assert H.shape == (1, 4, 4)
    assert np.array_equal(np.nonzero(H[0])[1], np.array([0, 1, 2, 3]))


def test_route_encoder_builds_paired_route_tensor_for_one_truncal():
    truncal = TruncalEnumerator().enumerate_for_state(3, 1)[0]
    routes = RouteEnumerator().enumerate_for_truncal(
        truncal=truncal,
        sample1_cn=(3, 1),
        sample2_cn=(3, 1),
        max_branch_events=3,
    )

    encoder = RouteEncoder()
    H = encoder.h_tensor(routes)

    assert H.shape[0] == len(routes)
    assert H.shape[2] == 4


def test_route_encoder_rejects_mixed_truncals():
    routes = RouteEnumerator().enumerate((3, 1), (3, 1), max_truncal_cn=4, max_branch_events=2)
    truncal_keys = {route.truncal_key for route in routes}
    assert len(truncal_keys) > 1

    try:
        RouteEncoder().h_tensor(routes)
    except ValueError as exc:
        assert "sharing one truncal" in str(exc)
    else:
        raise AssertionError("Expected ValueError for routes from mixed truncal trees")
