from spur_codex.route import Branch, Route, RouteType
from spur_codex.enumeration import TruncalEnumerator


def test_branch_properties():
    pure_gain = Branch(copy_fate_vector=(1, 2, 1), branch_index=0, sample_total_cn=4, sample_minor_cn=1)
    with_deletion = Branch(copy_fate_vector=(1, 0, 2), branch_index=1, sample_total_cn=3, sample_minor_cn=1)

    assert pure_gain.is_pure_gain
    assert not pure_gain.has_deletion
    assert pure_gain.n_events == 1

    assert with_deletion.has_deletion
    assert not with_deletion.is_pure_gain
    assert with_deletion.n_events == 2


def test_route_type_properties_are_pair_level():
    truncal = TruncalEnumerator().enumerate_for_state(2, 1)[0]
    pure_branch = Branch(copy_fate_vector=(1, 1, 1), branch_index=0, sample_total_cn=3, sample_minor_cn=1)
    deletion_branch = Branch(copy_fate_vector=(1, 0, 2), branch_index=1, sample_total_cn=3, sample_minor_cn=1)

    pure_route = Route(truncal=truncal, branch1=pure_branch, branch2=pure_branch)
    mixed_route = Route(truncal=truncal, branch1=pure_branch, branch2=deletion_branch)

    assert pure_route.is_pure_gain
    assert not pure_route.is_mixed_gain_loss
    assert pure_route.matches_type(RouteType.PURE_GAIN)

    assert mixed_route.is_mixed_gain_loss
    assert not mixed_route.is_pure_gain
    assert mixed_route.matches_type(RouteType.MIXED_GAIN_LOSS)
    assert mixed_route.matches_type(RouteType.ALL)
