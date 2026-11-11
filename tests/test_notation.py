from spur.notation import CopyNumberPair, CopyNumberState, canonical_cn_pair, parse_between_sample_cn, to_major_minor_tuple, total_minor_label, truncal_key


def test_copy_number_state_uses_total_minor_boundary_notation():
    state = CopyNumberState.from_label("4:2")

    assert state.total_cn == 4
    assert state.minor_cn == 2
    assert state.major_cn == 2
    assert state.cn_tuple == (4, 2)
    assert state.major_minor_tuple == (2, 2)
    assert state.label == "4:2"


def test_copy_number_pair_parses_between_sample_cn_and_canonicalizes():
    pair = CopyNumberPair.from_label("2:1_3:1").canonical()

    assert pair.label == "3:1_2:1"
    assert pair.sample1_cn_tuple == (3, 1)
    assert pair.sample2_cn_tuple == (2, 1)


def test_parse_between_sample_cn_returns_total_minor_tuples():
    assert parse_between_sample_cn("4:2_3:1") == ((4, 2), (3, 1))
    assert parse_between_sample_cn("3:1_4:2") == ((4, 2), (3, 1))
    assert parse_between_sample_cn("4:1_4:2") == ((4, 2), (4, 1))

def test_canonical_cn_pair_preserves_total_minor_tuple_output():
    assert canonical_cn_pair((2, 1), (3, 1)) == ((3, 1), (2, 1))


def test_total_minor_label_uses_total_colon_minor():
    assert total_minor_label((4, 2)) == "4:2"
    assert total_minor_label((3, 1)) == "3:1"


def test_total_minor_tuple_can_be_converted_at_old_spur_boundaries():
    assert to_major_minor_tuple((3, 1)) == (2, 1)


def test_truncal_key_uses_total_minor_cn_and_tree_shape():
    key = truncal_key(3, 1, "major=(L,L);minor=L")
    assert key == "truncal_3:1_major=(L,L);minor=L"
