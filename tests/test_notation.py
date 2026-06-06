from spur_codex.notation import CopyNumberPair, CopyNumberState, allele_cn_label, canonical_cn_pair, parse_between_sample_cn, total_minor_label, truncal_key


def test_copy_number_state_uses_total_minor_boundary_notation():
    state = CopyNumberState.from_label("4:2")

    assert state.total_cn == 4
    assert state.minor_cn == 2
    assert state.major_cn == 2
    assert state.allele_tuple == (2, 2)
    assert state.label == "4:2"


def test_copy_number_pair_parses_between_sample_cn_and_canonicalizes():
    pair = CopyNumberPair.from_label("2:1_3:1").canonical()

    assert pair.label == "3:1_2:1"
    assert pair.sample1_allele_tuple == (2, 1)
    assert pair.sample2_allele_tuple == (1, 1)


def test_parse_between_sample_cn_returns_internal_allele_tuples():
    assert parse_between_sample_cn("4:2_3:1") == ((2, 2), (2, 1))


def test_canonical_cn_pair_preserves_internal_allele_tuple_output():
    assert canonical_cn_pair((1, 1), (2, 1)) == ((2, 1), (1, 1))


def test_total_minor_label_uses_total_colon_minor():
    assert total_minor_label((2, 2)) == "4:2"
    assert total_minor_label((2, 1)) == "3:1"


def test_allele_cn_label_does_not_use_colon_notation():
    assert allele_cn_label(2, 1) == "major=2,minor=1"


def test_truncal_key_does_not_look_like_total_minor_label():
    key = truncal_key(2, 1, "major=(L,L);minor=L")
    assert key == "truncal_major2_minor1_major=(L,L);minor=L"
    assert "2:1" not in key
