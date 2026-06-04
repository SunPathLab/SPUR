from scripts.simulate_small_datasets import canonical_pair


def test_canonical_pair_puts_larger_total_cn_first():
    assert canonical_pair((2, 1), (2, 2)) == ((2, 2), (2, 1))
    assert canonical_pair((2, 2), (2, 1)) == ((2, 2), (2, 1))
