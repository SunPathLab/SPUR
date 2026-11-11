"""Build the predefined WGD-parsimony and non-WGD-parsimony routes."""

import numpy as np

from spur.enumeration import BranchEnumerator, TruncalEnumerator


# Non-WGD parsimony: the truncal state is produced by allele-specific gains.
NON_WGD_PARSIMONY_TRUNCAL_CN_BY_PAIR = {
    ((3, 1), (2, 1)): (2, 1),
    ((4, 2), (3, 1)): (3, 1),
    ((4, 1), (3, 1)): (3, 1),
    ((4, 1), (2, 1)): (2, 1),
    ((3, 1), (3, 1)): (3, 1),
    ((4, 2), (4, 2)): (4, 2),
}

# WGD parsimony: WGD occurs in the truncal history. Post-truncal branches are
# still selected using the same minimum-event rule as non-WGD parsimony.
WGD_PARSIMONY_TRUNCAL_CN_BY_PAIR = {
    ((3, 1), (2, 1)): (3, 1),
    ((4, 2), (3, 1)): (4, 2),
    ((5, 2), (4, 2)): (4, 2),
    ((5, 2), (3, 1)): (4, 2),
    ((4, 1), (3, 1)): (3, 1),
    ((3, 1), (3, 1)): (3, 1),
    ((4, 2), (4, 2)): (4, 2),
    ((4, 1), (4, 1)): (4, 1),
    ((5, 2), (5, 2)): (5, 2),
}

# Visible WGD states require equal major and minor root timings.
WGD_PARSIMONY_CONSTRAINT_BY_TRUNCAL_CN = {
    (4, 2): (1, 0, 0, -1, 0, 0),
    (5, 2): (1, 0, 0, 0, 0, -1, 0, 0),
}


# Shared helpers
def _get_truncal(truncal_cn_by_pair, nt1, nb1, nt2, nb2):
    cn_pair = ((nt1, nb1), (nt2, nb2))
    truncal_nt, truncal_nb = truncal_cn_by_pair[cn_pair]
    return TruncalEnumerator().enumerate_for_state(truncal_nt, truncal_nb)[0]


def _get_parsimony_branches(truncal, target_total_cn, target_minor_cn):
    minimum_events = abs(target_total_cn - truncal.truncal_total_cn)
    return BranchEnumerator().enumerate(
        truncal=truncal,
        target_total_cn=target_total_cn,
        target_minor_cn=target_minor_cn,
        max_events=minimum_events,
    )


# Non-WGD parsimony
def get_non_wgd_parsimony_route_constraint_matrix(nt1, nb1, nt2, nb2):
    truncal = _get_truncal(
        NON_WGD_PARSIMONY_TRUNCAL_CN_BY_PAIR, nt1, nb1, nt2, nb2
    )
    return truncal.constraint_matrix


def get_non_wgd_parsimony_branches(nt1, nb1, nt2, nb2):
    truncal = _get_truncal(
        NON_WGD_PARSIMONY_TRUNCAL_CN_BY_PAIR, nt1, nb1, nt2, nb2
    )
    sample1_branches = _get_parsimony_branches(truncal, nt1, nb1)
    sample2_branches = _get_parsimony_branches(truncal, nt2, nb2)
    return truncal, sample1_branches, sample2_branches


# WGD parsimony
def get_wgd_parsimony_route_constraint_matrix(nt1, nb1, nt2, nb2):
    truncal = _get_truncal(WGD_PARSIMONY_TRUNCAL_CN_BY_PAIR, nt1, nb1, nt2, nb2)
    truncal_cn = (truncal.truncal_total_cn, truncal.truncal_minor_cn)
    wgd_constraint = WGD_PARSIMONY_CONSTRAINT_BY_TRUNCAL_CN.get(truncal_cn)
    if wgd_constraint is None:
        return truncal.constraint_matrix
    return np.vstack((truncal.constraint_matrix, wgd_constraint))


def get_wgd_parsimony_branches(nt1, nb1, nt2, nb2):
    truncal = _get_truncal(WGD_PARSIMONY_TRUNCAL_CN_BY_PAIR, nt1, nb1, nt2, nb2)
    sample1_branches = _get_parsimony_branches(truncal, nt1, nb1)
    sample2_branches = _get_parsimony_branches(truncal, nt2, nb2)
    return truncal, sample1_branches, sample2_branches
