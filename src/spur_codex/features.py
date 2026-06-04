from __future__ import annotations

import numpy as np
from scipy.special import logsumexp


def score_features(pure_gain_scores: np.ndarray, mixed_gain_loss_scores: np.ndarray) -> dict[str, float | int]:
    pure_gain_scores = np.asarray(pure_gain_scores, dtype=float)
    mixed_gain_loss_scores = np.asarray(mixed_gain_loss_scores, dtype=float)
    all_scores = np.concatenate([pure_gain_scores, mixed_gain_loss_scores])

    best_pg = float(np.max(pure_gain_scores))
    best_mgl = float(np.max(mixed_gain_loss_scores))
    no_less_than_best_pg = all_scores >= best_pg

    return {
        "percentile_rank": float(1.0 - np.mean(no_less_than_best_pg)),
        "score_ratio": float(np.mean(all_scores[no_less_than_best_pg]) / best_pg),
        "best_pg": best_pg,
        "best_mixed_gain_loss": best_mgl,
        "delta_best_mixed_gain_loss_minus_pg": best_mgl - best_pg,
        "pg_evidence": float(logsumexp(pure_gain_scores) - np.log(len(pure_gain_scores))),
        "mixed_gain_loss_evidence": float(logsumexp(mixed_gain_loss_scores) - np.log(len(mixed_gain_loss_scores))),
        "n_pg_routes": int(len(pure_gain_scores)),
        "n_mixed_gain_loss_routes": int(len(mixed_gain_loss_scores)),
        "frac_mixed_gain_loss_ge_best_pg": float(np.mean(mixed_gain_loss_scores >= best_pg)),
    }
