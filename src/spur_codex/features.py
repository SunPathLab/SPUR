from __future__ import annotations
from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp

from spur_codex.route import Route

@dataclass(frozen=True)
class RouteScoreResult:
    routes: tuple[Route, ...]
    scores: np.ndarray

    def score_features(self) -> dict[str, float | int]:
        return score_features(self)

    
    def route_score_table(self, seed: int) -> dict[str, list[str] | np.ndarray]:
        return route_score_table(self, seed)
        

def score_features(route_score_result: RouteScoreResult) -> dict[str, float | int]:
    
    pure_gain_mask = np.array([route.is_pure_gain for route in route_score_result.routes])
    pure_gain_scores = route_score_result.scores[pure_gain_mask]
    mixed_gain_loss_scores = route_score_result.scores[~pure_gain_mask]
    all_scores = route_score_result.scores

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
    
    
def route_score_table(route_score_result: RouteScoreResult, seed: int) -> dict[str, list[str] | np.ndarray]:
    
    # route information
    route_score_dictionary = {
      "truncal": [r.truncal_key for r in route_score_result.routes],
      "branch1": [r.branch1.canonical_key for r in route_score_result.routes],
      "branch2": [r.branch2.canonical_key for r in route_score_result.routes],
    }

    # route score
    route_score_dictionary["score_seed_" + str(seed)] = route_score_result.scores
        
    return route_score_dictionary
    
    
