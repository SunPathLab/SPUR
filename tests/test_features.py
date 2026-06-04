import numpy as np

from spur_codex.features import score_features


def test_score_features_uses_mixed_gain_loss_names():
    features = score_features(
        pure_gain_scores=np.array([-10.0, -12.0]),
        mixed_gain_loss_scores=np.array([-9.0, -11.0, -13.0]),
    )

    assert features["best_pg"] == -10.0
    assert features["best_mixed_gain_loss"] == -9.0
    assert features["delta_best_mixed_gain_loss_minus_pg"] == 1.0
    assert features["n_pg_routes"] == 2
    assert features["n_mixed_gain_loss_routes"] == 3
