#!/usr/bin/env bash
set -e

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"
export PYTHONPATH="$REPO_ROOT/src"

OUT="run_outputs/permutation_comparison_4_seeds"
SEEDS="$OUT/seeds.txt"

conda activate spur

python -m spur_codex.main \
  --input_path /home/yunong-xia/Downloads/nonLOH_comparisons/A12_8_45223876_77358890_A12A_M1_21_A12D_M3_41_0.81397301703257_0.833884774.csv \
  --output_dir "$OUT/noncanonical_input_order" \
  --cn-pair 2:1_4:1 \
  --sample1_purity 0.81397301703257 \
  --sample2_purity 0.833884774 \
  --calculate-uncertainty \
  --seed-filepath "$SEEDS" \
  --float-dtype float32 \
  --route-filter total_branch_events \
  --route-filter-threshold 8 \
  --report-route-score

python -m spur_codex.main \
  --input_path "$OUT/canonical_input_order/input_4_1_then_2_1.csv" \
  --output_dir "$OUT/canonical_input_order" \
  --cn-pair 4:1_2:1 \
  --sample1_purity 0.833884774 \
  --sample2_purity 0.81397301703257 \
  --calculate-uncertainty \
  --seed-filepath "$SEEDS" \
  --float-dtype float32 \
  --route-filter total_branch_events \
  --route-filter-threshold 8 \
  --report-route-score
