Permutation comparison: canonical versus noncanonical input order.

noncanonical_input_order/: input columns and CN label are ordered 2:1 then 4:1. SPUR recognizes the CN pair as noncanonical and reorders reads and purities internally to canonical CN order.
canonical_input_order/: columns and CN label are ordered 4:1 then 2:1.

Both runs use seeds 75945, 85785, 7889, and 94073; float32; 4,000 Monte Carlo draws; scale 0.05; burn-in 200; skips 50; total_branch_events <= 8; and route-score reporting. The features and route scores match exactly between these input orderings.

To reproduce both runs from the repository root, execute:
  bash run_outputs/permutation_comparison_4_seeds/run_comparison.sh

The script reads the original patient CSV from /home/yunong-xia/Downloads/nonLOH_comparisons and the permuted CSV stored in canonical_input_order/.
