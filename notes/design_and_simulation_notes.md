# SPUR Codex Design And Simulation Notes

Date: 2026-06-04

This note documents the current `SPUR_codex` prototype design and the small simulation workflow.

## Design Goal

The package separates the route theory from the numeric scoring backend.

The main layers are:

```text
Object layer
    Defines theoretical route objects.

Enumeration layer
    Generates truncals, branches, and paired two-sample routes.

Encoding layer
    Converts route objects into matrices and J tensors.

Scoring / compatibility layer
    Computes T/J/B scores and pure-gain compatibility.

Simulation layer
    Samples timing vectors and ground-truth JMDs.
```

This is intentionally different from putting all behavior into one route-tree class. The goal is readability and easier changes to route definitions.

## Object Layer

Implemented in:

```text
src/spur_codex/objects.py
```

Core objects:

```text
GainHistoryTree
    major_tree
    minor_tree

Truncal
    major_cn
    minor_cn
    tree
    tree_hash

Branch
    copy_fate_vector
    branch_index

Route
    truncal
    branch1
    branch2
```

A `Branch.copy_fate_vector` records, for each MRCA copy, how many descendant copies remain in the final sample.

Examples:

```text
(1, 1, 1) = all MRCA copies retained
(2, 1, 1) = first MRCA copy gained
(0, 2, 1) = first MRCA copy deleted, second gained
```

## Route Definitions

The updated route definition is pair-level:

```python
route.is_mixed_gain_loss = route.branch1.has_deletion or route.branch2.has_deletion
route.is_pure_gain = not route.is_mixed_gain_loss
```

This avoids the old behavior where each branch was pre-filtered as complex and then cross-producted, which forced both branches to contain deletions.

## Enumeration Layer

Implemented in:

```text
src/spur_codex/enumeration.py
```

Enumerators:

```text
TruncalEnumerator
    Enumerates MRCA states and gain-history trees.

BranchEnumerator
    Enumerates copy-fate vectors for one sample under an event cap.

RouteEnumerator
    Combines branch1 and branch2 into paired routes and filters by RouteType.
```

For CN pairs, the simulation script canonicalizes `(sample1_cn, sample2_cn)` so reversed pairs are treated consistently. For example, if present, `3:1 vs 4:2` and `4:2 vs 3:1` are both represented as `4:2 vs 3:1`.

## Encoding Layer

Implemented in:

```text
src/spur_codex/encoding.py
```

The encoder converts routes into:

```text
copy-fate matrices
multiplicity matrices
paired J tensor
```

The J tensor shape is:

```text
n_routes x n_timing_stages x n_joint_multiplicity_groups
```

J uses sample1-major order:

```python
joint_idx = (sample1_multiplicity - 1) * sample2_major_cn + (sample2_multiplicity - 1)
```

For mixed-gain-loss routes, J is built from explicit route pairs. It is not a rectangular cross product of independently filtered branch lists.

## T/J/B Scoring

Implemented in:

```text
src/spur_codex/scoring.py
```

The scorer computes route log scores from:

```text
T: timing samples
J: timing contribution tensor
B: binomial likelihood tensor
```

A posterior helper normalizes route log scores:

```python
posterior = exp(log_scores - logsumexp(log_scores))
```

The posterior calculation is tested using a deterministic two-route example where route 0 has likelihood `0.9^2` and route 1 has likelihood `0.1^2`.

## Pure-Gain Compatibility

Implemented in:

```text
src/spur_codex/compatibility.py
```

The current compatibility checker preserves the old mathematical meaning:

Given a target JMD `q`, it checks whether any pure-gain route can generate `q` under some nonnegative timing vector `t` satisfying the truncal constraint system.

For each pure-gain route candidate:

```text
A t = 1
Operator t = q * sum(t)
```

The code uses the same feasibility construction as the old script:

```python
test_matrix = vstack((A, Operator - outer(q, ones(n_timings))))
test_b = hstack((ones(A.shape[0]), zeros_like(q)))
compatible = nnls(test_matrix, test_b).residual <= atol
```

Optimization currently used:

```text
1. Cache all pure-gain route feasibility candidates per CN pair.
2. Before NNLS, reject candidates whose support cannot cover q.
```

The support rejection is cheap and avoids many NNLS solves. It checks:

```python
if any(q > atol where candidate support is false):
    skip candidate
```

Potential future optimizations:

```text
- Group pure-gain candidates by support pattern.
- Use linear programming feasibility directly instead of NNLS residual checks.
- Precompute convex-hull or cone representations for each support class.
```

## Simulation Workflow

Implemented in:

```text
src/spur_codex/timing.py
src/spur_codex/simulation.py
scripts/simulate_small_datasets.py
```

The script reads CN pairs from:

```text
/home/yunong-xia/Downloads/May_27_real_data_indices.tsv
```

It selects a few small canonical CN pairs, then simulates examples for three classes:

```text
pure_gain
mixed_gain_loss
conclusive_mixed_gain_loss
```

Definitions:

```text
pure_gain
    Route has no deletion in either branch.

mixed_gain_loss
    Either branch has deletion. No compatibility restriction.

conclusive_mixed_gain_loss
    Mixed-gain-loss route whose sampled ground-truth JMD is not pure-gain compatible.
```

Default command:

```bash
PYTHONPATH=src python scripts/simulate_small_datasets.py \
  --cn-tsv /home/yunong-xia/Downloads/May_27_real_data_indices.tsv \
  --out-dir simulated_small \
  --n-cn-pairs 3 \
  --n-per-class 20
```

Current generated dataset:

```text
simulated_small/
```

It contains 180 examples:

```text
60 pure_gain
60 mixed_gain_loss
60 conclusive_mixed_gain_loss
```

## Lightweight Output Format

Each example is one compressed `.npz` file containing:

```text
jmd
    Ground-truth joint multiplicity density. Shape: n_joint_multiplicity_groups.

timing
    Sampled timing vector for the truncal tree.

branch1_copy_fate
    Copy-fate vector for sample 1.

branch2_copy_fate
    Copy-fate vector for sample 2.
```

The metadata table is:

```text
simulated_small/metadata.tsv
```

It records:

```text
example_id
class_label
npz_path
sample1_cn
sample2_cn
sample1_major_cn
sample1_minor_cn
sample2_major_cn
sample2_minor_cn
mrca_major_cn
mrca_minor_cn
truncal_tree_hash
branch1_index
branch2_index
branch1_events
branch2_events
pure_gain_compatible
```

This keeps the data compact while preserving enough information to reconstruct the biological route metadata and inspect the ground-truth JMD/timing.

## Tests

Tests are in:

```text
tests/
```

Because the active environment does not currently include `pytest`, a small runner is included:

```bash
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src:. python tests/run_tests.py
```

The tests cover:

```text
object properties
route enumeration and updated mixed-gain-loss definition
4:2 vs 3:1 route-count sanity checks
paired J tensor layout
T/J/B scoring
posterior normalization
pure-gain compatibility
simulation CN-pair canonicalization
```
