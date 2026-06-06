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

Route filter layer
    Filters, splits, and groups already-enumerated routes.

Encoding layer
    Converts route objects into matrices and H tensors.

Scoring / compatibility layer
    Computes T/H/B scores and pure-gain compatibility.

Simulation layer
    Samples timing vectors and ground-truth JMDs.
```

This is intentionally different from putting all behavior into one route-tree class. The goal is readability and easier changes to route definitions.

## Copy-Number Notation

Use `total:minor` only for observed CN labels. For example:

```text
4:2 = total copy number 4, minor copy number 2
```

The corresponding allele-specific counts are written as:

```text
major=2,minor=2
```

Do not use `major:minor` notation. In code, tuples such as `(2, 1)` still mean
`(major_cn, minor_cn)`, but display strings should use either `total:minor` or
explicit field names like `major=2,minor=1`.

CN-pair parsing is centralized in:

```text
src/spur_codex/notation.py
```

Boundary classes:

```text
CopyNumberState
    total_cn
    minor_cn
    major_cn = total_cn - minor_cn
    allele_tuple = (major_cn, minor_cn)

CopyNumberPair
    sample1
    sample2
    label = sample1_total:minor_sample2_total:minor
```

Scripts should read external CN labels as `total:minor` first, then convert to
internal allele-specific tuples only when calling route enumeration.

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
    canonical_key

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

For CN pairs, the simulation script canonicalizes `(sample1_cn, sample2_cn)` so reversed pairs are treated consistently. Displayed CN labels use `total:minor`, so if present, `3:1 vs 4:2` and `4:2 vs 3:1` are both represented as `4:2 vs 3:1`.

## Route Filter Layer

Implemented in:

```text
src/spur_codex/route_filter.py
```

The intended workflow is:

```python
all_routes = RouteEnumerator().enumerate(sample1_cn, sample2_cn, max_mrca_cn, max_branch_events)
split = split_routes(all_routes)
groups = group_by_truncal(split.mixed_gain_loss)
```

This keeps enumeration broad and moves route selection into a separate module.
Each `Truncal` has a canonical key that uses the readable branching shape
instead of the NetworkX hash:

```text
truncal_major{major_cn}_minor{minor_cn}_major={tree_shape};minor={tree_shape}
```

Tree shapes use:

```text
L = leaf
(L,L) = one split with two leaves
((L,L),L) = a three-leaf branching shape
```

For example:

```text
truncal_major3_minor1_major=((L,L),L);minor=L
```

Each `Route` exposes the same key through `route.truncal_key`.
Routes can also generate their own single-route H operator:

```python
H_route = route.h(sample1_major_cn, sample2_major_cn)
```

For scoring, routes should usually be grouped by truncal first. A `RouteGroup`
can collect the whole group H tensor:

```python
H_group = group.h_tensor(sample1_major_cn, sample2_major_cn)
```

## Encoding Layer

Implemented in:

```text
src/spur_codex/encoding.py
```

The encoder converts routes into:

```text
copy-fate matrices
multiplicity matrices
paired H tensor
```

The H tensor shape is:

```text
n_routes x n_timing_stages x n_joint_multiplicity_groups
```

H uses sample1-major order:

```python
joint_idx = (sample1_multiplicity - 1) * sample2_major_cn + (sample2_multiplicity - 1)
```

For mixed-gain-loss routes, H is built from explicit route pairs. It is not a rectangular cross product of independently filtered branch lists.

## T/H/B Scoring

Implemented in:

```text
src/spur_codex/scoring.py
```

The scorer computes route log scores from:

```text
T: timing samples
H: timing contribution tensor
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
  --n-per-class 20 \
  --max-branch-events 5
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
truncal_key
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
paired H tensor layout
T/H/B scoring
posterior normalization
pure-gain compatibility
simulation CN-pair canonicalization
```
