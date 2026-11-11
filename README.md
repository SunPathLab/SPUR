# SPUR Codex Prototype

This is a clean prototype package for SPUR route objects, route enumeration, route encoding, and T/H/B scoring.

The design separates these responsibilities:

- `tree.py`: generic tree node structure and tree-shape/path utilities only.
- `route.py`: biological route objects, including `Truncal.constraint_matrix`.
- `enumeration.py`: truncal, branch, and paired-route enumeration.
- `route_filter.py`: route filtering, pure/complex split, and truncal grouping.
- `encoding.py`: conversion from routes to numeric H arrays.
- `sampling.py`: timing and VAF sampling.
- `scoring.py`: B preparation/binomial likelihood helpers and T/H/B scoring.

The updated mixed-gain-loss definition is pair-level:

```python
route.is_mixed_gain_loss == route.branch1.has_deletion or route.branch2.has_deletion
```

## Notation

Use `total:minor` only for observed CN labels. For example, `4:2` means total
copy number 4 and minor copy number 2, so the allele-specific counts are
`major=2,minor=2`.

Public API tuples use `(total_cn, minor_cn)`. Major copy number is derived as
`total_cn - minor_cn` inside route objects and allele-tree math.

## Structure

```text
src/spur/route.py      theoretical route objects and labels
src/spur/tree.py         generic tree node utilities
src/spur/enumeration.py  truncal, branch, and paired-route enumeration
src/spur/route_filter.py route filtering, route splitting, and truncal grouping
src/spur/encoding.py     route-to-matrix and paired H tensor encoding
src/spur/hit_and_run.py  GRITIC-derived hit-and-run sampler
src/spur/sampling.py     timing/VAF sampling
src/spur/scoring.py      B preparation, T/H/B scoring, and binomial PMF helpers
src/spur/features.py     score-derived diagnostic features
scripts/simulate_small_datasets.py  small route/JMD simulation from real CN pairs
scripts/main.py  CLI smoke-test for route enumeration, grouping, timing, and H
```

## Tests

The active `spur` environment does not currently include `pytest`, so a small
test runner is included:

```bash
PYTHONPATH=src python tests/run_tests.py
```

## CLI Naming

Use `--max-truncal-cn` for the shared pre-branch copy-number cap. Commands and
output tables use truncal naming consistently.

## Small Simulation

Generate lightweight simulated route/JMD examples from the real-data CN-pair
TSV:

```bash
PYTHONPATH=src python scripts/simulate_small_datasets.py \
  --cn-tsv /home/yunong-xia/Downloads/May_27_real_data_indices.tsv \
  --out-dir simulated_small \
  --n-cn-pairs 3 \
  --n-per-class 20 \
  --max-branch-events 5
```

Each example is one compressed `.npz` file containing:

```text
jmd                 ground-truth joint multiplicity density
timing              sampled timing vector
branch1_copy_fate   sample 1 branch copy-fate vector
branch2_copy_fate   sample 2 branch copy-fate vector
```

`metadata.tsv` records the class label, observed CN pair in `total:minor`
notation, canonical `truncal_key`, total/minor CN columns, branch event counts,
and pure-gain compatibility flag.
