# SPUR Codex Prototype

This is a clean prototype package for SPUR route objects, route enumeration, route encoding, and T/J/B scoring.

The design separates three layers:

- `objects.py`: theoretical route objects and route labels.
- `enumeration.py`: truncal, branch, and paired-route enumeration.
- `route_filter.py`: route filtering, pure/complex split, and truncal grouping.
- `encoding.py` and `scoring.py`: conversion to numeric H/J arrays and T/H/B scoring.

The updated mixed-gain-loss definition is pair-level:

```python
route.is_mixed_gain_loss == route.branch1.has_deletion or route.branch2.has_deletion
```

## Notation

Use `total:minor` only for observed CN labels. For example, `4:2` means total
copy number 4 and minor copy number 2, so the allele-specific counts are
`major=2,minor=2`.

When referring directly to allele-specific major/minor counts, use
`major=...,minor=...` or field names like `sample1_major_cn`.

## Structure

```text
src/spur_codex/objects.py      theoretical route objects and labels
src/spur_codex/tree.py         gain-history tree utilities and constraint matrices
src/spur_codex/enumeration.py  truncal, branch, and paired-route enumeration
src/spur_codex/route_filter.py route filtering, route splitting, and truncal grouping
src/spur_codex/encoding.py     route-to-matrix and paired H/J tensor encoding
src/spur_codex/scoring.py      T/H/B scoring and binomial PMF tensor helpers
src/spur_codex/features.py     score-derived diagnostic features
scripts/simulate_small_datasets.py  small route/JMD simulation from real CN pairs
```

## Tests

The active `spur` environment does not currently include `pytest`, so a small
test runner is included:

```bash
PYTHONPATH=src python tests/run_tests.py
```

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
notation, truncal tree hash, allele-specific major/minor columns, branch
event counts, and pure-gain compatibility flag.
