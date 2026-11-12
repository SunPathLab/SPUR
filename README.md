# SPUR

SPUR evaluates whether a between-sample segmental gained aneuploidy (SGA) profile is better explained by a **unidirectional two-terminal SGA trajectory (UDT)** or a **bidirectional trajectory (BDT)**. It enumerates candidate evolutionary trajectories, calculates their read-data likelihood scores, and summarizes support for UDT and BDT trajectories.

The current software reports route scores and score-derived features. It does not apply a universal cutoff to emit a single binary UDT/BDT label; users can interpret the reported evidence or use the features in a downstream classifier.

Trajectory groups in the current implementation:

- **UDT:** routes with no deletion on either post-MRCA branch. The current implementation refers to this group as pure gain (`pg`) in some output fields.
- **BDT:** routes with a deletion on at least one post-MRCA branch. The current implementation refers to this group as mixed gain-loss in some output fields; `any_branch_del` is the default route definition.

Route scores are Monte Carlo log-likelihood estimates for the observed shared SSNV reads under each candidate route. They are not posterior probabilities over routes.

## Installation

SPUR requires Python 3.10 or later. Install the package and its dependencies from the repository root:

```bash
python -m pip install -e .
```

This installs the `spur-codex` command. The module can also be run from a source checkout with `PYTHONPATH=src python -m spur.main`.

## Input data

Provide one CSV file containing read counts for SSNVs shared by a pair of tumor samples from the same patient. Each row represents one shared SSNV. The **first four columns must be in this order**:

| Column | Meaning |
| --- | --- |
| Sample 1 alternate reads | Reads supporting the alternate allele in sample 1 |
| Sample 1 reference reads | Reads supporting the reference allele in sample 1 |
| Sample 2 alternate reads | Reads supporting the alternate allele in sample 2 |
| Sample 2 reference reads | Reads supporting the reference allele in sample 2 |

Example CSV with headers:

```csv
sample1altc,sample1refc,sample2altc,sample2refc
16,39,18,28
13,27,16,42
```

For example, the columns may be named `AOCS088_RAaltc`, `AOCS088_RArefc`, `AOCS088_PTaltc`, and `AOCS088_PTrefc`. SPUR reads the first four columns by position; use this order even if the column names differ. Rows are retained when both samples meet the alternate-read and total-depth thresholds (defaults: at least 3 alternate reads and depth 15).

The CN pair and purities are supplied as command-line arguments. CN states use `total:minor` notation, such as `4:2`; the example pair `4:2_3:1` describes sample 1 followed by sample 2. Purities and prevalence values must correspond to the same sample order. SPUR canonicalizes a noncanonical CN-pair ordering and swaps the read columns and sample-specific parameters with it.

A descriptive filename convention is:

```text
patient_chr_start_end_sample1_cn1_sample2_cn2_purity1_purity2.csv
```

For example:

```text
AOCS092_10_88827_38484114_AOCS092_RA_42_AOCS092_PT_42_0.631_0.58.csv
```

The filename is used to name output files; SPUR does not infer CN states or purities from it.

## Run SPUR

Example single-profile run with multiple Monte Carlo seeds:

```bash
spur-codex \
  --input_path path/to/shared_ssnv_reads.csv \
  --output_dir path/to/spur_results \
  --cn-pair "4:2_3:1" \
  --sample1_purity 0.8 \
  --sample2_purity 0.8 \
  --calculate-uncertainty \
  --seed-filepath path/to/seeds.txt \
  --report-route-score
```

The seed file contains one integer seed per line. Without `--calculate-uncertainty`, SPUR runs once using `--seed` (default `1`). For batch processing, provide either `--batch-filepath-list` (one input path per line) or `--batch-input-paths` together with `--batch-process`.

## Outputs

For each input, SPUR writes a feature table named `<input-basename>.tsv`. With uncertainty enabled, it contains one row per seed, including the seed, runtime, and score-derived features such as:

- `best_pg` and `best_mixed_gain_loss`: highest finite route scores for UDT and BDT, respectively. The `pg` and `mixed_gain_loss` names are retained in the current output schema.
- `delta_best_mixed_gain_loss_minus_pg`: difference between the best BDT and UDT scores.
- `percentile_rank`, `score_ratio`, and `frac_mixed_gain_loss_ge_best_pg`: compare route scores with the best UDT score. The last field is the fraction of finite BDT route scores at least as high as that score.
- `pg_evidence` and `mixed_gain_loss_evidence`: log-mean evidence summaries for UDT and BDT routes, respectively.
- `n_pg_routes` and `n_mixed_gain_loss_routes`: counts of finite-scoring UDT and BDT routes, respectively.

These features summarize the enumerated routes under the selected search settings; they are not calibrated probabilities. Non-finite route scores are omitted from feature calculations.

Add `--report-route-score` to also write `<input-basename>_route_scores.tsv`, containing each route's truncal and branch keys and one score column per seed. Set `--output_name_prefix` to choose a different output basename for a single-input run.

## Parameters to tune

The most consequential settings control read filtering, the route search space, or Monte Carlo accuracy/runtime.

| Option | Default | What it controls |
| --- | ---: | --- |
| `--alt-reads-threshold` | `3` | Minimum alternate reads required in each sample for an SSNV to be retained. |
| `--total-reads-threshold` | `15` | Minimum total depth required in each sample. |
| `--max-truncal-cn` | `min(max(sample total CN) + 2, 6)` | Maximum copy number allowed in the shared pre-branch truncal history. |
| `--max-branch-events` | `5` | Maximum events enumerated on each post-MRCA branch. |
| `--route-filter` | `none` | Optional route filter: `total_events`, `max_branch_events`, `total_branch_events`, or `single_sample_total_events`. |
| `--route-filter-threshold` | unset | Required when a route filter other than `none` is selected. |
| `--n_Monte_Carlo_samples` | `4000` | Number of timing/VAF Monte Carlo samples per run; increasing it can reduce sampling noise at higher runtime. |
| `--joint_vafs_scale` | `0.05` | Spread of sampled joint VAFs around their expected values. |
| `--burn-in` | `200` | Burn-in steps for timing-vector sampling. |
| `--skips` | `50` | Sampling steps between retained timing vectors. |
| `--float-dtype` | `float64` | Numeric precision for score calculations (`float32` or `float64`). |
| `--seed` / `--seed-filepath` | `1` / unset | Random seed for one run, or a list of seeds when `--calculate-uncertainty` is used. |

Sample purities (`--sample1_purity`, `--sample2_purity`) and prevalences (`--sample1_prevalence`, `--sample2_prevalence`) are also configurable; each defaults to `1`. The BDT route-definition option is named `--mixed-gain-loss-route-definition` in the CLI. Its default, `any_branch_del`, includes routes with a deletion on either branch; the experimental `both_branch_del` setting requires deletions on both branches.

## Manuscript

**Hidden bidirectional copy number trajectories shape the between-sample genomic divergence in human cancers**

Yunong Xia, Chenyu Wu, and Ruping Sun
Department of Laboratory Medicine and Pathology, University of Minnesota Twin Cities, Minneapolis, Minnesota, USA

The manuscript is in preparation for preprint upload. A public preprint link and citation will be added here when available.

## Development

Core implementation is in `src/spur/`: route enumeration and notation, route filtering, timing and VAF sampling, matrix encoding, scoring, and score-derived features. Run the test suite from the repository root with:

```bash
PYTHONPATH=src python tests/run_tests.py
```
