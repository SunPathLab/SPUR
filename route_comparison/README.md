# Route Comparison Tables

These files compare old SPUR complex-route generation against `SPUR_codex`
routes where both post-MRCA branches contain deletions.

## Kept Tables

```text
real_data_cn_pairs_old_spur_vs_codex_summary.tsv
```

Summary across canonical between-sample CN pairs observed in:

```text
/home/yunong-xia/Downloads/May_27_real_data_indices.tsv
```

This is the main validation table.
CN labels are parsed from the `between_sample_cn` column as `total:minor`.
For route enumeration, each label is converted internally as:

```text
major_cn = total_cn - minor_cn
minor_cn = minor_cn
```

```text
real_data_cn_pairs_one_branch_deletion_counts.tsv
```

Counts `SPUR_codex` complex routes for each real-data CN pair after splitting
them into:

```text
both_branch_deletion_routes
one_branch_deletion_routes
branch1_only_deletion_routes
branch2_only_deletion_routes
```

## Filters

The generated tables use:

```text
max_mrca_cn = min(max(sample total CN) + 2, 6)
max_branch_events = 5
```
