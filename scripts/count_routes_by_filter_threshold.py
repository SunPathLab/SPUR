from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

from spur_codex.enumeration import RouteEnumerator
from spur_codex.notation import parse_between_sample_cn
from spur_codex.route_filter import RouteFilter


FILTERS = {
    "total_events": RouteFilter.filter_by_total_events,
    "max_branch_events": RouteFilter.filter_by_max_branch_events,
    "total_branch_events": RouteFilter.filter_by_total_branch_events,
    "single_sample_total_events": RouteFilter.filter_by_single_sample_total_events,
}


def read_cn_pairs(path: Path) -> list[str]:
    seen = set()
    cn_pairs = []
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cn_pair = row["between_sample_cn"]
            if cn_pair in seen:
                continue
            seen.add(cn_pair)
            cn_pairs.append(cn_pair)
    return cn_pairs


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    cn_pair_path = repo_root / "route_counts" / "route_total_event_threshold8_counts.tsv"
    out_path = repo_root / "route_counts" / "route_counts_by_filter_threshold.tsv"

    cn_pairs = read_cn_pairs(cn_pair_path)
    enumerator = RouteEnumerator()
    rows = []

    for cn_pair in cn_pairs:
        sample1_cn, sample2_cn = parse_between_sample_cn(cn_pair)
        max_truncal_cn = min(max(sample1_cn[0], sample2_cn[0]) + 2, 6)
        routes = tuple(
            enumerator.enumerate(
                sample1_cn=sample1_cn,
                sample2_cn=sample2_cn,
                max_truncal_cn=max_truncal_cn,
                max_branch_events=5,
            )
        )

        for filter_name, filter_fn in FILTERS.items():
            for threshold in range(1, 11):
                filtered = filter_fn(routes, threshold)
                n_pure_gain = sum(route.is_pure_gain for route in filtered)
                n_mixed_gain_loss = len(filtered) - n_pure_gain
                rows.append(
                    {
                        "between_sample_cn": cn_pair,
                        "sample1_cn": f"{sample1_cn[0]}:{sample1_cn[1]}",
                        "sample2_cn": f"{sample2_cn[0]}:{sample2_cn[1]}",
                        "max_truncal_cn": max_truncal_cn,
                        "max_branch_events_for_enumeration": 5,
                        "route_filter": filter_name,
                        "threshold": threshold,
                        "n_pure_gain_routes": n_pure_gain,
                        "n_mixed_gain_loss_routes": n_mixed_gain_loss,
                        "n_total_routes": len(filtered),
                        "n_total_unfiltered_routes": len(routes),
                        "fraction_kept": len(filtered) / len(routes) if routes else 0.0,
                    }
                )

    df = pd.DataFrame(rows)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)

    print(f"Wrote {out_path}")
    print(f"CN pairs: {len(cn_pairs)}")
    print(f"Rows: {len(df)}")
    print()
    print("Summary by filter and threshold:")
    summary = (
        df.groupby(["route_filter", "threshold"], as_index=False)
        .agg(
            mean_total_routes=("n_total_routes", "mean"),
            median_total_routes=("n_total_routes", "median"),
            max_total_routes=("n_total_routes", "max"),
        )
        .sort_values(["route_filter", "threshold"])
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
