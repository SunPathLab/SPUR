from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from spur_codex.enumeration import RouteEnumerator
from spur_codex.route_filter import RouteFilter


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cn-profile-tsv",
        type=Path,
        default=Path(
            "/home/yunong-xia/Downloads/"
            "between_sample_cn_profiles_May_13_aneuploidy_thres_1.2e+07.tsv"
        ),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=repo_root
        / "route_counts"
        / "loh_loh_route_counts_by_total_event_threshold.tsv",
    )
    parser.add_argument("--max-observed-cn", type=int, default=6)
    parser.add_argument("--max-branch-events", type=int, default=5)
    parser.add_argument("--max-threshold", type=int, default=10)
    args = parser.parse_args()

    profiles = pd.read_csv(args.cn_profile_tsv, sep="\t")
    loh_loh = profiles[(profiles["minor1"] == 0) & (profiles["minor2"] == 0)].copy()
    loh_loh = loh_loh[
        (loh_loh["major1"] <= args.max_observed_cn)
        & (loh_loh["major2"] <= args.max_observed_cn)
    ]

    loh_loh["sample1_total_cn"] = loh_loh[["major1", "major2"]].max(axis=1)
    loh_loh["sample2_total_cn"] = loh_loh[["major1", "major2"]].min(axis=1)
    observed_pairs = (
        loh_loh.groupby(["sample1_total_cn", "sample2_total_cn"])
        .size()
        .reset_index(name="n_observed_segments")
    )

    enumerator = RouteEnumerator()
    rows = []
    for pair in observed_pairs.itertuples(index=False):
        sample1_total_cn = int(pair.sample1_total_cn)
        sample2_total_cn = int(pair.sample2_total_cn)
        max_truncal_cn = min(max(sample1_total_cn, sample2_total_cn) + 2, 6)
        if sample2_total_cn == 0:
            routes = tuple()
        else:
            routes = tuple(
                enumerator.enumerate(
                    sample1_cn=(sample1_total_cn, 0),
                    sample2_cn=(sample2_total_cn, 0),
                    max_truncal_cn=max_truncal_cn,
                    max_branch_events=args.max_branch_events,
                )
            )

        for threshold in range(1, args.max_threshold + 1):
            filtered_routes = RouteFilter.filter_by_total_events(routes, threshold)
            rows.append(
                {
                    "between_sample_cn": (
                        f"{sample1_total_cn}:0_{sample2_total_cn}:0"
                    ),
                    "sample1_cn": f"{sample1_total_cn}:0",
                    "sample2_cn": f"{sample2_total_cn}:0",
                    "n_observed_segments": int(pair.n_observed_segments),
                    "max_truncal_cn": max_truncal_cn,
                    "max_branch_events_for_enumeration": args.max_branch_events,
                    "total_event_threshold": threshold,
                    "n_routes": len(filtered_routes),
                    "n_routes_unfiltered": len(routes),
                }
            )

    result = pd.DataFrame(rows)
    result = result.sort_values(["between_sample_cn", "total_event_threshold"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, sep="\t", index=False)

    print(f"Wrote {args.out}")
    print(f"LOH-LOH CN pairs: {result['between_sample_cn'].nunique()}")
    print(
        result.pivot(
            index="between_sample_cn",
            columns="total_event_threshold",
            values="n_routes",
        ).to_string()
    )


if __name__ == "__main__":
    main()
