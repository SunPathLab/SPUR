from __future__ import annotations

import argparse
import csv
from pathlib import Path

from spur_codex.enumeration import RouteEnumerator
from spur_codex.notation import canonical_cn_pair, parse_between_sample_cn, total_minor_label
from spur_codex.route_filter import split_routes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cn-tsv", default="/home/yunong-xia/Downloads/May_27_real_data_indices.tsv")
    parser.add_argument("--max-truncal-cn", type=int, default=None)
    parser.add_argument("--max-branch-events", type=int, default=5)
    parser.add_argument("--out", default="route_comparison/real_data_cn_pairs_one_branch_deletion_counts.tsv")
    args = parser.parse_args()

    pairs = observed_cn_pairs_from_tsv(Path(args.cn_tsv))
    enumerator = RouteEnumerator()
    rows = []
    for sample1_cn, sample2_cn in pairs:
        max_truncal_cn = args.max_truncal_cn or pair_max_truncal_cn(sample1_cn, sample2_cn)
        routes = enumerator.enumerate(
            sample1_cn=sample1_cn,
            sample2_cn=sample2_cn,
            max_truncal_cn=max_truncal_cn,
            max_branch_events=args.max_branch_events,
        )
        split = split_routes(routes)
        branch1_only = 0
        branch2_only = 0
        both_branch = 0
        for route in split.mixed_gain_loss:
            branch1_has_deletion = route.branch1.has_deletion
            branch2_has_deletion = route.branch2.has_deletion
            if branch1_has_deletion and branch2_has_deletion:
                both_branch += 1
            elif branch1_has_deletion:
                branch1_only += 1
            elif branch2_has_deletion:
                branch2_only += 1

        one_branch = branch1_only + branch2_only
        rows.append(
            {
                "sample1_cn": total_minor_label(sample1_cn),
                "sample2_cn": total_minor_label(sample2_cn),
                "sample1_total_cn": sample1_cn[0],
                "sample1_minor_cn": sample1_cn[1],
                "sample2_total_cn": sample2_cn[0],
                "sample2_minor_cn": sample2_cn[1],
                "max_truncal_cn": max_truncal_cn,
                "max_branch_events": args.max_branch_events,
                "all_routes": len(routes),
                "pure_gain_routes": len(split.pure_gain),
                "complex_routes": len(split.mixed_gain_loss),
                "both_branch_deletion_routes": both_branch,
                "one_branch_deletion_routes": one_branch,
                "branch1_only_deletion_routes": branch1_only,
                "branch2_only_deletion_routes": branch2_only,
            }
        )
        print(
            f"{total_minor_label(sample1_cn)} vs {total_minor_label(sample2_cn)}\t"
            f"complex={len(split.mixed_gain_loss)}\t"
            f"both={both_branch}\t"
            f"one_branch={one_branch}"
        )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"counts_tsv\t{out_path}")


def observed_cn_pairs_from_tsv(path: Path) -> list[tuple[tuple[int, int], tuple[int, int]]]:
    pairs = set()
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample1_cn, sample2_cn = parse_between_sample_cn(row["between_sample_cn"])
            pairs.add(canonical_cn_pair(sample1_cn, sample2_cn))
    return sorted(pairs, key=lambda pair: (pair[0][0], pair[0][1], pair[1][0], pair[1][1]))


def pair_max_truncal_cn(sample1_cn: tuple[int, int], sample2_cn: tuple[int, int]) -> int:
    return min(max(sample1_cn[0], sample2_cn[0]) + 2, 6)


if __name__ == "__main__":
    main()
