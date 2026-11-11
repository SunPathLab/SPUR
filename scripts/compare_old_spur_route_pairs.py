from __future__ import annotations

import argparse
import csv
from pathlib import Path

from compare_old_spur_routes import (
    canonical_tree_multiset,
    codex_both_deletion_records,
    old_spur_complex_records,
)
from spur.notation import canonical_cn_pair, parse_between_sample_cn, total_minor_label


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--old-spur-src", default="../SPUR/src/spur")
    parser.add_argument("--cn-tsv", default="/home/yunong-xia/Downloads/May_27_real_data_indices.tsv")
    parser.add_argument("--max-total-cn", type=int, default=None)
    parser.add_argument("--max-truncal-cn", type=int, default=None)
    parser.add_argument("--max-branch-events", type=int, default=5)
    parser.add_argument("--out", default="route_comparison/real_data_cn_pairs_route_comparison_summary.tsv")
    args = parser.parse_args()

    pairs = observed_cn_pairs_from_tsv(Path(args.cn_tsv), args.max_total_cn)
    rows = []
    for sample1_cn, sample2_cn in pairs:
        max_truncal_cn = args.max_truncal_cn or pair_max_truncal_cn(sample1_cn, sample2_cn)
        old_records = old_spur_complex_records(
            old_spur_src=Path(args.old_spur_src),
            sample1_cn=sample1_cn,
            sample2_cn=sample2_cn,
            max_truncal_cn=max_truncal_cn,
            max_branch_events=args.max_branch_events,
        )
        codex_records = codex_both_deletion_records(
            sample1_cn=sample1_cn,
            sample2_cn=sample2_cn,
            max_truncal_cn=max_truncal_cn,
            max_branch_events=args.max_branch_events,
        )
        old_routes = canonical_tree_multiset(old_records)
        codex_routes = canonical_tree_multiset(codex_records)
        old_keys = set(old_routes)
        codex_keys = set(codex_routes)
        matched = sum(min(old_routes[key], codex_routes[key]) for key in old_keys & codex_keys)
        old_only = sum(old_routes[key] for key in old_keys - codex_keys)
        codex_only = sum(codex_routes[key] for key in codex_keys - old_keys)
        count_mismatch = sum(
            abs(old_routes[key] - codex_routes[key])
            for key in old_keys & codex_keys
            if old_routes[key] != codex_routes[key]
        )
        rows.append(
            {
                "sample1_cn": total_minor_label(sample1_cn),
                "sample2_cn": total_minor_label(sample2_cn),
                "sample1_total_cn": sample1_cn[0],
                "sample1_minor_cn": sample1_cn[1],
                "sample2_total_cn": sample2_cn[0],
                "sample2_minor_cn": sample2_cn[1],
                "cn_tsv": args.cn_tsv,
                "max_truncal_cn": max_truncal_cn,
                "max_branch_events": args.max_branch_events,
                "old_spur_complex_routes": len(old_records),
                "spur_both_deletion_routes": len(codex_records),
                "matched_routes": matched,
                "old_only_routes": old_only,
                "codex_only_routes": codex_only,
                "count_mismatch_routes": count_mismatch,
                "route_sets_equal": old_routes == codex_routes,
            }
        )
        print(
            f"{total_minor_label(sample1_cn)} vs {total_minor_label(sample2_cn)}\t"
            f"old={len(old_records)}\tcodex={len(codex_records)}\t"
            f"equal={old_routes == codex_routes}"
        )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"summary_tsv\t{out_path}")


def observed_cn_pairs_from_tsv(path: Path, max_total_cn: int | None) -> list[tuple[tuple[int, int], tuple[int, int]]]:
    pairs = set()
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample1_cn, sample2_cn = parse_between_sample_cn(row["between_sample_cn"])
            if max_total_cn is not None and (sample1_cn[0] > max_total_cn or sample2_cn[0] > max_total_cn):
                continue
            pairs.add(canonical_cn_pair(sample1_cn, sample2_cn))
    return sorted(pairs, key=lambda pair: (pair[0][0], pair[0][1], pair[1][0], pair[1][1]))


def synthetic_observed_cn_pairs(max_total_cn: int) -> list[tuple[tuple[int, int], tuple[int, int]]]:
    states = []
    for total_cn in range(2, max_total_cn + 1):
        for minor_cn in range(1, total_cn // 2 + 1):
            states.append((total_cn, minor_cn))

    pairs = []
    for i, sample1_cn in enumerate(states):
        for sample2_cn in states[: i + 1]:
            pairs.append(canonical_cn_pair(sample1_cn, sample2_cn))
    return pairs


def pair_max_truncal_cn(sample1_cn: tuple[int, int], sample2_cn: tuple[int, int]) -> int:
    return min(max(sample1_cn[0], sample2_cn[0]) + 2, 6)


if __name__ == "__main__":
    main()
