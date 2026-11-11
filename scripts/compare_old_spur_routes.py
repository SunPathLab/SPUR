from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from spur.enumeration import RouteEnumerator
from spur.notation import to_major_minor_tuple, total_minor_label
from spur.tree import union_tree_shape


@dataclass(frozen=True)
class RouteRecord:
    source: str
    truncal_total_cn: int
    truncal_minor_cn: int
    canonical_tree: str
    branch1_fate: tuple[int, ...]
    branch2_fate: tuple[int, ...]
    branch1_index: int
    branch2_index: int
    tree_index: int | None = None

    @property
    def key(self) -> tuple[int, int, str, tuple[int, ...], tuple[int, ...]]:
        return (
            self.truncal_total_cn,
            self.truncal_minor_cn,
            self.canonical_tree,
            self.branch1_fate,
            self.branch2_fate,
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--old-spur-src", default="../SPUR/src/spur")
    parser.add_argument("--sample1-total-cn", type=int, default=4)
    parser.add_argument("--sample1-minor-cn", type=int, default=2)
    parser.add_argument("--sample2-total-cn", type=int, default=3)
    parser.add_argument("--sample2-minor-cn", type=int, default=1)
    parser.add_argument("--max-truncal-cn", type=int, default=6)
    parser.add_argument("--max-branch-events", type=int, default=5)
    parser.add_argument("--out", default="route_comparison/4_2_vs_3_1_route_comparison.tsv")
    parser.add_argument("--fate-out", default="route_comparison/4_2_vs_3_1_route_comparison_fate_vectors.tsv")
    args = parser.parse_args()

    sample1_cn = (args.sample1_total_cn, args.sample1_minor_cn)
    sample2_cn = (args.sample2_total_cn, args.sample2_minor_cn)

    old_records = old_spur_complex_records(
        old_spur_src=Path(args.old_spur_src),
        sample1_cn=sample1_cn,
        sample2_cn=sample2_cn,
        max_truncal_cn=args.max_truncal_cn,
        max_branch_events=args.max_branch_events,
    )
    codex_records = codex_both_deletion_records(
        sample1_cn=sample1_cn,
        sample2_cn=sample2_cn,
        max_truncal_cn=args.max_truncal_cn,
        max_branch_events=args.max_branch_events,
    )

    old_by_key = {record.key: record for record in old_records}
    codex_by_key = {record.key: record for record in codex_records}
    all_keys = sorted(set(old_by_key) | set(codex_by_key), key=sort_key)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "comparison_status",
                "present_in_old_spur",
                "present_in_spur",
                "observed_sample1_cn",
                "observed_sample2_cn",
                "truncal_total_cn",
                "truncal_minor_cn",
                "canonical_tree",
                "old_tree_index",
                "old_branch1_index",
                "old_branch2_index",
                "codex_branch1_index",
                "codex_branch2_index",
                "branch1_fate",
                "branch2_fate",
                "branch1_events",
                "branch2_events",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for key in all_keys:
            old = old_by_key.get(key)
            codex = codex_by_key.get(key)
            record = old or codex
            assert record is not None
            present_old = old is not None
            present_codex = codex is not None
            if present_old and present_codex:
                status = "match"
            elif present_old:
                status = "old_only"
            else:
                status = "codex_only"
            writer.writerow(
                {
                    "comparison_status": status,
                    "present_in_old_spur": int(present_old),
                    "present_in_spur": int(present_codex),
                    "observed_sample1_cn": total_minor_label(sample1_cn),
                    "observed_sample2_cn": total_minor_label(sample2_cn),
                    "truncal_total_cn": record.truncal_total_cn,
                    "truncal_minor_cn": record.truncal_minor_cn,
                    "canonical_tree": record.canonical_tree,
                    "old_tree_index": "" if old is None else old.tree_index,
                    "old_branch1_index": "" if old is None else old.branch1_index,
                    "old_branch2_index": "" if old is None else old.branch2_index,
                    "codex_branch1_index": "" if codex is None else codex.branch1_index,
                    "codex_branch2_index": "" if codex is None else codex.branch2_index,
                    "branch1_fate": vector_label(record.branch1_fate),
                    "branch2_fate": vector_label(record.branch2_fate),
                    "branch1_events": event_count(record.branch1_fate),
                    "branch2_events": event_count(record.branch2_fate),
                }
            )

    fate_out_path = Path(args.fate_out)
    fate_out_path.parent.mkdir(parents=True, exist_ok=True)
    write_fate_comparison(
        fate_out_path,
        old_records=old_records,
        codex_records=codex_records,
        sample1_cn=sample1_cn,
        sample2_cn=sample2_cn,
    )

    print(f"old_spur_complex_routes\t{len(old_records)}")
    print(f"spur_both_deletion_routes\t{len(codex_records)}")
    print(f"matched_routes\t{len(set(old_by_key) & set(codex_by_key))}")
    print(f"old_only_routes\t{len(set(old_by_key) - set(codex_by_key))}")
    print(f"codex_only_routes\t{len(set(codex_by_key) - set(old_by_key))}")
    print(f"per_canonical_tree_route_set_equal\t{canonical_tree_multiset(old_records) == canonical_tree_multiset(codex_records)}")
    print(f"fate_multiset_equal_ignoring_tree\t{fate_multiset(old_records) == fate_multiset(codex_records)}")
    print(f"comparison_tsv\t{out_path}")
    print(f"fate_vector_comparison_tsv\t{fate_out_path}")


def old_spur_complex_records(
    old_spur_src: Path,
    sample1_cn: tuple[int, int],
    sample2_cn: tuple[int, int],
    max_truncal_cn: int,
    max_branch_events: int,
) -> list[RouteRecord]:
    old_src = old_spur_src.resolve()
    if str(old_src) not in sys.path:
        sys.path.insert(0, str(old_src))
    from mrca_tree import MRCARouteTreeCollection  # type: ignore[import-not-found]

    sample1_major_minor = to_major_minor_tuple(sample1_cn)
    sample2_major_minor = to_major_minor_tuple(sample2_cn)
    records = []
    for major_cn, minor_cn in truncal_states(max_truncal_cn):
        total_cn = major_cn + minor_cn
        for tree_index, old_tree in enumerate(MRCARouteTreeCollection(major_cn, minor_cn).all_route_tree_collection):
            branch1 = old_tree.get_complex_routes_multiplication_vector(sample1_major_minor[0], sample1_major_minor[1], max_branch_events)
            branch2 = old_tree.get_complex_routes_multiplication_vector(sample2_major_minor[0], sample2_major_minor[1], max_branch_events)
            if branch1.shape[0] == 0 or branch2.shape[0] == 0:
                continue
            for i, fate1 in enumerate(branch1):
                for j, fate2 in enumerate(branch2):
                    records.append(
                        RouteRecord(
                            source="old_spur",
                            truncal_total_cn=total_cn,
                            truncal_minor_cn=minor_cn,
                            canonical_tree=f"truncal_{total_cn}:{minor_cn}_{union_tree_shape(old_tree.major_tree, old_tree.minor_tree)}",
                            branch1_fate=array_to_tuple(fate1),
                            branch2_fate=array_to_tuple(fate2),
                            branch1_index=i,
                            branch2_index=j,
                            tree_index=tree_index,
                        )
                    )
    return records


def write_fate_comparison(
    path: Path,
    old_records: list[RouteRecord],
    codex_records: list[RouteRecord],
    sample1_cn: tuple[int, int],
    sample2_cn: tuple[int, int],
) -> None:
    old = canonical_tree_multiset(old_records)
    codex = canonical_tree_multiset(codex_records)
    keys = sorted(set(old) | set(codex), key=sort_key)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "comparison_status",
                "old_count",
                "codex_count",
                "observed_sample1_cn",
                "observed_sample2_cn",
                "truncal_total_cn",
                "truncal_minor_cn",
                "canonical_tree",
                "branch1_fate",
                "branch2_fate",
                "branch1_events",
                "branch2_events",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for key in keys:
            old_count = old.get(key, 0)
            codex_count = codex.get(key, 0)
            status = "match" if old_count == codex_count else "count_mismatch"
            total_cn, minor_cn, canonical_tree, branch1_fate, branch2_fate = key
            writer.writerow(
                {
                    "comparison_status": status,
                    "old_count": old_count,
                    "codex_count": codex_count,
                    "observed_sample1_cn": total_minor_label(sample1_cn),
                    "observed_sample2_cn": total_minor_label(sample2_cn),
                    "truncal_total_cn": total_cn,
                    "truncal_minor_cn": minor_cn,
                    "canonical_tree": canonical_tree,
                    "branch1_fate": vector_label(branch1_fate),
                    "branch2_fate": vector_label(branch2_fate),
                    "branch1_events": event_count(branch1_fate),
                    "branch2_events": event_count(branch2_fate),
                }
            )


def fate_multiset(records: list[RouteRecord]) -> dict[tuple[int, int, tuple[int, ...], tuple[int, ...]], int]:
    counts: dict[tuple[int, int, tuple[int, ...], tuple[int, ...]], int] = {}
    for record in records:
        key = (record.truncal_total_cn, record.truncal_minor_cn, record.branch1_fate, record.branch2_fate)
        counts[key] = counts.get(key, 0) + 1
    return counts


def canonical_tree_multiset(records: list[RouteRecord]) -> dict[tuple[int, int, str, tuple[int, ...], tuple[int, ...]], int]:
    counts: dict[tuple[int, int, str, tuple[int, ...], tuple[int, ...]], int] = {}
    for record in records:
        key = (
            record.truncal_total_cn,
            record.truncal_minor_cn,
            record.canonical_tree,
            record.branch1_fate,
            record.branch2_fate,
        )
        counts[key] = counts.get(key, 0) + 1
    return counts


def codex_both_deletion_records(
    sample1_cn: tuple[int, int],
    sample2_cn: tuple[int, int],
    max_truncal_cn: int,
    max_branch_events: int,
) -> list[RouteRecord]:
    routes = RouteEnumerator().enumerate(
        sample1_cn=sample1_cn,
        sample2_cn=sample2_cn,
        max_truncal_cn=max_truncal_cn,
        max_branch_events=max_branch_events,
    )
    records = []
    for route in routes:
        if not (route.branch1.has_deletion and route.branch2.has_deletion):
            continue
        records.append(
            RouteRecord(
                source="spur",
                truncal_total_cn=route.truncal.truncal_total_cn,
                truncal_minor_cn=route.truncal.truncal_minor_cn,
                canonical_tree=route.truncal.canonical_key,
                branch1_fate=route.branch1.copy_fate_vector,
                branch2_fate=route.branch2.copy_fate_vector,
                branch1_index=route.branch1.branch_index,
                branch2_index=route.branch2.branch_index,
            )
        )
    return records


def truncal_states(max_truncal_cn: int) -> list[tuple[int, int]]:
    return [
        (major_cn, minor_cn)
        for major_cn in range(1, max_truncal_cn)
        for minor_cn in range(1, major_cn + 1)
        if major_cn + minor_cn <= max_truncal_cn
    ]


def array_to_tuple(values: np.ndarray) -> tuple[int, ...]:
    return tuple(int(value) for value in values)


def vector_label(values: tuple[int, ...]) -> str:
    return ",".join(str(value) for value in values)


def event_count(values: tuple[int, ...]) -> int:
    return sum(abs(value - 1) for value in values)


def sort_key(key: tuple[int, int, str, tuple[int, ...], tuple[int, ...]]) -> tuple[object, ...]:
    total_cn, minor_cn, canonical_tree, branch1_fate, branch2_fate = key
    return (total_cn, minor_cn, canonical_tree, branch1_fate, branch2_fate)


def fate_sort_key(key: tuple[int, int, tuple[int, ...], tuple[int, ...]]) -> tuple[object, ...]:
    total_cn, minor_cn, branch1_fate, branch2_fate = key
    return (total_cn, minor_cn, branch1_fate, branch2_fate)


if __name__ == "__main__":
    main()
