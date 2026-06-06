from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import numpy as np

from spur_codex.compatibility import PureGainCompatibilityChecker
from spur_codex.enumeration import RouteEnumerator
from spur_codex.notation import parse_between_sample_cn, total_minor_label
from spur_codex.objects import Route
from spur_codex.route_filter import split_routes
from spur_codex.simulation import sample_route_jmd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cn-tsv", default="/home/yunong-xia/Downloads/May_27_real_data_indices.tsv")
    parser.add_argument("--out-dir", default="simulated_small")
    parser.add_argument("--n-cn-pairs", type=int, default=3)
    parser.add_argument("--n-per-class", type=int, default=20)
    parser.add_argument("--max-total-cn", type=int, default=5)
    parser.add_argument("--max-mrca-cn", type=int, default=6)
    parser.add_argument("--max-branch-events", type=int, default=5)
    parser.add_argument("--seed", type=int, default=20260604)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cn_pairs = select_cn_pairs(Path(args.cn_tsv), args.n_cn_pairs, args.max_total_cn)
    checker = PureGainCompatibilityChecker(args.max_mrca_cn, args.max_branch_events)
    enumerator = RouteEnumerator()

    rows = []
    for pair_index, (sample1_cn, sample2_cn) in enumerate(cn_pairs):
        all_routes = enumerator.enumerate(sample1_cn, sample2_cn, args.max_mrca_cn, args.max_branch_events)
        split = split_routes(all_routes)
        routes_by_class = {
            "pure_gain": list(split.pure_gain),
            "mixed_gain_loss": list(split.mixed_gain_loss),
        }
        if not routes_by_class["pure_gain"] or not routes_by_class["mixed_gain_loss"]:
            continue

        for label in ("pure_gain", "mixed_gain_loss"):
            rows.extend(
                simulate_class(
                    routes=routes_by_class[label],
                    label=label,
                    sample1_cn=sample1_cn,
                    sample2_cn=sample2_cn,
                    pair_index=pair_index,
                    n_examples=args.n_per_class,
                    out_dir=out_dir,
                    rng=rng,
                    checker=checker,
                    conclusive_required=False,
                )
            )

        rows.extend(
            simulate_class(
                routes=routes_by_class["mixed_gain_loss"],
                label="conclusive_mixed_gain_loss",
                sample1_cn=sample1_cn,
                sample2_cn=sample2_cn,
                pair_index=pair_index,
                n_examples=args.n_per_class,
                out_dir=out_dir,
                rng=rng,
                checker=checker,
                conclusive_required=True,
            )
        )

    write_metadata(out_dir / "metadata.tsv", rows)
    with (out_dir / "cn_pairs.json").open("w") as f:
        json.dump([{"sample1_cn": a, "sample2_cn": b} for a, b in cn_pairs], f, indent=2)


def select_cn_pairs(path: Path, n_pairs: int, max_total_cn: int) -> list[tuple[tuple[int, int], tuple[int, int]]]:
    counts: Counter[tuple[tuple[int, int], tuple[int, int]]] = Counter()
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample1, sample2 = parse_between_sample_cn(row["between_sample_cn"])
            if sum(sample1) > max_total_cn or sum(sample2) > max_total_cn:
                continue
            counts[canonical_pair(sample1, sample2)] += 1
    return [pair for pair, _ in counts.most_common(n_pairs)]


def canonical_pair(sample1_cn: tuple[int, int], sample2_cn: tuple[int, int]) -> tuple[tuple[int, int], tuple[int, int]]:
    return tuple(sorted((sample1_cn, sample2_cn), key=lambda x: (sum(x), x[1], x[0]), reverse=True))  # type: ignore[return-value]


def simulate_class(
    routes: list[Route],
    label: str,
    sample1_cn: tuple[int, int],
    sample2_cn: tuple[int, int],
    pair_index: int,
    n_examples: int,
    out_dir: Path,
    rng: np.random.Generator,
    checker: PureGainCompatibilityChecker,
    conclusive_required: bool,
) -> list[dict[str, object]]:
    rows = []
    attempts = 0
    max_attempts = max(1000, 100 * n_examples)
    while len(rows) < n_examples and attempts < max_attempts:
        attempts += 1
        route = routes[int(rng.integers(len(routes)))]
        timing, q = sample_route_jmd(route, sample1_cn, sample2_cn, rng)
        compatible = checker.is_compatible(q, sample1_cn, sample2_cn)
        if conclusive_required and compatible:
            continue

        example_id = f"pair{pair_index}_{label}_{len(rows):03d}"
        npz_name = f"{example_id}.npz"
        np.savez_compressed(
            out_dir / npz_name,
            jmd=q.astype(np.float32),
            timing=timing.astype(np.float32),
            branch1_copy_fate=np.asarray(route.branch1.copy_fate_vector, dtype=np.int16),
            branch2_copy_fate=np.asarray(route.branch2.copy_fate_vector, dtype=np.int16),
        )
        rows.append(
            {
                "example_id": example_id,
                "class_label": label,
                "npz_path": npz_name,
                "sample1_cn": total_minor_label(sample1_cn),
                "sample2_cn": total_minor_label(sample2_cn),
                "sample1_major_cn": sample1_cn[0],
                "sample1_minor_cn": sample1_cn[1],
                "sample2_major_cn": sample2_cn[0],
                "sample2_minor_cn": sample2_cn[1],
                "mrca_major_cn": route.truncal.truncal_major_cn,
                "mrca_minor_cn": route.truncal.truncal_minor_cn,
                "truncal_key": route.truncal_key,
                "branch1_index": route.branch1.branch_index,
                "branch2_index": route.branch2.branch_index,
                "branch1_events": route.branch1.n_events,
                "branch2_events": route.branch2.n_events,
                "pure_gain_compatible": compatible,
            }
        )
    return rows


def write_metadata(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
