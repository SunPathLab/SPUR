from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from pathlib import Path


ROUTE_TYPES = (
    "pure_gain",
    "mixed_gain_loss",
    "conclusive_mixed_gain_loss",
)


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cn-pair-tsv",
        default=str(repo_root / "route_counts" / "route_total_event_threshold8_counts.tsv"),
    )
    parser.add_argument("--out-dir", default=str(repo_root / "simulated_small"))
    parser.add_argument("--route-filter", default="total_branch_events")
    parser.add_argument("--route-filter-threshold", default="8")
    parser.add_argument("--n-samples", default="800")
    parser.add_argument("--python", default=sys.executable)
    args = parser.parse_args()

    cn_pairs = read_cn_pairs(Path(args.cn_pair_tsv))
    simulate_script = repo_root / "scripts" / "simulate_datasets.py"

    env = os.environ.copy()
    env["PYTHONPATH"] = str(repo_root / "src")

    print(f"CN pairs: {len(cn_pairs)}")
    print(f"Route types: {', '.join(ROUTE_TYPES)}")
    print(f"Output: {args.out_dir}")

    for cn_pair in cn_pairs:
        for route_type in ROUTE_TYPES:
            output_dir = Path(args.out_dir) / f"{cn_pair}_{route_type}"
            if output_dir.exists():
                print(f"Skipping existing output: {output_dir}")
                continue

            cmd = [
                args.python,
                str(simulate_script),
                "--cn-pair",
                cn_pair,
                "--out-dir",
                args.out_dir,
                "--route-filter",
                args.route_filter,
                "--route-filter-threshold",
                args.route_filter_threshold,
                "--n-samples",
                args.n_samples,
                "--route-type",
                route_type,
            ]

            print("\n" + "=" * 80)
            print(f"Simulating {cn_pair} / {route_type}")
            print(" ".join(cmd))
            subprocess.run(cmd, check=True, env=env)


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


if __name__ == "__main__":
    main()
