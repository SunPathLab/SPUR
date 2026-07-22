from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROUTE_TYPES = ("pure_gain", "mixed_gain_loss", "conclusive_mixed_gain_loss")


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
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--n-samples", type=int, default=800)
    parser.add_argument("--max-observed-cn", type=int, default=6)
    parser.add_argument("--total-event-threshold", type=int, default=8)
    parser.add_argument("--max-timing-samples-per-route", type=int, default=6400)
    parser.add_argument("--seed", type=int, default=1949)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument(
        "--route-types",
        nargs="+",
        choices=ROUTE_TYPES,
        default=list(ROUTE_TYPES),
    )
    args = parser.parse_args()

    profiles = pd.read_csv(args.cn_profile_tsv, sep="\t")
    loh_loh = profiles[(profiles["minor1"] == 0) & (profiles["minor2"] == 0)]
    loh_loh = loh_loh[
        (loh_loh["major1"] <= args.max_observed_cn)
        & (loh_loh["major2"] <= args.max_observed_cn)
        & (loh_loh["major1"] > 0)
        & (loh_loh["major2"] > 0)
    ]
    cn_pairs = sorted(
        {
            tuple(sorted((int(row.major1), int(row.major2)), reverse=True))
            for row in loh_loh.itertuples()
        }
    )

    simulate_script = repo_root / "scripts" / "simulate_datasets.py"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(repo_root / "src")

    print(f"LOH-LOH CN pairs: {len(cn_pairs)}", flush=True)
    print(f"Profiles per class: {args.n_samples}", flush=True)
    print(f"Output: {args.out_dir}", flush=True)

    for sample1_total_cn, sample2_total_cn in cn_pairs:
        cn_pair = f"{sample1_total_cn}:0_{sample2_total_cn}:0"
        for route_type in args.route_types:
            dataset_dir = args.out_dir / f"{cn_pair}_{route_type}"
            if dataset_dir.exists():
                print(f"Skipping existing output: {dataset_dir}")
                continue
            
            command = [
                args.python,
                "-u",
                str(simulate_script),
                "--cn-pair",
                cn_pair,
                "--out-dir",
                str(args.out_dir),
                "--route-filter",
                "total_events",
                "--route-filter-threshold",
                str(args.total_event_threshold),
                "--n-samples",
                str(args.n_samples),
                "--route-type",
                route_type,
                "--seed",
                str(args.seed),
                "--max-timing-samples-per-route",
                str(args.max_timing_samples_per_route),
            ]
            if route_type == "mixed_gain_loss" and sample1_total_cn >= 6:
                command.extend(["--n-timing-samples-per-route", "100"])
            print(f"Simulating {cn_pair} / {route_type}", flush=True)
            try:
                subprocess.run(command, check=True, env=env)
            except subprocess.CalledProcessError as exc:
                print(
                    f"Skipping failed simulation for {cn_pair} / {route_type}: exit code {exc.returncode}",
                    flush=True,
                )


if __name__ == "__main__":
    main()
