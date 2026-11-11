from __future__ import annotations

import argparse
import json
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

from simulate_datasets import clean_density_matrix, generate_alt_reads, generate_total_reads
from simulate_parsimony_route import (
    NON_WGD_PARSIMONY_TRUNCAL_CN_BY_PAIR,
    WGD_PARSIMONY_TRUNCAL_CN_BY_PAIR,
    get_non_wgd_parsimony_branches,
    get_non_wgd_parsimony_route_constraint_matrix,
    get_wgd_parsimony_branches,
    get_wgd_parsimony_route_constraint_matrix,
)
from spur.route import Route
from spur.sampling import JointVAFSampler, TimingSampler


PARSIMONY_TYPES = {
    "nonWGD_parsimony": (
        NON_WGD_PARSIMONY_TRUNCAL_CN_BY_PAIR,
        get_non_wgd_parsimony_branches,
        get_non_wgd_parsimony_route_constraint_matrix,
    ),
    "WGD_parsimony": (
        WGD_PARSIMONY_TRUNCAL_CN_BY_PAIR,
        get_wgd_parsimony_branches,
        get_wgd_parsimony_route_constraint_matrix,
    ),
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--n-samples", type=int, default=800)
    parser.add_argument("--n-timing-samples-per-route", type=int, default=10000)
    parser.add_argument("--n-ssnvs", type=int, default=80)
    parser.add_argument("--sample1-purity", type=float, default=0.8)
    parser.add_argument("--sample2-purity", type=float, default=0.8)
    parser.add_argument("--total-reads-mean", type=float, default=80)
    parser.add_argument("--total-reads-var", type=float, default=100)
    parser.add_argument("--density-thres", type=float, default=0.1)
    parser.add_argument("--density-thres-tol", type=float, default=1e-5)
    parser.add_argument("--scale", type=float, default=0.01)
    parser.add_argument("--seed", type=int, default=1949)
    parser.add_argument("--float-dtype", choices=["float32", "float64"], default="float64")
    parser.add_argument(
        "--route-type",
        choices=["all", "WGD_parsimony", "nonWGD_parsimony"],
        default="all",
    )
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    selected_types = PARSIMONY_TYPES
    if args.route_type != "all":
        selected_types = {args.route_type: PARSIMONY_TYPES[args.route_type]}

    for route_type, (cn_pairs, get_branches, get_constraint_matrix) in selected_types.items():
        for sample1_cn, sample2_cn in cn_pairs:
            simulate_cn_pair(
                sample1_cn,
                sample2_cn,
                route_type,
                get_branches,
                get_constraint_matrix,
                args,
                rng,
            )


def simulate_cn_pair(
    sample1_cn,
    sample2_cn,
    route_type,
    get_branches,
    get_constraint_matrix,
    args,
    rng,
):
    nt1, nb1 = sample1_cn
    nt2, nb2 = sample2_cn
    label = f"{nt1}:{nb1}_{nt2}:{nb2}_{route_type}"
    print(f"\nSimulating {label}")

    truncal, sample1_branches, sample2_branches = get_branches(nt1, nb1, nt2, nb2)
    routes = [
        Route(truncal=truncal, branch1=branch1, branch2=branch2)
        for branch1, branch2 in product(sample1_branches, sample2_branches)
    ]

    matrix = get_constraint_matrix(nt1, nb1, nt2, nb2)
    values = np.ones(matrix.shape[0])
    if matrix.shape[0] > truncal.constraint_matrix.shape[0]:
        values[-1] = 0

    q_rows = []
    for route in routes:
        timings = TimingSampler().sample_timing_vectors_from_constraints(
            matrix,
            values,
            rng,
            n_samples=args.n_timing_samples_per_route,
            burn_in=25,
            skips=15,
            dtype=np.dtype(args.float_dtype),
        )
        q = timings @ route.h(dtype=np.dtype(args.float_dtype))
        q_rows.append(q / q.sum(axis=1, keepdims=True))

    q = clean_density_matrix(np.vstack(q_rows), args.density_thres_tol)
    keep = ((q == 0) | (q >= args.density_thres - args.density_thres_tol)).all(axis=1)
    q = q[keep]
    print(f"  routes={len(routes)}, usable JMD rows={len(q)}")

    sampled_indices = rng.choice(len(q), size=args.n_samples, replace=False)
    sampled_q = q[sampled_indices]
    sampled_vafs = JointVAFSampler().sample_joint_vafs(
        sample1_total_cn=nt1,
        sample1_minor_cn=nb1,
        sample2_total_cn=nt2,
        sample2_minor_cn=nb2,
        sample1_pu=args.sample1_purity,
        sample2_pu=args.sample2_purity,
        n_samples=args.n_samples,
        scale=args.scale,
        rng=rng,
    )

    dataset_dir = args.out_dir / label
    dataset_dir.mkdir(parents=True, exist_ok=True)
    for i, (density, vafs) in enumerate(zip(sampled_q, sampled_vafs), start=1):
        total_reads1 = generate_total_reads(
            args.total_reads_mean, args.total_reads_var, args.n_ssnvs, rng
        )
        total_reads2 = generate_total_reads(
            args.total_reads_mean, args.total_reads_var, args.n_ssnvs, rng
        )
        alt_reads, _ = generate_alt_reads(
            total_reads1,
            total_reads2,
            density,
            vafs,
            rng,
        )
        pd.DataFrame(
            {
                "alt_reads_A": alt_reads[:, 0],
                "ref_reads_A": total_reads1 - alt_reads[:, 0],
                "alt_reads_B": alt_reads[:, 1],
                "ref_reads_B": total_reads2 - alt_reads[:, 1],
            }
        ).to_csv(dataset_dir / f"{i}.csv", index=False)

    np.savetxt(dataset_dir / "ground_truth_density.txt", sampled_q, delimiter=",")
    np.savetxt(dataset_dir / "ground_truth_expected_vafs.txt", sampled_vafs, delimiter=",")
    metadata = {
        "dataset_name": label,
        "route_type": route_type,
        "sample1_cn": sample1_cn,
        "sample2_cn": sample2_cn,
        "sample1_purity": args.sample1_purity,
        "sample2_purity": args.sample2_purity,
        "n_routes": len(routes),
        "n_samples": args.n_samples,
        "n_ssnvs": args.n_ssnvs,
        "seed": args.seed,
        "scale": args.scale,
    }
    (dataset_dir / "metadata.json").write_text(json.dumps(metadata, indent=2))
    print(f"  wrote {args.n_samples} datasets to {dataset_dir}")


if __name__ == "__main__":
    main()
