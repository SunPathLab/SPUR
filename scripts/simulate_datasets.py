from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from spur_codex.compatibility import PureGainCompatibilityChecker
from spur_codex.enumeration import RouteEnumerator
from spur_codex.notation import parse_between_sample_cn, total_minor_label
from spur_codex.route import Route
from spur_codex.route_filter import split_routes, RouteFilter
from spur_codex.simulation import sample_route_jmd
from spur_codex.sampling import TimingSampler, JointVAFSampler


def main() -> None:
    parser = argparse.ArgumentParser()
    # parser.add_argument("--cn-tsv", type = str)
    parser.add_argument("--cn-pair", default="4:2_3:1", help="Observed CN pair in total:minor notation, e.g. 4:2_3:1.")
    parser.add_argument("--sample1_purity", type=float, default=0.8)
    parser.add_argument("--sample2_purity", type=float, default=0.8)
    parser.add_argument("--out-dir", type = str)
    parser.add_argument("--seed", type=int, default=1949)
    
    # route limitation parameters
    parser.add_argument("--max-truncal-cn", type=int, default=None)
    parser.add_argument(
        "--route-filter",
        choices=[
            "none",
            "total_events",
            "max_branch_events",
            "total_branch_events",
            "single_sample_total_events",
        ],
        default="none",
    )
    parser.add_argument(
        "--route-filter-threshold",
        type=int,
        default=None,
    )
    # simulation parameters
    parser.add_argument("--route-type", type=str, default="pure_gain", choices=["mixed_gain_loss", "pure_gain", "conclusive_mixed_gain_loss"])
    parser.add_argument("--n-samples", type=int, default=800)
    parser.add_argument("--n-ssnvs", type=int, default=80)
    parser.add_argument("--total-reads-mean", type=float, default=80)
    parser.add_argument("--total-reads-var", type=float, default=100)
    parser.add_argument("--float-dtype", choices=["float32", "float64"], default="float64")
    # hit-and-run sampling parameters
    parser.add_argument("--n-timing-samples-per-route", type=int, default=None)
    parser.add_argument("--max-timing-samples-per-route", type=int, default=None)
    parser.add_argument("--hit-and-run-skip-steps", type=int, default=15)
    parser.add_argument("--hit-and-run-burn-in", type=int, default=25)
    # sample filtering parameters. Ensure sampled JMDs have significant density in non-zero density clusters.
    # Guarantee the existence of enough number of SSNVs belonging to the JM.
    parser.add_argument("--density-thres", type=float, default=0.1)
    parser.add_argument("--density-thres-tol", type=float, default=1e-5) 
    # parameters for generating synthetic datasets according to the sampled JMDs, i.e. how much noise to add to the true joint VAFs when generating alt reads according to the sampled densities, and the purity of the two samples which also affects the true joint VAFs. Note that these parameters do not affect the simulated JMDs, but only affect the synthetic datasets generated from the simulated JMDs.
    parser.add_argument("--scale", type=float, default=0.01)

    
    
    args = parser.parse_args()
    
    rng = np.random.default_rng(args.seed)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    print("Starting simulation")
    print(f"  cn_pair={args.cn_pair}")
    print(f"  route_type={args.route_type}")
    print(f"  route_filter={args.route_filter}, threshold={args.route_filter_threshold}")
    print(f"  n_samples={args.n_samples}, n_ssnvs={args.n_ssnvs}, seed={args.seed}")
    
    
    ########### route enumeration (with filter)
    
    # first determine the set of truncal_cn for enumeration
    sample1_cn, sample2_cn = parse_between_sample_cn(args.cn_pair)
    max_truncal_cn = args.max_truncal_cn
    if max_truncal_cn is None:
        max_truncal_cn = min(max(sample1_cn[0], sample2_cn[0]) + 2, 6)
    print(f"Resolved CN pair: sample1={sample1_cn}, sample2={sample2_cn}, max_truncal_cn={max_truncal_cn}")
    
    # second, enumerate routes
    print("Enumerating routes ...")
    routes = RouteEnumerator().enumerate(
        sample1_cn=sample1_cn,
        sample2_cn=sample2_cn,
        max_truncal_cn=max_truncal_cn
    )
    print(f"  enumerated routes: {len(routes)}")
    
    # third, filter the routes in terms of: 1. route types; 2. event numbers
    
    # 1. route types
    print("Splitting routes by route type ...")
    split = split_routes(routes)
    print(f"  pure_gain routes: {len(split.pure_gain)}")
    print(f"  mixed_gain_loss routes: {len(split.mixed_gain_loss)}")
    if args.route_type == "pure_gain":
        routes = split.pure_gain
    if args.route_type in ["mixed_gain_loss", "conclusive_mixed_gain_loss"]:
        routes = split.mixed_gain_loss
    print(f"  selected routes after route_type={args.route_type}: {len(routes)}")
    
    
    # 2. event numbers (according to different filter)
    if args.route_filter != "none":
        print(f"Applying route filter: {args.route_filter} <= {args.route_filter_threshold} ...")
        if args.route_filter_threshold is None:
            raise ValueError("--route-filter-threshold is required when --route-filter is not 'none'.")
        route_filter_map = {
            "total_events": RouteFilter.filter_by_total_events,
            "max_branch_events": RouteFilter.filter_by_max_branch_events,
            "total_branch_events": RouteFilter.filter_by_total_branch_events,
            "single_sample_total_events": RouteFilter.filter_by_single_sample_total_events,
        }
        routes = route_filter_map[args.route_filter](
            tuple(routes),
            args.route_filter_threshold,
        )
        print(f"  routes after event filter: {len(routes)}")
    
    
    ########### sample timings, calculate Q's by TH, filter Q
    
    # For pure-gain simulation, we sample more timings per route, since there are only a few pure-gain routes for each cn pair
    if args.n_timing_samples_per_route is not None:
        n_timing_samples_per_route = args.n_timing_samples_per_route
    elif args.route_type == "pure_gain":
        n_timing_samples_per_route = 2000
    else:
        if max(sample1_cn[0], sample2_cn[0]) >= 6:
            n_timing_samples_per_route = 20
        elif max(sample1_cn[0], sample2_cn[0]) >= 5:
            n_timing_samples_per_route = 50
        else:
            n_timing_samples_per_route = 100
    
    print(f"Sampling timings and computing Q for {len(routes)} routes ...")
    progress_every = max(1, len(routes) // 10) if routes else 1

    while True:
        print(f"  timing samples per route: {n_timing_samples_per_route}")
        q_rows = []
        for route_i, route in enumerate(routes, start=1):
            if route_i == 1 or route_i == len(routes) or route_i % progress_every == 0:
                print(f"  route {route_i}/{len(routes)}")
            T_route = TimingSampler().sample_timing_vectors(
                truncal=route.truncal,
                rng=rng,
                n_samples=n_timing_samples_per_route,
                burn_in=args.hit_and_run_burn_in,
                skips=args.hit_and_run_skip_steps,
                dtype=np.dtype(args.float_dtype),
            )
            Q_route = T_route @ route.h()
            row_sums = Q_route.sum(axis=1)
            Q_route = Q_route[row_sums > 0]
            Q_route = Q_route / Q_route.sum(axis=1, keepdims=True)
            q_rows.append(Q_route)

        Q = np.vstack(q_rows)
        Q = clean_density_matrix(Q, args.density_thres_tol)
        pass_mask = (
            (Q <= args.density_thres_tol)
            | (Q > args.density_thres - args.density_thres_tol)
        ).all(axis=1)
        Q = Q[pass_mask]
        Q = np.unique(np.round(Q, decimals=12), axis=0)

        if args.route_type == "conclusive_mixed_gain_loss":
            checker = PureGainCompatibilityChecker(
                max_truncal_cn=max_truncal_cn,
                route_filter="total_events",
                route_filter_threshold=8,
            )
            compatible = np.array(
                [checker.is_compatible(q, sample1_cn, sample2_cn) for q in Q],
                dtype=bool,
            )
            Q = Q[~compatible]

        print(f"  usable JMD rows: {len(Q)}")

        if len(Q) >= args.n_samples:
            break

        if (
            args.max_timing_samples_per_route is not None
            and n_timing_samples_per_route >= args.max_timing_samples_per_route
        ):
            raise ValueError(
                "Unable to generate enough usable JMD rows before reaching "
                f"--max-timing-samples-per-route={args.max_timing_samples_per_route}."
            )

        if len(Q) == 0:
            n_timing_samples_per_route *= 2
        else:
            n_timing_samples_per_route = int(
                np.ceil(
                    n_timing_samples_per_route
                    * args.n_samples
                    / len(Q)
                    * 1.2
                )
            )
        print("  insufficient usable JMD rows; sampling more timings")

    print(f"Sampling {args.n_samples} Q rows for read simulation ...")
    Q_sampled_indices = rng.choice(len(Q), size=args.n_samples, replace=False)
    Q_sampled = Q[Q_sampled_indices]
    
    ########## Reads simulation
    # sampled joint VAFs for each sampled Q
    print("Sampling expected joint VAFs ...")
    F_sampled = JointVAFSampler().sample_joint_vafs(
        sample1_total_cn=sample1_cn[0], sample1_minor_cn=sample1_cn[1],
        sample2_total_cn=sample2_cn[0], sample2_minor_cn=sample2_cn[1],
        sample1_pu=args.sample1_purity, sample2_pu=args.sample2_purity,
        n_samples=args.n_samples, scale=args.scale,
        rng=rng
    )
    print(f"  F_sampled shape: {F_sampled.shape}")
    
    # Simulate read counts given Q and F
    
    print("Generating read-count data frames ...")
    df_list = []
    for i in range(args.n_samples):
        if i == 0 or i + 1 == args.n_samples or (i + 1) % max(1, args.n_samples // 10) == 0:
            print(f"  dataset {i + 1}/{args.n_samples}")
        q = Q_sampled[i]
        f = F_sampled[i]
        total_reads_A = generate_total_reads(mean=args.total_reads_mean, var=args.total_reads_var, size=args.n_ssnvs, rng=rng)
        total_reads_B = generate_total_reads(mean=args.total_reads_mean, var=args.total_reads_var, size=args.n_ssnvs, rng=rng)
        alt_reads, cluster_assignments = generate_alt_reads(total_reads_A=total_reads_A,
                                        total_reads_B=total_reads_B,
                                        density=q,
                                        true_joint_vafs=f,
                                        rng=rng)
        df = {'alt_reads_A': alt_reads[:,0],
                'ref_reads_A': total_reads_A - alt_reads[:,0],
                'alt_reads_B': alt_reads[:,1],
                'ref_reads_B': total_reads_B - alt_reads[:,1]}
        
        df = pd.DataFrame(df)
        df_list.append(df)
        
    # write metadata
    dataset_dir = Path(out_dir) / f"{args.cn_pair}_{args.route_type}"
    dataset_dir.mkdir(parents=True, exist_ok=True)
    print(f"Writing outputs to {dataset_dir} ...")
    metadata = {
        "dataset_info": {
            "dataset_name": f"{args.cn_pair}_{args.route_type}".split("/")[-1],
            "date_generated": pd.Timestamp.now().isoformat(),
            "description": "Simulated shared mutation read count dataset",
            "version": "1.0"
        },

        "biological_parameters": {
            "cn_pair": args.cn_pair,
            "sample1_pu": float(args.sample1_purity),
            "sample2_pu": float(args.sample2_purity),
            "max_truncal_cn": max_truncal_cn,
        },

        "simulation_settings": {
            "seed": int(args.seed),
            "route_type": args.route_type,
            "n_timing_samples_per_route": int(n_timing_samples_per_route),
            "hit_and_run_skip_steps": int(args.hit_and_run_skip_steps),
            "n_samples": int(args.n_samples),
            "density_thres": float(args.density_thres),
            "density_thres_tol": float(args.density_thres_tol),
            "scale": float(args.scale),
            "total_reads_mean": float(args.total_reads_mean),
            "total_reads_var": float(args.total_reads_var),
            "n_ssnvs": int(args.n_ssnvs),
        },

        "ground_truth": {
            "route_type": args.route_type
        }
    }
    
    # write to dataset_dir
    for i, dir in enumerate(df_list):
        dir.to_csv(f"{str(dataset_dir)}/{i+1}.csv", index=False)
    # write metadata
    with open(f"{str(dataset_dir)}/metadata.json", "w") as f:
        json.dump(metadata, f)
    # write ground truth density
    ground_truth_density = np.array(Q_sampled)
    np.savetxt(f"{str(dataset_dir)}/ground_truth_density.txt", ground_truth_density, delimiter=',')
    # write ground truth expected VAFs
    np.savetxt(f"{str(dataset_dir)}/ground_truth_expected_vafs.txt", F_sampled, delimiter=',')
    print("Simulation complete")
    
    return
    
    
def select_cn_pairs(path: Path, n_pairs: int, max_total_cn: int) -> list[tuple[tuple[int, int], tuple[int, int]]]:
    counts: Counter[tuple[tuple[int, int], tuple[int, int]]] = Counter()
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sample1, sample2 = parse_between_sample_cn(row["between_sample_cn"])
            if sample1[0] > max_total_cn or sample2[0] > max_total_cn:
                continue
            counts[canonical_pair(sample1, sample2)] += 1
    return [pair for pair, _ in counts.most_common(n_pairs)]

def canonical_pair(sample1_cn: tuple[int, int], sample2_cn: tuple[int, int]) -> tuple[tuple[int, int], tuple[int, int]]:
    return tuple(sorted((sample1_cn, sample2_cn), key=lambda x: (x[0], x[1]), reverse=True))  # type: ignore[return-value]


# from chatGPT
# try X = np.array([1e-8, 0.2, -1e-7, 0.5])
def clean_density_matrix(X, tol=1e-12):
    X = np.where(np.abs(X) < tol, 0.0, X)
    X = np.maximum(X, 0.0)
    row_sums = X.sum(axis=1, keepdims=True)
    return X / row_sums


# generate total reads by negative inomial distr
def generate_total_reads(mean, var, size, rng):
    if var < mean:
        print("Error: mean is less than var")
        return None
    p = mean/var
    n = mean**2/(var-mean)
    return rng.negative_binomial(n, p, size)

# generate alt reads by binomial distr, given total reads, Q, and true joint VAFs
def generate_alt_reads(total_reads_A, total_reads_B, density, true_joint_vafs, rng):
    n_snvs = len(total_reads_A)
    density = np.asarray(density, dtype=float)
    density = np.maximum(density, 0.0)
    density_sum = density.sum()
    if density_sum <= 0:
        raise ValueError("density must have positive total mass.")
    density = density / density_sum
    true_joint_vafs = np.asarray(true_joint_vafs, dtype=float)
    n_joint_groups = density.size
    if true_joint_vafs.ndim == 1:
        if true_joint_vafs.size != 2 * n_joint_groups:
            raise ValueError(
                f"true_joint_vafs has length {true_joint_vafs.size}, expected {2 * n_joint_groups} for density length {n_joint_groups}."
            )
        sample1_vafs = true_joint_vafs[:n_joint_groups]
        sample2_vafs = true_joint_vafs[n_joint_groups:]
    elif true_joint_vafs.shape == (2, n_joint_groups):
        sample1_vafs = true_joint_vafs[0]
        sample2_vafs = true_joint_vafs[1]
    else:
        raise ValueError(
            f"true_joint_vafs must have shape ({2 * n_joint_groups},) or (2, {n_joint_groups}); got {true_joint_vafs.shape}."
        )

    # https://stackoverflow.com/questions/4265988/generate-random-numbers-with-a-given-numerical-distribution
    # https://cs.stackexchange.com/questions/3227/uniform-sampling-from-a-simplex
    cluster_assignments = rng.choice(np.arange(n_joint_groups), size=n_snvs, p=density)
    

    alt_reads = np.zeros((n_snvs, 2), dtype=int) # reads in sample A and B
    for i in range(len(alt_reads)):
        # sample A reads
        try:
            #print(r[i]*true_joint_vafs.shape[2]+c[i] == cluster_assignments[i])
            alt_reads[i, 0] = rng.binomial(total_reads_A[i], sample1_vafs[cluster_assignments[i]])
        except:
            print(f"total reads A: {total_reads_A[i]}, p: {sample1_vafs[cluster_assignments[i]]}")
        # sample B reads
        try:
            alt_reads[i, 1] = rng.binomial(total_reads_B[i], sample2_vafs[cluster_assignments[i]])
        except:
            print(f"total reads B: {total_reads_B[i]}, p: {sample2_vafs[cluster_assignments[i]]}")
    return alt_reads, cluster_assignments


if __name__ == "__main__":
    main()
