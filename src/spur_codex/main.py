from __future__ import annotations

import time

from pathlib import Path

from .arg_parser import build_parser
from .enumeration import RouteEnumerator
from .notation import parse_between_sample_cn
from .scoring import Scorer
from .route_filter import RouteFilter, filter_routes, split_routes



import numpy as np
import pandas as pd

def process_data(df, alt_thres = 3, depth_thres = 15):
    """Process the original data frame, where there are 4 cols of alt and ref reads. Filter the SNV based on threshold
    """
    r_a = df.iloc[:,0]
    d_a = df.iloc[:,0]+df.iloc[:,1]
    r_b = df.iloc[:,2]
    d_b = df.iloc[:,2]+df.iloc[:,3]
    
    cleaned_idices = (d_a >= depth_thres) & (d_b >= depth_thres) & (r_a >= alt_thres) & (r_b >= alt_thres)
    
    return r_a[cleaned_idices].to_numpy(), d_a[cleaned_idices].to_numpy(), r_b[cleaned_idices].to_numpy(), d_b[cleaned_idices].to_numpy()


def load_input_records(args) -> list[dict[str, str | np.ndarray]]:
    if not args.batch_process:
        filepaths = [args.input_path]
    elif args.batch_input_paths is not None:
        filepaths = args.batch_input_paths
    else:
        path_file = Path(args.batch_filepath_list)
        filepaths = [
            line.strip()
            for line in path_file.read_text().splitlines()
            if line.strip()
        ]

    inputs = []
    for filepath in filepaths:
        shared_ssnvs_reads_df = pd.read_csv(filepath)
        alt1, dep1, alt2, dep2 = process_data(
            shared_ssnvs_reads_df,
            alt_thres=args.alt_reads_threshold,
            depth_thres=args.total_reads_threshold,
        )
        inputs.append(
            {
                "input_path": filepath,
                "alt1": alt1,
                "dep1": dep1,
                "alt2": alt2,
                "dep2": dep2,
            }
        )
    return inputs


def print_run_parameters(args, inputs, seeds, routes, max_truncal_cn) -> None:
    print("SPUR Codex run parameters")
    print(f"  cn_pair: {args.cn_pair}")
    print(f"  sample1_purity: {args.sample1_purity}")
    print(f"  sample2_purity: {args.sample2_purity}")
    print(f"  sample1_prevalence: {args.sample1_prevalence}")
    print(f"  sample2_prevalence: {args.sample2_prevalence}")
    print(f"  max_truncal_cn: {max_truncal_cn}")
    print(f"  max_branch_events: {args.max_branch_events}")
    print(f"  route_filter: {args.route_filter}")
    print(f"  route_filter_threshold: {args.route_filter_threshold}")
    print(f"  mixed_gain_loss_route_definition: {args.mixed_gain_loss_route_definition}")
    print(f"  n_routes_scored: {len(routes)}")
    print(f"  n_inputs: {len(inputs)}")
    print(f"  seeds: {','.join(map(str, seeds))}")
    print(f"  n_Monte_Carlo_samples: {args.n_Monte_Carlo_samples}")
    print(f"  joint_vafs_scale: {args.joint_vafs_scale}")
    print(f"  burn_in: {args.burn_in}")
    print(f"  skips: {args.skips}")
    print(f"  float_dtype: {args.float_dtype}")
    print(f"  alt_reads_threshold: {args.alt_reads_threshold}")
    print(f"  total_reads_threshold: {args.total_reads_threshold}")
    print(f"  output_dir: {args.output_dir}")
    print(f"  report_route_score: {args.report_route_score}", flush=True)
    


def main() -> int:
    
    ########## parse arguments
    
    # see details in arg_parser
    args = build_parser().parse_args()
    
    # error handling
    if args.batch_process is True:
        if args.batch_filepath_list is None and args.batch_input_paths is None:
            raise ValueError("--batch-filepath-list or --batch-input-paths is required when --batch-process is set.")
        if args.batch_filepath_list is not None and args.batch_input_paths is not None:
            raise ValueError("Use only one of --batch-filepath-list or --batch-input-paths.")
    else:
        if args.input_path is None:
            raise ValueError("--input_path is required.")
        
    if args.route_filter != "none":
        if args.route_filter_threshold is None:
            raise ValueError("--route-filter-threshold is required when --route-filter is not 'none'.")

    if args.output_dir is None:
        raise ValueError("--output_dir is required.")

    if args.calculate_uncertainty and args.seed_filepath is None:
        raise ValueError("--seed-filepath is required when --calculate-uncertainty is set.")
    
    
    ########## route enumeration
    # first determine the set of truncal_cn for enumeration
    sample1_cn, sample2_cn = parse_between_sample_cn(args.cn_pair)
    max_truncal_cn = args.max_truncal_cn
    if max_truncal_cn is None:
        max_truncal_cn = min(max(sample1_cn[0], sample2_cn[0]) + 2, 6)
    # second, enumerate routes
    routes = RouteEnumerator().enumerate(
        sample1_cn=sample1_cn,
        sample2_cn=sample2_cn,
        max_truncal_cn=max_truncal_cn,
        max_branch_events=args.max_branch_events,
    )
    
    
    if args.route_filter != "none":
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
    
    # this is only for test. use old complex definition
    if args.mixed_gain_loss_route_definition == "both_branch_del":
        # split routes into pure-gain and mixed-gain-loss
        split = split_routes(routes)
        pure_gain_routes = split.pure_gain
        mixed_gain_loss_routes = split.mixed_gain_loss
        # then get routes with both branches with DEL
        both_branch_del_routes = filter_routes(
            mixed_gain_loss_routes,
            predicate=lambda r: r.branch1.has_deletion and r.branch2.has_deletion,
        )
        # combine pure-gain routes and both-branch-del routes
        routes = pure_gain_routes + both_branch_del_routes
    if not routes:
        raise ValueError("No routes remain after enumeration/filtering.")
    
    
    ##### route scoring
    
    # single input, or batch input
    inputs = load_input_records(args)
    
    # single seed, or multiple seeds.
    if args.calculate_uncertainty:
        seeds = [int(line) for line in Path(args.seed_filepath).read_text().splitlines() if line.strip()]
    else:
        seeds = [args.seed]
        
    print_run_parameters(args, inputs, seeds, routes, max_truncal_cn)
        
    ##### Write output

    # get the output dir, create if necessary
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for input_record in inputs:
        # initialize result data frame
        feature_df = None
        if args.report_route_score:
            route_score_df = None

        for seed in seeds:
            # create rng
            rng = np.random.default_rng(seed)

            # score enumerated routes
            start_time = time.perf_counter()
            score_results = Scorer().score(
                alt1=input_record["alt1"],
                dep1=input_record["dep1"],
                alt2=input_record["alt2"],
                dep2=input_record["dep2"],
                routes=routes,
                sample1_pu=args.sample1_purity, sample2_pu=args.sample2_purity,
                sample1_prevalence=args.sample1_prevalence,
                sample2_prevalence=args.sample2_prevalence,
                burn_in=args.burn_in, skips=args.skips,
                n_samples=args.n_Monte_Carlo_samples,
                scale=args.joint_vafs_scale,
                rng=rng,
                dtype=args.float_dtype,
            )
            runtime_seconds = time.perf_counter() - start_time
            print(
                f"Scored {input_record['input_path']} "
                f"seed={seed} runtime_seconds={runtime_seconds:.3f}",
                flush=True,
            )

            feature_row = {
                "input_path": input_record["input_path"],
                "seed": seed,
                "runtime_seconds": runtime_seconds,
                **score_results.score_features(),
            }

            # create feature_df if it does not exist
            if feature_df is None:
                feature_df = pd.DataFrame([feature_row])

            # concatenate the new row to the feature_df
            else:
                feature_df = pd.concat([feature_df,
                                        pd.DataFrame([feature_row])], axis=0, ignore_index=True)

            if args.report_route_score:
                # create new score table if not created
                if route_score_df is None:
                    route_score_df = pd.DataFrame(score_results.route_score_table(seed))
                # merge to the existing table if exsits
                else:
                    route_score_df = route_score_df.merge(
                        pd.DataFrame(score_results.route_score_table(seed)),
                        on=["truncal", "branch1", "branch2"],
                        how="outer",
                    )

        # get the name of the input
        input_path = Path(str(input_record["input_path"]))
        input_filename = input_path.stem # no extension (e.g. no ".tsv" or ".csv")

        # feature table (supporting multiple seeds)
        feature_df.to_csv(output_dir / f"{input_filename}.tsv", sep="\t", index=False)

        # route score table (supporting multiple seeds)
        if args.report_route_score:
            route_score_df.to_csv(output_dir / f"{input_filename}_route_scores.tsv", sep="\t", index=False)
    
    return 0


if __name__ == "__main__":
    main()
