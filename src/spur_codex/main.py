from __future__ import annotations

import time

from pathlib import Path

from .arg_parser import build_parser
from .enumeration import RouteEnumerator
from .notation import parse_between_sample_cn
from .scoring import Scorer
from .route_filter import RouteFilter

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

def main() -> int:
    
    ########## parse arguments
    
    # see details in arg_parser
    args = build_parser().parse_args()

    

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
    
    # then filter the routes according to the 
    if args.route_filter != "none":
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
    
    
    ##### route scoring
    shared_ssnvs_reads_df = pd.read_csv(args.input_path)
    
    # get read counts
    alt1, dep1, alt2, dep2 = process_data(shared_ssnvs_reads_df, 
                                          alt_thres=args.alt_reads_threshold, 
                                          depth_thres=args.total_reads_threshold)
    
    # single seed, or multiple seeds.
    if args.calculate_uncertainty:
        seeds = [int(line) for line in Path(args.seed_filepath).read_text().splitlines() if line.strip()]
    else:
        seeds = [args.seed]
        
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
            alt1=alt1, dep1=dep1, alt2=alt2, dep2=dep2,
            routes=routes,
            sample1_pu=args.sample1_purity, sample2_pu=args.sample2_purity,
            sample1_prevalence=args.sample1_prevalence,
            sample2_prevalence=args.sample2_prevalence,
            burn_in=args.burn_in, skips=args.skips,
            n_samples=args.n_Monte_Carlo_samples,
            rng=rng,
            dtype=args.float_dtype,
        )
        runtime_seconds = time.perf_counter() - start_time
        
        feature_row = {"seed": seed, "runtime_seconds": runtime_seconds, **score_results.score_features()}

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
        
        
        
    ##### Write output
    
    # get the name of the input
    input_path = Path(args.input_path)
    input_filename = input_path.stem # no extension (e.g. no ".tsv" or ".csv")
    
    # get the output dir, create if necessary
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # feature table (supporting multiple seeds)
    feature_df.to_csv(output_dir / f"{input_filename}.tsv", sep="\t", index=False)
    
    # route score table (supporting multiple seeds)
    if args.report_route_score:
        route_score_df.to_csv(output_dir / f"{input_filename}_route_scores.tsv", sep="\t", index=False)
    
    return 0


if __name__ == "__main__":
    main()
