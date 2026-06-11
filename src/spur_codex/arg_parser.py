from __future__ import annotations

import argparse


def build_parser() -> argparse.ArgumentParser:
    
    
    parser = argparse.ArgumentParser(description="SPUR arguments.")
    
    # input data
    parser.add_argument("--input_path", type=str)
    parser.add_argument("--cn-pair", default="4:2_3:1", help="Observed CN pair in total:minor notation, e.g. 4:2_3:1.")
    parser.add_argument("--sample1_purity", type=float, default=1)
    parser.add_argument("--sample2_purity", type=float, default=1)
    parser.add_argument("--sample1_prevalence", type=float, default=1)
    parser.add_argument("--sample2_prevalence", type=float, default=1)
    
    # batch input
    parser.add_argument(
        "--batch-process",
        action="store_true",
        help="Run in batch-processing mode.",
    )
    parser.add_argument("--batch-filepath-list", type=str, help="A file where each row is a read count filepath")
    parser.add_argument(
        "--batch-input-paths",
        nargs="+",
        help="Read-count file paths for batch scoring without a filepath-list file.",
    )
    
    # data process
    parser.add_argument("--alt-reads-threshold", type=int, default=3)
    parser.add_argument("--total-reads-threshold", type=int, default=15)
    
    # Route options
    parser.add_argument(
        "--max-truncal-cn",
        type=int,
        default=None,
        help="Defaults to min(max(total CN) + 2, 6).",
    )
    parser.add_argument("--max-branch-events", type=int, default=5)
    # parser.add_argument(
    #     "--route-set",
    #     choices=["all", "pure_gain", "complex", "both_del", "one_branch_del", "sample1_branch_del", "sample2_branch_del"],
    #     default="all",
    #     help="Subset to group and optionally encode.",
    # )
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
    
    # Computation options
    parser.add_argument("--n_Monte_Carlo_samples", type=int, default=4000)
    parser.add_argument("--joint_vafs_scale", type=float, default=0.1)
    parser.add_argument("--burn-in", type=int, default=200)
    parser.add_argument("--skips", type=int, default=50)
    parser.add_argument("--float-dtype", choices=["float32", "float64"], default="float64")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument(
        "--calculate-uncertainty",
        action="store_true",
        help="Calculate scores and features of score distribution from multiple seeds.",
    )
    parser.add_argument(
        "--seed-filepath",
        type=str,
        help="Path of the file containing seeds when calculate-uncertainty is turning on.",
    )
    
    
        
    #output
    parser.add_argument("--output_dir", type=str)
    parser.add_argument(
        "--report-route-score",
        action="store_true",
        help="Generate tables of scores of enumerated route score",
    )

    return parser
