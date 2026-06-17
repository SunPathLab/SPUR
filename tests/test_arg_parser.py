from spur_codex.arg_parser import build_parser


def test_build_parser_preserves_main_cli_options():
    # The parser keeps the CLI boundary stable while main.py handles orchestration.
    args = build_parser().parse_args(
        [
            "--cn-pair",
            "3:1_2:1",
            "--float-dtype",
            "float32",
            "--route-filter",
            "max_branch_events",
            "--route-filter-threshold",
            "5",
            "--max-truncal-cn",
            "6",
            "--calculate-uncertainty",
            "--seed-filepath",
            "seeds.txt",
            "--batch-input-paths",
            "a.csv",
            "b.csv",
        ]
    )

    assert args.cn_pair == "3:1_2:1"
    assert args.float_dtype == "float32"
    assert args.route_filter == "max_branch_events"
    assert args.route_filter_threshold == 5
    assert args.max_truncal_cn == 6
    assert args.calculate_uncertainty is True
    assert args.seed_filepath == "seeds.txt"
    assert args.batch_input_paths == ["a.csv", "b.csv"]


def test_build_parser_defaults_to_old_spur_vaf_scale():
    args = build_parser().parse_args([])

    assert args.joint_vafs_scale == 0.05
    assert args.mixed_gain_loss_route_definition == "any_branch_del"
