import csv
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

from spur.main import main


def test_main_runs_example_reads_csv_and_writes_feature_tsv():
    with TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        reads_csv = tmp_path / "example_reads.csv"
        output_dir = tmp_path / "scores"
        reads_csv.write_text(
            "alt1,ref1,alt2,ref2\n"
            "5,20,6,20\n"
            "7,20,8,20\n"
        )

        old_argv = sys.argv
        sys.argv = [
            "spur-codex",
            "--input_path",
            str(reads_csv),
            "--output_dir",
            str(output_dir),
            "--cn-pair",
            "3:1_3:1",
            "--max-truncal-cn",
            "3",
            "--max-branch-events",
            "2",
            "--n_Monte_Carlo_samples",
            "4",
            "--burn-in",
            "2",
            "--skips",
            "1",
            "--seed",
            "1",
            "--sample1_purity",
            "0.8",
            "--sample2_purity",
            "0.8",
        ]
        try:
            assert main() == 0
        finally:
            sys.argv = old_argv

        output_tsv = output_dir / "example_reads.tsv"
        assert output_tsv.exists()

        with output_tsv.open() as f:
            rows = list(csv.DictReader(f, delimiter="\t"))

        assert len(rows) == 1
        assert int(rows[0]["n_pg_routes"]) > 0
        assert int(rows[0]["n_mixed_gain_loss_routes"]) > 0
        assert "best_pg" in rows[0]
        assert "best_mixed_gain_loss" in rows[0]


def test_main_runs_multiple_seeds_in_single_program_run():
    with TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        reads_csv = tmp_path / "example_reads.csv"
        seed_file = tmp_path / "seeds.txt"
        output_dir = tmp_path / "scores"
        reads_csv.write_text(
            "alt1,ref1,alt2,ref2\n"
            "5,20,6,20\n"
            "7,20,8,20\n"
        )
        seed_file.write_text("1\n2\n")

        old_argv = sys.argv
        sys.argv = [
            "spur-codex",
            "--input_path",
            str(reads_csv),
            "--output_dir",
            str(output_dir),
            "--cn-pair",
            "3:1_3:1",
            "--max-truncal-cn",
            "3",
            "--max-branch-events",
            "2",
            "--n_Monte_Carlo_samples",
            "4",
            "--burn-in",
            "2",
            "--skips",
            "1",
            "--sample1_purity",
            "0.8",
            "--sample2_purity",
            "0.8",
            "--calculate-uncertainty",
            "--seed-filepath",
            str(seed_file),
            "--report-route-score",
        ]
        try:
            assert main() == 0
        finally:
            sys.argv = old_argv

        feature_tsv = output_dir / "example_reads.tsv"
        route_score_tsv = output_dir / "example_reads_route_scores.tsv"
        assert feature_tsv.exists()
        assert route_score_tsv.exists()

        with feature_tsv.open() as f:
            feature_rows = list(csv.DictReader(f, delimiter="\t"))
        assert [row["seed"] for row in feature_rows] == ["1", "2"]

        with route_score_tsv.open() as f:
            route_score_rows = list(csv.DictReader(f, delimiter="\t"))
        assert route_score_rows
        assert "score_seed_1" in route_score_rows[0]
        assert "score_seed_2" in route_score_rows[0]


def test_main_runs_batch_reads_file_list_and_writes_one_output_per_input():
    with TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        reads_csv_1 = tmp_path / "example_reads_1.csv"
        reads_csv_2 = tmp_path / "example_reads_2.csv"
        filepath_list = tmp_path / "inputs.txt"
        output_dir = tmp_path / "scores"
        reads_csv_1.write_text(
            "alt1,ref1,alt2,ref2\n"
            "5,20,6,20\n"
            "7,20,8,20\n"
        )
        reads_csv_2.write_text(
            "alt1,ref1,alt2,ref2\n"
            "6,20,5,20\n"
            "8,20,7,20\n"
        )
        filepath_list.write_text(f"{reads_csv_1}\n{reads_csv_2}\n")

        old_argv = sys.argv
        sys.argv = [
            "spur-codex",
            "--batch-process",
            "--batch-filepath-list",
            str(filepath_list),
            "--output_dir",
            str(output_dir),
            "--cn-pair",
            "3:1_3:1",
            "--max-truncal-cn",
            "3",
            "--max-branch-events",
            "2",
            "--n_Monte_Carlo_samples",
            "4",
            "--burn-in",
            "2",
            "--skips",
            "1",
            "--seed",
            "1",
            "--sample1_purity",
            "0.8",
            "--sample2_purity",
            "0.8",
        ]
        try:
            assert main() == 0
        finally:
            sys.argv = old_argv

        output_tsv_1 = output_dir / "example_reads_1.tsv"
        output_tsv_2 = output_dir / "example_reads_2.tsv"
        assert output_tsv_1.exists()
        assert output_tsv_2.exists()

        with output_tsv_1.open() as f:
            rows_1 = list(csv.DictReader(f, delimiter="\t"))
        with output_tsv_2.open() as f:
            rows_2 = list(csv.DictReader(f, delimiter="\t"))

        assert rows_1[0]["input_path"] == str(reads_csv_1)
        assert rows_2[0]["input_path"] == str(reads_csv_2)


def test_main_runs_batch_direct_input_paths_without_file_list():
    with TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        reads_csv_1 = tmp_path / "direct_reads_1.csv"
        reads_csv_2 = tmp_path / "direct_reads_2.csv"
        output_dir = tmp_path / "scores"
        reads_csv_1.write_text(
            "alt1,ref1,alt2,ref2\n"
            "5,20,6,20\n"
            "7,20,8,20\n"
        )
        reads_csv_2.write_text(
            "alt1,ref1,alt2,ref2\n"
            "6,20,5,20\n"
            "8,20,7,20\n"
        )

        old_argv = sys.argv
        sys.argv = [
            "spur-codex",
            "--batch-process",
            "--batch-input-paths",
            str(reads_csv_1),
            str(reads_csv_2),
            "--output_dir",
            str(output_dir),
            "--cn-pair",
            "3:1_3:1",
            "--max-truncal-cn",
            "3",
            "--max-branch-events",
            "2",
            "--n_Monte_Carlo_samples",
            "4",
            "--burn-in",
            "2",
            "--skips",
            "1",
            "--seed",
            "1",
            "--sample1_purity",
            "0.8",
            "--sample2_purity",
            "0.8",
        ]
        try:
            assert main() == 0
        finally:
            sys.argv = old_argv

        output_tsv_1 = output_dir / "direct_reads_1.tsv"
        output_tsv_2 = output_dir / "direct_reads_2.tsv"
        assert output_tsv_1.exists()
        assert output_tsv_2.exists()

        with output_tsv_1.open() as f:
            rows_1 = list(csv.DictReader(f, delimiter="\t"))
        with output_tsv_2.open() as f:
            rows_2 = list(csv.DictReader(f, delimiter="\t"))

        assert rows_1[0]["input_path"] == str(reads_csv_1)
        assert rows_2[0]["input_path"] == str(reads_csv_2)


@pytest.mark.parametrize("prefix", ["run", ""])
@pytest.mark.parametrize(
    "batch_inputs",
    [["--batch-input-paths", "a.csv", "b.csv"], ["--batch-filepath-list", "inputs.txt"]],
)
def test_main_rejects_output_name_prefix_in_batch_mode(monkeypatch, tmp_path, prefix, batch_inputs):
    output_dir = tmp_path / "scores"
    monkeypatch.setattr(sys, "argv", [
        "spur-codex", "--batch-process", *batch_inputs,
        "--output_name_prefix", prefix, "--output_dir", str(output_dir),
    ])

    with pytest.raises(ValueError, match="--output_name_prefix cannot be used with --batch-process"):
        main()

    assert not output_dir.exists()


def test_main_requires_output_dir():
    with TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        reads_csv = tmp_path / "example_reads.csv"
        reads_csv.write_text(
            "alt1,ref1,alt2,ref2\n"
            "5,20,6,20\n"
            "7,20,8,20\n"
        )

        old_argv = sys.argv
        sys.argv = [
            "spur-codex",
            "--input_path",
            str(reads_csv),
            "--cn-pair",
            "3:1_3:1",
            "--sample1_purity",
            "0.8",
            "--sample2_purity",
            "0.8",
        ]
        try:
            with pytest.raises(ValueError, match="--output_dir"):
                main()
        finally:
            sys.argv = old_argv


def test_main_requires_seed_filepath_for_uncertainty():
    with TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        reads_csv = tmp_path / "example_reads.csv"
        output_dir = tmp_path / "scores"
        reads_csv.write_text(
            "alt1,ref1,alt2,ref2\n"
            "5,20,6,20\n"
            "7,20,8,20\n"
        )

        old_argv = sys.argv
        sys.argv = [
            "spur-codex",
            "--input_path",
            str(reads_csv),
            "--output_dir",
            str(output_dir),
            "--cn-pair",
            "3:1_3:1",
            "--sample1_purity",
            "0.8",
            "--sample2_purity",
            "0.8",
            "--calculate-uncertainty",
        ]
        try:
            with pytest.raises(ValueError, match="--seed-filepath"):
                main()
        finally:
            sys.argv = old_argv
