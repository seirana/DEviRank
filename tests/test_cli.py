import sys

from scr.DEviRank import DEFAULT_SEED
from scr import run_comparison, run_devirank


def test_devirank_cli_default_seed(monkeypatch):
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_devirank.py",
            "--disease_file",
            "disease.csv",
            "--output_folder",
            "out",
        ],
    )

    args = run_devirank.parse_args()

    assert args.seed == DEFAULT_SEED
    assert args.sampling_size == 100000


def test_comparison_cli_default_seed(monkeypatch):
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_comparison.py",
            "--disease_file",
            "disease.csv",
            "--output_folder",
            "out",
        ],
    )

    args = run_comparison.parse_args()

    assert args.seed == DEFAULT_SEED
    assert args.sampling_size == 100000
