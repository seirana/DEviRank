#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CLI runner for the full DEviRank pipeline (proximity + drug scoring + disease-gene scoring).

Designed to work with Docker commands like:
  docker run --rm -v "$REPO_DIR:/app" devirank:latest \
    python scr/run_devirank.py \
      --disease_file /app/data/disease_target_genes.csv \
      --sampling_size 1000 \
      --output_folder /app/experiments/results_quick_test
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]      # /app
SRC_DIR = REPO_ROOT / "src"                          # /app/src

sys.path.insert(0, str(REPO_ROOT))                   # allow imports from repo root
sys.path.insert(0, str(SRC_DIR))                     # allow importing DEviRank.py from src/

def _ensure_repo_on_syspath() -> Path:
    """
    Ensure the repo root is on sys.path so 'src.DEviRank' imports work
    when running as: python experiments/run_devirank.py
    """
    repo_root = Path(__file__).resolve().parents[1]  # <repo>/
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
    return repo_root


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Run DEviRank (full pipeline: proximity + scoring).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--disease_file",
        required=True,
        help="Path to disease gene CSV containing column 'ENSEMBL ID'. Can be anywhere (must be accessible inside container).",
    )
    p.add_argument(
        "--sampling_size",
        type=int,
        default=None,
        help="Monte Carlo sample size for DEviRank proximity. If omitted, DEviRank defaults to 100000.",
    )
    p.add_argument(
        "--output_folder",
        required=True,
        help="Directory where results will be written (e.g., /app/experiments/results/run1).",
    )
    p.add_argument(
        "--p_value",
        type=float,
        default=0.05,
        help="P-value threshold used during scoring.",
    )
    p.add_argument(
        "--z_score",
        type=float,
        default=-1.96,
        help="Z-score threshold used during scoring.",
    )
    return p.parse_args()


def main() -> int:
    _ensure_repo_on_syspath()

    # Import after sys.path fix
    from src.DEviRank import suggested_drugs_DEviRank  # noqa: E402

    args = parse_args()

    disease_path = Path(args.disease_file).expanduser().resolve()
    out_dir = Path(args.output_folder).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if not disease_path.exists():
        raise FileNotFoundError(f"disease_file not found: {disease_path}")

    suggested_drugs_DEviRank(
        disease_file=disease_path,
        out_dir=out_dir,
        sampling_size=args.sampling_size,
        p_value=args.p_value,
        z_score=args.z_score,
    )

    print(f" DEviRank finished. Results in: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
