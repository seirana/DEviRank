#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CLI runner to compare DEviRank vs Nbisdes baseline.

Designed to work with Docker commands like:
  docker run --rm -v "$REPO_DIR:/app" devirank:latest \
    python experiments/run_comparison.py \
      --disease_file /app/data/disease_target_genes.csv \
      --sampling_size 1000 \
      --output_folder /app/experiments/results/compare_run
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _ensure_repo_on_syspath() -> Path:
    repo_root = Path(__file__).resolve().parents[1]  # <repo>/
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))
    return repo_root


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compare DEviRank vs Nbisdes baseline (proximity + comparison table).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--disease_file",
        required=True,
        help="Path to disease gene CSV containing column 'ENSEMBL ID'.",
    )
    p.add_argument(
        "--sampling_size",
        type=int,
        default=None,
        help="Monte Carlo sample size. Note: Nbisdes sampling is fixed to 1000 inside DEviRank.py.",
    )
    p.add_argument(
        "--output_folder",
        required=True,
        help="Directory where results will be written (e.g., /app/experiments/results/compare_run).",
    )
    return p.parse_args()


def main() -> int:
    _ensure_repo_on_syspath()

    from scr.DEviRank import compare_DEviRank_Nbisdes  # noqa: E402

    args = parse_args()

    disease_path = Path(args.disease_file).expanduser().resolve()
    out_dir = Path(args.output_folder).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if not disease_path.exists():
        raise FileNotFoundError(f"disease_file not found: {disease_path}")

    compare_DEviRank_Nbisdes(
        disease_file=disease_path,
        out_dir=out_dir,
        sampling_size=args.sampling_size,
    )

    print(f" Comparison finished. Results in: {out_dir}")
    print("   - drug_scores_DEviRank.csv")
    print("   - drug_scores_Nbisdes.csv")
    print("   - DEviRankVSNbisdes.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
