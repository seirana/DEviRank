#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CLI runner for DEviRank (Docker-friendly).

Example (as you run it):
sudo docker run --rm -v "$REPO_DIR:/app" -w /app devirank:latest \
  python /app/scr/run_devirank.py \
    --disease_file /app/data/disease_target_genes.csv \
    --sampling_size 1000 \
    --output_folder /app/experiments/results_quick_test
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from scr.DEviRank import suggested_drugs_DEviRank  # adjust if your core file name differs


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run DEviRank pipeline.")
    p.add_argument("--disease_file", required=True, type=str, help="Path to disease_target_genes.csv")
    p.add_argument("--sampling_size", required=True, type=int, help="Random sampling size")
    p.add_argument("--output_folder", required=True, type=str, help="Output directory for results")
    p.add_argument("--chunk_size", default=1000, type=int, help="Chunk size for proximity collecting (default: 1000)")
    p.add_argument("--max_drugs", default=None, type=int, help="Run only first N drugs (quick test)")
    p.add_argument("--p_value", default=0.05, type=float, help="P-value threshold for scoring (default: 0.05)")
    p.add_argument("--z_score", default=-1.96, type=float, help="Z-score threshold for scoring (default: -1.96)")
    return p.parse_args()


def main() -> int:
    args = _parse_args()

    # In Docker we mount repo at /app, so set REPO_DIR accordingly unless already set
    os.environ.setdefault("REPO_DIR", "/app")

    suggested_drugs_DEviRank(
        disease_file=args.disease_file,
        out_dir=args.output_folder,
        sampling_size=args.sampling_size,
        chunk_size=args.chunk_size,
        max_drugs=args.max_drugs,
        p_value=args.p_value,
        z_score=args.z_score,
    )

    print(f"Done. Results in: {Path(args.output_folder).resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
