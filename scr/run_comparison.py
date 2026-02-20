#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: seirana
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

REPO_ROOT = Path(os.environ.get("REPO_DIR", Path(__file__).resolve().parents[1])).resolve()
sys.path.insert(0, str(REPO_ROOT))

from scr.DEviRank import compare_DEviRank_Nbisdes  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description="Run DEviRank vs Nbisdes comparison.")
    p.add_argument("--disease_file", required=True, help="Path to disease gene CSV")
    p.add_argument("--sampling_size",  default=100000, type=int, help="Sampling size for DEviRank (Nbisdes uses its own setting)")
    p.add_argument("--output_folder", required=True, help="Output directory")
    p.add_argument("--chunk_size", default=10, type=int, help="Chunk size for proximity computation")
    p.add_argument("--max_drugs", default=None, type=int, help="Quick test: limit to first N drugs")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    # If you set REPO_DIR, respect it. Otherwise infer from this file.
    os.environ.setdefault("REPO_DIR", str(REPO_ROOT))

    compare_DEviRank_Nbisdes(
        disease_file=args.disease_file,
        out_dir=args.output_folder,
        sampling_size=args.sampling_size,
        chunk_size=args.chunk_size,
        max_drugs=args.max_drugs,
    )

    print(f"Done. Results in: {Path(args.output_folder).resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
