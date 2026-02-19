#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# Ensure repo root (/app in Docker) is on sys.path
REPO_ROOT = Path(os.environ.get("REPO_DIR", "/app")).resolve()
sys.path.insert(0, str(REPO_ROOT))

# Now import your core module file (adjust filename if needed)
# If your core file is /app/scr/DEviRank.py, do this:
from scr.DEviRank import suggested_drugs_DEviRank  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--disease_file", required=True)
    p.add_argument("--sampling_size", required=True, type=int)
    p.add_argument("--output_folder", required=True)
    p.add_argument("--chunk_size", default=1000, type=int)
    p.add_argument("--max_drugs", default=None, type=int)
    p.add_argument("--p_value", default=0.05, type=float)
    p.add_argument("--z_score", default=-1.96, type=float)
    return p.parse_args()


def main() -> int:
    args = parse_args()

    # In your docker command you mount repo to /app
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
