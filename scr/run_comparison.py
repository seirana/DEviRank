#!/usr/bin/env python3
import os
from pathlib import Path

from src.devirank_core import compare_DEviRank_Nbisdes

if __name__ == "__main__":
    repo_dir = os.environ.get("REPO_DIR", None)
    if repo_dir is None:
        repo_dir = str(Path(__file__).resolve().parents[1])
        os.environ["REPO_DIR"] = repo_dir

    disease_file = os.environ.get("DISEASE_FILE", f"{repo_dir}/data/disease_target_genes.csv")
    out_dir      = os.environ.get("OUT_DIR",      f"{repo_dir}/experiments/results")

    sampling_size = int(os.environ.get("SAMPLING_SIZE", "100000"))
    chunk_size    = int(os.environ.get("CHUNK_SIZE", "1000"))
    max_drugs_env = os.environ.get("MAX_DRUGS", "")
    max_drugs     = int(max_drugs_env) if max_drugs_env else None

    compare_DEviRank_Nbisdes(
        disease_file=disease_file,
        out_dir=out_dir,
        sampling_size=sampling_size,
        chunk_size=chunk_size,
        max_drugs=max_drugs,
    )
