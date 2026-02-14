#!/usr/bin/env bash
set -euo pipefail

IMAGE_NAME="${IMAGE_NAME:-devirank:latest}"
PY_SCRIPT="${PY_SCRIPT:-scripts/devirank_pipeline.py}"  # adjust to your file path in repo

DISEASE_FILE="${1:?Usage: $0 <disease_file.csv> [sampling] [p_value] [z_score]}"
SAMPLING="${2:-100000}"
P_VALUE="${3:-0.05}"
Z_SCORE="${4:--1.96}"

: "${DATA_DIR:?export DATA_DIR=/absolute/path/to/data_root}"
: "${OUT_DIR:?export OUT_DIR=/absolute/path/to/output_root}"

# Build image if not present (or always build if you prefer)
docker build -t "$IMAGE_NAME" .

docker run --rm \
  -e REPO_DIR=/app \
  -e DATA_DIR=/data \
  -e OUT_DIR=/out \
  -e PY_SCRIPT="/app/$PY_SCRIPT" \
  -e DISEASE_FILE="/data/$(basename "$DISEASE_FILE")" \
  -e SAMPLING="$SAMPLING" \
  -e P_VALUE="$P_VALUE" \
  -e Z_SCORE="$Z_SCORE" \
  -v "$DATA_DIR":/data:ro \
  -v "$OUT_DIR":/out \
  "$IMAGE_NAME" \
  python3 - <<'PY'
import os
from pathlib import Path
import importlib.util

py_script = os.environ["PY_SCRIPT"]
disease_file = os.environ["DISEASE_FILE"]
sampling = int(os.environ["SAMPLING"])
p_value = float(os.environ["P_VALUE"])
z_score = float(os.environ["Z_SCORE"])

spec = importlib.util.spec_from_file_location("devirank_pipeline", py_script)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

mod.suggested_drugs_DEviRank(
    disease_file=Path(disease_file),
    sampling_size=sampling,
    p_value=p_value,
    z_score=z_score,
)
PY

