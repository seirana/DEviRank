from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Mapping


TRACKED_PACKAGES = (
    "numpy",
    "pandas",
    "networkx",
    "scipy",
)


def _package_versions() -> dict[str, str]:
    versions: dict[str, str] = {}

    for package in TRACKED_PACKAGES:
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = "not-installed"

    return versions


def write_run_metadata(
    output_dir: str | Path,
    *,
    command: str,
    parameters: Mapping[str, Any],
) -> Path:
    """Write machine-readable provenance next to experiment outputs."""

    destination = Path(output_dir).expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)

    metadata = {
        "command": command,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "parameters": dict(parameters),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": _package_versions(),
    }

    path = destination / "run_metadata.json"
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(
        json.dumps(metadata, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    tmp.replace(path)

    return path
