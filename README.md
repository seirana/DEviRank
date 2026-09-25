# DEviRank

**Evidence-Weighted Drug Ranking via Network-Based Proximity Analysis**

DEviRank is a network-based drug-prioritization framework that ranks candidate drugs according to their network proximity to disease-associated genes in a protein-protein interaction (PPI) network. It combines bounded path-based network evidence with curated drug-gene interaction confidence scores and evaluates network proximity using degree-matched random sampling.

This repository accompanies the DEviRank method and is organized so the computational workflow can be run from the command line, in Docker, or imported as Python code.

## Method overview

The implementation follows three main ideas:

1. **Disease-gene weighting** — disease-associated genes receive weights based on local PPI connectivity.
2. **Bounded path-based target scoring** — drug targets are connected to disease genes through bounded simple paths, with PPI confidence contributing to path weights.
3. **Drug-level evidence aggregation** — target-level evidence is combined with curated drug-gene interaction confidence scores.

The code also contains a comparison workflow for the Nbisdes network-proximity baseline.

## Repository layout

```text
DEviRank/
├── data/                       repository input tables
├── experiments/                generated run outputs
├── scr/                        implementation and CLI runners
│   ├── DEviRank.py
│   ├── reproducibility.py
│   ├── run_devirank.py
│   └── run_comparison.py
├── supplementary/              paper/supplementary artifacts
├── tests/                      automated unit tests
├── .github/workflows/ci.yml    continuous integration
├── pyproject.toml              package + development configuration
├── requirements.txt
├── Dockerfile
├── LICENSE
└── README.md
```

The historical directory name `scr/` is retained for compatibility with the accompanying research materials.

## Requirements

- Python 3.10+
- NumPy
- pandas
- NetworkX

The package metadata and dependency ranges are defined in `pyproject.toml`.

## Installation

Clone the repository and create an isolated environment:

```bash
git clone https://github.com/seirana/DEviRank.git
cd DEviRank

python -m venv .venv
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -e .
```

For development and testing:

```bash
python -m pip install -e ".[dev]"
```

## Inputs

The default example inputs are stored under `data/`. The main workflow expects:

- a disease-gene CSV with an `ENSEMBL ID` column;
- a PPI table containing `gene1`, `gene2`, and the confidence values used by the scoring code;
- the drug-target matrix;
- the corresponding drug-gene interaction confidence matrix;
- a drug-name table;
- a protein-coding-gene table;
- the repeated-row mapping used to avoid recomputing duplicate target profiles.

See `data/README.md` for file-level roles and schema notes.

## Quick test

A small run is useful for checking installation and pipeline integrity before a publication-scale experiment:

```bash
devirank \
  --disease_file data/disease_target_genes.csv \
  --sampling_size 10 \
  --max_drugs 2 \
  --output_folder experiments/quick_test
```

The default random seed is `452456`. It can be changed explicitly:

```bash
devirank \
  --disease_file data/disease_target_genes.csv \
  --sampling_size 10 \
  --max_drugs 2 \
  --seed 12345 \
  --output_folder experiments/quick_test
```

Each CLI run writes `run_metadata.json` next to the results. It records the command, parameters, random seed, Python version, platform, and core package versions.

## Full DEviRank run

The default DEviRank sampling size is 100,000:

```bash
devirank \
  --disease_file data/disease_target_genes.csv \
  --output_folder experiments/results_DEviRank
```

A full run can be computationally expensive. The quick-test settings above should be used for installation checks and CI-style validation.

## Compare DEviRank with Nbisdes

```bash
devirank-compare \
  --disease_file data/disease_target_genes.csv \
  --sampling_size 1000 \
  --output_folder experiments/results_DEviRank_vs_Nbisdes
```

The comparison command uses the same reproducible seed mechanism and also writes run metadata.

## Running the original scripts directly

The script entry points remain available:

```bash
python scr/run_devirank.py \
  --disease_file data/disease_target_genes.csv \
  --sampling_size 10 \
  --max_drugs 2 \
  --output_folder experiments/quick_test

python scr/run_comparison.py \
  --disease_file data/disease_target_genes.csv \
  --sampling_size 10 \
  --max_drugs 2 \
  --output_folder experiments/comparison_quick_test
```

No hard-coded user home directory is required. Repository-local data are resolved relative to the cloned project, or the root can be supplied through `REPO_DIR` for IDE/HPC workflows.

## Docker

Build:

```bash
docker build -t devirank:latest .
```

Run a quick test while preserving outputs on the host:

```bash
mkdir -p experiments

docker run --rm \
  -v "$PWD/experiments:/app/experiments" \
  devirank:latest \
  --disease_file /app/data/disease_target_genes.csv \
  --sampling_size 10 \
  --max_drugs 2 \
  --output_folder /app/experiments/quick_test
```

Run the comparison entry point:

```bash
docker run --rm \
  --entrypoint devirank-compare \
  -v "$PWD/experiments:/app/experiments" \
  devirank:latest \
  --disease_file /app/data/disease_target_genes.csv \
  --sampling_size 10 \
  --max_drugs 2 \
  --output_folder /app/experiments/comparison_quick_test
```

## Statistical evaluation

DEviRank uses degree-matched random sampling to estimate a null distribution of network proximity scores. The implementation reports:

- observed proximity;
- mean and standard deviation of the sampled null distribution;
- z-score;
- empirical p-value;
- target and disease-gene counts;
- observed shortest-distance values.

The current empirical p-value formula is preserved from the research implementation. This modernization does **not** silently change the published/statistical method.

The sampling seed is now exposed as an explicit parameter so runs can be repeated exactly under the same software and input conditions.

## Reproducibility

Reproducibility is supported through:

- configurable and explicit random seed;
- machine-readable `run_metadata.json`;
- bounded dependency ranges;
- Docker execution;
- automated tests;
- CI across Python 3.10, 3.11, and 3.12;
- repository-relative data paths;
- generated outputs kept separate from source code.

See `REPRODUCIBILITY.md` for the recommended workflow.

## Testing and CI

Run the test suite:

```bash
python -m pytest
```

Run correctness-oriented lint checks:

```bash
python -m ruff check scr tests
```

The automated tests cover portable path resolution, graph construction validation, shortest-path edge cases, deterministic random sampling, CSV I/O, CLI defaults, and experiment metadata.

GitHub Actions runs the checks on Python 3.10, 3.11, and 3.12.

The CI suite intentionally uses small synthetic inputs. It does not attempt the full 100k-sample research experiment on every commit.

## Computational complexity

DEviRank uses bounded path enumeration around drug targets and disease genes. Runtime depends strongly on target-set size, local node degree, network density, and sampling size. Drug evaluations are independent at the workflow level, although the current reference implementation executes them serially to keep the implementation transparent.

## Research-software scope

This repository is research software. The quality upgrade focuses on portability, testing, provenance, reproducibility, and maintainability while preserving the scientific scoring logic and default seed.

The repository does not claim that a ranked drug is clinically effective. Ranking results require biological interpretation and appropriate downstream validation.

## Citation

If you use DEviRank in research, cite the associated publication when the final citation is available.

## License

MIT License. See `LICENSE`.

## Contact

For questions or reproducibility issues, open a GitHub issue.
