# Reproducibility Guide

This document describes the minimum information needed to reproduce a DEviRank run.

## 1. Record the exact code version

Use a tagged release or record the Git commit SHA used for the analysis.

## 2. Keep input data fixed

The numerical results depend on the exact contents and row ordering of the input tables under `data/` and on the disease-gene file supplied to the CLI.

Do not replace input files without recording their provenance and version.

## 3. Record parameters

At minimum record:

- disease-gene file;
- sampling size;
- random seed;
- chunk size;
- maximum number of drugs, if a subset was used;
- p-value threshold;
- z-score threshold;
- output directory.

The CLI writes these settings to `run_metadata.json`.

## 4. Record the software environment

`run_metadata.json` also records:

- Python version;
- operating-system/platform information;
- NumPy version;
- pandas version;
- NetworkX version.

For stronger isolation, use the Docker image built from the same repository commit.

## 5. Use an explicit seed

The default is:

```text
452456
```

For a paper, benchmark, or comparison, pass the seed explicitly rather than relying only on the default:

```bash
devirank \
  --disease_file data/disease_target_genes.csv \
  --sampling_size 100000 \
  --seed 452456 \
  --output_folder experiments/results_DEviRank
```

## 6. Separate quick tests from research-scale runs

A quick test verifies that installation and data plumbing work. It is not a replacement for the sampling regime used in a scientific analysis.

Example:

```bash
devirank \
  --disease_file data/disease_target_genes.csv \
  --sampling_size 10 \
  --max_drugs 2 \
  --seed 452456 \
  --output_folder experiments/quick_test
```

## 7. Preserve raw outputs

Generated files under `experiments/` are intentionally ignored by Git. For a publication or benchmark, archive the exact output directory together with:

- `run_metadata.json`;
- the commit SHA;
- the input-data version/provenance;
- any downstream analysis scripts used on those outputs.

## Statistical-method note

The modernization work intentionally preserves the existing proximity, z-score, and empirical p-value calculations. Methodological changes should be introduced separately, documented explicitly, and compared against the reference implementation rather than being mixed into software-maintenance changes.
