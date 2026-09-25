# Experiments

This directory is reserved for generated DEviRank outputs.

Run outputs are intentionally excluded from version control so source code and generated results remain separate. Each CLI run writes a `run_metadata.json` file alongside the result CSV files.

Recommended structure:

```text
experiments/
├── quick_test/
├── results_DEviRank/
└── results_DEviRank_vs_Nbisdes/
```

For publication or benchmark archiving, store the output directory together with the Git commit SHA and exact input-data provenance.
