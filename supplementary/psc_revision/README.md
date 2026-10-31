# PSC revision results

These are the complete final tables accompanying the PSC-only manuscript revision.
See [PSC_REVISION.md](../../PSC_REVISION.md) for the methods, interpretation and file map.
Run the following from the repository root to regenerate them:

```bash
python -m pip install -e ".[analysis]"
python -m scr.psc_revision --samples 100000 --seed 452456 --output experiments/psc_revision
```

Saved-screen correction and fresh sensitivity calculations are separate analyses.
At 100,000 samples the saved screen gives 10 BH and six BY discoveries at 5%.
The fresh reconstructed analysis gives 11 BH and seven BY discoveries at each
penalty. The original 36-drug list remains exploratory. The explicit reconstruction
reproduces all 36 normalization denominators and 2,839 of 2,868 saved distances;
the remaining 29 are retained in the input audit. No additional disease benchmark
was performed. Input checksums, software versions, seed and generator are in
`analysis_metadata.json`.

These files reproduce the final results delivered in the manuscript revision archive;
that archive's code overlay describes how to update the preceding source snapshot.
