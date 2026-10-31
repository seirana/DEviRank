# PSC manuscript reanalysis

This code supports **DEviRank: Evidence-Weighted Network-Based Drug Ranking with a Case Study in Primary Sclerosing Cholangitis**. PSC is the only disease analysed.

## Run

From this repository checkout, run:

```bash
python -m pip install -e ".[analysis,dev]"
python -m pytest
python -m ruff check scr tests
python -m scr.psc_revision --samples 100000 --seed 452456 --output experiments/psc_revision
```

`--max-profiles` is a pilot-only option. Never use partial-family FDR results in a manuscript.

The final complete tables are also available in [`supplementary/psc_revision/`](supplementary/psc_revision/). The reanalysis input snapshot is based on source commit `738c2af027605f4d8c2ea1c12106405429037aaf`; input checksums and the reconstructed gene mapping are recorded explicitly. The deposited results are the same final analysis delivered with the revised manuscript.

## Statistical correction

The saved p-values were `b/B`. Their tail counts reconstruct as integers at the recorded sample size. The correction is `(b+1)/(B+1)`, including tied distances. Six zero p-values at 100,000 draws become `1/100001`.

There are 2,868 evaluable saved tests among 2,871 rows. Atovaquone, glucarpidase and kanamycin A have 20,000 sentinels; the saved output does not explain why. Those p-values are unavailable, not valid numerical tests. BH and BY are applied before filtering. At 100,000 samples, 10 pass BH and six pass BY at 5%; seven and five belong to the original 36-drug list. Including the unavailable rows as non-rejections in a 2,871-test family changes adjusted values slightly but not those counts.

The 36-drug ranking remains exploratory. A further z-score filter does not automatically preserve FDR control for a subset. BH requires independent or suitably positively dependent valid tests; BY permits arbitrary dependence. Neither adjustment validates the biological null model or predicts clinical effectiveness. Tazemetostat's saved BH-adjusted p-value is about 0.488, and BY-adjusted p is 1.

The sampling-convergence table is recomputed from the workbook because its previous manuscript values did not reproduce exactly. The corrected z-score SDD is 0.204 for 1,000 versus 100,000 draws and 0.011 for 10,000 versus 100,000 draws. Small-tail FDR decisions still differ across sampling sizes.

## Sensitivity and input provenance

The unrevised default disease file has 71 identifiers, whereas Supplementary Table 1 lists 54 symbols. The revision provides `data/psc_manuscript_targets.csv` and a complete gene mapping and provenance record. All 54 supplied symbols map uniquely on primary-assembly chromosome 6 using Ensembl Genes 116 annotations distributed in annotables. The protein universe is the supplied protein-coding list union the PPI identifiers: 21,857 nodes and 231,641 edges. Eligible isolated proteins remain; out-of-universe targets are excluded. This reconstruction matches 2,839 of 2,868 evaluable saved distances. The remaining 29 differences prevent a claim of exact historical replication.

Fresh calculations use 100,000 sets per distinct eligible profile, seed 452456 and NumPy PCG64. Each target draws from its degree bin excluding itself; collisions are resampled until unique target-set size is preserved. Atovaquone, glucarpidase and kanamycin A have no target in the declared protein universe and are unavailable tests in this reconstruction. The duplicate isatuximab target is collapsed for proximity. These draws are explicitly different from the saved Python-random stream.

Common random sets are reused across penalties 100, 1000 and 20000. All exceed the maximum finite closest distance of 10. Observed and random distances are recalculated together, using integer sums to preserve ties. All 2,868 evaluable fresh tests are adjusted together. At every penalty, the nominal rule selects the same 35 drugs, BH the same 11 and BY the same seven. These fresh counts do not replace the corrected saved counts of 36 nominal, 10 BH and six BY.

Path length counts edges. The deposited baseline is two edges; one, three and four are also tested. The 36 saved names and each drug's eligible-target denominator remain fixed across limits. The denominator includes eligible targets with zero path score and reproduces all 36 saved normalization denominators exactly. Ties receive average ranks; top-ten ties are resolved alphabetically. Correlations with two edges are 0.826, 0.932 and 0.910 for one, three and four edges; top-ten overlap is six, nine and nine. The reconstructed two-edge order agrees with the saved order at 0.99994, with the residual difference attributable to tie handling. Raw-score magnitudes and 29 observed distances are not claimed to reproduce exactly.

Path scoring retains first-match undirected PPI confidence and reference first-endpoint disease-weight aggregation. Self loops remain in degree matching but are excluded from simple paths. Shared targets contribute a zero-edge path. Sparse non-backtracking walk sums are corrected for cycles and are tested against explicit weighted simple-path enumeration through four edges.

## Outputs

| File | Contents |
|---|---|
| `saved_PSC_1000_FDR.csv`, `saved_PSC_10000_FDR.csv`, `saved_PSC_100000_FDR.csv` | Complete saved families with tail counts, original and corrected p-values and BH/BY adjustment |
| `saved_36_candidates_FDR.csv` | Saved ranks and scores with FDR annotations |
| `saved_FDR_summary.csv` | Counts at each sampling size |
| `saved_convergence_recomputed.csv` | Revised Table 1 values |
| `saved_denominator_audit.csv` | Independent verification of all 36 saved normalization denominators |
| `input_reproduction_audit.csv` | Observed-distance comparison for every drug |
| `repository_penalty_results.csv` | Every fresh drug test at every penalty |
| `penalty_sensitivity_summary.csv` | Counts and Jaccard overlaps |
| `repository_path_scores_all_drugs.csv` | Raw and normalized scores at each edge limit |
| `repository_path_scores_fixed_36.csv` | Fixed-list path sensitivity scores and ranks |
| `path_sensitivity_summary.csv` | Correlations, overlaps and rank changes |
| `analysis_metadata.json` | Seed, generator, software versions, source commit and input checksums |

Exact reproduction of the historical screen still requires reconciliation of its original identifiers, PPI snapshot or target data. The 29 remaining distance discrepancies are named in `input_reproduction_audit.csv`. The revised manuscript states this limitation and includes a localized Word comment. No additional disease benchmark is introduced. The gene reconstruction was defined from the supplied symbols and independent annotations, without adding genes to improve agreement.
