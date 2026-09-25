# DEviRank input data

The files in this directory are the repository-local inputs used by the reference implementation.

| File | Role in the code |
|---|---|
| `disease_target_genes.csv` | Default disease-associated genes. Expected column: `ENSEMBL ID`. |
| `drugs(filtered).csv` | Drug names aligned by row with the drug-target matrices. Expected column: `drug_name`. |
| `DtoGI_ENSEMBL(filtered).csv` | Wide drug-to-gene target matrix. Non-zero entries represent gene identifiers used as targets. |
| `DtoGI_scores(filtered).csv` | Wide matrix of drug-gene interaction confidence values aligned with `DtoGI_ENSEMBL(filtered).csv`. |
| `gene_gene_PPI700_ENSEMBL.csv` | PPI network used for graph construction and PPI weighting. The implementation expects `gene1`, `gene2`, and uses `max_ppi` when scoring paths. |
| `protein_coding_genes_ENSEMBL.csv` | Protein-coding genes. Expected column: `Gene stable ID`. |
| `repeated(filtered).csv` | Row-index mapping used to reuse results for repeated drug-target profiles. |
| `drugs_links.csv` | Auxiliary drug-link table included with the research data. |
| `proteins.csv` | Auxiliary protein table included with the research data. |

## Alignment requirement

The drug rows in the drug-name, target, interaction-score, and repeated-index tables are positionally aligned by the current implementation. Reordering one table independently can invalidate results.

## Provenance

The code repository does not contain enough machine-readable metadata to reconstruct the external source/version of every included biological table. Do not infer provenance from filenames alone.

For a new publication-scale analysis, record the source database, release/version, retrieval date, filtering rules, identifier mapping procedure, and checksum for each regenerated input table.
