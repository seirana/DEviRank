# DEviRank

**Evidence-Weighted Drug Ranking via Network-Based Proximity Analysis**

DEviRank is a network-based drug efficacy screening framework that ranks candidate drugs based on their proximity to disease-associated genes in a protein–protein interaction (PPI) network. It extends classical interactome-based proximity methods by integrating **path-level network evidence** with **curated drug–gene interaction confidence**, enabling interpretable and statistically robust drug prioritization.

This repository accompanies the DEviRank method described in our LNCS publication and provides a fully reproducible implementation, supplementary analyses, and documentation.
---

---

## Table of Contents
- [DEviRank](#devirank)
- [Method Overview](#method-overview)
- [Key Features](#-key-features)
- [Quick Start (Drug Ranking)](#quick-start-drug-ranking)
  - [0) Clone the repository](#0-clone-the-repository)
  - [1) Fast sanity check (minutes)](#1-fast-sanity-check-minutes)
  - [2) Full drug ranking (hours-to-days)](#2-full-drug-ranking-hours-to-days)
  - [3) Compare against Guney-style proximity baseline](#3-compare-against-guney-style-proximity-baseline)
- [Repository Structure](#repository-structure)
- [Installation](#️-installation)
- [Usage](#usage)
- [Run with Docker](#run-with-docker)
- [Statistical Evaluation](#statistical-evaluation)
- [Computational Complexity](#️-computational-complexity)
- [Reproducibility](#reproducibility)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## 🔬 Method Overview

DEviRank builds on the observation that effective drugs typically interact with **a subset of disease-associated proteins located within a local network neighborhood**, rather than targeting the entire disease module.

The method consists of three main steps:

1. **Disease gene weighting**
   Disease-associated genes are weighted by their local PPI connectivity to reflect their network influence.

2. **Path-based drug target scoring**
   For each drug target gene, DEviRank aggregates evidence from **all simple paths of bounded length (≤ 3)** connecting the drug target to disease-associated genes.
   Path contributions are weighted by:

   * PPI interaction confidence
   * disease-gene importance

3. **Drug-level aggregation**
   Drug target scores are combined using **curated drug–gene interaction (DGI) confidence**, yielding a final, interpretable drug score.

Unlike end-to-end learning approaches, DEviRank emphasizes **interpretability, biological grounding, and statistical stability**, making it suitable for settings with limited labeled data.

---

## ✨ Key Features

* Network-based drug ranking with **explicit scoring formulation**
* Integration of **PPI confidence** and **curated DGI evidence**
* Bounded **simple-path enumeration** (max length = 3)
* Robust statistical evaluation via **degree-preserving random sampling**
* Designed for **parallel execution** across drugs
* Fully reproducible and transparent implementation

---

## 📂 Repository Structure

```
DEviRank/
│
├── src/
│   ├── DEviRank.py  
│
├── data/
│   ├── disease_target_genes.csv/
│   ├── drugs(filtered).csv/ 
│   ├── drugs_links.csv/ 
│   ├── DtoDGI_ENSEMBL(filtered).csve/ 
│   ├── DtoGI_scores(filtered).csv/
│   ├── gene_gene_PPI700_ENSEMBL.csv/
│   ├── protein_coding_genes_ENSEMBL.csv        
│   └── repeated(filtered).csv /      
│
├── experiments/
│   ├── run_devirank.py     # Main experiment runner
│   ├── sampling_analysis/  # Random sampling robustness experiments
│   └── results/            # Output rankings and statistics
│
├── supplementary/
│   ├── time_complexity.pdf # Step-by-step complexity derivation
│   └── additional_tables/  # Supplementary tables and figures
│
├── requirements.txt
├── LICENSE
├── Dockerfile
└── README.md
```

---

## Quick Start (Drug Ranking)

DEviRank can take **hours to days** depending on network size, number of drugs, and the random sampling size.
So: **run a fast sanity check first**, then scale up.

### 1. Clone the repository
```bash
cd ~
git clone https://github.com/seirana/DEviRank.git
cd DEviRank
```

---

### ⚙️ 2. Installation

Requirements

* Python ≥ 3.9
* NumPy
* NetworkX
* pandas

DEviRank supports two installation methods.

Option A — Native Python

```bash
pip install -r requirements.txt
```

Requires Python ≥ 3.9.

Option B — Docker (Recommended for Reproducibility)

```bash
docker build -t devirank:latest .
```

---

### 🚀 3. Usage

3.1. Prepare Inputs

You will need the following input data:

* A set of disease-associated genes (replace it with your desired genes)

* A weighted protein–protein interaction (PPI) network,

* Curated drug–gene interaction data with confidence scores,

* A gene–protein mapping table (retrieved from Ensembl BioMart), and

* All required input files are available in the ./data/ directory.

---

🐳 3.2 Build the Docker Image

All experiments are executed inside Docker to ensure reproducibility and consistent environments.

```bash
docker build -t devirank:latest .
```

3.3 Run DEviRank

Option A — Quick Test (Sanity Check, Minutes)

Runs a small random sampling to verify installation and pipeline integrity.

```bash
docker run --rm \
  -v "$(pwd)/data:/app/data" \
  -v "$(pwd)/experiments/results:/app/experiments/results" \
  devirank:latest \
  python experiments/run_devirank.py \
    --ppi data/ppi/ppi_network.tsv \
    --disease data/disease_genes/disease_genes.txt \
    --drug_gene data/drug_gene/drug_gene.tsv \
    --max_path_length 2 \
    --n_random 1000
```

Option B — Full Drug Ranking (Hours to Days)

High-precision Monte Carlo estimation.

```bash
docker run --rm \
  -v "$(pwd)/data:/app/data" \
  -v "$(pwd)/experiments/results:/app/experiments/results" \
  devirank:latest \
  python experiments/run_devirank.py \
    --ppi data/ppi/ppi_network.tsv \
    --disease data/disease_genes/disease_genes.txt \
    --drug_gene data/drug_gene/drug_gene.tsv \
    --max_path_length 3 \
    --n_random 100000
 ```
   
Output:

* Ranked list of drugs
* z-scores and p-values from random sampling
* Intermediate statistics

Option C — Comparison with Network Proximity–Based Baseline

To compare DEviRank against the shortest-path proximity baseline:

```bash
docker run --rm \
  -v "$(pwd)/data:/app/data" \
  -v "$(pwd)/experiments/results:/app/experiments/results" \
  devirank:latest \
  python experiments/compare_devirank_vs_baseline.py \
    --ppi data/ppi/ppi_network.tsv \
    --disease data/disease_genes/disease_genes.txt \
    --drug_gene data/drug_gene/drug_gene.tsv \
    --n_random 10000  
```

Output:

* z-scores and p-values from random sampling
* Intermediate statistics

Option D — If Docker is not available:

```bash
pip install -r requirements.txt
python experiments/run_devirank.py ...
```

---

## 📊 Statistical Evaluation

DEviRank evaluates drug–disease proximity using **degree-preserving random sampling**, following established interactome-based methods.

Random sampling is used to:

* estimate null distributions
* compute z-scores and p-values
* assess robustness and convergence

### Sampling sizes

We provide results for:

* 1k samples (baseline, following prior work)
* 10k samples (stable and efficient)
* 100k samples (high-precision reference)

A detailed robustness and convergence analysis is included in the **supplementary material**.

---

## ⏱️ Computational Complexity

DEviRank enumerates **simple paths of bounded length (L ≤ 3)**.

Per drug, the runtime is bounded by:

[
O!\left(\sum_{t \in T} \deg(t) + \sum_{s \in S} \deg(s)\Delta^2\right)
]

and by (O(|S|\Delta^3)) in the worst case, where:

* (S) = drug target genes
* (T) = disease genes
* (\Delta) = maximum PPI degree

Since drugs are evaluated independently, DEviRank is **embarrassingly parallel**.

A step-by-step complexity derivation is provided in:

```
supplementary/time_complexity.pdf
```

---

## 🔁 Reproducibility

All experiments reported in the paper can be reproduced using this repository.

To ensure reproducibility:

* random seeds are fixed
* sampling sizes are configurable
* intermediate results are logged

---

## 📄 Citation

If you use DEviRank in your research, please cite:

```
[Will be added after publication]
```

---

## 📜 License

This project is released under the **MIT License**.
See the `LICENSE` file for details.

---

## 🤝 Contact

For questions, issues, or collaboration requests, please open a GitHub issue or contact:

**Seirana Hashemi**
GitHub: [https://github.com/seirana](https://github.com/seirana)

---

## 🧠 Notes for Reviewers

* This repository serves as **supplementary material** for the associated LNCS paper.
* All design choices are explicitly documented.
* The implementation favors **clarity and interpretability** over black-box optimization.

---
