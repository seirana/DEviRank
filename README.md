# DEviRank

**Evidence-Weighted Drug Ranking via Network-Based Proximity Analysis**

DEviRank is a network-based drug prioritization framework that ranks candidate drugs according to their network proximity to disease-associated genes in a protein–protein interaction (PPI) network. Building on classical interactome-based proximity measures, DEviRank aggregates evidence from bounded simple paths between drug targets and disease genes and integrates this with curated drug–gene interaction confidence scores. The resulting formulation provides an interpretable, evidence-weighted ranking while enabling statistical assessment through degree-preserving random sampling.

This repository accompanies the DEviRank method described in my LNCS publication and provides a fully reproducible implementation, supplementary analyses, and documentation.

---

## Table of Contents

- [DEviRank](#devirank)
- [🔬 Method Overview](#-method-overview)
- [✨ Key Features](#-key-features)
- [📂 Repository Structure](#-repository-structure)
- [🚀 Quick Start](#-quick-start)
  - [1. ⬇️ Download System Requirements](#1-⬇️-download-system-requirements)
  - [2. 📥 Clone the Repository](#2-📥-clone-the-repository)
  - [3. ⚙️ Installation](#3-⚙️-installation)
  - [4. 🧬 Usage](#4-🧬-usage)
    - [4.1. 📑 Prepare Inputs](#41-📑-prepare-inputs)
    - [4.2. ▶️ Run DEviRank](#42-▶️-run-devirank)
- [📊 Statistical Evaluation](#-statistical-evaluation)
- [⏱️ Computational Complexity](#️-computational-complexity)
- [🔁 Reproducibility](#-reproducibility)
- [📄 Citation](#-citation)
- [📜 License](#-license)
- [🤝 Contact](#-contact)
- [🧠 Notes for Reviewers](#-notes-for-reviewers)

---

## 🔬 Method Overview

DEviRank builds on the observation that therapeutically relevant drugs often modulate a subset of disease-associated proteins located within a local network neighborhood, rather than targeting the entire disease module.

The method consists of three main components:

1. **Disease Gene Weighting**  
   Each disease-associated gene is assigned a weight based on its local connectivity within the PPI network, reflecting its relative network influence.

2. **Bounded Path-Based Target Scoring**  
   For each drug target gene, DEviRank aggregates evidence from all simple paths of bounded length (≤ 3) connecting the drug target to disease-associated genes.  
   Each path contribution is weighted by:
   - PPI interaction confidence (edge weights)
   - the precomputed importance of the corresponding disease gene

3. **Drug-Level Evidence Aggregation**  
   Drug-level scores are obtained by combining target-level scores using curated drug–gene interaction (DGI) confidence, producing an interpretable, evidence-weighted ranking.

Unlike end-to-end learning approaches, DEviRank is explicitly model-driven and emphasizes interpretability, biological transparency, and statistical validation through degree-preserving random sampling. This makes the framework particularly suitable for settings with limited labeled data or where methodological transparency is required.

---

## ✨ Key Features

- Explicit evidence-weighted scoring formulation combining network topology and interaction confidence  
- Integration of weighted PPI edges and curated drug–gene interaction (DGI) confidence  
- Aggregation over bounded simple paths (maximum length = 3) to capture local network effects  
- Degree-preserving random sampling for statistically grounded proximity assessment  
- Embarrassingly parallel per-drug evaluation for scalable computation  
- Fully reproducible and transparent research implementation

---

## 📂 Repository Structure

```
DEviRank/
│
├── data/
│   ├── disease_target_genes.csv
│   ├── drugs(filtered).csv
│   ├── drugs_links.csv
│   ├── DtoGI_ENSEMBL(filtered).csv
│   ├── DtoGI_scores(filtered).csv
│   ├── gene_gene_PPI700_ENSEMBL.csv
│   ├── protein_coding_genes_ENSEMBL.csv
│   ├── proteins.csv      
│   └── repeated(filtered).csv
│
├── docs/
│   ├── DEviRank Overview.png  
│   └── DeviRank Psudocode.txt 
│
├── experiments/
│   ├── results_devirank/  # DEviRank output
│   ├── results_quick_test/ # Quick-test output (DEviRank)
│   └── results_DEviRank_vs_Nbisdes/ # Comparison output: DEviRank vs Nbisdes
│
├── scr/
│   ├── DEviRank.py
│   ├── __init__.py
│   ├── run_devirank.py      
│   └── run_comparision.py
│
├── supplementary/
│   ├── DeviRank_TimeComplexityProof.pdf
│   └── Supplementary tables.xlsx
│
├── requirements.txt
├── LICENSE
├── Dockerfile
└── README.md
```

---

## 🚀 Quick Start

Requirements

* Python ≥ 3.9
* numpy ≥ 1.24
* pandas ≥ 2.0
* networkx ≥ 3.1

### 1. ⬇️ Download system requirements 

```bash
sudo apt update
sudo apt install -y \
    git \
    build-essential \
    wget curl unzip \
    python3.12 python3.12-venv
```
    
### 2. 📥 Clone the repository
```bash
cd ~
git clone https://github.com/seirana/DEviRank.git
cd DEviRank
```

### ⚙️ 3. Installation

DEviRank supports two installation methods.

*Option A — Native setup (for HPC / no-Docker environments)*
```bash
cd ~
REPO_DIR="$(find . -maxdepth 5 -type f -name  requirements.txt -path '*/DEviRank/*' -print -quit | xargs -r dirname)"
echo "Using repo: $REPO_DIR"
cd "$REPO_DIR"

conda create -n devirank python=3.12
conda activate devirank
python -m pip install --upgrade pip
pip install -r requirements.txt
```

*Option B — 🐳 Docker (Recommended for Reproducibility)*

```bash
cd ~
REPO_DIR="$(find "$HOME" -maxdepth 5 -type f -name Dockerfile -path '*/DEviRank/*' -print -quit | xargs -r dirname)"
echo "Using repo: $REPO_DIR"
cd "$REPO_DIR"

sudo docker build -t devirank:latest .
```

### 4. 🧬 Usage

#### 4.1. 📑 Prepare Inputs

You will need the following input data:

* A set of disease-associated genes (replace it with your desired genes)

* A weighted protein–protein interaction (PPI) network,

* Curated drug–gene interaction data with confidence scores,

* A gene–protein mapping table (retrieved from Ensembl BioMart), and

* All required input files are available in the ./data/ directory.
  

All experiments are executed inside Docker to ensure reproducibility and consistent environments.

#### 4.2 ▶️ Run DEviRank

#### Option A — Quick Test (Sanity Check, Minutes)

Runs a small random sampling to verify installation and pipeline integrity.

```bash
cd ~
REPO_DIR="$(find "$HOME" -maxdepth 5 -type f -name Dockerfile -path '*/DEviRank/*' -print -quit | xargs -r dirname)"
echo "Using repository at: $REPO_DIR"

sudo docker run --rm -v "$REPO_DIR:/app" -w /app devirank:latest \
  python /app/scr/run_devirank.py \
    --disease_file /app/data/disease_target_genes.csv \
    --sampling_size 1000 \
    --output_folder /app/experiments/results_quick_test
```

without Docker:

```bash
cd ~
REPO_DIR="$(find "$HOME" -maxdepth 5 -type f -name Dockerfile -path '*/DEviRank/*' -print -quit | xargs -r dirname)"
echo "Using repository at: $REPO_DIR"

python scr/run_devirank.py\
  --disease_file data/disease_target_genes.csv \
  --sampling_size 1000 \
  --output_folder experiments/results_quick_test
```

#### Option B — Full Drug Ranking (Hours to Days)

High-precision Monte Carlo estimation.

```bash
cd ~
REPO_DIR="$(find "$HOME" -maxdepth 5 -type f -name Dockerfile -path '*/DEviRank/*' -print -quit | xargs -r dirname)"
echo "Using repository at: $REPO_DIR"

sudo docker run --rm \
  -v "$REPO_DIR:/app" \
  devirank:latest \
  python scr/run_devirank.py \
    --disease_file /app/data/disease_target_genes.csv \  
    --output_folder /app/experiments/result_DEviRank \ 
 ```
 without Docker:

```bash
cd ~
REPO_DIR="$(find "$HOME" -maxdepth 5 -type f -name Dockerfile -path '*/DEviRank/*' -print -quit | xargs -r dirname)"
echo "Using repository at: $REPO_DIR"

python scr/run_devirank.py \
  --disease_file /app/data/disease_target_genes.csv \  
  --output_folder /app/experiments/result_DEviRank \ 
```
  
Output:

* Ranked list of drugs
* z-scores and p-values from random sampling
* Intermediate statistics

#### Option C — Comparison with Network Proximity–Based Baseline

To compare DEviRank against the shortest-path proximity baseline:

```bash
cd ~
REPO_DIR="$(find "$HOME" -maxdepth 5 -type f -name Dockerfile -path '*/DEviRank/*' -print -quit | xargs -r dirname)"
echo "Using repository at: $REPO_DIR"

sudo docker run --rm \
  -v "$REPO_DIR:/app" \
  devirank:latest \
  python scr/run_comparison.py\
    --disease_file /app/data/disease_target_genes.csv \  
    --output_folder /app/experiments/results_DEviRank_vs_Nbisdes
```

 without Docker:

```bash
cd ~
REPO_DIR="$(find "$HOME" -maxdepth 5 -type f -name Dockerfile -path '*/DEviRank/*' -print -quit | xargs -r dirname)"
echo "Using repository at: $REPO_DIR"

python scr/run_comparison.py \
  --disease_file /app/data/disease_target_genes.csv \  
  --output_folder /app/experiments/results_DEviRank_vs_Nbisdes
```

Output:

* z-scores and p-values from random sampling
* Intermediate statistics

---

## 📊 Statistical Evaluation

DEviRank assesses drug–disease proximity using degree-preserving random sampling in the PPI network. For each drug, random protein sets are generated that match the size and degree distribution of the original drug target set, ensuring a biologically meaningful null model.

Random sampling is used to:
- estimate null distributions of network proximity scores,
- compute z-scores and empirical p-values,
- evaluate stability and convergence of the ranking results.

### Sampling Sizes

To assess robustness, we report results using three sampling regimes:

- **1k samples** – baseline configuration consistent with prior interactome-based studies  
- **10k samples** – computationally efficient and empirically stable  
- **100k samples** – high-resolution reference for convergence analysis  

A detailed robustness and convergence analysis, including variability metrics and theoretical considerations, is provided in the supplementary material.

---

## ⏱️ Computational Complexity

DEviRank enumerates simple paths of bounded length (L ≤ 3) between drug target genes and disease-associated genes in the PPI network.

Let:
- |S| denote the number of drug target genes,
- |T| denote the number of disease-associated genes,
- Δ denote the maximum node degree in the PPI network.

The computation of disease gene weights scales with the local neighborhood size and requires:

O(∑_{t ∈ T} deg(t)).

Because path enumeration is restricted to simple paths of maximum length 3, the number of candidate paths originating from a drug target s is bounded by O(deg(s) Δ²). Summing over all s ∈ S yields a per-drug runtime of:

O(∑_{t ∈ T} deg(t) + ∑_{s ∈ S} deg(s) Δ²),

and in the worst case:

O(|S| Δ³).

Since drug evaluations are independent, DEviRank is embarrassingly parallel and scales linearly with the number of drugs.

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
