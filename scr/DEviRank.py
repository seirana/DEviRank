#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

"""
DEviRank / Nbisdes proximity + scoring pipeline.

Design (as requested):
  - disease_file: can be anywhere (absolute/relative path)
  - static project data: always read from ./DEviRank/data/ (repo-local)
  - results: written to an output directory chosen by the user

Optional environment variable:
  - REPO_DIR : repo root (auto-inferred if not set)
"""

import math
import os
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Optional, Sequence, Union

import networkx as nx
import numpy as np
import pandas as pd


# ======================================================================================
# Portable paths
# ======================================================================================

PathLike = Union[str, Path]


def _resolve_repo_root() -> Path:
    """Resolve repo root (folder that contains 'data/' and 'experiments/')."""
    repo_env = os.environ.get("REPO_DIR")
    if repo_env:
        return Path(repo_env).expanduser().resolve()
    # If this file is in src/, repo root is parents[1]
    return Path(__file__).resolve().parents[1]


def _resolve_out_dir(out_dir: Optional[PathLike]) -> Path:
    """
    Resolve output directory.
    If None, default to <repo>/experiments/results
    """
    if out_dir is None:
        out = _resolve_repo_root() / "experiments" / "results"
    else:
        out = Path(out_dir).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    return out


@dataclass(frozen=True)
class ProjectPaths:
    repo_dir: Path
    data_dir: Path


def get_project_paths() -> ProjectPaths:
    repo_dir = _resolve_repo_root()
    data_dir = repo_dir / "data"
    if not data_dir.exists():
        raise FileNotFoundError(
            f"Repo-local data directory not found: {data_dir}\n"
            "Expected: <repo>/data/ (i.e., ./DEviRank/data/ after cloning)."
        )
    return ProjectPaths(repo_dir=repo_dir, data_dir=data_dir)


PATHS = get_project_paths()
DEVI_DATA = PATHS.data_dir


# ======================================================================================
# CSV read & write (robust)
# ======================================================================================

def _ensure_csv_suffix(p: Path) -> Path:
    # Only append .csv if there is no suffix at all.
    return p if p.suffix else p.with_suffix(".csv")


def read_csv(path_like: PathLike, *, must_exist: bool = True, **read_kwargs: Any) -> pd.DataFrame:
    p = _ensure_csv_suffix(Path(path_like)).expanduser().resolve()
    if must_exist and not p.exists():
        raise FileNotFoundError(f"CSV not found: {p}")
    try:
        return pd.read_csv(p, **read_kwargs)
    except Exception as e:
        raise RuntimeError(f"Failed to read CSV: {p}\nReason: {e}") from e


def write_csv(
    df: pd.DataFrame,
    path_like: PathLike,
    *,
    index: bool = False,
    mkdir: bool = True,
    atomic: bool = True,
    **to_csv_kwargs: Any,
) -> Path:
    p = _ensure_csv_suffix(Path(path_like)).expanduser().resolve()
    if mkdir:
        p.parent.mkdir(parents=True, exist_ok=True)

    try:
        if atomic:
            tmp = p.with_suffix(p.suffix + ".tmp")
            df.to_csv(tmp, index=index, **to_csv_kwargs)
            tmp.replace(p)
        else:
            df.to_csv(p, index=index, **to_csv_kwargs)
    except Exception as e:
        raise RuntimeError(f"Failed to write CSV: {p}\nReason: {e}") from e

    return p


# ======================================================================================
# Core constants
# ======================================================================================

FINITE_INFINITY = 20000  # > longest path in gene network
DEFAULT_SEED = 452456
DEFAULT_MIN_BIN_SIZE = 100
DEFAULT_CHUNK_SIZE = 1000


# ======================================================================================
# Network builders
# ======================================================================================

def _build_graph_from_ppi(ppi: pd.DataFrame, *, include_all_genes: Optional[Iterable[str]] = None) -> nx.Graph:
    """
    Build an undirected graph from a PPI dataframe with columns: gene1, gene2.
    If include_all_genes is provided, all those genes are added as isolated nodes too.
    """
    required = {"gene1", "gene2"}
    missing = required - set(ppi.columns)
    if missing:
        raise ValueError(f"PPI dataframe missing columns: {sorted(missing)}")

    g = nx.Graph()

    if include_all_genes is not None:
        for gene in include_all_genes:
            g.add_node(gene)

    for gene1, gene2 in zip(ppi["gene1"], ppi["gene2"]):
        g.add_edge(gene1, gene2)

    return g


def matrix_to_network_Nbisdes() -> nx.Graph:
    """Nbisdes: PPI network with genes present in the PPI."""
    ppi = read_csv(DEVI_DATA / "gene_gene_PPI700_ENSEMBL")
    all_genes = pd.concat([ppi["gene1"], ppi["gene2"]]).drop_duplicates()
    return _build_graph_from_ppi(ppi, include_all_genes=all_genes)


def matrix_to_network_DEviRank() -> nx.Graph:
    """DEviRank: PPI network plus all protein-coding genes (isolated nodes allowed)."""
    pcg = read_csv(DEVI_DATA / "protein_coding_genes_ENSEMBL")
    ppi = read_csv(DEVI_DATA / "gene_gene_PPI700_ENSEMBL")

    if "Gene stable ID" not in pcg.columns:
        raise ValueError("protein_coding_genes_ENSEMBL must include column 'Gene stable ID'")

    all_genes = pd.concat([pcg["Gene stable ID"], ppi["gene1"], ppi["gene2"]]).drop_duplicates()
    return _build_graph_from_ppi(ppi, include_all_genes=all_genes)


# ======================================================================================
# Shortest path utilities (finite infinity)
# ======================================================================================

def _bidirectional_pred_succ(G: nx.Graph, source: str, target: str):
    if target == source:
        return ({target: None}, {source: None}, source)

    Gpred = G.adj
    Gsucc = G.adj

    pred = {source: None}
    succ = {target: None}

    forward_fringe = [source]
    reverse_fringe = [target]
    while forward_fringe and reverse_fringe:
        if len(forward_fringe) <= len(reverse_fringe):
            this_level = forward_fringe
            forward_fringe = []
            for v in this_level:
                for w in Gsucc[v]:
                    if w not in pred:
                        forward_fringe.append(w)
                        pred[w] = v
                    if w in succ:
                        return pred, succ, w
        else:
            this_level = reverse_fringe
            reverse_fringe = []
            for v in this_level:
                for w in Gpred[v]:
                    if w not in succ:
                        succ[w] = v
                        reverse_fringe.append(w)
                    if w in pred:
                        return pred, succ, w

    return FINITE_INFINITY


def bidirectional_shortest_path(G: nx.Graph, source: str, target: str):
    if source not in G or target not in G:
        raise nx.NodeNotFound(f"Either source {source} or target {target} is not in G")

    results = _bidirectional_pred_succ(G, source, target)
    if results == FINITE_INFINITY:
        return FINITE_INFINITY

    pred, succ, w = results

    path = []
    while w is not None:
        path.append(w)
        w = pred[w]
    path.reverse()

    w = succ[path[-1]]
    while w is not None:
        path.append(w)
        w = succ[w]
    return path


def shortest_path_length(G: nx.Graph, source: str, target: str) -> int:
    p = bidirectional_shortest_path(G, source, target)
    if p == FINITE_INFINITY:
        return FINITE_INFINITY
    return len(p) - 1


# ======================================================================================
# Degree-aware random sampling
# ======================================================================================

def get_degree_binning(g: nx.Graph, bin_size: int, lengths=None):
    degree_to_nodes = {}
    for node, degree in g.degree():
        if lengths is not None and node not in lengths:
            continue
        degree_to_nodes.setdefault(degree, []).append(node)

    values = sorted(degree_to_nodes.keys())

    bins = []
    i = 0
    while i < len(values):
        low = values[i]
        val = degree_to_nodes[values[i]]
        while len(val) < bin_size:
            i += 1
            if i == len(values):
                break
            val.extend(degree_to_nodes[values[i]])
        if i == len(values):
            i -= 1
        high = values[i]
        i += 1

        if len(val) < bin_size and bins:
            low_, high_, val_ = bins[-1]
            bins[-1] = (low_, high, val_ + val)
        else:
            bins.append((low, high, val))
    return bins


def get_degree_equivalents(seeds: Sequence[str], bins, g: nx.Graph):
    seed_to_nodes = {}
    for seed in seeds:
        d = g.degree(seed)
        for l, h, nodes in bins:
            if l <= d <= h:
                mod_nodes = list(nodes)
                if seed in mod_nodes:
                    mod_nodes.remove(seed)
                seed_to_nodes[seed] = mod_nodes
                break
    return seed_to_nodes


def pick_random_nodes_matching_selected(
    network: nx.Graph,
    bins,
    nodes_selected: Sequence[str],
    n_random: int,
    *,
    degree_aware: bool = True,
    connected: bool = False,
    seed: Optional[int] = None,
):
    if seed is not None:
        random.seed(seed)

    nodes = list(network.nodes())
    values = []

    for _ in range(n_random):
        if degree_aware:
            if connected:
                raise ValueError("connected=True not implemented for degree_aware=True")
            nodes_random = set()
            node_to_equiv = get_degree_equivalents(nodes_selected, bins, network)
            for _, equiv_nodes in node_to_equiv.items():
                if not equiv_nodes:
                    continue
                chosen = random.choice(equiv_nodes)
                for _k in range(20):
                    if chosen in nodes_random:
                        chosen = random.choice(equiv_nodes)
                nodes_random.add(chosen)
            values.append(list(nodes_random))
        else:
            if connected:
                nodes_random = [random.choice(nodes)]
                k = 1
                while k < len(nodes_selected):
                    node_random = random.choice(nodes_random)
                    neighbor = random.choice(list(network.neighbors(node_random)))
                    if neighbor in nodes_random:
                        continue
                    nodes_random.append(neighbor)
                    k += 1
                values.append(nodes_random)
            else:
                values.append(random.sample(sorted(nodes), len(nodes_selected)))

    return values


def get_random_nodes(
    nodes: Sequence[str],
    network: nx.Graph,
    *,
    bins=None,
    n_random: int = 1000,
    min_bin_size: int = DEFAULT_MIN_BIN_SIZE,
    degree_aware: bool = True,
    seed: Optional[int] = DEFAULT_SEED,
):
    if bins is None:
        bins = get_degree_binning(network, min_bin_size)
    return pick_random_nodes_matching_selected(network, bins, nodes, n_random, degree_aware=degree_aware, seed=seed)


# ======================================================================================
# Proximity
# ======================================================================================

def calculate_closest_distance(network: nx.Graph, nodes_from: Iterable[str], nodes_to: Iterable[str]) -> list[int]:
    nodes_to = list(nodes_to)
    out = []
    for node_from in nodes_from:
        dists = [shortest_path_length(network, node_from, node_to) for node_to in nodes_to]
        out.append(min(dists))
    return out


def calculate_proximity(network: nx.Graph, nodes_from: Sequence[str], nodes_to: Sequence[str], n_random: int):
    nodes_network = set(network.nodes())

    nodes_from_pruned = set(nodes_from) & nodes_network
    nodes_to_pruned = set(nodes_to) & nodes_network

    if not nodes_from_pruned or not nodes_to_pruned:
        if not (nodes_from_pruned & nodes_to_pruned):
            return (
                FINITE_INFINITY, FINITE_INFINITY, FINITE_INFINITY, FINITE_INFINITY, FINITE_INFINITY,
                len(nodes_from_pruned), len(nodes_to_pruned), FINITE_INFINITY
            )

    shortest_distances = calculate_closest_distance(network, nodes_from_pruned, nodes_to_pruned)
    d = float(np.mean(shortest_distances))
    all_distances = [list(shortest_distances)]

    bins = get_degree_binning(network, DEFAULT_MIN_BIN_SIZE, lengths=None)
    nodes_from_random = get_random_nodes(
        list(nodes_from_pruned),
        network,
        bins=bins,
        n_random=n_random,
        min_bin_size=DEFAULT_MIN_BIN_SIZE,
        seed=DEFAULT_SEED,
    )

    nodes_to_fixed = list(nodes_to_pruned)
    values = np.empty(n_random, dtype=float)

    for i in range(n_random):
        frm = nodes_from_random[i] if not isinstance(nodes_from_random[i], int) else [nodes_from_random[i]]
        dists = calculate_closest_distance(network, frm, nodes_to_fixed)
        values[i] = float(np.mean(dists))

    pval = float(np.sum(values <= d)) / len(values)
    m, s = float(np.mean(values)), float(np.std(values))
    z = 0.0 if s == 0 else (d - m) / s

    return d, z, pval, m, s, len(nodes_from_pruned), len(nodes_to_pruned), all_distances


def calculate_proximity_collecting(
    disease_file: PathLike,
    *,
    out_dir: PathLike,
    sampling: Optional[int] = None,
    which_method: str = "DEviRank",
    chunk_size: int = DEFAULT_CHUNK_SIZE,
):
    """
    Compute drug proximity scores against a disease gene set.

    Writes outputs to: out_dir/

    Rules:
      - DEviRank: sampling defaults to 100k if not provided; chunked by chunk_size drugs per part.
      - Nbisdes: sampling fixed to 1000; runs as a single chunk.
    """
    out_dir = _resolve_out_dir(out_dir)

    disease_to_genes = read_csv(disease_file)
    drug_to_targets = read_csv(DEVI_DATA / "DtoDGI_ENSEMBL(filtered)")
    drugs = read_csv(DEVI_DATA / "drugs(filtered)")
    check = read_csv(DEVI_DATA / "repeated(filtered)")

    if "ENSEMBL ID" not in disease_to_genes.columns:
        raise ValueError("disease_file must contain column 'ENSEMBL ID'")
    if "drug_name" not in drugs.columns:
        raise ValueError("drugs(filtered) must contain column 'drug_name'")

    l = len(drug_to_targets)

    which_method = str(which_method)
    if which_method == "DEviRank":
        network = matrix_to_network_DEviRank()
        sampling_eff = 100000 if sampling is None else int(sampling)
        parts = int(math.ceil(l / chunk_size))
        chunk = int(chunk_size)
    elif which_method == "Nbisdes":
        network = matrix_to_network_Nbisdes()
        sampling_eff = 1000
        parts = 1
        chunk = l
    else:
        raise ValueError("which_method must be 'DEviRank' or 'Nbisdes'")

    cols = [
        "drug", "shortest distance", "z score", "p-value", "mean", "SD",
        "number of drug target genes with PPI", "number of disease target genes with PPI",
        "distances", "number of drug target genes", "number of disease target genes"
    ]

    nodes_to = list(disease_to_genes["ENSEMBL ID"])

    for part in range(parts):
        b = part * chunk
        e = min((part + 1) * chunk, l)

        output = pd.DataFrame(index=range(e - b), columns=cols)

        for i in range(b, e):
            ind = [j for j, x in enumerate(drug_to_targets.loc[i, :]) if x > 0]
            nodes_from = list(drug_to_targets.iloc[i, ind])
            drug_name = drugs.loc[i, "drug_name"]

            if int(check.iloc[i, 0]) == i:
                d, z, pval, m, s, len_from_ppi, len_to_ppi, all_distances = calculate_proximity(
                    network, nodes_from, nodes_to, sampling_eff
                )
                output.loc[i - b, :] = [
                    drug_name, d, z, pval, m, s,
                    len_from_ppi, len_to_ppi,
                    all_distances, len(nodes_from), len(nodes_to)
                ]
            else:
                output.loc[i - b, :] = [drug_name, "", "", "", "", "", "", "", "", "", ""]

        write_csv(output, out_dir / f"drug_scores_{which_method}_{part}")

    if which_method == "DEviRank":
        merged = pd.concat(
            [read_csv(out_dir / f"drug_scores_DEviRank_{part}") for part in range(parts)],
            axis=0,
            ignore_index=True,
        )
        write_csv(merged, out_dir / "drug_scores_DEviRank")

        for part in range(parts):
            p = _ensure_csv_suffix(out_dir / f"drug_scores_DEviRank_{part}")
            if p.exists():
                p.unlink()


# ======================================================================================
# Compare DEviRank vs Nbisdes
# ======================================================================================

def intersection(lst1: Sequence[int], lst2: Sequence[int]) -> list[int]:
    s2 = set(lst2)
    return [v for v in lst1 if v in s2]


def DEviRankVSNbisdes(*, out_dir: PathLike):
    out_dir = _resolve_out_dir(out_dir)

    nbisdes = read_csv(out_dir / "drug_scores_Nbisdes")
    devirank = read_csv(out_dir / "drug_scores_DEviRank")

    p_value = 0.05
    z_score = -1.96
    ind_z = [j for j, x in enumerate(pd.to_numeric(devirank["z score"], errors="coerce")) if x <= z_score]
    ind_p = [j for j, x in enumerate(pd.to_numeric(devirank["p-value"], errors="coerce")) if x <= p_value]
    ind_devirank = intersection(ind_z, ind_p)

    p_value = 1
    z_score = -0.15
    ind_z = [j for j, x in enumerate(pd.to_numeric(nbisdes["z score"], errors="coerce")) if x <= z_score]
    ind_p = [j for j, x in enumerate(pd.to_numeric(nbisdes["p-value"], errors="coerce")) if x <= p_value]
    ind_nbisdes = intersection(ind_z, ind_p)

    ind_both = intersection(ind_devirank, ind_nbisdes)

    db = pd.DataFrame(
        {
            "DEviRank_method": [len(ind_devirank)],
            "Nbisdes_original method": [len(ind_nbisdes)],
            "mutually suggested drugs": [len(ind_both)],
        }
    )
    print(db)
    write_csv(db, out_dir / "DEviRankVSNbisdes")


# ======================================================================================
# Drug scoring
# ======================================================================================

def matrix_to_network(
    drug_genes: pd.DataFrame,
    ppi: pd.DataFrame,
    disease_genes: pd.DataFrame,
    pc_genes: pd.DataFrame,
) -> nx.Graph:
    network = nx.Graph()

    all_drug_genes = pd.DataFrame(drug_genes.values.ravel(), columns=["merged"])
    genes = pd.concat(
        [
            pc_genes["Gene stable ID"],
            ppi["gene1"],
            ppi["gene2"],
            disease_genes["ENSEMBL ID"],
            all_drug_genes["merged"],
        ],
        axis=0,
    ).drop_duplicates().sort_values()

    for g in genes:
        network.add_node(g)

    for g1, g2 in zip(ppi["gene1"], ppi["gene2"]):
        network.add_edge(g1, g2)

    return network


def weight_a_disease_gene(network: nx.Graph, disease_gene: str, ppi: pd.DataFrame) -> float:
    weight = 1.0
    if network.has_node(disease_gene) and any(True for _ in network.neighbors(disease_gene)):
        inds = [k for k, x in enumerate(ppi["gene1"]) if x == disease_gene]
        for j in inds:
            weight += float(ppi.loc[j, "max_ppi"])
    return weight


def PPI_weight(gene1: str, gene2: str, ppi: pd.DataFrame) -> float:
    m1 = (ppi["gene1"] == gene1) & (ppi["gene2"] == gene2)
    m2 = (ppi["gene1"] == gene2) & (ppi["gene2"] == gene1)
    hits = ppi.loc[m1 | m2, "max_ppi"]
    if hits.empty:
        return 0.0
    return float(hits.iloc[0])


def sum_weight_all_paths_between_a_disease_gene_and_a_drug_gene(
    network: nx.Graph,
    drug_gene: str,
    disease_gene: str,
    ppi: pd.DataFrame,
) -> float:
    weight = 0.0
    for path in nx.all_simple_paths(network, source=drug_gene, target=disease_gene, cutoff=2):
        weight_path = 1.0
        for k in range(len(path) - 1):
            weight_path *= PPI_weight(path[k], path[k + 1], ppi)
        weight += weight_path * weight_a_disease_gene(network, disease_gene, ppi)
    return weight


def weight_a_drug(
    network: nx.Graph,
    drug_idx: int,
    dgi: pd.DataFrame,
    drug_genes: pd.DataFrame,
    ppi: pd.DataFrame,
    disease_genes: pd.DataFrame,
) -> float:
    total = 0.0
    clmns = [k for k, x in enumerate(drug_genes.iloc[drug_idx, :]) if x > 0]
    for col_idx in clmns:
        drug_gene = drug_genes.iloc[drug_idx, col_idx]
        gene_weight = 0.0
        for j in range(len(disease_genes)):
            gene_weight += sum_weight_all_paths_between_a_disease_gene_and_a_drug_gene(
                network,
                str(drug_gene),
                str(disease_genes.iloc[j, 0]),
                ppi,
            )
        total += gene_weight * float(dgi.iloc[drug_idx, col_idx])
    return total


def score_drugs(
    network: nx.Graph,
    dgi: pd.DataFrame,
    drug_genes: pd.DataFrame,
    drugs: pd.DataFrame,
    ppi: pd.DataFrame,
    disease_genes: pd.DataFrame,
    candidate_drugs: pd.DataFrame,
    *,
    p_value: float = 0.05,
    z_score: float = -1.96,
) -> np.ndarray:
    weights = np.zeros((len(drugs), 1), dtype=float)

    z_vals = pd.to_numeric(candidate_drugs["z score"], errors="coerce")
    p_vals = pd.to_numeric(candidate_drugs["p-value"], errors="coerce")
    selected = (z_vals <= z_score) & (p_vals <= p_value)

    c = 0
    for drug_idx in range(len(drugs)):
        if selected.iloc[drug_idx]:
            c += 1
            print(c)
            weights[drug_idx, 0] = weight_a_drug(network, drug_idx, dgi, drug_genes, ppi, disease_genes)

    return weights


def drug_scoring(
    disease_file: PathLike,
    *,
    out_dir: PathLike,
    p_value: float = 0.05,
    z_score: float = -1.96,
):
    out_dir = _resolve_out_dir(out_dir)

    dgi = read_csv(DEVI_DATA / "DtoGI_scores(filtered)")
    drug_genes = read_csv(DEVI_DATA / "DtoGI_ENSEMBL(filtered)")
    drugs = read_csv(DEVI_DATA / "drugs(filtered)")
    pc_genes = read_csv(DEVI_DATA / "protein_coding_genes_ENSEMBL")

    ppi = read_csv(DEVI_DATA / "gene_gene_PPI700_ENSEMBL").fillna(0)
    if "max_ppi" in ppi.columns:
        ppi["max_ppi"] = ppi["max_ppi"] / 1000

    disease_genes = read_csv(disease_file)
    candidate_drugs = read_csv(out_dir / "drug_scores_DEviRank")

    network = matrix_to_network(drug_genes, ppi, disease_genes, pc_genes)

    scores = score_drugs(
        network, dgi, drug_genes, drugs, ppi, disease_genes, candidate_drugs,
        p_value=p_value, z_score=z_score
    )
    scores_df = pd.DataFrame(scores, columns=["DEviRank_Drug_score"])

    out = pd.concat([candidate_drugs, scores_df], axis=1)
    write_csv(out, out_dir / "scoring_drugs")


# ======================================================================================
# Disease target gene scoring
# ======================================================================================

def sum_weight_all_paths_between_a_drug_gene_and_a_disease_gene(
    network: nx.Graph,
    drug_gene: str,
    disease_gene: str,
    ppi: pd.DataFrame,
) -> float:
    weight = 0.0
    for path in nx.all_simple_paths(network, source=drug_gene, target=disease_gene, cutoff=2):
        weight_path = 1.0
        for k in range(len(path) - 1):
            weight_path *= PPI_weight(path[k], path[k + 1], ppi)
        weight += weight_path
    return weight


def weight_a_drug_for_disease_gene(
    network: nx.Graph,
    drug_idx: int,
    disease_idx: int,
    ppi: pd.DataFrame,
    drug_genes: pd.DataFrame,
    disease_genes: pd.DataFrame,
    dgi: pd.DataFrame,
) -> float:
    total = 0.0
    clmns = [k for k, x in enumerate(drug_genes.iloc[drug_idx, :]) if x > 0]
    disease_gene = str(disease_genes.iloc[disease_idx, 0])

    for col_idx in clmns:
        drug_gene = str(drug_genes.iloc[drug_idx, col_idx])
        w = sum_weight_all_paths_between_a_drug_gene_and_a_disease_gene(network, drug_gene, disease_gene, ppi)
        total += w * float(dgi.iloc[drug_idx, col_idx])

    return total


def score_disease_genes(
    network: nx.Graph,
    ppi: pd.DataFrame,
    drug_genes: pd.DataFrame,
    disease_genes: pd.DataFrame,
    drugs: pd.DataFrame,
    dgi: pd.DataFrame,
    candidate_drugs: pd.DataFrame,
    *,
    p_value: float = 0.05,
    z_score: float = -1.96,
) -> np.ndarray:
    weights = np.zeros((len(drugs), len(disease_genes)), dtype=float)

    z_vals = pd.to_numeric(candidate_drugs["z score"], errors="coerce")
    p_vals = pd.to_numeric(candidate_drugs["p-value"], errors="coerce")
    selected = (z_vals <= z_score) & (p_vals <= p_value)

    c = 0
    for drug_idx in range(len(drugs)):
        if selected.iloc[drug_idx]:
            c += 1
            print(c)
            for disease_idx in range(len(disease_genes)):
                weights[drug_idx, disease_idx] = weight_a_drug_for_disease_gene(
                    network, drug_idx, disease_idx, ppi, drug_genes, disease_genes, dgi
                )

    return weights


def disease_target_gene_scoring(
    disease_file: PathLike,
    *,
    out_dir: PathLike,
    p_value: float = 0.05,
    z_score: float = -1.96,
):
    out_dir = _resolve_out_dir(out_dir)

    dgi = read_csv(DEVI_DATA / "DtoGI_scores(filtered)")
    drug_genes = read_csv(DEVI_DATA / "DtoGI_ENSEMBL(filtered)")
    drugs = read_csv(DEVI_DATA / "drugs(filtered)")
    pc_genes = read_csv(DEVI_DATA / "protein_coding_genes_ENSEMBL")

    ppi = read_csv(DEVI_DATA / "gene_gene_PPI700_ENSEMBL").fillna(0)
    if "max_ppi" in ppi.columns:
        ppi["max_ppi"] = ppi["max_ppi"] / 1000

    disease_genes = read_csv(disease_file)
    candidate_drugs = read_csv(out_dir / "drug_scores_DEviRank")

    network = matrix_to_network(drug_genes, ppi, disease_genes, pc_genes)

    scores = score_disease_genes(
        network, ppi, drug_genes, disease_genes, drugs, dgi, candidate_drugs,
        p_value=p_value, z_score=z_score
    )

    scores_df = pd.DataFrame(scores)
    if "ENSEMBL ID" in disease_genes.columns:
        scores_df.columns = list(disease_genes["ENSEMBL ID"])

    out = pd.concat([candidate_drugs.reset_index(drop=True), scores_df.reset_index(drop=True)], axis=1)
    write_csv(out, out_dir / "drug_scores_DEviRank_gene_scores")


# ======================================================================================
# Public entry points
# ======================================================================================

def suggested_drugs_DEviRank(
    disease_file: PathLike,
    *,
    out_dir: PathLike,
    sampling_size: Optional[int] = None,
    p_value: float = 0.05,
    z_score: float = -1.96,
):
    """
    Full DEviRank pipeline:
      1) proximity + random sampling  -> out_dir/drug_scores_DEviRank.csv
      2) drug scoring                 -> out_dir/scoring_drugs.csv
      3) disease-gene scoring         -> out_dir/drug_scores_DEviRank_gene_scores.csv
    """
    calculate_proximity_collecting(
        disease_file,
        out_dir=out_dir,
        sampling=sampling_size,
        which_method="DEviRank",
    )
    drug_scoring(disease_file, out_dir=out_dir, p_value=p_value, z_score=z_score)
    disease_target_gene_scoring(disease_file, out_dir=out_dir, p_value=p_value, z_score=z_score)


def compare_DEviRank_Nbisdes(
    disease_file: PathLike,
    *,
    out_dir: PathLike,
    sampling_size: Optional[int] = None,
):
    """
    Runs DEviRank + Nbisdes and writes comparison table to out_dir/DEviRankVSNbisdes.csv
    """
    calculate_proximity_collecting(
        disease_file,
        out_dir=out_dir,
        sampling=sampling_size,
        which_method="DEviRank",
    )
    calculate_proximity_collecting(
        disease_file,
        out_dir=out_dir,
        sampling=sampling_size,
        which_method="Nbisdes",
    )
    DEviRankVSNbisdes(out_dir=out_dir)


if __name__ == "__main__":
    # Intentionally empty. Use experiments/run_devirank.py as CLI entry point.
    pass
