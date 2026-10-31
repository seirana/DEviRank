"""Scientific checks for Monte Carlo adjustment and simple-path enumeration."""

import networkx as nx
import numpy as np
import pytest
from scipy import sparse

from scr.psc_revision import adjusted_pvalues, reference_null_distances, weighted_simple_path_scores
from scr.DEviRank import (
    add_fdr_columns,
    calculate_closest_distance,
    calculate_proximity,
    sum_weight_all_paths_between_a_disease_gene_and_a_drug_gene,
)
import pandas as pd
from scipy.stats import false_discovery_control


def test_core_adjustment_matches_scipy_and_preserves_missing_tests():
    frame = pd.DataFrame({"p-value": [0.001, 0.01, 0.7, float("nan"), 20000]})
    observed = add_fdr_columns(frame)
    for method in ["bh", "by"]:
        np.testing.assert_allclose(
            observed["p_" + method.upper()].iloc[:3],
            false_discovery_control(frame["p-value"].iloc[:3], method=method),
        )
        assert observed["p_" + method.upper()].iloc[3:].isna().all()


def test_numeric_gene_ids_and_configured_edge_cutoff():
    graph = nx.path_graph([1, 2, 3])
    ppi = pd.DataFrame({"gene1": [1, 2], "gene2": [2, 3], "max_ppi": [0.8, 0.9]})
    f = sum_weight_all_paths_between_a_disease_gene_and_a_drug_gene
    assert f(graph, 1, 3, ppi, max_path_length=1) == 0
    assert f(graph, 1, 3, ppi, max_path_length=2) > 0
    assert f(graph, 3, 3, ppi, max_path_length=1) == 1


def test_missing_test_has_no_numerical_pvalue():
    result = calculate_proximity(nx.path_graph(3), [], [0], 10)
    assert np.isnan(result[2])


def test_disconnected_penalty_changes_unreachable_distances_only():
    graph = nx.path_graph(3)
    graph.add_node(9)
    assert calculate_closest_distance(graph, [0, 9], [2], disconnected_penalty=100) == [2, 100]
    assert calculate_closest_distance(graph, [0, 9], [2], disconnected_penalty=1000) == [2, 1000]


def test_eligible_isolates_are_retained_but_out_of_universe_targets_are_not_tested():
    graph = nx.path_graph(3)
    absent = calculate_proximity(graph, [0, 9], [2], 10, disconnected_penalty=100)
    assert absent[0] == 2 and absent[5] == 1
    graph.add_node(9)
    present = calculate_proximity(graph, [0, 9], [2], 10, disconnected_penalty=100)
    assert present[0] == 51 and present[5] == 2


@pytest.mark.parametrize("seed", range(8))
def test_weighted_simple_paths_equal_explicit_enumeration(seed):
    rng = np.random.default_rng(seed)
    graph = nx.gnp_random_graph(8, 0.45, seed=seed)
    for a, b in graph.edges:
        graph[a][b]["weight"] = float(rng.uniform(0.1, 1))
    W = nx.to_scipy_sparse_array(graph, nodelist=list(range(8)), weight="weight", format="csr")
    endpoints = np.zeros(8)
    endpoints[[0, 3, 7]] = [1, 3, 2]
    actual = weighted_simple_path_scores(W, endpoints, 4)
    for L in range(1, 5):
        expected = np.zeros(8)
        for s in graph:
            for t in np.flatnonzero(endpoints):
                if s == t:
                    expected[s] += endpoints[t]
                    continue
                for path in nx.all_simple_paths(graph, s, int(t), cutoff=L):
                    weight = 1
                    for a, b in zip(path, path[1:]):
                        weight *= graph[a][b]["weight"]
                    expected[s] += weight * endpoints[t]
        np.testing.assert_allclose(actual[L], expected, rtol=1e-11, atol=1e-11)


def test_fdr_adjusts_full_valid_family_not_selected_subset():
    p = np.array([0.001, 0.01, 0.7, 20000, np.nan])
    bh, by = adjusted_pvalues(p)
    np.testing.assert_allclose(bh[:3], [0.003, 0.015, 0.7])
    assert np.isnan(bh[3:]).all()
    assert np.all(by[:3] >= bh[:3])


def test_zero_monte_carlo_tail_is_positive():
    B = 100000
    p = (0 + 1) / (B + 1)
    assert 0 < p < 0.000011


def test_null_target_count_preserved_in_crowded_degree_bin():
    graph = nx.cycle_graph(8)
    bins = [(2, 2, list(graph))]
    nearest = {n: 10 for n in graph}
    a, b = reference_null_distances(graph, bins, [0, 1, 2, 3, 4], nearest, 1000, 11)
    np.testing.assert_array_equal(a, np.full(1000, 10.0))
    np.testing.assert_array_equal(b, np.zeros(1000))


def test_unreachable_affine_penalty_and_shared_target():
    W = sparse.csr_matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]], dtype=float)
    scores = weighted_simple_path_scores(W, [1, 0, 0], 4)
    for L in scores:
        np.testing.assert_allclose(scores[L], [1, 1, 0])
