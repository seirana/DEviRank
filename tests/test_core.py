import json
import random
import networkx as nx
import pandas as pd
import pytest

from scr.DEviRank import (
    FINITE_INFINITY,
    _build_graph_from_ppi,
    _find_repo_root,
    calculate_proximity,
    pick_random_nodes_matching_selected,
    read_csv,
    shortest_path_length,
    write_csv,
)
from scr.reproducibility import write_run_metadata


def test_find_repo_root_does_not_require_experiments_directory(tmp_path):
    repo = tmp_path / "DEviRank"
    nested = repo / "scr" / "nested"
    (repo / "data").mkdir(parents=True)
    nested.mkdir(parents=True)

    assert _find_repo_root(nested) == repo


def test_build_graph_requires_ppi_columns():
    ppi = pd.DataFrame({"gene1": ["A"]})

    with pytest.raises(ValueError, match="gene2"):
        _build_graph_from_ppi(ppi)


def test_shortest_path_handles_missing_nodes_without_mutating_graph():
    graph = nx.Graph()
    graph.add_edge("A", "B")
    before = set(graph.nodes())

    assert shortest_path_length(graph, "A", "B") == 1
    assert shortest_path_length(graph, "MISSING", "MISSING") == 0
    assert shortest_path_length(graph, "A", "MISSING") == FINITE_INFINITY
    assert set(graph.nodes()) == before


def test_empty_proximity_input_returns_finite_infinity_sentinel():
    graph = nx.Graph()
    graph.add_edge("A", "B")

    result = calculate_proximity(
        graph,
        [],
        ["B"],
        n_random=10,
    )

    assert result[0] == FINITE_INFINITY
    assert result[5] == 0
    assert result[6] == 1


def test_proximity_rejects_non_positive_sampling():
    graph = nx.Graph()
    graph.add_edge("A", "B")

    with pytest.raises(ValueError, match="n_random"):
        calculate_proximity(
            graph,
            ["A"],
            ["B"],
            n_random=0,
        )


def test_degree_matched_sampling_is_seeded_and_does_not_change_global_rng():
    graph = nx.cycle_graph(["A", "B", "C", "D", "E", "F"])
    bins = [(2, 2, list(graph.nodes()))]

    random.seed(123)
    first = random.random()

    sample_1 = pick_random_nodes_matching_selected(
        graph,
        bins,
        ["A", "B"],
        n_random=5,
        seed=7,
    )
    after = random.random()

    random.seed(123)
    assert first == random.random()
    assert after == random.random()

    sample_2 = pick_random_nodes_matching_selected(
        graph,
        bins,
        ["A", "B"],
        n_random=5,
        seed=7,
    )

    def normalize(samples):
        return [sorted(sample) for sample in samples]

    assert normalize(sample_1) == normalize(sample_2)


def test_csv_round_trip_is_portable(tmp_path):
    expected = pd.DataFrame(
        {
            "drug": ["A", "B"],
            "score": [1.5, 2.5],
        }
    )

    path = write_csv(expected, tmp_path / "results")
    observed = read_csv(path)

    assert path == (tmp_path / "results.csv").resolve()
    pd.testing.assert_frame_equal(observed, expected)


def test_run_metadata_records_parameters(tmp_path):
    path = write_run_metadata(
        tmp_path,
        command="unit-test",
        parameters={"seed": 42, "sampling_size": 100},
    )

    payload = json.loads(path.read_text(encoding="utf-8"))

    assert payload["command"] == "unit-test"
    assert payload["parameters"]["seed"] == 42
    assert payload["parameters"]["sampling_size"] == 100
    assert "python" in payload
    assert "packages" in payload
