"""Reanalyse the saved PSC results and test the available repository inputs.

The saved supplement and repository disease modules differ. They are analysed
separately; this script never substitutes fresh results for the historical table.
Run: python -m scr.psc_revision --output experiments/psc_revision
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time

import networkx as nx
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import false_discovery_control, spearmanr

from scr import DEviRank as core


def adjusted_pvalues(pvalues):
    """Adjust the complete valid test family, leaving unavailable tests as NaN."""
    p = np.asarray(pvalues, dtype=float)
    valid = np.isfinite(p) & (p >= 0) & (p <= 1)
    bh = np.full(p.shape, np.nan)
    by = np.full(p.shape, np.nan)
    if valid.any():
        bh[valid] = false_discovery_control(p[valid], method="bh")
        by[valid] = false_discovery_control(p[valid], method="by")
    return bh, by


def reanalyse_saved(workbook: Path, output: Path):
    summaries = []
    final = None
    for table, B in [(2, 1000), (3, 10000), (4, 100000)]:
        df = pd.read_excel(workbook, sheet_name=f"Supplementary Table {table}", header=2)
        df = df.rename(columns={df.columns[0]: "drug"})
        raw = pd.to_numeric(df["p-value"], errors="coerce").to_numpy()
        valid = np.isfinite(raw) & (raw >= 0) & (raw <= 1)
        counts = np.rint(raw[valid] * B)
        if not np.allclose(raw[valid] * B, counts, atol=1e-7, rtol=0):
            raise ValueError("Saved p-values do not reconstruct integer tail counts")
        df["p_original"] = raw
        df["null_tail_count"] = np.nan
        df.loc[valid, "null_tail_count"] = counts
        df["p_MC"] = np.nan
        df.loc[valid, "p_MC"] = (counts + 1) / (B + 1)
        df["p_BH"], df["p_BY"] = adjusted_pvalues(df["p_MC"])
        df["test_status"] = np.where(valid, "evaluable", "unavailable_in_saved_output")
        df["original_exploratory"] = valid & (raw <= 0.05) & (df["z score"] <= -1.96)
        df["corrected_exploratory"] = valid & (df["p_MC"] <= 0.05) & (df["z score"] <= -1.96)
        df["BH_0.05"] = df["p_BH"] <= 0.05
        df["BY_0.05"] = df["p_BY"] <= 0.05
        # Invalid sentinels are retained in p_original for audit, not in p-value.
        df["p-value"] = df["p_MC"]
        df.to_csv(output / f"saved_PSC_{B}_FDR.csv", index=False)
        summaries.append(
            {
                "samples": B,
                "screened": len(df),
                "evaluable": int(valid.sum()),
                "unavailable": int((~valid).sum()),
                "original_exploratory": int(df.original_exploratory.sum()),
                "corrected_exploratory": int(df.corrected_exploratory.sum()),
                "BH": int(df["BH_0.05"].sum()),
                "BY": int(df["BY_0.05"].sum()),
                "original_and_BH": int((df.original_exploratory & df["BH_0.05"]).sum()),
                "original_and_BY": int((df.original_exploratory & df["BY_0.05"]).sum()),
                "zero_original_p": int((raw[valid] == 0).sum()),
            }
        )
        if B == 100000:
            final = df
    pd.DataFrame(summaries).to_csv(output / "saved_FDR_summary.csv", index=False)
    # Preserve the historical ranking, adding evidence rather than redefining it.
    ranking = pd.read_excel(workbook, sheet_name="Supplementary Table 6", header=3)
    ranking = ranking.iloc[:, :5]
    ranking.columns = [
        "drug",
        "raw_score_saved",
        "raw_rank_saved",
        "normalized_score_saved",
        "normalized_rank_saved",
    ]
    ranking = ranking[ranking.drug.notna()].copy()
    for column in ranking.columns[1:]:
        ranking[column] = pd.to_numeric(ranking[column], errors="raise")
    ranking = ranking.merge(
        final[["drug", "null_tail_count", "p_original", "p_MC", "p_BH", "p_BY"]],
        on="drug",
        validate="one_to_one",
    )
    ranking.to_csv(output / "saved_36_candidates_FDR.csv", index=False)
    convergence = []
    loaded = {B: pd.read_csv(output / f"saved_PSC_{B}_FDR.csv") for B in (1000, 10000, 100000)}
    for first, second in [(1000, 100000), (10000, 100000), (1000, 10000)]:
        a, b = loaded[first], loaded[second]
        if a.drug.tolist() != b.drug.tolist():
            raise ValueError("Drug order differs between saved sampling tables")
        valid = a.p_MC.notna() & b.p_MC.notna()
        for column in ["z score", "p_MC", "mean", " Standard Deviation"]:
            diff = a.loc[valid, column].to_numpy() - b.loc[valid, column].to_numpy()
            convergence.append(
                {
                    "first_samples": first,
                    "second_samples": second,
                    "metric": column,
                    "evaluable": int(valid.sum()),
                    "SDD": float(diff.std(ddof=1)),
                    "signed_mean_difference": float(diff.mean()),
                    "range_of_differences": float(np.ptp(diff)),
                }
            )
    pd.DataFrame(convergence).to_csv(output / "saved_convergence_recomputed.csv", index=False)
    return final, ranking, summaries


def weighted_simple_path_scores(W, endpoint_weights, max_length=4):
    """Exact weighted simple-path sums through four edges for every start node.

    The matrix has zero diagonal, symmetric nonnegative weights, and no parallel
    edges. Non-backtracking edge messages count open walks; triangle and four-
    cycle corrections remove every remaining repeated-vertex walk. Length-zero
    paths are included for shared drug/disease targets, as in NetworkX 3.7.
    Tested against explicit all_simple_paths on small weighted graphs.
    """
    if max_length not in (1, 2, 3, 4):
        raise ValueError("max_length must be 1, 2, 3, or 4")
    W = sparse.csr_matrix(W, dtype=float)
    if np.any(W.diagonal() != 0):
        raise ValueError("Self loops must be removed before simple-path scoring")
    u = np.asarray(endpoint_weights, dtype=float)
    coo = W.tocoo()
    src, dst, w = coo.row, coo.col, coo.data
    n = W.shape[0]
    # CSR/COO entries have a unique reverse entry because W is symmetric.
    keys = src.astype(np.int64) * n + dst
    sorter = np.argsort(keys)
    reverse = sorter[np.searchsorted(keys[sorter], dst.astype(np.int64) * n + src)]
    message = w * u[dst]
    cumulative = u.copy()
    out = {}
    triangle = None
    repeated_triangle = None
    four_cycle = None
    if max_length >= 3:
        W2 = W @ W
        triangle = np.asarray(W2.multiply(W).sum(axis=1)).ravel()
        if max_length >= 4:
            Wsquare = W.multiply(W)
            q = np.asarray(Wsquare.sum(axis=1)).ravel()
            r4 = np.asarray(Wsquare.multiply(Wsquare).sum(axis=1)).ravel()
            four_cycle = np.asarray(W2.multiply(W2).sum(axis=1)).ravel() - q * q - Wsquare @ q + r4
            repeated_triangle = W2.multiply(Wsquare)
        del W2
    for k in range(1, max_length + 1):
        nb = np.bincount(src, weights=message, minlength=n)
        exact = nb.copy()
        if k == 3:
            exact -= triangle * u
        elif k == 4:
            exact -= triangle * (W @ u) + W @ (triangle * u)
            exact += 3 * (repeated_triangle @ u)
            exact -= four_cycle * u
        # Roundoff at nodes whose exact sum is zero must not create contributors.
        tol = np.finfo(float).eps * np.maximum(1, np.abs(nb)) * 100
        exact[np.abs(exact) <= tol] = 0
        if exact.min() < -1e-6:
            raise ArithmeticError("Negative path sum exceeds floating-point tolerance")
        exact = np.maximum(exact, 0)
        cumulative = cumulative + exact
        out[k] = cumulative.copy()
        if k < max_length:
            message = w * (nb[dst] - message[reverse])
    return out


def reference_weighted_matrix(graph, ppi):
    """Use the first undirected confidence match, as in core.PPI_weight."""
    nodes = list(graph.nodes())
    index = {node: i for i, node in enumerate(nodes)}
    pair = ppi.copy()
    pair["a"] = np.minimum(pair.gene1, pair.gene2)
    pair["b"] = np.maximum(pair.gene1, pair.gene2)
    pair = pair.drop_duplicates(["a", "b"], keep="first")
    pair = pair[pair.a != pair.b]
    a = np.array([index[int(x)] for x in pair.a], dtype=int)
    b = np.array([index[int(x)] for x in pair.b], dtype=int)
    w = pair.max_ppi.to_numpy(dtype=float) / 1000
    W = sparse.csr_matrix((np.r_[w, w], (np.r_[a, b], np.r_[b, a])), shape=(len(nodes), len(nodes)))
    # Match reference disease weights, including its first-endpoint aggregation.
    local = ppi.groupby("gene1").max_ppi.sum().to_dict()
    return W, nodes, index, local


def reference_null_distances(graph, bins, selected, nearest, B, seed):
    """Vectorised reference matching with exact unique-target-set sizes.

    Each target draws from its degree bin excluding itself; collisions with
    earlier draws in that bin are resampled. NumPy PCG64 is used explicitly.
    Common random sets are reused at every disconnected penalty.
    """
    rng = np.random.default_rng(seed)
    selected = sorted(set(int(x) for x in selected))
    if not selected:
        raise ValueError("Empty target set")
    by_degree = {}
    for low, high, nodes in bins:
        pool = np.array(sorted(nodes), dtype=np.int64)
        for degree in range(low, high + 1):
            by_degree[degree] = pool
    columns = []
    previous = {}
    for node in selected:
        pool = by_degree[graph.degree(node)]
        pos = int(np.searchsorted(pool, node))
        if len(pool) <= 1:
            raise ValueError("No degree-matched alternative")
        sampled = rng.integers(0, len(pool) - 1, size=B)
        sampled += sampled >= pos
        ids = pool[sampled]
        key = int(pool[0])
        earlier = previous.setdefault(key, [])
        for attempt in range(1000):
            collision = np.zeros(B, dtype=bool)
            for old in earlier:
                collision |= ids == old
            if not collision.any():
                break
            redrawn = rng.integers(0, len(pool) - 1, size=int(collision.sum()))
            redrawn += redrawn >= pos
            ids[collision] = pool[redrawn]
            sampled[collision] = redrawn
        else:
            raise RuntimeError("Could not preserve unique target count")
        earlier.append(ids)
        # Convert only the small eligible pool, not B Python lookups per target.
        # Collision resampling updates the eligible-pool indices as well.
        distance_by_pool = np.array([nearest.get(int(x), -1) for x in pool], dtype=np.int16)
        columns.append(distance_by_pool[sampled])
    distances = np.column_stack(columns)
    finite_sum = np.maximum(distances, 0).sum(axis=1).astype(float)
    missing = (distances < 0).sum(axis=1)
    return finite_sum / len(selected), missing / len(selected)


def jaccard(a, b):
    return len(a & b) / len(a | b) if a | b else 1.0


def run(args):
    start = time.perf_counter()
    root = Path(__file__).resolve().parents[1]
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    workbook = root / "supplementary" / "DEviRank-Supplementary tables.xlsx"
    saved, ranking, saved_summary = reanalyse_saved(workbook, output)
    print("Saved FDR analysis:", saved_summary, flush=True)
    graph = core.matrix_to_network_DEviRank()
    ppi = pd.read_csv(root / "data/gene_gene_PPI700_ENSEMBL.csv")
    targets = pd.read_csv(root / "data/DtoGI_ENSEMBL(filtered).csv").to_numpy()
    dgi = pd.read_csv(root / "data/DtoGI_scores(filtered).csv").to_numpy()
    names = pd.read_csv(root / "data/drugs(filtered).csv").drug_name.to_numpy()
    disease_path = Path(args.disease_file)
    if not disease_path.is_absolute():
        disease_path = root / disease_path
    disease = pd.read_csv(disease_path).iloc[:, 0].astype(int).tolist()
    if not set(disease) <= set(graph):
        raise ValueError("Disease IDs must belong to the predefined protein universe")
    nearest = nx.multi_source_dijkstra_path_length(graph, disease, weight=None)
    max_finite = max(nearest.values())
    if min(args.penalties) <= max_finite:
        raise ValueError("Every penalty must exceed the largest finite nearest distance")
    profiles = [
        tuple(sorted(set(int(v) for v in row if v > 0 and int(v) in graph))) for row in targets
    ]
    observed = np.array(
        [sum(nearest.get(v, 20000) for v in p) / len(p) if p else np.nan for p in profiles]
    )
    raw_counts = (targets > 0).sum(axis=1)
    eligible_counts = np.array([sum(v > 0 and int(v) in graph for v in row) for row in targets])
    denom_audit = ranking[["drug", "raw_score_saved", "normalized_score_saved"]].copy()
    denom_audit["saved_denominator"] = np.rint(
        denom_audit.raw_score_saved / denom_audit.normalized_score_saved
    ).astype(int)
    lookup = {name: i for i, name in enumerate(names)}
    denom_audit["raw_target_entries"] = [raw_counts[lookup[name]] for name in denom_audit.drug]
    denom_audit["eligible_target_entries"] = [
        eligible_counts[lookup[name]] for name in denom_audit.drug
    ]
    denominator_mismatches = int(
        (denom_audit.saved_denominator != denom_audit.eligible_target_entries).sum()
    )
    denom_audit.to_csv(output / "saved_denominator_audit.csv", index=False)
    if denominator_mismatches:
        raise ValueError("Eligible-target denominator does not reproduce the saved normalization")
    old_dist = saved["Averag Drug–Disease Distance"].to_numpy(dtype=float)
    available = saved.test_status.eq("evaluable").to_numpy()
    mismatch = available & ~np.isclose(observed, old_dist, atol=1e-7, rtol=0)
    pd.DataFrame(
        {
            "drug": names,
            "distance_saved": old_dist,
            "distance_repository": observed,
            "saved_evaluable": available,
            "distance_mismatch": mismatch,
        }
    ).to_csv(output / "input_reproduction_audit.csv", index=False)
    W, nodes, index, local = reference_weighted_matrix(graph, ppi)
    u = np.zeros(len(nodes))
    for gene in disease:
        u[index[gene]] = 1 + float(local.get(gene, 0)) / 1000
    stamp = time.perf_counter()
    path_scores = weighted_simple_path_scores(W, u, max_length=max(args.path_lengths))
    print(f"All-node simple-path scores computed in {time.perf_counter() - stamp:.2f}s", flush=True)
    score_rows = []
    for L in args.path_lengths:
        for i, name in enumerate(names):
            active = [j for j in np.flatnonzero(targets[i] > 0) if int(targets[i, j]) in index]
            gene_score = np.array([path_scores[L][index[int(targets[i, j])]] for j in active])
            contributing = gene_score > 0
            raw = float(np.dot(gene_score, dgi[i, active]))
            count = len(active)
            score_rows.append(
                {
                    "drug": name,
                    "path_length_edges": L,
                    "raw_score": raw,
                    "normalization_targets": count,
                    "path_connected_targets": int(contributing.sum()),
                    "normalized_score": raw / count if count else 0,
                }
            )
    scores = pd.DataFrame(score_rows)
    scores.to_csv(output / "repository_path_scores_all_drugs.csv", index=False)
    # Path-length sensitivity holds the historical 36-drug list fixed.
    focus = scores[scores.drug.isin(ranking.drug)].copy()
    focus["rank"] = focus.groupby("path_length_edges").normalized_score.rank(
        method="average", ascending=False
    )
    focus.to_csv(output / "repository_path_scores_fixed_36.csv", index=False)
    baseline = focus[focus.path_length_edges == 2].set_index("drug")
    path_summary = []
    for L in args.path_lengths:
        cur = focus[focus.path_length_edges == L].set_index("drug").reindex(baseline.index)
        changes = (cur["rank"] - baseline["rank"]).abs()
        top0 = set(
            baseline.sort_values(["normalized_score", "drug"], ascending=[False, True])
            .head(10)
            .index
        )
        top = set(
            cur.sort_values(["normalized_score", "drug"], ascending=[False, True]).head(10).index
        )
        path_summary.append(
            {
                "path_length_edges": L,
                "n_drugs": len(cur),
                "spearman_vs_L2": float(spearmanr(baseline["rank"], cur["rank"]).statistic),
                "top10_overlap_vs_L2": len(top & top0),
                "median_abs_rank_change": float(changes.median()),
                "maximum_abs_rank_change": float(changes.max()),
                "positive_scores": int((cur.normalized_score > 0).sum()),
            }
        )
    pd.DataFrame(path_summary).to_csv(output / "path_sensitivity_summary.csv", index=False)
    old = ranking.set_index("drug").reindex(baseline.index)
    old_agreement = {
        str(L): float(
            spearmanr(
                old.normalized_rank_saved,
                focus[focus.path_length_edges == L].set_index("drug").reindex(old.index)["rank"],
            ).statistic
        )
        for L in args.path_lengths
    }
    bins = core.get_degree_binning(graph, 100)
    cached = {}
    null_rows = []
    for i, profile in enumerate(profiles):
        if profile not in cached:
            if not profile:
                cached[profile] = [
                    {
                        "penalty": C,
                        "observed_distance": np.nan,
                        "null_mean": np.nan,
                        "null_sd": np.nan,
                        "z": np.nan,
                        "tail_count": np.nan,
                        "p_MC": np.nan,
                    }
                    for C in args.penalties
                ]
                for stats in cached[profile]:
                    null_rows.append({"drug": names[i], **stats})
                continue
            a, b = reference_null_distances(graph, bins, profile, nearest, args.samples, args.seed)
            observed_dist = np.array([nearest.get(v, -1) for v in profile])
            cached[profile] = []
            finite_sum = np.rint(a * len(profile)).astype(np.int64)
            missing_count = np.rint(b * len(profile)).astype(np.int64)
            for C in args.penalties:
                null_sums = finite_sum + missing_count * C
                observed_sum = (
                    int(np.maximum(observed_dist, 0).sum()) + int((observed_dist < 0).sum()) * C
                )
                values = null_sums / len(profile)
                d = observed_sum / len(profile)
                sd = float(values.std())
                mean = float(values.mean())
                # Integer sums make tied distances exact at every penalty.
                tail = int(np.count_nonzero(null_sums <= observed_sum))
                cached[profile].append(
                    {
                        "penalty": C,
                        "observed_distance": d,
                        "null_mean": mean,
                        "null_sd": sd,
                        "z": (d - mean) / sd if sd > 0 else 0.0,
                        "tail_count": tail,
                        "p_MC": (tail + 1) / (args.samples + 1),
                    }
                )
        for stats in cached[profile]:
            null_rows.append({"drug": names[i], **stats})
        if i % 100 == 0:
            print(
                f"Null sets {i + 1}/{len(profiles)}, unique profiles {len(cached)}, elapsed {time.perf_counter() - start:.1f}s",
                flush=True,
            )
        if args.max_profiles and len(cached) >= args.max_profiles:
            break
    null = pd.DataFrame(null_rows)
    for C, idx in null.groupby("penalty").groups.items():
        null.loc[idx, "p_BH"], null.loc[idx, "p_BY"] = adjusted_pvalues(null.loc[idx, "p_MC"])
    null["nominal_exploratory"] = (null.z <= -1.96) & (null.p_MC <= 0.05)
    null.to_csv(output / "repository_penalty_results.csv", index=False)
    ref = null[null.penalty == 20000]
    pen_summary = []
    for C, cur in null.groupby("penalty"):
        stats = {"penalty": C, "evaluable": int(cur.p_MC.notna().sum())}
        for key, mask in [
            ("nominal_exploratory", cur.nominal_exploratory),
            ("BH", cur.p_BH <= 0.05),
            ("BY", cur.p_BY <= 0.05),
        ]:
            sel = set(cur.loc[mask, "drug"])
            refmask = (
                ref.nominal_exploratory if key == "nominal_exploratory" else ref["p_" + key] <= 0.05
            )
            sel0 = set(ref.loc[refmask, "drug"])
            stats[key + "_count"] = len(sel)
            stats[key + "_jaccard_vs_20000"] = jaccard(sel, sel0)
        pen_summary.append(stats)
    pd.DataFrame(pen_summary).to_csv(output / "penalty_sensitivity_summary.csv", index=False)
    files = [p for p in (root / "data").iterdir() if p.is_file()] + [workbook]
    metadata = {
        "source_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True
        ).strip(),
        "samples": args.samples,
        "seed": args.seed,
        "rng": "numpy.random.Generator PCG64",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "networkx": nx.__version__,
        "penalties": args.penalties,
        "path_lengths_edges": args.path_lengths,
        "baseline_path_edges": 2,
        "graph_nodes": len(graph),
        "graph_edges": graph.number_of_edges(),
        "self_loops": nx.number_of_selfloops(graph),
        "repository_disease_genes": len(disease),
        "disease_file": str(disease_path.relative_to(root))
        if disease_path.is_relative_to(root)
        else str(disease_path),
        "protein_universe": "provided protein-coding list union PPI IDs; retain eligible isolates",
        "normalization_denominator": "all eligible target entries, fixed across path limits",
        "saved_denominator_mismatches": denominator_mismatches,
        "supplement_disease_symbols": 54,
        "max_finite_nearest_distance": max_finite,
        "unreachable_nodes": len(graph) - len(nearest),
        "observed_distance_mismatches": int(mismatch.sum()),
        "unavailable_saved_drugs": saved.loc[~available, "drug"].tolist(),
        "old_rank_vs_repository_rank_spearman": old_agreement,
        "saved_summary": saved_summary,
        "penalty_summary": pen_summary,
        "path_summary": path_summary,
        "complete_drug_universe": not bool(args.max_profiles),
        "unique_target_profiles": sum(bool(profile) for profile in cached),
        "elapsed_seconds": time.perf_counter() - start,
        "input_sha256": {
            str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files
        },
    }
    (output / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="experiments/psc_revision")
    parser.add_argument("--disease-file", default="data/psc_manuscript_targets.csv")
    parser.add_argument("--samples", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=452456)
    parser.add_argument("--penalties", type=int, nargs="+", default=[100, 1000, 20000])
    parser.add_argument("--path-lengths", type=int, nargs="+", default=[1, 2, 3, 4])
    parser.add_argument(
        "--max-profiles",
        type=int,
        default=None,
        help="Pilot only: never use partial FDR results in the manuscript",
    )
    args = parser.parse_args()
    if args.samples < 1:
        parser.error("--samples must be positive")
    run(args)


if __name__ == "__main__":
    main()
