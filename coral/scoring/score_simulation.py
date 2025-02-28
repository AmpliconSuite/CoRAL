"""
A script for computing the performance of a reconstruction algorithm on simulated
amplicon structures.
"""

from __future__ import annotations

import io
import os
import pathlib
import re
import sys
import warnings
from itertools import combinations, product
from typing import Optional

from coral.breakpoint import parse_graph
from coral.breakpoint.parse_graph import parse_breakpoint_graph

warnings.simplefilter(action="ignore", category=FutureWarning)

import networkx as nx
import numpy as np
import pandas as pd
import pyranges
import scipy
import seaborn as sns
import tqdm
import typer

from coral.scoring import io_utils, scoring_utils

AMPLICON_COLS = [
    "test_id",
    "algorithm",
    "amplicon",
    "coverage",
    "n_breakpoints_simulated",
    "n_amplified_intervals",
    "n_sequence_edges",
    "n_breakpoint_edges",
    "amplified_intervals_covered",
    "sequence_edges_covered",
    "breakpoint_edges_covered",
    "total_reconstructed_sequence_edges",
    "total_reconstructed_breakpoint_edges",
    "fragment_overlap_unweighted",
    "reconstruction_length_ratio",
    "cycle_triplets_correct",
    "best_copy_number_ratio",
    "top_three_copy_number_ratio",
    "maximum_lcs_length",
    "maximum_normalized_lcs_length",
]


def score_reconstruction(
    true_cycles_file: io.TextIOWrapper,
    true_graph_file: io.TextIOWrapper,
    bp_graph_dir: pathlib.Path,
    cnv_seeds_file: io.TextIOWrapper,
    tolerance: int = 100,
    decoil: bool = False,
) -> None:
    # define statistics
    n_amplified_intervals = 0
    n_sequence_edges = 0
    n_breakpoint_edges = 0
    n_amplified_intervals_covered = 0
    n_sequence_edges_covered = 0
    n_breakpoint_edges_covered = 0
    n_reconstructed_sequence_edges_total = 0
    n_reconstructed_breakpoint_edges_total = 0
    fragment_overlap = 0
    reconstruction_length_ratio = 0
    cycle_triplets_correct = np.nan
    best_copy_number_ratio = 0
    top_three_copy_number_ratio = 0
    overall_max_lcs = 0
    overall_max_normalized_lcs = 0

    true_intervals, _ = io_utils.read_cycles_intervals_to_bed(true_cycles_file)
    true_cycles_bed = io_utils.read_cycles_file_to_bed(true_cycles_file)
    cnv_seeds = pyranges.read_bed(cnv_seeds_file.name)
    cnv_seeds_relaxed = cnv_seeds.extend(int(tolerance))

    # compute overlapping CNV intervals
    overlapping_cnv_seeds = true_intervals.intersect(cnv_seeds_relaxed)
    n_amplified_intervals = len(true_intervals)
    n_amplified_intervals_covered = len(overlapping_cnv_seeds)

    # compute overlapping sequence & discordant edges
    (
        simulated_seq_edges,
        simulated_bp_edges,
    ) = io_utils.read_breakpoint_graph(true_graph_file.name)
    true_graph = parse_graph.parse_breakpoint_graph(true_graph_file)

    reconstruction_to_statistics: dict[str, tuple[int, int, int, int]] = {}
    for bp_filepath in bp_graph_dir.glob("*_graph.txt"):
        with bp_filepath.open("r") as f:
            reconstructed_graph = parse_graph.parse_breakpoint_graph(f)

        (num_overlapping_sequence, _) = (
            scoring_utils.find_overlapping_sequence_edges(
                true_graph,
                reconstructed_graph,
                tolerance=tolerance,
            )
        )

        (num_overlapping_bp, _) = scoring_utils.find_overlapping_bp_edges(
            true_graph,
            reconstructed_graph,
            tolerance=tolerance,
        )

        reconstruction_to_statistics[bp_filepath.name] = (
            num_overlapping_sequence,
            num_overlapping_bp,
            len(reconstructed_graph.sequence_edges),
            len(reconstructed_graph.breakpoint_edges),
        )

        # compute fragment overlap
        reconstructed_cycles_file_name = re.sub(
            r"graphs/amplicon(\d+)_graph.txt",
            r"cycles/amplicon\1_cycles.txt",
            bp_filepath.name,
        )
        reconstructed_cycles_filepath = pathlib.Path(
            reconstructed_cycles_file_name
        )
        if reconstructed_cycles_filepath.exists():
            with reconstructed_cycles_filepath.open("r") as f:
                reconstructed_cycles_bed = io_utils.read_cycles_file_to_bed(f)

            binned_genome, _ = io_utils.bin_genome(
                true_cycles_bed.df, reconstructed_cycles_bed.df
            )
            (
                _fragment_overlap,
                _reconstruction_length_ratio,
            ) = scoring_utils.get_fragments_similarity_unweighted(
                true_cycles_bed,
                reconstructed_cycles_bed,
                pyranges.PyRanges(binned_genome),
            )

            _best_copy_number_ratio = scoring_utils.get_cycle_copy_number_ratio(
                reconstructed_cycles_bed.df,
                reconstructed_graph.sequence_edges,
                n=1,
            )
            _top_three_copy_number_ratio = (
                scoring_utils.get_cycle_copy_number_ratio(
                    reconstructed_cycles_bed.df,
                    reconstructed_graph.sequence_edges,
                    n=3,
                )
            )

            if _fragment_overlap > fragment_overlap:
                fragment_overlap = _fragment_overlap
                reconstruction_length_ratio = _reconstruction_length_ratio

                best_copy_number_ratio = _best_copy_number_ratio
                top_three_copy_number_ratio = _top_three_copy_number_ratio

            if len(true_cycles_bed.df) < 3:
                cycle_triplets_correct = np.nan
            else:
                _cycle_triplets_correct = scoring_utils.score_triplets_correct(
                    true_cycles_bed.df, reconstructed_cycles_bed.df
                )

                if not np.isnan(_cycle_triplets_correct) and np.isnan(
                    cycle_triplets_correct
                ):
                    cycle_triplets_correct = _cycle_triplets_correct

                elif _cycle_triplets_correct > cycle_triplets_correct:
                    cycle_triplets_correct = _cycle_triplets_correct

            (
                max_lcs,
                max_normalized_lcs,
            ) = scoring_utils.compute_normalized_longest_cycle_subsequence(
                true_cycles_bed.df, reconstructed_cycles_bed.df, binned_genome
            )
            if max_lcs > overall_max_lcs:
                overall_max_lcs, overall_max_normalized_lcs = (
                    max_lcs,
                    max_normalized_lcs,
                )

    for reconstruction_name in reconstruction_to_statistics:
        if (
            reconstruction_to_statistics[reconstruction_name][1]
            > n_breakpoint_edges_covered
        ):
            (
                n_sequence_edges_covered,
                n_breakpoint_edges_covered,
                n_reconstructed_sequence_edges_total,
                n_reconstructed_breakpoint_edges_total,
            ) = reconstruction_to_statistics[reconstruction_name]

    return (
        n_amplified_intervals,
        n_sequence_edges,
        n_breakpoint_edges,
        n_amplified_intervals_covered,
        n_sequence_edges_covered,
        n_breakpoint_edges_covered,
        n_reconstructed_sequence_edges_total,
        n_reconstructed_breakpoint_edges_total,
        fragment_overlap,
        reconstruction_length_ratio,
        cycle_triplets_correct,
        best_copy_number_ratio,
        top_three_copy_number_ratio,
        overall_max_lcs,
        overall_max_normalized_lcs,
    )


def score_simulations(
    simulation_dir_path: pathlib.Path,
    output_dir: pathlib.Path,
    tolerance: int,
) -> None:
    amplicon_statistics = pd.DataFrame(columns=AMPLICON_COLS)

    to_skip: list[str] = []

    simulation_directory = str(simulation_dir_path)
    solver, date, threads = "scip", "20241211", "1"
    # solver_path = f"{solver}/{threads}/{date}"
    solver_path = f"{solver}/{threads}"  # /{date}"
    datasets = [x for x in os.listdir(simulation_directory)]
    for dataset in tqdm.tqdm(datasets, desc="Iterating through test datasets"):
        if dataset in to_skip:
            print(f"Skipping {dataset}!")
            continue

        simulated_amplicons = np.unique(
            [
                amplicon.split("_graph")[0]
                for amplicon in os.listdir(
                    f"{simulation_directory}/{dataset}/ref/graphs"
                )
            ]
        )

        amplicon_summaries = pd.read_csv(
            f"{simulation_directory}/{dataset}/simulated_amplicons.txt",
            sep="\t",
            header=None,
        )
        amplicon_summaries.columns = [
            "type",
            "breakpoints",
            "coverage",
            "replicate",
        ]
        amplicon_summaries.index = amplicon_summaries.apply(
            lambda x: f"{x.type}_{x.breakpoints}_{x.replicate}", axis=1
        )

        reconstructed_amplicons = np.unique(
            [
                amplicon.split("_graph")[0]
                for amplicon in os.listdir(
                    f"{simulation_directory}/{dataset}/{solver_path}/graphs"
                )
            ]
        )
        for amplicon in tqdm.tqdm(
            simulated_amplicons, desc="Iterating through simulated amplicons"
        ):
            coverage = amplicon_summaries.loc[amplicon, "coverage"]

            (
                _,
                simulated_bp_edges,
            ) = io_utils.read_breakpoint_graph(
                f"{simulation_directory}/{dataset}/ref/graphs/{amplicon}_graph.txt"
            )
            n_breakpoints = len(simulated_bp_edges)

            results = score_reconstruction(
                f"{simulation_directory}/{dataset}/ref/cycles/{amplicon}_cycles.txt",
                f"{simulation_directory}/{dataset}/ref/graphs/{amplicon}_graph.txt",
                [
                    f"{simulation_directory}/{dataset}/{solver_path}/graphs/{reconstructed_amplicon}_graph.txt"
                    for reconstructed_amplicon in reconstructed_amplicons
                ],
                f"{simulation_directory}/{dataset}/CNV_SEEDS.bed",
                tolerance,
            )
            results_df = pd.DataFrame(
                [
                    [dataset, "CoRAL", amplicon, coverage, n_breakpoints]
                    + list(results)
                ],
                columns=amplicon_statistics.columns,
            )

            amplicon_statistics = pd.concat(
                [
                    amplicon_statistics,
                    results_df,
                ]
            )

    amplicon_statistics.index = range(len(amplicon_statistics))
    amplicon_statistics["breakpoint_accuracy"] = amplicon_statistics.apply(
        lambda x: x["breakpoint_edges_covered"] / x["n_breakpoint_edges"],
        axis=1,
    )
    amplicon_statistics["breakpoint_class"] = amplicon_statistics.apply(
        lambda x: x.amplicon.split("_")[0], axis=1
    )

    amplicon_statistics.to_csv(f"{output_dir}/amplicon_stats.csv", sep="\t")
