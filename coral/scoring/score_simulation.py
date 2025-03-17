"""
A script for computing the performance of a reconstruction algorithm on simulated
amplicon structures.
"""

from __future__ import annotations

import concurrent.futures
import dataclasses
import io
import os
import pathlib
import re
import sys
import warnings
from itertools import combinations, product
from typing import Optional

import pandera as pa
import pandera.typing as pat

from coral.breakpoint import parse_graph
from coral.breakpoint.parse_graph import parse_breakpoint_graph
from coral.scoring import scoring_types
from coral.scoring.scoring_types import (
    AmpliconReconstructionStats,
    TrueAmpliconStats,
)

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
    "has_cycle_match",
]


def score_reconstruction(
    ground_truth_dir: pathlib.Path,
    reconstruction_dir: pathlib.Path,
    cnv_seeds_file: pathlib.Path,
    amplicon: str,
    tolerance: int = 100,
    decoil: bool = False,
) -> TrueAmpliconStats:
    # define statistics

    # compute overlapping sequence & discordant edges
    with (ground_truth_dir / f"{amplicon}_graph.txt").open("r") as f:
        true_graph = parse_graph.parse_breakpoint_graph(f)
    cycle_filepath = ground_truth_dir / f"{amplicon}_cycles.txt"
    true_intervals = io_utils.read_cycles_intervals_to_bed(cycle_filepath)
    true_cycles_bed = io_utils.read_cycles_file_to_bed(cycle_filepath)
    cnv_seeds = pyranges.read_bed(str(cnv_seeds_file))
    cnv_seeds_relaxed = cnv_seeds.extend(int(tolerance))

    # compute overlapping CNV intervals
    overlapping_cnv_seeds = true_intervals.intersect(cnv_seeds_relaxed)

    recon_stats = scoring_types.TrueAmpliconStats(
        n_amplified_intervals=len(true_intervals),
        n_amplified_intervals_covered=len(overlapping_cnv_seeds),
        n_sequence_edges=len(true_graph.sequence_edges),
        n_breakpoint_edges=len(true_graph.discordant_edges),
    )

    per_amplicon_stats: dict[str, AmpliconReconstructionStats] = {}

    for bp_graph_path in reconstruction_dir.glob("*_graph.txt"):
        with bp_graph_path.open("r") as f:
            reconstructed_graph = parse_graph.parse_breakpoint_graph(f)

        num_overlapping_seq_edges = (
            scoring_utils.find_overlapping_sequence_edges(
                true_graph,
                reconstructed_graph,
                tolerance=tolerance,
            )
        )

        num_overlapping_bp_edges = scoring_utils.find_overlapping_bp_edges(
            true_graph,
            reconstructed_graph,
            tolerance=tolerance,
        )

        per_amplicon_stats[bp_graph_path.name] = AmpliconReconstructionStats(
            num_overlapping_seq_edges,
            num_overlapping_bp_edges,
            len(reconstructed_graph.sequence_edges),
            len(reconstructed_graph.breakpoint_edges),
        )

        reconstructed_cycles_file_name = re.sub(
            r"amplicon(\d+)_graph.txt",
            r"amplicon\1_cycles.txt",
            bp_graph_path.name,
        )
        reconstructed_cycles_filepath = pathlib.Path(
            reconstruction_dir / f"{reconstructed_cycles_file_name}"
        )
        if reconstructed_cycles_filepath.exists():
            reconstructed_cycles_bed = io_utils.read_cycles_file_to_bed(
                reconstructed_cycles_filepath
            )
            per_amplicon_stats[bp_graph_path.name].has_cycle_match = (
                len(reconstructed_cycles_bed.df) > 0
            )

            # compute fragment overlap
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

            if _fragment_overlap > recon_stats.fragment_overlap:
                recon_stats.fragment_overlap = _fragment_overlap
                recon_stats.reconstruction_length_ratio = (
                    _reconstruction_length_ratio
                )

                recon_stats.best_copy_number_ratio = _best_copy_number_ratio
                recon_stats.top_three_copy_number_ratio = (
                    _top_three_copy_number_ratio
                )

            if len(true_cycles_bed.df) < 3:
                recon_stats.cycle_triplets_correct = np.nan
            else:
                _cycle_triplets_correct = scoring_utils.score_triplets_correct(
                    true_cycles_bed.df, reconstructed_cycles_bed.df
                )

                if (
                    not np.isnan(_cycle_triplets_correct)
                    and np.isnan(recon_stats.cycle_triplets_correct)
                ) or (
                    _cycle_triplets_correct > recon_stats.cycle_triplets_correct
                ):
                    recon_stats.cycle_triplets_correct = _cycle_triplets_correct

            (
                max_lcs,
                max_normalized_lcs,
            ) = scoring_utils.compute_normalized_longest_cycle_subsequence(
                true_cycles_bed.df, reconstructed_cycles_bed.df, binned_genome
            )
            if max_lcs > recon_stats.overall_max_lcs:
                recon_stats.overall_max_lcs = max_lcs
                recon_stats.overall_max_normalized_lcs = max_normalized_lcs

    for name, amplicon_stats in per_amplicon_stats.items():
        if (
            amplicon_stats.n_overlapping_bp_edges
            > recon_stats.n_breakpoint_edges_covered
        ):
            recon_stats.n_breakpoint_edges_covered = (
                amplicon_stats.n_overlapping_bp_edges
            )
            recon_stats.n_sequence_edges_covered = (
                amplicon_stats.n_overlapping_seq_edges
            )
            recon_stats.n_reconstructed_sequence_edges_total = (
                amplicon_stats.n_reconstructed_sequence_edges_total
            )
            recon_stats.n_reconstructed_breakpoint_edges_total = (
                amplicon_stats.n_reconstructed_breakpoint_edges_total
            )
            recon_stats.has_cycle_match = amplicon_stats.has_cycle_match
            recon_stats.matched_amplicon = name
    return recon_stats


def score_simulations(
    ground_truth: pathlib.Path,
    reconstruction_dir: pathlib.Path,
    output_dir: pathlib.Path,
    tolerance: int,
    to_skip: list[str],
    num_threads: int = 1,
) -> None:
    amplicon_statistics: pat.DataFrame[
        scoring_types.ReconstructionScoreModel
    ] = pd.DataFrame(
        columns=scoring_types.ReconstructionScoreSchema.columns.keys()
    )

    to_skip: list[str] = []

    simulation_directory = str(ground_truth)
    solver, date, threads = "scip", "20241211", "1"
    # solver_path = f"{solver}/{threads}/{date}"
    solver_path = f"{solver}/{threads}"  # /{date}"
    dataset_paths = list(ground_truth.iterdir())

    for dataset_path in tqdm.tqdm(
        dataset_paths, desc="Iterating through test datasets"
    ):
        if dataset_path.name in to_skip:
            print(f"Skipping {dataset_path.name}!")
            continue
        simulated_amplicons: np.ndarray[str, np.dtype[np.str_]] = np.unique(
            [
                amplicon.split("_graph")[0]
                for amplicon in os.listdir(
                    ground_truth / dataset_path / "ref" / "graphs"
                )
            ]
        )
        raw_df = pd.read_csv(
            ground_truth / dataset_path / "simulated_amplicons.txt",
            sep="\t",
            names=[
                "type",
                "breakpoints",
                "coverage",
                "replicate",
            ],
        )
        amplicon_summaries: pat.DataFrame[
            scoring_types.SimulatedAmpliconSchema
        ] = scoring_types.SimulatedAmpliconSchema.validate(raw_df)
        amplicon_summaries.set_index(
            amplicon_summaries.apply(
                lambda x: f"{x.type}_bp{x.breakpoints}_amplicon{x.replicate}",
                axis=1,
            ),
            inplace=True,
        )

        for amplicon in tqdm.tqdm(
            simulated_amplicons,
            desc=f"Iterating through simulated amplicons for {dataset_path.name}",
        ):
            coverage = amplicon_summaries.loc[amplicon, "coverage"]

            (
                _,
                simulated_bp_edges,
            ) = io_utils.read_breakpoint_graph(
                ground_truth / "ref" / dataset_path / f"{amplicon}_graph.txt"
            )
            n_breakpoints = len(simulated_bp_edges)
            recon_stats = score_reconstruction(
                ground_truth_dir=dataset_path / "ref",
                reconstruction_dir=reconstruction_dir / dataset_path.name,
                cnv_seeds_file=simulation_directory
                / dataset_path
                / "CNV_SEEDS.bed",
                amplicon=amplicon,
                tolerance=tolerance,
            )
            results_df = pd.DataFrame(
                [
                    [
                        dataset_path.name,
                        "CoRAL",
                        amplicon,
                        coverage,
                        n_breakpoints,
                    ]
                    + list(dataclasses.astuple(recon_stats))
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
        lambda x: x["n_breakpoint_edges_covered"] / x["n_breakpoint_edges"],
        axis=1,
    )
    amplicon_statistics["breakpoint_class"] = amplicon_statistics.apply(
        lambda x: x.amplicon.split("_")[0], axis=1
    )

    amplicon_statistics.to_csv(f"{output_dir}/amplicon_scores.tsv", sep="\t")
