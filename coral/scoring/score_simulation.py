"""
A script for computing the performance of a reconstruction algorithm on simulated
amplicon structures.
"""
import os
import sys

from itertools import combinations, product
from typing import Optional
import warnings

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


app = typer.Typer()


def score_reconstruction(
    true_cycles_file,
    true_graph_file,
    reconstructed_graph_files,
    cnv_seeds_file,
    tolerance=100,
    decoil=False,
):
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
    
    true_intervals, _ = io_utils.read_cycles_intervals_to_bed(
        true_cycles_file
    )
    true_cycles_bed = io_utils.read_cycles_file_to_bed(true_cycles_file)
    cnv_seeds = pyranges.read_bed(cnv_seeds_file)
    cnv_seeds_relaxed = cnv_seeds.extend(int(tolerance))

    # compute overlapping CNV intervals
    overlapping_cnv_seeds = true_intervals.intersect(cnv_seeds_relaxed)
    n_amplified_intervals = len(true_intervals)
    n_amplified_intervals_covered = len(overlapping_cnv_seeds)

    # compute overlapping sequence & discordant edges
    (
        simulated_seq_edges,
        simulated_bp_edges,
    ) = io_utils.read_breakpoint_graph(true_graph_file)

    n_sequence_edges = len(simulated_seq_edges)
    n_breakpoint_edges = len(simulated_bp_edges)

    reconstruction_to_statistics = {}
    for reconstructed_graph_file in reconstructed_graph_files:

        (
            reconstructed_seq_edges,
            reconstructed_bp_edges,
        ) = io_utils.read_breakpoint_graph(reconstructed_graph_file)

        (
            num_overlapping_sequence,
            _,
        ) = scoring_utils.find_overlapping_sequence_edges(
            simulated_seq_edges,
            reconstructed_seq_edges,
            tolerance=tolerance,
        )

        (num_overlapping_bp, _) = scoring_utils.find_overlapping_bp_edges(
            simulated_bp_edges,
            reconstructed_bp_edges,
            tolerance=tolerance,
        )

        reconstruction_to_statistics[reconstructed_graph_file] = (
            num_overlapping_sequence,
            num_overlapping_bp,
            len(reconstructed_seq_edges),
            len(reconstructed_bp_edges),
        )

        # compute fragment overlap
        reconstructed_cycles_file = reconstructed_graph_file.replace(
            "graph.txt", "cycles.txt"
        )
        if os.path.exists(reconstructed_cycles_file):
            reconstructed_cycles_bed = io_utils.read_cycles_file_to_bed(
                reconstructed_cycles_file
            )

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

            _best_copy_number_ratio = (
                scoring_utils.get_cycle_copy_number_ratio(
                    reconstructed_cycles_bed.df, reconstructed_seq_edges, n=1
                )
            )
            _top_three_copy_number_ratio = (
                scoring_utils.get_cycle_copy_number_ratio(
                    reconstructed_cycles_bed.df, reconstructed_seq_edges, n=3
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
                _cycle_triplets_correct = (
                    scoring_utils.score_triplets_correct(
                        true_cycles_bed.df, reconstructed_cycles_bed.df
                    )
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

    for reconstruction in reconstruction_to_statistics.keys():

        if (
            reconstruction_to_statistics[reconstruction][1]
            > n_breakpoint_edges_covered
        ):
            (
                n_sequence_edges_covered,
                n_breakpoint_edges_covered,
                n_reconstructed_sequence_edges_total,
                n_reconstructed_breakpoint_edges_total,
            ) = reconstruction_to_statistics[reconstruction]

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

def score_decoil_reconstruction(
    true_cycles_file,
    true_graph_file,
    reconstruction_bed_file,
    tolerance=100,
    base_coverage=13.0,
):
    # define statistics
    n_amplified_intervals = np.nan
    n_sequence_edges = 0
    n_breakpoint_edges = 0
    n_amplified_intervals_covered = np.nan
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
    
    true_intervals, _ = io_utils.read_cycles_intervals_to_bed(
        true_cycles_file
    )
    true_cycles_bed = io_utils.read_cycles_file_to_bed(true_cycles_file)

    # compute overlapping sequence & discordant edges
    (
        simulated_seq_edges,
        simulated_bp_edges,
    ) = io_utils.read_breakpoint_graph(true_graph_file)

    n_sequence_edges = len(simulated_seq_edges)
    n_breakpoint_edges = len(simulated_bp_edges)

    # read in decoil file
    decoil_reconstruction = pd.read_csv(reconstruction_bed_file, sep='\t')
    decoil_reconstruction = decoil_reconstruction[['#chr', 'start', 'end', 'strand', 'circ_id', 'fragment_id', "coverage", "estimated_proportions"]]
    decoil_reconstruction.rename(columns={'#chr': "Chromosome", 'start': 'Start', 'end': 'End', 'strand': 'orientation', 'circ_id': 'cycle_id', 'fragment_id': 'segment_label'}, inplace=True)
    decoil_reconstruction['iscyclic'] = True
    total_proportion = decoil_reconstruction["estimated_proportions"].sum()
    decoil_reconstruction['weight'] = decoil_reconstruction.apply(lambda x: x.estimated_proportions / total_proportion, axis=1)
    # decoil_reconstruction['weight'] = decoil_reconstruction['coverage']/base_coverage

    reconstruction_to_statistics = {}

    (
        reconstructed_bp_edges,
        reconstructed_seq_edges
    ) = io_utils.get_seq_edges_from_decoil(decoil_reconstruction, base_coverage=base_coverage)

    (   
    num_overlapping_sequence,
        _,
    ) = scoring_utils.find_overlapping_sequence_edges(
        simulated_seq_edges,
        reconstructed_seq_edges,
        tolerance=tolerance,
    )

    (num_overlapping_bp, _) = scoring_utils.find_overlapping_bp_edges(
        simulated_bp_edges,
        reconstructed_bp_edges,
        tolerance=tolerance,
    )

    (
        n_sequence_edges_covered,
        n_breakpoint_edges_covered,
        n_reconstructed_sequence_edges_total,
        n_reconstructed_breakpoint_edges_total,
    ) = (num_overlapping_sequence, num_overlapping_bp, len(reconstructed_seq_edges), len(reconstructed_bp_edges))

    for ecdna, reconstructed_cycles_bed in decoil_reconstruction.groupby('cycle_id'):

        # reconstruction_to_statistics[ecdna] = (
        #     num_overlapping_sequence,
        #     num_overlapping_bp,
        #     len(reconstructed_seq_edges),
        #     len(reconstructed_bp_edges),
        # )

        # compute fragment overlap
        binned_genome, _ = io_utils.bin_genome(
            true_cycles_bed.df, reconstructed_cycles_bed
        )
        (
            _fragment_overlap,
            _reconstruction_length_ratio,
        ) = scoring_utils.get_fragments_similarity_unweighted(
            true_cycles_bed,
            pyranges.PyRanges(reconstructed_cycles_bed),
            pyranges.PyRanges(binned_genome),
        )

        _best_copy_number_ratio = np.nan
        _top_three_copy_number_ratio = np.nan

        # _best_copy_number_ratio = (
        #     scoring_utils.get_cycle_copy_number_ratio(
        #         reconstructed_cycles_bed, reconstructed_seq_edges, n=1
        #     )
        # )
        # _top_three_copy_number_ratio = (
        #     scoring_utils.get_cycle_copy_number_ratio(
        #         reconstructed_cycles_bed, reconstructed_seq_edges, n=3
        #     )
        # )

        if _fragment_overlap > fragment_overlap:
            fragment_overlap = _fragment_overlap
            reconstruction_length_ratio = _reconstruction_length_ratio

            best_copy_number_ratio = _best_copy_number_ratio
            top_three_copy_number_ratio = _top_three_copy_number_ratio

        if len(true_cycles_bed.df) < 3:
            cycle_triplets_correct = np.nan
        else:
            _cycle_triplets_correct = (
                scoring_utils.score_triplets_correct(
                    true_cycles_bed.df, reconstructed_cycles_bed
                )
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
            true_cycles_bed.df, reconstructed_cycles_bed, binned_genome
        )
        if max_lcs > overall_max_lcs:
            overall_max_lcs, overall_max_normalized_lcs = (
                max_lcs,
                max_normalized_lcs,
            )

    # for reconstruction in reconstruction_to_statistics.keys():

    #     if (
    #         reconstruction_to_statistics[reconstruction][1]
    #         > n_breakpoint_edges_covered
    #     ):
    #         (
    #             n_sequence_edges_covered,
    #             n_breakpoint_edges_covered,
    #             n_reconstructed_sequence_edges_total,
    #             n_reconstructed_breakpoint_edges_total,
    #         ) = reconstruction_to_statistics[reconstruction]

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


@app.command()
def score_simulations(
    simulation_directory: str = typer.Argument(
        ..., help="Path to directory containing simulations."
    ),
    output_file: str = typer.Argument(..., help="Output directory path."),
    tolerance: int = typer.Option(
        100, help="Tolerance in basepairs for matching regions."
    ),
    base_coverage: float = 13.0,
    debug: Optional[str] = typer.Option(None, help="Debug specific test."),
):

    amplicon_statistics = pd.DataFrame(
        columns=[
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
    )

    to_skip = []

    datasets = [x for x in os.listdir(simulation_directory) if os.path.isdir(x)]
    for dataset in tqdm.tqdm(datasets, desc="Iterating through test datasets"):
        if dataset in to_skip:
            print(f"Skipping {dataset}!")
            continue

        if debug and dataset != debug:
            continue

        simulated_amplicons = np.unique(
            [
                amplicon.split("_graph")[0]
                for amplicon in os.listdir(
                    f"{simulation_directory}/{dataset}/ground_truth"
                )
                if "graph" in amplicon
            ]
        )

        amplicon_summaries = pd.read_csv(
            f"{simulation_directory}/{dataset}/ground_truth/simulated_amplicons.txt",
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

        nanopore_reconstructed_amplicons = np.unique(
            [
                amplicon.split("_graph")[0]
                for amplicon in os.listdir(
                    f"{simulation_directory}/{dataset}/reconstruction_020524/"
                )
                if "graph" in amplicon
            ]
        )
        illumina_reconstructed_amplicons = np.unique(
            [
                amplicon.split("_graph")[0]
                for amplicon in os.listdir(
                    f"{simulation_directory}/{dataset}/illumina/"
                )
                if "graph" in amplicon
            ]
        )
        hybrid_reconstructed_amplicons = np.unique(
            [
                amplicon.split("_graph")[0]
                for amplicon in os.listdir(
                    f"{simulation_directory}/{dataset}/hybrid_aa/"
                )
                if "graph" in amplicon
            ]
        )
        for amplicon in tqdm.tqdm(
            simulated_amplicons, desc="Iterating through simulated amplicons"
        ):

            coverage = amplicon_summaries.loc[amplicon, "coverage"]

            (_, simulated_bp_edges,) = io_utils.read_breakpoint_graph(
                f"{simulation_directory}/{dataset}/ground_truth/{amplicon}_graph.txt"
            )
            n_breakpoints = len(simulated_bp_edges)

            nanopore_results = score_reconstruction(
                f"{simulation_directory}/{dataset}/ground_truth/{amplicon}_cycles.txt",
                f"{simulation_directory}/{dataset}/ground_truth/{amplicon}_graph.txt",
                [
                    f"{simulation_directory}/{dataset}/reconstruction_020524/{reconstructed_amplicon}_graph.txt"
                    for reconstructed_amplicon in nanopore_reconstructed_amplicons
                ],
                f"{simulation_directory}/{dataset}/reconstruction_020524/{dataset}_sorted_CNV_SEEDS.bed",
                tolerance,
            )
            nanopore_results_df = pd.DataFrame(
                [
                    [dataset, "CoRAL", amplicon, coverage, n_breakpoints]
                    + list(nanopore_results)
                ],
                columns=amplicon_statistics.columns,
            )

            illumina_results = score_reconstruction(
                f"{simulation_directory}/{dataset}/ground_truth/{amplicon}_cycles.txt",
                f"{simulation_directory}/{dataset}/ground_truth/{amplicon}_graph.txt",
                [
                    f"{simulation_directory}/{dataset}/illumina/{reconstructed_amplicon}_graph.txt"
                    for reconstructed_amplicon in illumina_reconstructed_amplicons
                ],
                f"{simulation_directory}/{dataset}/illumina/{dataset}_AA_CNV_SEEDS.bed",
                tolerance,
            )
            illumina_results_df = pd.DataFrame(
                [
                    [
                        dataset,
                        "AmpliconArchitect",
                        amplicon,
                        coverage,
                        n_breakpoints,
                    ]
                    + list(illumina_results)
                ],
                columns=amplicon_statistics.columns,
            )

            hybrid_results = score_reconstruction(
                f"{simulation_directory}/{dataset}/ground_truth/{amplicon}_cycles.txt",
                f"{simulation_directory}/{dataset}/ground_truth/{amplicon}_graph.txt",
                [
                    f"{simulation_directory}/{dataset}/hybrid_aa/{reconstructed_amplicon}_graph.txt"
                    for reconstructed_amplicon in hybrid_reconstructed_amplicons
                ],
                f"{simulation_directory}/{dataset}/illumina/{dataset}_AA_CNV_SEEDS.bed",
                tolerance,
            )
            hybrid_results_df = pd.DataFrame(
                [
                    [dataset, "CoRAL-Hybrid", amplicon, coverage, n_breakpoints]
                    + list(hybrid_results)
                ],
                columns=amplicon_statistics.columns,
            )

            decoil_results = score_decoil_reconstruction(
                f"{simulation_directory}/{dataset}/ground_truth/{amplicon}_cycles.txt",
                f"{simulation_directory}/{dataset}/ground_truth/{amplicon}_graph.txt",
                f"{simulation_directory}/{dataset}/decoil/reconstruct.bed",
                tolerance,
                base_coverage,
            )
            decoil_results_df = pd.DataFrame(
                [
                    [dataset, "Decoil", amplicon, coverage, n_breakpoints]
                    + list(decoil_results)
                ],
                columns=amplicon_statistics.columns,
            )

            amplicon_statistics = pd.concat(
                [
                    amplicon_statistics,
                    nanopore_results_df,
                    hybrid_results_df,
                    illumina_results_df,
                    decoil_results_df,
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

    amplicon_statistics.to_csv(output_file, sep="\t")


if __name__ == "__main__":
    app()
