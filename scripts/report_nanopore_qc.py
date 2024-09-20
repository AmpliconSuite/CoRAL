"""Script for reporting QC on Nanopore reads.
"""

import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns


def summarize_quality_control(fastq_file: str, output_dir: str, verbose: bool = False):
    """Reports QC from input BAM file.

    Args:
        fastq_file: input FASTQ file
        output_dir: where to output plots
        verbose: Report progress

    Returns:
        None. Produces plots in output_dir.

    """
    _iter = 0

    mean_qualities = []
    mean_lengths = []
    mapping_qualities = []

    with pysam.FastxFile(fastq_file) as fastq:
        for record in fastq:
            sequence, quality = (
                record.sequence,
                np.array([convert_character_to_quality(char) for char in record.quality]),
            )

            if sequence:
                mean_lengths.append(len(sequence))
                mean_qualities.append(np.mean(quality))

            _iter += 1
            if verbose and _iter % 1e6 == 0:
                print(f"Processed {_iter} records.")

    h = plt.figure(figsize=(10, 5))
    sns.histplot(mean_lengths)
    plt.xlabel("Mean Sequence Length")
    plt.ylabel("Frequency")
    plt.title(f"Mean Length of Nanopore Sequences (mean = {np.mean(mean_lengths)})")
    plt.savefig(f"{output_dir}/mean_length_histogram.png", dpi=300)
    plt.close()

    h = plt.figure(figsize=(10, 5))
    sns.histplot(mean_qualities)
    plt.xlabel("Mean Sequence Quality")
    plt.ylabel("Frequency")
    plt.title(f"Mean Quality of Nanopore Sequences (mean = {np.mean(mean_qualities)})")
    plt.savefig(f"{output_dir}/mean_sequence_quality_histogram.png", dpi=300)
    plt.close()

    summary_data_frame = pd.DataFrame(columns=["Q25", "Q50", "Q75"])
    summary_data_frame.loc["mean_length"] = [
        np.percentile(mean_lengths, 25),
        np.percentile(mean_lengths, 50),
        np.percentile(mean_lengths, 75),
    ]
    summary_data_frame.loc["mean_sequence_quality"] = [
        np.percentile(mean_qualities, 25),
        np.percentile(mean_qualities, 50),
        np.percentile(mean_qualities, 75),
    ]

    summary_data_frame.to_csv(f"{output_dir}/quality_control_summary.tsv", sep="\t")


def convert_character_to_quality(character: str) -> int:
    """Converts a ASCII character from FASTQ to Q score."""
    ascii = ord(character)
    return ascii - 33


def main():
    parser = argparse.ArgumentParser(description="Quality filter reads.")
    parser.add_argument("fastq_file", type=str, help="Path to FASTQ file.")
    parser.add_argument("output_directory", type=str, help="Where to write plots.")
    parser.add_argument("--verbose", action="store_true")

    args = parser.parse_args()

    fastq_file = args.fastq_file
    output_directory = args.output_directory
    verbose = args.verbose

    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    summarize_quality_control(fastq_file, output_directory, verbose=verbose)


if __name__ == "__main__":
    main()
