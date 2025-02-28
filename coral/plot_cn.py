import pathlib

import matplotlib
import matplotlib.axes
import matplotlib.pyplot as plt
import typer
from matplotlib import rcParams, ticker

from coral.cnv_seed import parse_centromere_arms

rcParams["pdf.fonttype"] = 42
font = {"family": "Arial", "size": 17}
plt.rc("font", **font)

CHR_NAME_TO_IDX = {"X": 23, "Y": 24}
CHR_IDX_TO_NAME = {23: "X", 24: "Y"}


def parse_cnr_file(file: typer.FileText) -> list:
    """
    Parses a CNR (Copy Number Ratio) file from CNVkit output.

    Parameters:
    filename (str): The path to the CNR file to be parsed.

    Returns:
    list of dict: A list where each item represents a line from the CNR file as
    a dictionary.
    """
    cnr_data = []
    # Skip the header line
    next(file)
    for line in file:
        # Split each line into its components
        parts = line.strip().split("\t")
        # Assuming the file structure is:
        # chromosome, start, end, gene (optional), log2, depth, probes
        chr_tag = parts[0].split("chr")[1]

        record = {
            "chromosome": int(CHR_NAME_TO_IDX.get(chr_tag, chr_tag)),
            "start": int(parts[1]),
            "end": int(parts[2]),
            "gene": parts[3],
            "log2": 2 * (2 ** float(parts[5])),
            "depth": float(parts[4]) if parts[4] else 0,
            "weight": float(parts[6]),
        }
        if file.name.endswith("cns"):
            record["log2"] = 2 * (2 ** float(parts[4]))
            record["depth"] = float(parts[5])
        cnr_data.append(record)
    return cnr_data


def plot_cnr_data(cnr_data: list[dict], output_dir: str, name: str) -> None:
    chr_arms = parse_centromere_arms()
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(20, 30))
    fig.suptitle(f"{name} CN Plots", fontsize=30)
    fig.supylabel("Copy Number")
    fig.supxlabel("Position")
    # rows = 6
    # column = 4
    # grid = plt.GridSpec(rows, column, wspace=0.25, hspace=0.25)
    # for i in range(24):
    i = 9
    row, col = divmod(i, 4)
    selected_x = []
    selected_y = []
    chr_idx = i + 1  # Chromosomes are 1-indexed
    for record in cnr_data:
        if record["chromosome"] == chr_idx:
            selected_x.append((record["start"] + record["end"]) / 2)
            selected_y.append(record["log2"])

    subplot: matplotlib.axes.Axes = ax
    subplot.set_ylim((0, 10))
    subplot.xaxis.set_major_formatter(ticker.EngFormatter())
    subplot.scatter(x=selected_x, y=selected_y, alpha=0.5, s=10, c="#0072b2")

    chr_name = "chr" + str(CHR_IDX_TO_NAME.get(chr_idx, chr_idx))
    centromere_intv = chr_arms[chr_name].interval
    centromere_midpt = (centromere_intv.start + centromere_intv.end) / 2
    subplot.axvline(
        x=centromere_midpt, color="red", linestyle="--", label="Centromere"
    )
    if chr_name == "chr8":
        subplot.axvline(
            x=128643133, color="purple", linestyle="--", label="MYC"
        )
        handles, labels = subplot.get_legend_handles_labels()
        fig.legend(
            handles=handles,
            labels=labels,
            loc="upper right",
            # bbox_to_anchor=(1, -0.1),
            ncol=len(labels),
            bbox_transform=fig.transFigure,
        )

    subplot.set_title(chr_name)

    # labels = ["Centromere", "MYC"]

    fig.tight_layout()

    plt.savefig(f"{output_dir}/{name}.pdf", dpi=30, bbox_inches="tight")


def plot_cnr(
    cnr_file: typer.FileText, output_dir: pathlib.Path, name: str
) -> None:
    parsed_data = parse_cnr_file(cnr_file)
    plot_cnr_data(parsed_data, output_dir, name)
