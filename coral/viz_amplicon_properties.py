#!/usr/bin/env python3

import logging
import pathlib
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from coral.breakpoint import breakpoint_graph, parse_graph
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.classify import classify_types, classify_utils

# Setup logging
logger = logging.getLogger()
logger.setLevel(logging.ERROR)
warnings.filterwarnings("ignore")

# Define directories
ecDNA_dir = pathlib.Path("/mnt/gq/bfb/ecDNA")
bfb_dir = pathlib.Path("/mnt/gq/bfb/BFB")
bfb_sim_dir = pathlib.Path("/mnt/gq/bfb/BFB/simulation")
ecDNA_sim_dir = pathlib.Path(
    "/mnt/gq/bfb/ecDNA/simulation"
)  # New ecDNA simulation directory


def load_graphs(dir: pathlib.Path) -> list[BreakpointGraph]:
    graphs = []
    for sample_dir in dir.iterdir():
        if not sample_dir.is_dir():
            continue
        for bp_graph in sample_dir.glob("*_graph.txt"):
            graphs.append(parse_graph.parse_breakpoint_graph(bp_graph.open()))
    return graphs


def main():
    # Load graph data
    ecdna_props = [
        classify_types.AmpliconProperties.from_breakpoint_graph(bp_graph)
        for bp_graph in load_graphs(ecDNA_dir)
    ]
    bfb_props = [
        classify_types.AmpliconProperties.from_breakpoint_graph(bp_graph)
        for bp_graph in load_graphs(bfb_dir)
    ]
    bfb_sim_props = [
        classify_types.AmpliconProperties.from_breakpoint_graph(bp_graph)
        for bp_graph in load_graphs(bfb_sim_dir)
    ]
    # Add ecDNA simulation properties
    ecdna_sim_props = [
        classify_types.AmpliconProperties.from_breakpoint_graph(bp_graph)
        for bp_graph in load_graphs(ecDNA_sim_dir)
    ]

    # Create a dataframe with all the properties
    data = []

    # Add ecDNA properties
    for prop in ecdna_props:
        data.append(
            {
                "type": "ecDNA",
                "rearranged_edges": prop.rearranged_edges,
                "fb_edges": prop.fb_edges,
                "fb_prop": prop.fb_prop,
                "max_cn": prop.max_cn,
                "tot_amplified_seq_len": prop.tot_amplified_seq_len,
            }
        )

    # Add BFB properties
    for prop in bfb_props:
        data.append(
            {
                "type": "BFB",
                "rearranged_edges": prop.rearranged_edges,
                "fb_edges": prop.fb_edges,
                "fb_prop": prop.fb_prop,
                "max_cn": prop.max_cn,
                "tot_amplified_seq_len": prop.tot_amplified_seq_len,
            }
        )

    # Add BFB simulation properties
    for prop in bfb_sim_props:
        data.append(
            {
                "type": "BFB Sim",
                "rearranged_edges": prop.rearranged_edges,
                "fb_edges": prop.fb_edges,
                "fb_prop": prop.fb_prop,
                "max_cn": prop.max_cn,
                "tot_amplified_seq_len": prop.tot_amplified_seq_len,
            }
        )

    # Add ecDNA simulation properties
    for prop in ecdna_sim_props:
        data.append(
            {
                "type": "ecDNA Sim",
                "rearranged_edges": prop.rearranged_edges,
                "fb_edges": prop.fb_edges,
                "fb_prop": prop.fb_prop,
                "max_cn": prop.max_cn,
                "tot_amplified_seq_len": prop.tot_amplified_seq_len,
            }
        )

    df = pd.DataFrame(data)

    # Convert dataframe to long format for easier plotting
    df_long = pd.melt(
        df, id_vars=["type"], var_name="metric", value_name="value"
    )

    # Set the style
    sns.set(style="whitegrid")

    # Plot each metric in a separate subplot
    metrics = [
        "rearranged_edges",
        "fb_edges",
        "fb_prop",
        "max_cn",
        "tot_amplified_seq_len",
    ]
    fig, axes = plt.subplots(3, 2, figsize=(16, 18))
    axes = axes.flatten()

    # Color palette for different types
    palette = {
        "ecDNA": "blue",
        "BFB": "red",
        "BFB Sim": "green",
        "ecDNA Sim": "purple",
    }

    # Plot each metric
    for i, metric in enumerate(metrics):
        # Filter data for the current metric
        metric_data = df_long[df_long["metric"] == metric]

        # Create boxplot with individual points
        ax = axes[i]
        sns.boxplot(
            x="type", y="value", data=metric_data, palette=palette, ax=ax
        )
        sns.stripplot(
            x="type",
            y="value",
            data=metric_data,
            jitter=True,
            dodge=False,
            alpha=0.7,
            palette=palette,
            ax=ax,
        )

        # Set titles and labels
        ax.set_title(metric.replace("_", " ").title(), fontsize=14)
        ax.set_xlabel("Amplicon Type", fontsize=12)
        ax.set_ylabel("Value", fontsize=12)

        # For very large values, use log scale
        if metric == "tot_amplified_seq_len":
            ax.set_yscale("log")
            ax.set_title(
                metric.replace("_", " ").title() + " (log scale)", fontsize=14
            )

    # Add a histogram for fb_prop in the last subplot slot
    ax = axes[5]
    for amplicon_type, color in palette.items():
        fb_prop_values = df[df["type"] == amplicon_type]["fb_prop"]
        sns.kdeplot(
            fb_prop_values,
            ax=ax,
            label=amplicon_type,
            color=color,
            fill=True,
            alpha=0.3,
        )

    ax.set_title("FB Prop Distribution", fontsize=14)
    ax.set_xlabel("FB Prop Value", fontsize=12)
    ax.set_ylabel("Density", fontsize=12)
    ax.legend()

    # Adjust layout and show the plot
    plt.tight_layout()
    plt.savefig(
        "amplicon_property_comparison.png", dpi=300, bbox_inches="tight"
    )
    plt.show()

    print("Visualization complete. Saved to 'amplicon_property_comparison.png'")

    # Create a new scatter plot showing FB prop vs max CN colored by type
    plt.figure(figsize=(12, 8))

    # Create scatter plot with fb_prop on y-axis and max_cn on x-axis
    scatter = sns.scatterplot(
        data=df,
        x="max_cn",
        y="fb_prop",
        hue="type",
        palette=palette,
        s=100,
        alpha=0.7,
        style="type",
        edgecolor="black",
    )

    # Add title and labels
    plt.title("FB Proportion vs Max Copy Number by Amplicon Type", fontsize=16)
    plt.xlabel("Max Copy Number", fontsize=14)
    plt.ylabel("FB Proportion", fontsize=14)

    # Add a legend
    plt.legend(title="Amplicon Type", fontsize=12, title_fontsize=14)

    # Add a grid
    plt.grid(True, linestyle="--", alpha=0.7)

    # Improve aesthetics
    plt.tight_layout()

    # Save the plot
    plt.savefig("fb_prop_vs_max_cn.png", dpi=300, bbox_inches="tight")
    plt.show()

    # Print summary statistics for each metric by type
    print("\nSummary Statistics:")
    for metric in metrics:
        print(f"\n{metric.replace('_', ' ').title()}:")
        print(df.groupby("type")[metric].describe())


if __name__ == "__main__":
    main()
