import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def plot_amplicon_properties(ecdna_props, bfb_props, bfb_sim_props):
    """
    Create plots for each metric in AmpliconProperties, with different colors for ecDNA, BFB, and BFB simulation data.

    Args:
        ecdna_props: List of AmpliconProperties from ecDNA samples
        bfb_props: List of AmpliconProperties from BFB samples
        bfb_sim_props: List of AmpliconProperties from BFB simulation
    """
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
    palette = {"ecDNA": "blue", "BFB": "red", "BFB Sim": "green"}

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
    return plt
