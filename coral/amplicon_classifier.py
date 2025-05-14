#!/usr/bin/env python3

import dataclasses
import logging
import os
import pathlib
import traceback
import warnings
from io import StringIO

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pydotplus
import scipy.stats as stats
import seaborn as sns
from matplotlib.patches import Ellipse
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import mutual_info_classif
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import (
    StratifiedKFold,
    cross_val_predict,
    cross_val_score,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier, export_graphviz, plot_tree

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
# ecDNA_sim_dir = pathlib.Path(
#     "/mnt/gq/globus/sims/reconstructions/default/20250311"
# )
ecDNA_sim_dir = pathlib.Path("/mnt/gq/globus/sims/ref")


def load_graphs(
    dir: pathlib.Path, suffix: str | None = None
) -> list[tuple[BreakpointGraph, str]]:
    """Load breakpoint graphs and return a list of tuples (graph, filepath)"""
    graphs = []
    for sample_dir in dir.iterdir():
        if not sample_dir.is_dir():
            continue
        sample_dir = sample_dir / suffix if suffix is not None else sample_dir
        for bp_graph in sample_dir.glob("*_graph.txt"):
            try:
                graph = parse_graph.parse_breakpoint_graph(bp_graph.open())
                graphs.append((graph, str(bp_graph)))
            except Exception as e:
                print(f"Error loading graph {bp_graph}: {e}")
    return graphs


def get_properties(
    graphs: list[tuple[BreakpointGraph, str]],
) -> list[tuple[classify_types.AmpliconProperties, str]]:
    """Extract AmpliconProperties from a list of breakpoint graphs and return with filepaths"""
    properties = []
    low_complexity_region_map = classify_utils.build_low_complexity_region_map()
    for graph, filepath in graphs:
        try:
            props = classify_types.AmpliconProperties.from_breakpoint_graph(
                graph, low_complexity_region_map
            )
            properties.append((props, filepath))
        except Exception as e:
            print(f"Error processing graph: {e}")
    return properties


def plot_bfb_sim_fb_prop_distribution(properties_list):
    """
    Plot the distribution of foldback proportion in BFB simulation samples.

    Args:
        properties_list: List of tuples containing (AmpliconProperties, filepath)
    """
    plots_dir = "classify_plots"
    print(f"\nStarting BFB simulation foldback proportion analysis...")
    print(f"Total properties to analyze: {len(properties_list)}")

    try:
        # Ensure the plots directory exists
        os.makedirs(plots_dir, exist_ok=True)

        # Filter for BFB simulation samples
        bfb_sim_props = [
            (props, filepath)
            for props, filepath in properties_list
            if "BFB/simulation" in filepath
        ]

        print(f"Found {len(bfb_sim_props)} BFB simulation samples for analysis")

        if not bfb_sim_props:
            print("No BFB simulation samples found. Skipping plot generation.")
            return

        # Extract the foldback proportions
        fb_props = [props.fb_prop for props, _ in bfb_sim_props]

        # Print sample info
        for i, (props, filepath) in enumerate(bfb_sim_props):
            print(
                f"Sample {i+1}: {filepath}, Foldback proportion: {props.fb_prop:.4f}"
            )

        # Create plot
        plt.figure(figsize=(10, 6))

        # Plot histogram
        plt.hist(fb_props, bins=15, alpha=0.7, color="blue", edgecolor="black")

        # Add details
        plt.xlabel("Foldback Proportion")
        plt.ylabel("Frequency")
        plt.title(
            "Distribution of Foldback Proportion in BFB Simulation Samples"
        )

        # Add statistics to plot
        mean_value = np.mean(fb_props)
        median_value = np.median(fb_props)

        plt.axvline(
            mean_value,
            color="red",
            linestyle="dashed",
            linewidth=1,
            label=f"Mean: {mean_value:.4f}",
        )
        plt.axvline(
            median_value,
            color="green",
            linestyle="dashed",
            linewidth=1,
            label=f"Median: {median_value:.4f}",
        )

        plt.legend()
        plt.grid(True, alpha=0.3)

        # Save plot
        output_file = os.path.join(
            plots_dir, "bfb_sim_fb_prop_distribution.png"
        )
        print(f"Saving plot to: {os.path.abspath(output_file)}")
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close()

        if os.path.exists(output_file):
            print(f"Successfully saved plot to {output_file}")
        else:
            print(f"ERROR: Failed to save plot to {output_file}")

        # Create summary statistics
        stats = {
            "mean": mean_value,
            "median": median_value,
            "min": min(fb_props),
            "max": max(fb_props),
            "std": np.std(fb_props),
            "count": len(fb_props),
        }

        # Save statistics to file
        stats_file = os.path.join(plots_dir, "bfb_sim_fb_prop_stats.txt")
        print(f"Saving statistics to: {os.path.abspath(stats_file)}")
        with open(stats_file, "w") as f:
            f.write("BFB Simulation Foldback Proportion Statistics\n")
            f.write("-" * 45 + "\n")
            f.write(f"Number of samples: {stats['count']}\n")
            f.write(f"Mean: {stats['mean']:.4f}\n")
            f.write(f"Median: {stats['median']:.4f}\n")
            f.write(f"Standard Deviation: {stats['std']:.4f}\n")
            f.write(f"Minimum: {stats['min']:.4f}\n")
            f.write(f"Maximum: {stats['max']:.4f}\n")

        if os.path.exists(stats_file):
            print(f"Successfully saved statistics to {stats_file}")
        else:
            print(f"ERROR: Failed to save statistics to {stats_file}")

        print(f"BFB simulation foldback proportion analysis completed")

    except Exception as e:
        print(f"Error in plot_bfb_sim_fb_prop_distribution: {e}")
        traceback.print_exc()


def generate_plots_without_bfb_sim(df):
    """
    Generate a new set of plots excluding BFB simulation data.

    Args:
        df: DataFrame containing all amplicon data
    """
    # Create directory for plots without BFB sim
    plots_dir = os.path.join("classify_plots", "no_bfb_sim")
    os.makedirs(plots_dir, exist_ok=True)
    print(f"\nGenerating plots without BFB simulation data in {plots_dir}")

    try:
        print("Step 1: Filtering out BFB Sim samples")
        # Filter out BFB Sim samples
        no_bfb_sim_df = df[df["original_type"] != "BFB Sim"].copy()
        print(f"Dataset without BFB Sim: {len(no_bfb_sim_df)} samples")

        # Count by type after filtering
        type_counts = no_bfb_sim_df["original_type"].value_counts()
        print("\nSample counts after filtering:")
        for idx, count in type_counts.items():
            print(f"{idx}: {count}")

        print("\nStep 2: Preparing data for classification")
        try:
            # Use the same pipeline as in the main analysis
            X = no_bfb_sim_df.drop(
                ["type", "original_type", "filepath"], axis=1
            )
            y = no_bfb_sim_df["type"]

            # Run classification analysis and get predictions
            results_df = run_classification_analysis(
                X, y, no_bfb_sim_df, plots_dir
            )
        except Exception as e:
            print(f"Error in classification steps: {e}")
            traceback.print_exc()
            # Continue with the rest of the visualization without relying on predictions
            no_bfb_sim_df["predicted"] = no_bfb_sim_df[
                "type"
            ]  # Default to true labels
            results_df = no_bfb_sim_df

        print("\nStep 3: Creating detailed confusion matrix by original type")
        try:
            # Create detailed confusion matrix by original type
            original_types = results_df["original_type"].unique()
            predicted_types = ["BFB", "ecDNA"]

            # Initialize matrix with zeros
            detailed_cm = np.zeros((len(original_types), len(predicted_types)))

            # Fill the matrix
            for i, orig_type in enumerate(original_types):
                for j, pred_type in enumerate(predicted_types):
                    detailed_cm[i, j] = (
                        (results_df["original_type"] == orig_type)
                        & (results_df["predicted"] == pred_type)
                    ).sum()

            # Plot detailed confusion matrix
            plt.figure(figsize=(10, 8))
            sns.heatmap(
                detailed_cm,
                annot=True,
                fmt=".0f",
                cmap="Blues",
                xticklabels=predicted_types,
                yticklabels=original_types,
            )
            plt.xlabel("Predicted Type")
            plt.ylabel("Original Type")
            plt.title("Detailed Confusion Matrix by Original Amplicon Type")
            plt.savefig(
                os.path.join(plots_dir, "detailed_confusion_matrix.png")
            )
            plt.close()

            # Calculate accuracy for each original type
            print("\nAccuracy by original type:")
            accuracies = {}
            for orig_type in original_types:
                type_subset = results_df[
                    results_df["original_type"] == orig_type
                ]
                if len(type_subset) > 0:
                    if orig_type.startswith("ecDNA"):
                        correct = (type_subset["predicted"] == "ecDNA").sum()
                    else:
                        correct = (type_subset["predicted"] == "BFB").sum()
                    accuracies[orig_type] = correct / len(type_subset)
                    print(f"{orig_type}: {accuracies[orig_type]:.4f}")

        except Exception as e:
            print(f"Error creating confusion matrices: {e}")
            traceback.print_exc()

        print("\nStep 4: Handling outliers for visualization")
        try:
            # Handle outliers for visualization
            Q1 = results_df["max_cn"].quantile(0.25)
            Q3 = results_df["max_cn"].quantile(0.75)
            IQR = Q3 - Q1
            outlier_threshold = Q3 + 1.5 * IQR

            viz_df = results_df[results_df["max_cn"] <= outlier_threshold]
            print(
                f"\nRemoved {len(results_df) - len(viz_df)} max_cn outliers for visualization"
            )
        except Exception as e:
            print(f"Error handling outliers: {e}")
            traceback.print_exc()
            viz_df = results_df  # Use all data if outlier removal fails

        print("\nStep 5: Creating amplicon distribution plot")
        try:
            # Use VisualizationUtils to create the scatter plot
            VisualizationUtils.plot_scatter(
                x="max_cn",
                y="fb_prop",
                data=viz_df,
                output_path=os.path.join(
                    plots_dir, "amplicon_distribution.png"
                ),
                title="Distribution of Amplicon Types without BFB Simulation",
                xlabel="Maximum Copy Number (max_cn)",
                ylabel="Foldback Proportion (fb_prop)",
                hue="original_type",
                style="original_type",
            )

            # Find a meaningful decision boundary and add to plot if needed
            bfb_samples = viz_df[viz_df["original_type"] == "BFB"]
            if len(bfb_samples) > 0:
                bfb_median_fb_prop = bfb_samples["fb_prop"].median()
                print(
                    f"BFB median foldback proportion: {bfb_median_fb_prop:.4f}"
                )
        except Exception as e:
            print(f"Error creating amplicon distribution plot: {e}")
            traceback.print_exc()

        print("\nStep 6: Creating ecDNA comparison plot")
        try:
            # Create ecDNA vs ecDNA Sim plot
            ecdna_df = viz_df[
                viz_df["original_type"].isin(["ecDNA", "ecDNA Sim"])
            ]

            if len(ecdna_df) > 0:
                VisualizationUtils.plot_scatter(
                    x="max_cn",
                    y="fb_prop",
                    data=ecdna_df,
                    output_path=os.path.join(plots_dir, "ecdna_comparison.png"),
                    title="ecDNA vs ecDNA Simulation Comparison",
                    xlabel="Maximum Copy Number (max_cn)",
                    ylabel="Foldback Proportion (fb_prop)",
                    hue="original_type",
                    style="original_type",
                )
            else:
                print("No ecDNA samples found for comparison plot")
        except Exception as e:
            print(f"Error creating ecDNA comparison plot: {e}")
            traceback.print_exc()

        print("\nStep 7: Creating feature comparison grid")
        try:
            # Create feature comparison plots (2x2 grid for different features)
            feature_pairs = [
                ("fb_prop", "max_cn"),
                ("fb_edges", "rearranged_edges"),
                ("fb_prop", "fb_edges"),
                ("max_cn", "rearranged_edges"),
            ]

            plt.figure(figsize=(14, 12))

            for i, (x_feature, y_feature) in enumerate(feature_pairs):
                plt.subplot(2, 2, i + 1)

                sns.scatterplot(
                    data=viz_df,
                    x=x_feature,
                    y=y_feature,
                    hue="original_type",
                    style="original_type",
                    alpha=0.7,
                    s=80,
                    legend=False if i > 0 else True,
                )

                plt.xlabel(x_feature)
                plt.ylabel(y_feature)

                if i == 0:
                    plt.legend(title="Sample Types", loc="best")

            plt.suptitle(
                "Feature Comparison without BFB Simulation", fontsize=16
            )
            plt.tight_layout(rect=[0, 0, 1, 0.96])
            plt.savefig(os.path.join(plots_dir, "feature_comparison.png"))
            plt.close()
        except Exception as e:
            print(f"Error creating feature comparison grid: {e}")
            traceback.print_exc()

        # Generate the sequence feature plots for the no_bfb_sim analysis
        print("\nGenerating sequence feature plots for no_bfb_sim analysis...")
        plot_sequence_feature_distributions(results_df, plots_dir)

        print("\nStep 8: Performing PCA analysis")
        try:
            # Use VisualizationUtils for PCA analysis
            VisualizationUtils.create_pca_plots(
                X, results_df, plots_dir, title_suffix=" (without BFB Sim)"
            )
        except Exception as e:
            print(f"Error performing PCA analysis: {e}")
            traceback.print_exc()

        print(f"\nGenerated plots without BFB simulation data in {plots_dir}:")
        print(
            "1. binary_confusion_matrix.png - Classification performance without BFB Sim"
        )
        print(
            "2. detailed_confusion_matrix.png - How each original type was classified"
        )
        print("3. feature_importance.png - Top features for classification")
        print(
            "4. amplicon_distribution.png - Distribution of samples without BFB Sim"
        )
        print(
            "5. ecdna_comparison.png - Comparison of ecDNA vs ecDNA Sim samples"
        )
        print("6. feature_comparison.png - Grid of feature relationships")
        print(
            "7. pca_scree_plot.png - Explained variance by each principal component"
        )
        print(
            "8. pca_loadings_heatmap.png - Loading coefficients showing feature contributions to PCs"
        )
        print(
            "9. pca_top_contributors.png - Top features contributing to each principal component"
        )
        print(
            "10. feature_correlation_matrix.png - Correlation matrix of original features"
        )
        print("11. pca_analysis.png - PCA with feature contribution vectors")
        print("12. pca_clusters.png - PCA with cluster visualization")
        print("13. pca_biplot.png - PCA biplot highlighting top 5 features")
        print(
            "\nAdditional sequence feature plots in 'sequence_features' subdirectory:"
        )
        print("- Feature boxplots, scatterplots, and distributions")
        print("- Sequence feature correlations heatmap")
        print("\nAdditional files:")
        print(
            "pca_loadings.csv - Full matrix of loading coefficients for all principal components"
        )

        # Generate feature entropy plot only (skip mutual information and decision tree)
        try:
            title_suffix = " (without BFB Sim)"

            # Create entropy directory if it doesn't exist
            entropy_dir = os.path.join(plots_dir, "entropy_analysis")
            os.makedirs(entropy_dir, exist_ok=True)

            # Calculate feature entropy directly
            print("\nCalculating feature entropy for no_bfb_sim dataset...")
            X_entropy = results_df.drop(
                ["type", "original_type", "filepath"], axis=1
            )
            feature_names = X_entropy.columns.tolist()

            # Calculate entropy values
            entropy_values = {}
            for feature in feature_names:
                values = X_entropy[feature].values

                # For continuous features
                if X_entropy[feature].dtype in [
                    "float64",
                    "float32",
                    "int64",
                    "int32",
                ]:
                    # Use histogram to bin the data
                    hist, _ = np.histogram(values, bins=20)
                    probabilities = hist / len(values)
                    # Calculate entropy (only on non-zero probabilities)
                    probabilities = probabilities[probabilities > 0]
                    entropy = -np.sum(probabilities * np.log2(probabilities))
                else:
                    # For categorical variables
                    value_counts = X_entropy[feature].value_counts(
                        normalize=True
                    )
                    entropy = stats.entropy(value_counts, base=2)

                entropy_values[feature] = entropy

            # Create entropy plot
            fig, ax = plt.subplots(figsize=(14, 8))

            # Sort features by entropy
            sorted_features = sorted(
                entropy_values.items(), key=lambda x: x[1], reverse=True
            )
            features = [x[0] for x in sorted_features]
            entropies = [x[1] for x in sorted_features]

            # Plot bars
            bars = ax.bar(range(len(features)), entropies, color="skyblue")

            # Add color gradient
            norm = mcolors.Normalize(min(entropies), max(entropies))
            for i, bar in enumerate(bars):
                bar.set_color(plt.cm.viridis(norm(entropies[i])))

            # Add labels and title
            ax.set_xlabel("Feature")
            ax.set_ylabel("Entropy (bits)")
            ax.set_title(f"Feature Entropy{title_suffix}")
            ax.set_xticks(range(len(features)))
            ax.set_xticklabels(features, rotation=90)
            ax.grid(axis="y", alpha=0.3)

            # Save plot
            output_file = os.path.join(
                entropy_dir,
                f'feature_entropy{title_suffix.replace(" ", "_").lower()}.png',
            )
            fig.tight_layout()
            fig.savefig(output_file, dpi=300, bbox_inches="tight")
            plt.close(fig)

            print(f"Created and saved entropy visualization to {output_file}")

        except Exception as e:
            print(f"Error in entropy analysis for no_bfb_sim dataset: {e}")
            traceback.print_exc()

    except Exception as e:
        print(f"Error generating plots without BFB Sim: {e}")
        traceback.print_exc()


def generate_plots_without_simulations(df):
    """
    Generate a new set of plots excluding all simulation data (both BFB and ecDNA simulations).
    Plots will be saved to 'classify_plots/no_sims' directory.

    Args:
        df: DataFrame containing all amplicon data
    """
    # Create directory for plots without simulations
    plots_dir = os.path.join("classify_plots", "no_sims")
    os.makedirs(plots_dir, exist_ok=True)
    print(f"\nGenerating plots without simulation data in {plots_dir}")

    try:
        print("Step 1: Filtering out simulation samples")
        # Filter out all simulation samples
        no_sims_df = df[
            ~df["original_type"].isin(["BFB Sim", "ecDNA Sim"])
        ].copy()
        print(f"Dataset without simulations: {len(no_sims_df)} samples")

        # Count by type after filtering
        type_counts = no_sims_df["original_type"].value_counts()
        print("\nSample counts after filtering:")
        for idx, count in type_counts.items():
            print(f"{idx}: {count}")

        print("\nStep 2: Preparing data for classification")
        try:
            # Prepare data for classification
            X = no_sims_df.drop(["type", "original_type", "filepath"], axis=1)
            y = no_sims_df["type"]

            # Run classification analysis and get results
            results_df = run_classification_analysis(
                X, y, no_sims_df, plots_dir
            )
        except Exception as e:
            print(f"Error in classification steps: {e}")
            traceback.print_exc()
            # Continue with the rest of the visualization without relying on predictions
            no_sims_df["predicted"] = no_sims_df[
                "type"
            ]  # Default to true labels
            results_df = no_sims_df

        print("\nStep 3: Handling outliers for visualization")
        try:
            # Handle outliers for visualization
            Q1 = results_df["max_cn"].quantile(0.25)
            Q3 = results_df["max_cn"].quantile(0.75)
            IQR = Q3 - Q1
            outlier_threshold = Q3 + 1.5 * IQR

            viz_df = results_df[results_df["max_cn"] <= outlier_threshold]
            print(
                f"\nRemoved {len(results_df) - len(viz_df)} max_cn outliers for visualization"
            )
        except Exception as e:
            print(f"Error handling outliers: {e}")
            traceback.print_exc()
            viz_df = results_df  # Use all data if outlier removal fails

        print("\nStep 4: Creating amplicon distribution plot")
        try:
            # Use VisualizationUtils to create the scatter plot
            VisualizationUtils.plot_scatter(
                x="max_cn",
                y="fb_prop",
                data=viz_df,
                output_path=os.path.join(
                    plots_dir, "amplicon_distribution.png"
                ),
                title="Distribution of Real Amplicon Types",
                xlabel="Maximum Copy Number (max_cn)",
                ylabel="Foldback Proportion (fb_prop)",
                hue="original_type",
                style="original_type",
            )

            # Find a meaningful decision boundary and add to plot if needed
            bfb_samples = viz_df[viz_df["original_type"] == "BFB"]
            ecdna_samples = viz_df[viz_df["original_type"] == "ecDNA"]

            if len(bfb_samples) > 0 and len(ecdna_samples) > 0:
                bfb_median_fb_prop = bfb_samples["fb_prop"].median()
                print(
                    f"BFB median foldback proportion: {bfb_median_fb_prop:.4f}"
                )
        except Exception as e:
            print(f"Error creating amplicon distribution plot: {e}")
            traceback.print_exc()

        print("\nStep 5: Creating feature comparison grid")
        try:
            # Create feature comparison plots (2x2 grid for different features)
            feature_pairs = [
                ("fb_prop", "max_cn"),
                ("fb_edges", "rearranged_edges"),
                ("fb_prop", "fb_edges"),
                ("max_cn", "rearranged_edges"),
            ]

            plt.figure(figsize=(14, 12))

            for i, (x_feature, y_feature) in enumerate(feature_pairs):
                plt.subplot(2, 2, i + 1)

                sns.scatterplot(
                    data=viz_df,
                    x=x_feature,
                    y=y_feature,
                    hue="original_type",
                    style="original_type",
                    alpha=0.7,
                    s=80,
                    legend=False if i > 0 else True,
                )

                plt.xlabel(x_feature)
                plt.ylabel(y_feature)

                if i == 0:
                    plt.legend(title="Sample Types", loc="best")

            plt.suptitle("Feature Comparison of Real Amplicons", fontsize=16)
            plt.tight_layout(rect=[0, 0, 1, 0.96])
            plt.savefig(os.path.join(plots_dir, "feature_comparison.png"))
            plt.close()
        except Exception as e:
            print(f"Error creating feature comparison grid: {e}")
            traceback.print_exc()

        # Generate the sequence feature plots
        print("\nGenerating sequence feature plots for real data...")
        plot_sequence_feature_distributions(results_df, plots_dir)

        print("\nStep 6: Performing PCA analysis")
        try:
            # Use VisualizationUtils for PCA analysis
            VisualizationUtils.create_pca_plots(
                X, results_df, plots_dir, title_suffix=" (Real Data Only)"
            )
        except Exception as e:
            print(f"Error performing PCA analysis: {e}")
            traceback.print_exc()

        print(f"\nGenerated plots without simulation data in {plots_dir}:")
        print(
            "1. binary_confusion_matrix.png - Classification performance on real data"
        )
        print(
            "2. feature_importance.png - Top features for classification of real data"
        )
        print("3. amplicon_distribution.png - Distribution of real samples")
        print("4. feature_comparison.png - Grid of feature relationships")
        print(
            "5. pca_scree_plot.png - Explained variance by each principal component"
        )
        print(
            "6. pca_loadings_heatmap.png - Loading coefficients showing feature contributions"
        )
        print("7. pca_analysis.png - PCA with feature contribution vectors")
        print("8. pca_clusters.png - PCA with cluster visualization")
        print(
            "\nAdditional sequence feature plots in 'sequence_features' subdirectory"
        )
        print("\nAdditional files:")
        print(
            "pca_loadings.csv - Full matrix of loading coefficients for all principal components"
        )

        # Generate feature entropy plot only (skip mutual information and decision tree)
        try:
            title_suffix = " (Real Data Only)"

            # Create entropy directory if it doesn't exist
            entropy_dir = os.path.join(plots_dir, "entropy_analysis")
            os.makedirs(entropy_dir, exist_ok=True)

            # Calculate feature entropy directly
            print("\nCalculating feature entropy for real data only dataset...")
            X_entropy = results_df.drop(
                ["type", "original_type", "filepath"], axis=1
            )
            feature_names = X_entropy.columns.tolist()

            # Calculate entropy values
            entropy_values = {}
            for feature in feature_names:
                values = X_entropy[feature].values

                # For continuous features
                if X_entropy[feature].dtype in [
                    "float64",
                    "float32",
                    "int64",
                    "int32",
                ]:
                    # Use histogram to bin the data
                    hist, _ = np.histogram(values, bins=20)
                    probabilities = hist / len(values)
                    # Calculate entropy (only on non-zero probabilities)
                    probabilities = probabilities[probabilities > 0]
                    entropy = -np.sum(probabilities * np.log2(probabilities))
                else:
                    # For categorical variables
                    value_counts = X_entropy[feature].value_counts(
                        normalize=True
                    )
                    entropy = stats.entropy(value_counts, base=2)

                entropy_values[feature] = entropy

            # Create entropy plot
            fig, ax = plt.subplots(figsize=(14, 8))

            # Sort features by entropy
            sorted_features = sorted(
                entropy_values.items(), key=lambda x: x[1], reverse=True
            )
            features = [x[0] for x in sorted_features]
            entropies = [x[1] for x in sorted_features]

            # Plot bars
            bars = ax.bar(range(len(features)), entropies, color="skyblue")

            # Add color gradient
            norm = mcolors.Normalize(min(entropies), max(entropies))
            for i, bar in enumerate(bars):
                bar.set_color(plt.cm.viridis(norm(entropies[i])))

            # Add labels and title
            ax.set_xlabel("Feature")
            ax.set_ylabel("Entropy (bits)")
            ax.set_title(f"Feature Entropy{title_suffix}")
            ax.set_xticks(range(len(features)))
            ax.set_xticklabels(features, rotation=90)
            ax.grid(axis="y", alpha=0.3)

            # Save plot
            output_file = os.path.join(
                entropy_dir,
                f'feature_entropy{title_suffix.replace(" ", "_").lower()}.png',
            )
            fig.tight_layout()
            fig.savefig(output_file, dpi=300, bbox_inches="tight")
            plt.close(fig)

            print(f"Created and saved entropy visualization to {output_file}")

        except Exception as e:
            print(f"Error in entropy analysis for real data only dataset: {e}")
            traceback.print_exc()

    except Exception as e:
        print(f"Error generating plots without simulations: {e}")
        traceback.print_exc()


def main():
    # Create directory for plots if it doesn't exist
    plots_dir = "classify_plots"
    os.makedirs(plots_dir, exist_ok=True)
    print(f"Saving plots to {plots_dir} directory")

    # Load and process data
    data_loader = AmpliconDataLoader(
        ecDNA_dir=ecDNA_dir,
        bfb_dir=bfb_dir,
        bfb_sim_dir=bfb_sim_dir,
        ecDNA_sim_dir=ecDNA_sim_dir,
    )
    df = data_loader.load_and_process()

    # Plot the distribution of foldback proportion in BFB simulation samples
    plot_bfb_sim_fb_prop_distribution(data_loader.bfb_sim_props)

    # Check for and remove NaN values
    nan_counts = df.isna().sum()
    if nan_counts.sum() > 0:
        print("\nNaN counts per column:")
        for col, count in nan_counts.items():
            if count > 0:
                print(f"{col}: {count}")

        original_rows = len(df)
        df = df.dropna()
        print(f"Dropped {original_rows - len(df)} rows with NaN values")

    print(f"Final DataFrame shape: {df.shape}")

    # Prepare data for training
    X = df.drop(["type", "original_type", "filepath"], axis=1)
    y = df["type"]  # Only predict binary classes: ecDNA or BFB

    # Create and run the classifier
    run_classification_analysis(X, y, df, plots_dir)

    # Generate the sequence feature plots
    plot_sequence_feature_distributions(df, plots_dir)

    # Generate entropy and decision tree plots
    print("\nGenerating entropy and decision tree plots...")
    generate_entropy_and_tree_plots(df, plots_dir)

    # Generate additional plot sets
    print("\nGenerating additional plot sets...")
    generate_plots_without_bfb_sim(df)
    generate_plots_without_simulations(df)

    print_summary()


class AmpliconDataLoader:
    """Class for loading and preprocessing amplicon data."""

    def __init__(self, ecDNA_dir, bfb_dir, bfb_sim_dir, ecDNA_sim_dir):
        """Initialize with directories containing amplicon data."""
        self.ecDNA_dir = ecDNA_dir
        self.bfb_dir = bfb_dir
        self.bfb_sim_dir = bfb_sim_dir
        self.ecDNA_sim_dir = ecDNA_sim_dir
        self.ecdna_graphs = None
        self.bfb_graphs = None
        self.bfb_sim_graphs = None
        self.ecdna_sim_graphs = None
        self.ecdna_props = None
        self.bfb_props = None
        self.bfb_sim_props = None
        self.ecdna_sim_props = None

    def load_graphs(self):
        """Load all breakpoint graphs from directories."""
        print("Loading breakpoint graphs...")
        try:
            print("Loading ecDNA graphs...")
            self.ecdna_graphs = load_graphs(self.ecDNA_dir)
            print("Loading BFB graphs...")
            self.bfb_graphs = load_graphs(self.bfb_dir)
            print("Loading BFB simulation graphs...")
            self.bfb_sim_graphs = load_graphs(self.bfb_sim_dir, suffix="coral")
            print("Loading ecDNA simulation graphs...")
            self.ecdna_sim_graphs = load_graphs(self.ecDNA_sim_dir)
        except Exception as e:
            print(f"Error loading graphs: {e}")
            traceback.print_exc()
            return False

        print(f"Loaded {len(self.ecdna_graphs)} ecDNA graphs")
        print(f"Loaded {len(self.bfb_graphs)} BFB graphs")
        print(f"Loaded {len(self.bfb_sim_graphs)} BFB simulation graphs")
        print(f"Loaded {len(self.ecdna_sim_graphs)} ecDNA simulation graphs")
        return True

    def extract_properties(self):
        """Extract AmpliconProperties from loaded graphs."""
        print("Extracting properties...")
        self.ecdna_props = get_properties(self.ecdna_graphs)
        self.bfb_props = get_properties(self.bfb_graphs)
        self.bfb_sim_props = get_properties(self.bfb_sim_graphs)
        self.ecdna_sim_props = get_properties(self.ecdna_sim_graphs)

        print(f"Extracted {len(self.ecdna_props)} ecDNA properties")
        print(f"Extracted {len(self.bfb_props)} BFB properties")
        print(f"Extracted {len(self.bfb_sim_props)} BFB simulation properties")
        print(
            f"Extracted {len(self.ecdna_sim_props)} ecDNA simulation properties"
        )

    def create_dataframe(self):
        """Create a dataframe with all the extracted properties."""
        print("Creating dataset...")
        data = []

        # Add ecDNA samples
        for props, filepath in self.ecdna_props:
            data.append(
                {
                    **dataclasses.asdict(props),
                    "type": "ecDNA",
                    "original_type": "ecDNA",
                    "filepath": filepath,
                }
            )
        print(f"Added {len(self.ecdna_props)} ecDNA samples to dataset")

        # Add ecDNA simulation samples
        for props, filepath in self.ecdna_sim_props:
            data.append(
                {
                    **dataclasses.asdict(props),
                    "type": "ecDNA",
                    "original_type": "ecDNA Sim",
                    "filepath": filepath,
                }
            )
        print(f"Added {len(self.ecdna_sim_props)} ecDNA Sim samples to dataset")

        # Add BFB samples
        for props, filepath in self.bfb_props:
            data.append(
                {
                    **dataclasses.asdict(props),
                    "type": "BFB",
                    "original_type": "BFB",
                    "filepath": filepath,
                }
            )
        print(f"Added {len(self.bfb_props)} BFB samples to dataset")

        # Add BFB simulation samples
        for props, filepath in self.bfb_sim_props:
            data.append(
                {
                    **dataclasses.asdict(props),
                    "type": "BFB",
                    "original_type": "BFB Sim",
                    "filepath": filepath,
                }
            )
        print(f"Added {len(self.bfb_sim_props)} BFB Sim samples to dataset")

        df = pd.DataFrame(data)
        print(f"DataFrame shape after creation: {df.shape}")

        # Count samples by original type
        original_type_counts = df["original_type"].value_counts()
        print("\nCounts by original type:")
        for idx, count in original_type_counts.items():
            print(f"{idx}: {count}")

        return df

    def process_list_fields(self, df):
        """Process list fields into statistical features."""
        print("Processing list fields for modeling...")

        # Handle bp_read_supports field
        df["mean_bp_read_support"] = df["bp_read_supports"].apply(
            lambda x: np.mean(x) if len(x) > 0 else 0
        )
        df["max_bp_read_support"] = df["bp_read_supports"].apply(
            lambda x: max(x) if len(x) > 0 else 0
        )
        df["min_bp_read_support"] = df["bp_read_supports"].apply(
            lambda x: min(x) if len(x) > 0 else 0
        )
        df["median_bp_read_support"] = df["bp_read_supports"].apply(
            lambda x: np.median(x) if len(x) > 0 else 0
        )
        df["std_bp_read_support"] = df["bp_read_supports"].apply(
            lambda x: np.std(x) if len(x) > 0 else 0
        )
        df["count_bp_read_supports"] = df["bp_read_supports"].apply(len)

        # Handle bp_cns field
        df["mean_bp_cn"] = df["bp_cns"].apply(
            lambda x: np.mean(x) if len(x) > 0 else 0
        )
        df["max_bp_cn"] = df["bp_cns"].apply(
            lambda x: max(x) if len(x) > 0 else 0
        )
        df["min_bp_cn"] = df["bp_cns"].apply(
            lambda x: min(x) if len(x) > 0 else 0
        )
        df["median_bp_cn"] = df["bp_cns"].apply(
            lambda x: np.median(x) if len(x) > 0 else 0
        )
        df["std_bp_cn"] = df["bp_cns"].apply(
            lambda x: np.std(x) if len(x) > 0 else 0
        )

        # Handle seq_cns field
        df["mean_seq_cn"] = df["seq_cns"].apply(
            lambda x: np.mean(x) if len(x) > 0 else 0
        )
        df["max_seq_cn"] = df["seq_cns"].apply(
            lambda x: max(x) if len(x) > 0 else 0
        )
        df["min_seq_cn"] = df["seq_cns"].apply(
            lambda x: min(x) if len(x) > 0 else 0
        )
        df["median_seq_cn"] = df["seq_cns"].apply(
            lambda x: np.median(x) if len(x) > 0 else 0
        )
        df["std_seq_cn"] = df["seq_cns"].apply(
            lambda x: np.std(x) if len(x) > 0 else 0
        )

        # Handle seq_coverages field
        df["mean_seq_coverage"] = df["seq_coverages"].apply(
            lambda x: np.mean(x) if len(x) > 0 else 0
        )
        df["max_seq_coverage"] = df["seq_coverages"].apply(
            lambda x: max(x) if len(x) > 0 else 0
        )
        df["min_seq_coverage"] = df["seq_coverages"].apply(
            lambda x: min(x) if len(x) > 0 else 0
        )
        df["median_seq_coverage"] = df["seq_coverages"].apply(
            lambda x: np.median(x) if len(x) > 0 else 0
        )
        df["std_seq_coverage"] = df["seq_coverages"].apply(
            lambda x: np.std(x) if len(x) > 0 else 0
        )

        # Calculate additional ratios that might be helpful
        df["cn_coverage_ratio"] = df.apply(
            lambda row: row["mean_seq_cn"] / row["mean_seq_coverage"]
            if row["mean_seq_coverage"] > 0
            else 0,
            axis=1,
        )
        df["bp_seq_cn_ratio"] = df.apply(
            lambda row: row["mean_bp_cn"] / row["mean_seq_cn"]
            if row["mean_seq_cn"] > 0
            else 0,
            axis=1,
        )

        # Remove the original list columns
        print("Removing original list columns...")
        df = df.drop(
            ["bp_read_supports", "bp_cns", "seq_cns", "seq_coverages"], axis=1
        )

        return df

    def display_sequence_stats(self, df):
        """Display summary statistics for sequence features by original type."""
        # Display summary statistics for sequence-related features by original type
        print("\nSequence features statistics by original type:")
        seq_features = [
            "mean_seq_cn",
            "max_seq_cn",
            "std_seq_cn",
            "mean_seq_coverage",
            "max_seq_coverage",
            "std_seq_coverage",
            "mean_bp_cn",
            "max_bp_cn",
            "cn_coverage_ratio",
            "bp_seq_cn_ratio",
        ]

        original_types = ["ecDNA", "ecDNA Sim", "BFB", "BFB Sim"]
        for orig_type in original_types:
            print(
                f"\n{orig_type} samples (n={len(df[df['original_type'] == orig_type])}):"
            )
            subset = df[df["original_type"] == orig_type]
            for feature in seq_features:
                mean_value = subset[feature].mean()
                std_value = subset[feature].std()
                min_value = subset[feature].min()
                max_value = subset[feature].max()
                print(
                    f"  {feature}: mean={mean_value:.4f}, std={std_value:.4f}, min={min_value:.4f}, max={max_value:.4f}"
                )

    def load_and_process(self):
        """Load, extract properties, and process data into a DataFrame."""
        if not self.load_graphs():
            return None

        self.extract_properties()
        df = self.create_dataframe()
        df = self.process_list_fields(df)
        self.display_sequence_stats(df)

        return df


def run_classification_analysis(X, y, df, plots_dir):
    """Run the classification analysis using the prepared data."""
    # Create pipeline for classification
    pipeline = Pipeline(
        [
            ("scaler", StandardScaler()),
            (
                "classifier",
                RandomForestClassifier(
                    n_estimators=100,
                    random_state=42,
                    class_weight="balanced",  # Add class weight to handle class imbalance
                ),
            ),
        ]
    )

    # Use stratified cross-validation for evaluation
    print("\nPerforming 5-fold stratified cross-validation...")
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    # Get cross-validated predictions for each sample
    y_pred = cross_val_predict(pipeline, X, y, cv=cv)

    # Create a copy of the dataframe with predictions
    results_df = df.copy()
    results_df["predicted"] = y_pred

    # Overall classification report
    print("\nBinary Classification Report (ecDNA vs BFB):")
    print(classification_report(y, y_pred))

    # Confusion matrix
    print("\nConfusion Matrix:")
    cm = confusion_matrix(y, y_pred)
    plt.figure(figsize=(8, 6))
    sns.heatmap(
        cm,
        annot=True,
        fmt="d",
        cmap="Blues",
        xticklabels=["BFB", "ecDNA"],
        yticklabels=["BFB", "ecDNA"],
    )
    plt.xlabel("Predicted")
    plt.ylabel("True")
    plt.title("Binary Confusion Matrix (ecDNA vs BFB)")
    plt.savefig(os.path.join(plots_dir, "binary_confusion_matrix.png"))
    plt.close()

    # Train the model on the full dataset for feature importance
    pipeline.fit(X, y)

    # Feature importance
    print("\nFeature Importance:")
    feature_importance = pipeline.named_steps["classifier"].feature_importances_
    feature_names = X.columns

    importance_df = pd.DataFrame(
        {"Feature": feature_names, "Importance": feature_importance}
    )
    importance_df = importance_df.sort_values("Importance", ascending=False)

    plt.figure(figsize=(12, 8))
    sns.barplot(x="Importance", y="Feature", data=importance_df.head(10))
    plt.title("Top 10 Feature Importance")
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "feature_importance.png"))
    plt.close()

    return results_df


def print_summary():
    """Print a summary of all generated plots and files."""
    print(
        "\nAnalysis complete. Check the generated plots in the 'classify_plots' directory:"
    )
    print(
        "1. binary_confusion_matrix.png - Overall binary classification performance"
    )
    print("2. feature_importance.png - Top features for classification")
    print(
        "3. detailed_confusion_matrix.png - How each original type was classified"
    )
    print(
        "4. fb_prop_vs_max_cn.png - Scatter plot with all samples by original type"
    )
    print(
        "5. amplicon_distribution.png - All samples with misclassifications highlighted"
    )
    print(
        "6. bfb_classification_analysis.png - Feature space analysis of BFB classifications"
    )
    print(
        "7. high_fb_ecdna_sim.png - Visualization of ecDNA Sim samples with high foldback proportion"
    )
    print(
        "8. bfb_sim_fb_prop_distribution.png - Distribution of foldback proportion in BFB simulation samples"
    )
    print("9. pca_analysis.png - PCA with feature contribution vectors")
    print("10. pca_clusters.png - PCA with cluster visualization")
    print("11. pca_biplot.png - PCA biplot highlighting top 5 features")
    print("12. sequence_features/* - Various sequence feature analyses")
    print(
        "13. entropy_analysis/* - Feature entropy and decision tree visualizations"
    )

    print(
        "\nAdditional plots without BFB simulation in 'classify_plots/no_bfb_sim' directory:"
    )
    print(
        "1. binary_confusion_matrix.png - Classification performance without BFB Sim"
    )
    print(
        "2. detailed_confusion_matrix.png - How each original type was classified"
    )
    print("3. feature_importance.png - Top features for classification")
    print(
        "4. amplicon_distribution.png - Distribution of samples without BFB Sim"
    )
    print("5. ecdna_comparison.png - Comparison of ecDNA vs ecDNA Sim samples")
    print("6. feature_comparison.png - Grid of feature relationships")
    print(
        "7. pca_scree_plot.png - Explained variance by each principal component"
    )
    print(
        "8. pca_loadings_heatmap.png - Loading coefficients showing feature contributions to PCs"
    )
    print(
        "9. pca_top_contributors.png - Top features contributing to each principal component"
    )
    print(
        "10. feature_correlation_matrix.png - Correlation matrix of original features"
    )
    print("11. pca_analysis.png - PCA with feature contribution vectors")
    print("12. pca_clusters.png - PCA with cluster visualization")
    print("13. pca_biplot.png - PCA biplot highlighting top 5 features")
    print(
        "\nAdditional sequence feature plots in 'sequence_features' subdirectory:"
    )
    print("- Feature boxplots, scatterplots, and distributions")
    print("- Sequence feature correlations heatmap")
    print("\nAdditional files:")
    print(
        "pca_loadings.csv - Full matrix of loading coefficients for all principal components"
    )

    # Remove the invalid call to generate_entropy_and_tree_plots with undefined variables
    print("\nGeneration of all plots is complete.")


class VisualizationUtils:
    """Utility class for creating visualizations and plots."""

    @staticmethod
    def create_pca_plots(X, df, plots_dir, title_suffix=""):
        """Create PCA analysis plots."""
        # Standardize the data
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)

        # Initialize PCA with multiple components for scree plot
        n_components = min(len(X.columns), 10)
        pca_full = PCA(n_components=n_components)
        pca_full.fit(X_scaled)

        # Create scree plot
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

        # Individual explained variance
        ax1.bar(
            range(1, len(pca_full.explained_variance_ratio_) + 1),
            pca_full.explained_variance_ratio_,
            alpha=0.7,
            color="skyblue",
        )
        ax1.step(
            range(1, len(pca_full.explained_variance_ratio_) + 1),
            pca_full.explained_variance_ratio_.cumsum(),
            where="mid",
            color="red",
            marker="o",
        )
        ax1.axhline(
            y=0.95,
            color="green",
            linestyle="--",
            alpha=0.7,
            label="95% Explained Variance",
        )
        ax1.axhline(
            y=0.8,
            color="orange",
            linestyle="--",
            alpha=0.7,
            label="80% Explained Variance",
        )

        # Add annotations for cumulative explained variance
        for i, ratio in enumerate(pca_full.explained_variance_ratio_.cumsum()):
            if i % 2 == 0 or ratio > 0.9:
                ax1.text(i + 1.2, ratio, f"{ratio:.2%}", va="center")

        ax1.set_xlabel("Principal Component")
        ax1.set_ylabel("Explained Variance Ratio")
        ax1.set_title("Scree Plot - Individual Component Contribution")
        ax1.set_xticks(range(1, len(pca_full.explained_variance_ratio_) + 1))
        ax1.legend()
        ax1.grid(alpha=0.3)

        # Cumulative explained variance
        cumulative = np.cumsum(pca_full.explained_variance_ratio_)
        ax2.plot(
            range(1, len(cumulative) + 1),
            cumulative,
            "o-",
            linewidth=2,
            color="blue",
        )

        # Add threshold lines
        ax2.axhline(
            y=0.95,
            color="green",
            linestyle="--",
            alpha=0.7,
            label="95% Explained Variance",
        )
        ax2.axhline(
            y=0.8,
            color="orange",
            linestyle="--",
            alpha=0.7,
            label="80% Explained Variance",
        )

        # Get the number of components needed for these thresholds
        components_80 = (
            np.where(cumulative >= 0.8)[0][0] + 1
            if any(cumulative >= 0.8)
            else n_components
        )
        components_95 = (
            np.where(cumulative >= 0.95)[0][0] + 1
            if any(cumulative >= 0.95)
            else n_components
        )

        # Add vertical lines to show number of components needed
        ax2.axvline(x=components_80, color="orange", linestyle="--", alpha=0.7)
        ax2.axvline(x=components_95, color="green", linestyle="--", alpha=0.7)

        # Add annotations
        for i, percentage in enumerate(cumulative):
            ax2.annotate(
                f"{percentage:.2%}",
                xy=(i + 1, percentage),
                xytext=(3, 0),
                textcoords="offset points",
                ha="left",
                va="center",
                fontsize=8,
            )

        ax2.text(
            components_80,
            0.4,
            f"{components_80} components\nfor 80% variance",
            ha="center",
            va="center",
            bbox=dict(facecolor="white", alpha=0.7),
        )
        ax2.text(
            components_95,
            0.6,
            f"{components_95} components\nfor 95% variance",
            ha="center",
            va="center",
            bbox=dict(facecolor="white", alpha=0.7),
        )

        ax2.set_xlabel("Number of Components")
        ax2.set_ylabel("Cumulative Explained Variance")
        ax2.set_title("Cumulative Explained Variance")
        ax2.set_xticks(range(1, len(cumulative) + 1))
        ax2.grid(alpha=0.3)

        fig.tight_layout()
        fig.savefig(os.path.join(plots_dir, "pca_scree_plot.png"))
        plt.close(fig)

        print(
            f"Created scree plot showing variance explained by {n_components} components"
        )
        print(f"Components needed for 80% variance: {components_80}")
        print(f"Components needed for 95% variance: {components_95}")

        # Create loadings heatmap
        loadings = pca_full.components_
        loading_df = pd.DataFrame(
            loadings.T,
            columns=[f"PC{i+1}" for i in range(loadings.shape[0])],
            index=X.columns,
        )

        # Plot heatmap of loadings for the first few components
        fig, ax = plt.subplots(figsize=(12, 10))
        num_comp_to_show = min(5, loadings.shape[0])
        sns.heatmap(
            loading_df.iloc[:, :num_comp_to_show],
            cmap="coolwarm",
            annot=True,
            fmt=".2f",
            linewidths=0.5,
            center=0,
            cbar_kws={"label": "Loading Coefficient"},
            ax=ax,
        )
        plt.title("PCA Loading Coefficients Heatmap")
        fig.tight_layout()
        fig.savefig(os.path.join(plots_dir, "pca_loadings_heatmap.png"))
        plt.close(fig)

        # Save loading coefficients to a CSV file
        loading_csv = os.path.join(plots_dir, "pca_loadings.csv")
        loading_df.to_csv(loading_csv)

        # Create 2D PCA plot
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(X_scaled)

        # Create a DataFrame with the PCA results
        pca_df = pd.DataFrame(
            data=pca_result,
            columns=["Principal Component 1", "Principal Component 2"],
        )

        # Add the original type information for coloring
        pca_df["original_type"] = df["original_type"].values

        # Calculate the explained variance
        explained_variance = pca.explained_variance_ratio_

        # Create PCA plot with clusters
        fig, ax = plt.subplots(figsize=(12, 10))
        scatter = sns.scatterplot(
            data=pca_df,
            x="Principal Component 1",
            y="Principal Component 2",
            hue="original_type",
            style="original_type",
            s=100,
            alpha=0.8,
            ax=ax,
        )

        # Add ellipses to show clusters
        for i, orig_type in enumerate(pca_df["original_type"].unique()):
            subset = pca_df[pca_df["original_type"] == orig_type]

            # Calculate the mean and covariance
            if (
                len(subset) >= 2
            ):  # Need at least 2 points to calculate covariance
                mean_x = subset["Principal Component 1"].mean()
                mean_y = subset["Principal Component 2"].mean()
                cov = np.cov(
                    subset["Principal Component 1"],
                    subset["Principal Component 2"],
                )

                # Calculate eigenvalues and eigenvectors
                evals, evecs = np.linalg.eig(cov)

                # Sort eigenvalues in decreasing order
                sort_indices = np.argsort(evals)[::-1]
                evals = evals[sort_indices]
                evecs = evecs[:, sort_indices]

                # Angle in degrees
                angle = np.degrees(np.arctan2(evecs[1, 0], evecs[0, 0]))

                # 95% confidence ellipse
                width, height = 2 * np.sqrt(
                    5.991 * evals
                )  # 5.991 is chi-square value for 95% confidence and 2 DOF

                # Get color from the scatter plot legend
                try:
                    # First try to get color from the scatter plot legend
                    color = (
                        scatter.get_legend().legendHandles[i].get_facecolor()
                    )
                except:
                    # Fallback to default colors
                    color = plt.cm.tab10(i)

                # Create ellipse
                ellipse = Ellipse(
                    xy=(mean_x, mean_y),
                    width=width,
                    height=height,
                    angle=angle,
                    edgecolor=color,
                    fc="none",
                    lw=2,
                    alpha=0.3,
                )
                ax.add_patch(ellipse)

                # Add text label in the center of each cluster
                ax.text(
                    mean_x,
                    mean_y,
                    orig_type,
                    horizontalalignment="center",
                    verticalalignment="center",
                    size=12,
                    weight="bold",
                    alpha=0.7,
                )

        # Add title and labels
        ax.set_title(
            f"PCA of Amplicon Features with Cluster Visualization{title_suffix}"
        )
        ax.set_xlabel(
            f"Principal Component 1 ({explained_variance[0]:.2%} variance)"
        )
        ax.set_ylabel(
            f"Principal Component 2 ({explained_variance[1]:.2%} variance)"
        )
        ax.axhline(y=0, color="gray", linestyle="--", alpha=0.3)
        ax.axvline(x=0, color="gray", linestyle="--", alpha=0.3)
        ax.grid(True, alpha=0.3)

        # Save the plot
        fig.tight_layout()
        fig.savefig(os.path.join(plots_dir, "pca_clusters.png"))
        plt.close(fig)

        return pca, pca_result, loading_df, components_80, components_95

    @staticmethod
    def plot_confusion_matrix(
        y_true, y_pred, output_path, labels=None, title="Confusion Matrix"
    ):
        """Create and save a confusion matrix plot."""
        if labels is None:
            labels = ["BFB", "ecDNA"]

        cm = confusion_matrix(y_true, y_pred)
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.heatmap(
            cm,
            annot=True,
            fmt="d",
            cmap="Blues",
            xticklabels=labels,
            yticklabels=labels,
            ax=ax,
        )
        ax.set_xlabel("Predicted")
        ax.set_ylabel("True")
        ax.set_title(title)
        fig.tight_layout()
        fig.savefig(output_path)
        plt.close(fig)

        return cm

    @staticmethod
    def plot_feature_importance(
        importance_df, output_path, title="Feature Importance", top_n=10
    ):
        """Create and save a feature importance bar plot."""
        sorted_df = importance_df.sort_values(
            "Importance", ascending=False
        ).head(top_n)

        fig, ax = plt.subplots(figsize=(12, 8))
        sns.barplot(x="Importance", y="Feature", data=sorted_df, ax=ax)
        ax.set_title(title)
        fig.tight_layout()
        fig.savefig(output_path)
        plt.close(fig)

        return sorted_df

    @staticmethod
    def plot_scatter(
        x,
        y,
        data,
        output_path,
        title="Scatter Plot",
        xlabel="X",
        ylabel="Y",
        hue=None,
        style=None,
    ):
        """Create and save a scatter plot."""
        fig, ax = plt.subplots(figsize=(12, 10))
        scatter = sns.scatterplot(
            data=data, x=x, y=y, hue=hue, style=style, s=100, alpha=0.7, ax=ax
        )

        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend(title=(hue if hue else "Legend"), loc="best")
        fig.tight_layout()
        fig.savefig(output_path)
        plt.close(fig)

        return scatter


def plot_sequence_feature_distributions(df, plots_dir):
    """
    Generate visualizations for sequence-related features to better understand
    the distribution and relationships between these properties in different amplicon types.

    Args:
        df: DataFrame containing the amplicon data with sequence features
        plots_dir: Directory to save the plots to
    """
    print(f"Generating sequence feature plots...")

    # Create a subdirectory for sequence feature plots
    seq_plots_dir = os.path.join(plots_dir, "sequence_features")
    os.makedirs(seq_plots_dir, exist_ok=True)

    try:
        # Filter sequence-related features for analysis
        seq_features = [
            "mean_seq_cn",
            "max_seq_cn",
            "std_seq_cn",
            "mean_seq_coverage",
            "max_seq_coverage",
            "std_seq_coverage",
            "mean_bp_cn",
            "max_bp_cn",
            "std_bp_cn",
            "cn_coverage_ratio",
            "bp_seq_cn_ratio",
        ]

        # Create boxplots comparing feature distributions by amplicon type
        plt.figure(figsize=(15, 10))
        for i, feature in enumerate(seq_features):
            plt.subplot(3, 4, i + 1)
            sns.boxplot(x="original_type", y=feature, data=df)
            plt.title(feature)
            plt.xticks(rotation=90)

        plt.tight_layout()
        plt.savefig(os.path.join(seq_plots_dir, "seq_features_boxplot.png"))
        plt.close()

        # Create correlation heatmap for sequence features
        plt.figure(figsize=(12, 10))
        seq_corr = df[seq_features].corr()
        sns.heatmap(seq_corr, annot=True, cmap="coolwarm", fmt=".2f")
        plt.title("Correlation Matrix of Sequence Features")
        plt.tight_layout()
        plt.savefig(os.path.join(seq_plots_dir, "seq_features_correlation.png"))
        plt.close()

        # Create scatter plots comparing important sequence features
        plt.figure(figsize=(15, 15))
        feature_pairs = [
            ("mean_seq_cn", "mean_seq_coverage"),
            ("max_seq_cn", "max_seq_coverage"),
            ("mean_bp_cn", "mean_seq_cn"),
            ("cn_coverage_ratio", "bp_seq_cn_ratio"),
        ]

        for i, (x_feat, y_feat) in enumerate(feature_pairs):
            plt.subplot(2, 2, i + 1)
            sns.scatterplot(x=x_feat, y=y_feat, hue="original_type", data=df)
            plt.title(f"{x_feat} vs {y_feat}")
            if i == 0:
                plt.legend(
                    title="Sample Type",
                    bbox_to_anchor=(1.05, 1),
                    loc="upper left",
                )
            else:
                plt.legend([], [], frameon=False)

        plt.tight_layout()
        plt.savefig(os.path.join(seq_plots_dir, "seq_features_scatter.png"))
        plt.close()

        # Create distribution plots for key sequence features
        plt.figure(figsize=(15, 12))
        for i, feature in enumerate(
            seq_features[:6]
        ):  # Just show first 6 features
            plt.subplot(2, 3, i + 1)
            for amplicon_type in df["original_type"].unique():
                subset = df[df["original_type"] == amplicon_type]
                sns.kdeplot(subset[feature], label=amplicon_type)
            plt.title(f"Distribution of {feature}")
            if i == 0:
                plt.legend(title="Sample Type")

        plt.tight_layout()
        plt.savefig(
            os.path.join(seq_plots_dir, "seq_features_distribution.png")
        )
        plt.close()

        print(f"Created sequence feature plots in {seq_plots_dir}")

    except Exception as e:
        print(f"Error generating sequence feature plots: {e}")
        traceback.print_exc()


# Add the EntropyAnalyzer class
class EntropyAnalyzer:
    """Class for performing entropy analysis and creating decision trees."""

    def __init__(self, df, plots_dir, title_suffix=""):
        """Initialize with the data and output directories."""
        self.df = df
        self.plots_dir = plots_dir
        self.title_suffix = title_suffix
        self.entropy_dir = os.path.join(plots_dir, "entropy_analysis")
        os.makedirs(self.entropy_dir, exist_ok=True)

        # Prepare data
        self.X = df.drop(["type", "original_type", "filepath"], axis=1)
        self.y = df["type"]
        self.feature_names = self.X.columns.tolist()
        # Ensure we have a numeric target for mutual info calculation
        # (0 for BFB, 1 for ecDNA)
        self.y_numeric = np.array([0 if y == "BFB" else 1 for y in self.y])

    def calculate_feature_entropy(self):
        """Calculate the Shannon entropy for each feature."""
        print("Calculating feature entropy...")
        entropy_values = {}

        for feature in self.feature_names:
            # Get feature values
            values = self.X[feature].values

            # Bin the values for continuous features
            if self.X[feature].dtype in [
                "float64",
                "float32",
                "int64",
                "int32",
            ]:
                # Use Freedman-Diaconis rule for bin width
                q75, q25 = np.percentile(values, [75, 25])
                iqr = q75 - q25
                bin_width = (
                    2 * iqr * (len(values) ** (-1 / 3)) if iqr > 0 else 0.5
                )
                n_bins = (
                    int(np.ceil((np.max(values) - np.min(values)) / bin_width))
                    if bin_width > 0
                    else 10
                )
                n_bins = min(
                    max(10, n_bins), 50
                )  # Limit bins between 10 and 50

                # Use histogram to bin the data
                hist, bin_edges = np.histogram(values, bins=n_bins)
                probabilities = hist / len(values)

                # Calculate entropy (only on non-zero probabilities)
                probabilities = probabilities[probabilities > 0]
                entropy = -np.sum(probabilities * np.log2(probabilities))
            else:
                # For categorical variables, use value counts directly
                value_counts = self.X[feature].value_counts(normalize=True)
                entropy = stats.entropy(value_counts, base=2)

            entropy_values[feature] = entropy

        return entropy_values

    def plot_feature_entropy(self, entropy_values):
        """Create a bar plot of feature entropy values."""
        print("Creating feature entropy plot...")
        # Sort features by entropy
        sorted_features = sorted(
            entropy_values.items(), key=lambda x: x[1], reverse=True
        )
        features = [x[0] for x in sorted_features]
        entropies = [x[1] for x in sorted_features]

        # Create the plot
        fig, ax = plt.subplots(figsize=(14, 8))
        bars = ax.bar(range(len(features)), entropies, color="skyblue")

        # Add a color gradient based on entropy values
        norm = mcolors.Normalize(min(entropies), max(entropies))
        sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=norm)
        sm.set_array([])

        for i, bar in enumerate(bars):
            bar.set_color(plt.cm.viridis(norm(entropies[i])))

        # Add the colorbar with the correct axes reference
        fig.colorbar(sm, ax=ax, label="Entropy (bits)")

        # Add labels and title
        ax.set_xlabel("Feature")
        ax.set_ylabel("Entropy (bits)")
        ax.set_title(f"Feature Entropy{self.title_suffix}")
        ax.set_xticks(range(len(features)))
        ax.set_xticklabels(features, rotation=90)
        ax.grid(axis="y", alpha=0.3)

        # Add value annotations
        for i, v in enumerate(entropies):
            ax.text(
                i,
                v + 0.05,
                f"{v:.2f}",
                ha="center",
                va="bottom",
                fontsize=8,
                rotation=90,
            )

        fig.tight_layout()

        # Save the plot
        output_file = os.path.join(
            self.entropy_dir,
            f'feature_entropy{self.title_suffix.replace(" ", "_").lower()}.png',
        )
        fig.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close(fig)

        print(f"Saved feature entropy plot to {output_file}")
        return output_file

    def calculate_mutual_information(self):
        """Calculate mutual information between features and target."""
        print("Calculating mutual information...")
        mutual_info = mutual_info_classif(self.X, self.y_numeric)

        # Create a DataFrame for sorting
        mi_df = pd.DataFrame(
            {"Feature": self.feature_names, "Mutual Information": mutual_info}
        )
        mi_df = mi_df.sort_values("Mutual Information", ascending=False)

        return mi_df

    def plot_mutual_information(self, mi_df):
        """Create a bar plot of mutual information between features and target."""
        # Create the plot
        fig, ax = plt.subplots(figsize=(14, 8))
        bars = ax.bar(
            range(len(mi_df)), mi_df["Mutual Information"], color="skyblue"
        )

        # Add a color gradient based on MI values
        norm = mcolors.Normalize(
            min(mi_df["Mutual Information"]), max(mi_df["Mutual Information"])
        )
        sm = plt.cm.ScalarMappable(cmap=plt.cm.plasma, norm=norm)
        sm.set_array([])

        for i, bar in enumerate(bars):
            bar.set_color(
                plt.cm.plasma(norm(mi_df["Mutual Information"].iloc[i]))
            )

        fig.colorbar(sm, ax=ax, label="Mutual Information")

        # Add labels and title
        ax.set_xlabel("Feature")
        ax.set_ylabel("Mutual Information with Target")
        ax.set_title(
            f"Feature Mutual Information with Target{self.title_suffix}"
        )
        ax.set_xticks(range(len(mi_df)))
        ax.set_xticklabels(mi_df["Feature"], rotation=90)
        ax.grid(axis="y", alpha=0.3)

        # Add value annotations
        for i, v in enumerate(mi_df["Mutual Information"]):
            ax.text(
                i,
                v + 0.002,
                f"{v:.4f}",
                ha="center",
                va="bottom",
                fontsize=8,
                rotation=90,
            )

        fig.tight_layout()

        # Save the plot
        output_file = os.path.join(
            self.entropy_dir,
            f'mutual_information{self.title_suffix.replace(" ", "_").lower()}.png',
        )
        fig.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close(fig)

        print(f"Saved mutual information plot to {output_file}")
        return output_file

    def plot_entropy_decision_tree(
        self, X_top, top_features, entropy_values, max_depth=4
    ):
        """Create a decision tree visualization with node colors representing entropy."""
        print("Creating decision tree with entropy visualization...")
        # Train a decision tree
        dtree = DecisionTreeClassifier(max_depth=max_depth, random_state=42)
        dtree.fit(X_top, self.y)

        # Calculate class probabilities at each node
        n_nodes = dtree.tree_.node_count
        children_left = dtree.tree_.children_left
        children_right = dtree.tree_.children_right
        feature = dtree.tree_.feature
        threshold = dtree.tree_.threshold

        # Calculate node entropy
        node_entropy = np.zeros(n_nodes)
        node_samples = dtree.tree_.weighted_n_node_samples
        node_impurity = dtree.tree_.impurity  # Gini impurity by default

        # For leaf nodes, use class probabilities to calculate entropy
        for i in range(n_nodes):
            if children_left[i] == children_right[i]:  # Leaf node
                # Convert Gini impurity to entropy for visualization consistency
                node_entropy[i] = 1 - node_impurity[i]  # Approximation
            else:
                # For internal nodes, use the feature's entropy
                selected_feature = feature[i]
                if selected_feature >= 0:  # Check if it's a valid feature index
                    feature_name = top_features[selected_feature]
                    node_entropy[i] = entropy_values.get(feature_name, 0)

        # Normalize entropy for coloring
        max_entropy = max(max(node_entropy), 0.001)  # Avoid division by zero
        normalized_entropy = node_entropy / max_entropy

        # Create color map for nodes
        node_colors = []
        for e in normalized_entropy:
            if e <= 0.001:  # Leaf nodes or very low entropy
                node_colors.append("#f0f0f0")  # Light gray
            else:
                # Use a color gradient from blue (low entropy) to red (high entropy)
                color = plt.cm.coolwarm(e)
                node_colors.append(mcolors.rgb2hex(color))

        # Calculate information gain for edges
        edge_info_gain = {}
        for i in range(n_nodes):
            if children_left[i] != children_right[i]:  # Internal node
                # Information gain = parent impurity - weighted average of children impurity
                parent_impurity = node_impurity[i]
                left_child = children_left[i]
                right_child = children_right[i]
                left_impurity = node_impurity[left_child]
                right_impurity = node_impurity[right_child]
                left_weight = node_samples[left_child] / node_samples[i]
                right_weight = node_samples[right_child] / node_samples[i]
                info_gain = parent_impurity - (
                    left_weight * left_impurity + right_weight * right_impurity
                )

                edge_info_gain[(i, left_child)] = info_gain
                edge_info_gain[(i, right_child)] = info_gain

        # Normalize information gain for edge thickness
        if edge_info_gain:
            max_info_gain = max(edge_info_gain.values())
            min_edge_width = 1.0
            max_edge_width = 4.0
            for edge, gain in edge_info_gain.items():
                edge_info_gain[edge] = min_edge_width + (
                    gain / max_info_gain
                ) * (max_edge_width - min_edge_width)

        # Create figure for matplotlib plot
        fig, ax = plt.subplots(figsize=(20, 10))

        # Plot the tree
        plot_tree(
            dtree,
            feature_names=top_features,
            class_names=["BFB", "ecDNA"],
            filled=True,
            rounded=True,
            fontsize=10,
            node_ids=True,
            precision=2,
            ax=ax,
        )

        ax.set_title(f"Decision Tree with Entropy Coloring{self.title_suffix}")

        # Save the matplotlib version
        output_file_mpl = os.path.join(
            self.entropy_dir,
            f'decision_tree_entropy_mpl{self.title_suffix.replace(" ", "_").lower()}.png',
        )
        fig.savefig(output_file_mpl, dpi=300, bbox_inches="tight")
        plt.close(fig)

        # Try to create GraphViz version if available
        output_file_dot = None
        try:
            # Create a more detailed GraphViz version
            dot_data = StringIO()
            export_graphviz(
                dtree,
                out_file=dot_data,
                feature_names=top_features,
                class_names=["BFB", "ecDNA"],
                filled=True,
                rounded=True,
                special_characters=True,
                leaves_parallel=True,
                impurity=True,
                proportion=True,
                node_ids=True,
            )

            # Get the GraphViz representation
            graph = pydotplus.graph_from_dot_data(dot_data.getvalue())

            # Color the nodes based on entropy
            for i, node in enumerate(graph.get_nodes()):
                if node.get_name() not in ("node", "edge"):
                    try:
                        node_id = int(node.get_name().split()[0])
                        if node_id < len(node_colors):
                            node.set_fillcolor(node_colors[node_id])
                    except ValueError:
                        # Skip nodes that don't have a numeric ID
                        continue

            # Set edge thickness based on information gain
            for edge in graph.get_edges():
                try:
                    source = int(edge.get_source().split()[0])
                    target = int(edge.get_destination().split()[0])
                    edge_key = (source, target)
                    if edge_key in edge_info_gain:
                        edge.set_penwidth(str(edge_info_gain[edge_key]))
                except ValueError:
                    # Skip edges that don't have numeric source/target
                    continue

            # Save the GraphViz version
            output_file_dot = os.path.join(
                self.entropy_dir,
                f'decision_tree_entropy{self.title_suffix.replace(" ", "_").lower()}.png',
            )
            graph.write_png(output_file_dot)
            print(
                f"Created GraphViz decision tree visualization: {output_file_dot}"
            )
        except Exception as e:
            print(f"GraphViz visualization failed: {e}")
            print("Using matplotlib visualization only.")
            output_file_dot = output_file_mpl

        return output_file_mpl, output_file_dot


# Add function to handle entropy and decision tree plots generation
def generate_entropy_and_tree_plots(df, plots_dir, title_suffix=""):
    """
    Generate entropy and decision tree plots for the given dataset.

    Args:
        df: DataFrame containing the data
        plots_dir: Directory to save the plots to
        title_suffix: Optional suffix for plot titles
    """
    print(f"\nGenerating entropy and decision tree plots{title_suffix}...")

    # Initialize return values to None
    entropy_plot = None
    mi_plot = None
    tree_plot_dot = None
    entropy_mi_plot = None

    try:
        # Create an instance of the EntropyAnalyzer class
        analyzer = EntropyAnalyzer(df, plots_dir, title_suffix)

        # Explicitly ensure y_numeric is set for mutual information calculation
        analyzer.y_numeric = np.array(
            [0 if y == "BFB" else 1 for y in analyzer.y]
        )

        # Calculate feature entropy
        print("Calculating feature entropy...")
        entropy_values = analyzer.calculate_feature_entropy()

        # Plot feature entropy
        print("Creating feature entropy plot...")
        entropy_plot = analyzer.plot_feature_entropy(entropy_values)
        print(f"Saved feature entropy plot to {entropy_plot}")

        # Calculate mutual information - handle potential errors
        print("Calculating mutual information...")
        try:
            mi_df = analyzer.calculate_mutual_information()
            mi_plot = analyzer.plot_mutual_information(mi_df)
            print(f"Saved mutual information plot to {mi_plot}")

            # Select top features for decision tree
            top_features = mi_df.head(15)["Feature"].tolist()
            X_top = analyzer.X[top_features]

            # Plot decision tree
            print("Creating decision tree with entropy visualization...")
            tree_plot_mpl, tree_plot_dot = analyzer.plot_entropy_decision_tree(
                X_top, top_features, entropy_values, max_depth=4
            )
            print(f"Saved decision tree plots to:\n{tree_plot_mpl}")
            if (
                tree_plot_dot != tree_plot_mpl
            ):  # If GraphViz version was created
                print(f"{tree_plot_dot}")
        except Exception as e:
            print(f"Error in mutual information or decision tree: {e}")
            traceback.print_exc()
            # Continue with entropy vs MI plot even if MI calculation failed

        # Create entropy vs. mutual information scatter plot
        try:
            if "mi_df" in locals() and mi_df is not None:
                fig, ax = plt.subplots(figsize=(12, 8))

                # Get mutual information for all features
                mi_values = dict(
                    zip(mi_df["Feature"], mi_df["Mutual Information"])
                )

                # Prepare data for scatter plot
                scatter_data = []
                for feature, entropy in entropy_values.items():
                    if feature in mi_values:
                        scatter_data.append(
                            {
                                "Feature": feature,
                                "Entropy": entropy,
                                "Mutual Information": mi_values[feature],
                            }
                        )

                scatter_df = pd.DataFrame(scatter_data)

                # Sort by mutual information for coloring
                scatter_df = scatter_df.sort_values(
                    "Mutual Information", ascending=False
                )

                # Create scatter plot
                scatter = ax.scatter(
                    scatter_df["Entropy"],
                    scatter_df["Mutual Information"],
                    c=plt.cm.viridis(np.linspace(0, 1, len(scatter_df))),
                    s=80,
                    alpha=0.7,
                )

                # Add feature labels
                for i, row in scatter_df.iterrows():
                    ax.annotate(
                        row["Feature"],
                        xy=(row["Entropy"], row["Mutual Information"]),
                        xytext=(5, 0),
                        textcoords="offset points",
                        fontsize=8,
                        alpha=0.8,
                    )

                ax.set_xlabel("Feature Entropy (bits)")
                ax.set_ylabel("Mutual Information with Target")
                ax.set_title(
                    f"Feature Entropy vs. Mutual Information{title_suffix}"
                )
                ax.grid(alpha=0.3)

                # Add trend line
                if len(scatter_df) > 1:
                    z = np.polyfit(
                        scatter_df["Entropy"],
                        scatter_df["Mutual Information"],
                        1,
                    )
                    p = np.poly1d(z)
                    ax.plot(
                        np.sort(scatter_df["Entropy"]),
                        p(np.sort(scatter_df["Entropy"])),
                        "r--",
                        alpha=0.8,
                        label=f"Trend: y={z[0]:.4f}x+{z[1]:.4f}",
                    )
                    ax.legend()

                fig.tight_layout()

                # Save the plot
                entropy_mi_plot = os.path.join(
                    analyzer.entropy_dir,
                    f'entropy_vs_mi{title_suffix.replace(" ", "_").lower()}.png',
                )
                fig.savefig(entropy_mi_plot, dpi=300, bbox_inches="tight")
                plt.close(fig)
                print(
                    f"Saved entropy vs. mutual information plot to {entropy_mi_plot}"
                )
        except Exception as e:
            print(f"Error creating entropy vs. mutual information plot: {e}")
            traceback.print_exc()

        return entropy_plot, mi_plot, tree_plot_dot, entropy_mi_plot

    except Exception as e:
        print(f"Error generating entropy and decision tree plots: {e}")
        traceback.print_exc()
        return entropy_plot, mi_plot, tree_plot_dot, entropy_mi_plot


# Add a special EntropyAnalyzer class for subsets that explicitly handles string target values
class EntropyAnalyzerForSubsets(EntropyAnalyzer):
    """Special version of EntropyAnalyzer that properly handles string target values for subset analyses."""

    def __init__(self, df, plots_dir, title_suffix=""):
        """Initialize with the data and output directories."""
        # Initialize all the same variables as the parent
        self.df = df
        self.plots_dir = plots_dir
        self.title_suffix = title_suffix
        self.entropy_dir = os.path.join(plots_dir, "entropy_analysis")
        os.makedirs(self.entropy_dir, exist_ok=True)

        # Prepare data
        self.X = df.drop(["type", "original_type", "filepath"], axis=1)
        self.y = df["type"]
        self.feature_names = self.X.columns.tolist()

        # Explicitly convert string labels to numeric values (0 for BFB, 1 for ecDNA)
        # This ensures y_numeric is correctly set and of numeric type
        self.y_numeric = np.zeros(len(self.y), dtype=int)
        for i, val in enumerate(self.y):
            if str(val).strip() != "BFB":
                self.y_numeric[i] = 1

    def calculate_mutual_information(self):
        """Calculate mutual information between features and target."""
        print("Calculating mutual information...")
        try:
            # Use explicit integer array for mutual information calculation
            y_values = np.zeros(len(self.y), dtype=int)
            for i, val in enumerate(self.y):
                if str(val).strip() != "BFB":
                    y_values[i] = 1

            # Now calculate mutual information with the integer array
            mutual_info = mutual_info_classif(self.X, y_values)

            # Create a DataFrame for sorting
            mi_df = pd.DataFrame(
                {
                    "Feature": self.feature_names,
                    "Mutual Information": mutual_info,
                }
            )
            mi_df = mi_df.sort_values("Mutual Information", ascending=False)

            return mi_df
        except Exception as e:
            print(f"Error in mutual information calculation: {e}")
            traceback.print_exc()
            # Return empty DataFrame in case of error
            return pd.DataFrame({"Feature": [], "Mutual Information": []})


if __name__ == "__main__":
    main()
