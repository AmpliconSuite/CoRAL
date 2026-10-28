#!/usr/bin/env python3
from __future__ import annotations

import functools
import importlib.resources
import io
import logging
import math
import os
import pathlib
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, DefaultDict, Literal, Optional

import colorama
import intervaltree
import matplotlib as mpl
import matplotlib.colors as mcolors
import numpy as np
import typer

from coral import datatypes
from coral.breakpoint import (
    breakpoint_utilities,  # type: ignore[import-untyped]
)
from coral.breakpoint.breakpoint_graph import BreakpointGraph
from coral.breakpoint.parse_graph import parse_breakpoint_graph
from coral.datatypes import Interval
from coral.summary.parsing import parse_cycle_file

mpl.use("Agg")
import matplotlib.pyplot as plt
import pysam
from matplotlib import gridspec, ticker
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import Arc, FancyArrowPatch, Patch, Rectangle
from matplotlib.text import Text
from matplotlib.textpath import TextPath
from pylab import rcParams  # type: ignore[import-untyped]

from coral import (
    bam_types,
    core_types,
    core_utils,
    supplemental_data,
)

rcParams["pdf.fonttype"] = 42


logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class CyclePlotStyle:
    """Visual defaults shared by standalone and combined cycle plots."""

    segment_line_width: float
    connection_line_width: float
    cycle_color: str
    path_color: str
    segment_face_color: str
    segment_edge_color: str


DISCORDANT_EDGE_COLORS = {
    "+-": "red",
    "++": "magenta",
    "-+": (139 / 256.0, 69 / 256.0, 19 / 256.0),
    "--": "teal",
    "interchromosomal": "blue",
}
AA_DISCORDANT_EDGE_MIN_MAX_READ_COUNT = 4.0
DISCORDANT_EDGE_MIN_LINE_WIDTH = 0.30
DISCORDANT_EDGE_MAX_LINE_WIDTH = 4.00
DEFAULT_COVERAGE_HEADROOM = 1.25
DEFAULT_ARC_TOP_PADDING = 0.05
MIN_ARC_APEX_CN_SCALE = 1.25
MAX_ARC_APEX_CN_SCALE = 2.8
DEFAULT_GENE_FONT_SIZE = 12.0
DEFAULT_PLOT_FONT_SIZE = 18.0
DEFAULT_LEGEND_FONT_SIZE = 10.0
MAX_TEXT_LAYOUT_SCALE = 2.0
DEFAULT_GRAPH_WIDTH = 12.0
DEFAULT_GRAPH_HEIGHT = 5.0
DEFAULT_CYCLE_WIDTH = 12.0
DEFAULT_INTERVAL_OFFSET = 0.10
DEFAULT_MIN_COORD_WIDTH = 0.0
DEFAULT_MIN_INTERVAL_DISPLAY_WIDTH = 0.005
DEFAULT_CYCLE_PLOT_STYLE = CyclePlotStyle(
    segment_line_width=1.8,
    connection_line_width=1.5,
    cycle_color="#5F78B5",
    path_color="#D55E00",
    segment_face_color="#F2D2A2",
    segment_edge_color="#3F3F3F",
)
CYCLE_SEGMENT_LINE_WIDTH = DEFAULT_CYCLE_PLOT_STYLE.segment_line_width
CYCLE_CONNECTION_LINE_WIDTH = DEFAULT_CYCLE_PLOT_STYLE.connection_line_width
CYCLE_COLOR = DEFAULT_CYCLE_PLOT_STYLE.cycle_color
PATH_COLOR = DEFAULT_CYCLE_PLOT_STYLE.path_color
CYCLE_SEGMENT_FACE_COLOR = DEFAULT_CYCLE_PLOT_STYLE.segment_face_color
CYCLE_SEGMENT_EDGE_COLOR = DEFAULT_CYCLE_PLOT_STYLE.segment_edge_color
COMBINED_AXIS_LEFT = 0.23
COMBINED_AXIS_RIGHT = 0.88
PLOT_DATA_HEIGHT_RATIO = 8.0
GENE_TRACK_HEIGHT_RATIO = 3.0
CYCLE_GENE_TRACK_BASE_HEIGHT_INCHES = 0.8
CYCLE_GENE_TRACK_HEIGHT_PER_LANE_INCHES = 0.4
CYCLE_GENE_TRACK_MAX_HEIGHT_INCHES = 3.0
CYCLE_GENE_TRACK_MAX_FIGURE_FRACTION = 0.45
GRAPH_MAIN_AXIS_BOTTOM_WITH_GENES = 0.30 + (0.88 - 0.30) * (
    GENE_TRACK_HEIGHT_RATIO / (PLOT_DATA_HEIGHT_RATIO + GENE_TRACK_HEIGHT_RATIO)
)
GRAPH_MAIN_AXIS_BOTTOM_WITHOUT_GENES = 0.12
GRAPH_BASELINE_LABEL_MASK_HEIGHT = 0.05
GENE_LABEL_HORIZONTAL_PADDING_POINTS = 6.0
GENE_LABEL_LANE_GAP_FRACTION = 0.004
GENOMIC_COORDINATE_FONT_SCALE = 0.68
CHROMOSOME_GUIDE_LINE_WIDTH = 1.25
INTERVAL_GUIDE_LINE_WIDTH = 1.0
GUIDE_ALPHA = 0.75


@dataclass(frozen=True)
class GraphAxisLimits:
    """Final limits for the overlaid coverage and CN axes."""

    coverage_ymax: float
    cn_ymax: float
    expansion_factor: float


@dataclass(frozen=True)
class CycleSelection:
    cycle_id: int
    expected_is_cyclic: bool


@dataclass(frozen=True)
class GenomicAxisLayout:
    """One reference-coordinate transform shared by aligned plot panels."""

    total_genomic_length: float
    interval_starts: dict[str, list[float]]
    axis_max: float
    chromosome_guides: tuple[float, ...]
    interval_guides: tuple[float, ...]


def align_cn_axis_ticks(
    coverage_axis: Axes,
    cn_axis: Axes,
    axis_limits: GraphAxisLimits,
) -> None:
    """Place CN ticks at the same heights as coverage-axis ticks."""
    coverage_ticks = [
        float(tick)
        for tick in coverage_axis.get_yticks()
        if 0 <= tick <= axis_limits.coverage_ymax
    ]
    cn_per_coverage = axis_limits.cn_ymax / axis_limits.coverage_ymax
    cn_axis.set_yticks([tick * cn_per_coverage for tick in coverage_ticks])


def build_genomic_axis_layout(
    intervals_by_chrom: dict[str, list[datatypes.Interval]],
    interval_offset: float,
) -> GenomicAxisLayout:
    """Map genomic intervals to a common display axis."""
    sorted_chromosomes = breakpoint_utilities.sort_chrom_names(
        intervals_by_chrom.keys()
    )
    interval_count = sum(
        len(intervals_by_chrom[chrom]) for chrom in sorted_chromosomes
    )
    total_genomic_length = float(
        sum(
            len(interval)
            for chrom in sorted_chromosomes
            for interval in intervals_by_chrom[chrom]
        )
    )
    if interval_count == 0 or total_genomic_length <= 0:
        raise ValueError("cannot build a genomic axis without intervals")

    margin = get_interval_margin(interval_count, interval_offset)
    interval_starts: dict[str, list[float]] = {}
    chromosome_guides: list[float] = []
    interval_guides: list[float] = []
    x = margin
    for chromosome_index, chrom in enumerate(sorted_chromosomes):
        interval_starts[chrom] = []
        for interval_index, interval in enumerate(intervals_by_chrom[chrom]):
            interval_starts[chrom].append(x)
            if interval_index > 0:
                interval_guides.append(x - margin * 0.5)
            elif chromosome_index > 0:
                chromosome_guides.append(x - margin * 0.5)
            x += (interval.end - interval.start) * 100.0 / total_genomic_length
            x += margin

    return GenomicAxisLayout(
        total_genomic_length=total_genomic_length,
        interval_starts=interval_starts,
        axis_max=100.0 + (interval_count + 1) * margin,
        chromosome_guides=tuple(chromosome_guides),
        interval_guides=tuple(interval_guides),
    )


def genomic_position_to_axis(
    position: int,
    chrom: str,
    intervals_by_chrom: dict[str, list[datatypes.Interval]],
    layout: GenomicAxisLayout,
) -> float:
    """Project one reported genomic coordinate onto a plot axis."""
    for interval_index, interval in enumerate(intervals_by_chrom[chrom]):
        if interval.start <= position <= interval.end:
            return layout.interval_starts[chrom][interval_index] + (
                (position - interval.start)
                * 100.0
                / layout.total_genomic_length
            )
    raise ValueError(f"coordinate {chrom}:{position} is outside the plot axis")


def draw_interval_guides(
    data_axis: Axes,
    coordinate_axis: Axes,
    layout: GenomicAxisLayout,
) -> None:
    """Draw interval guides through the data and bottom coordinate axis."""
    guide_groups = (
        (
            layout.chromosome_guides,
            "--",
            "black",
            CHROMOSOME_GUIDE_LINE_WIDTH,
        ),
        (
            layout.interval_guides,
            ":",
            "0.25",
            INTERVAL_GUIDE_LINE_WIDTH,
        ),
    )
    for positions, linestyle, color, linewidth in guide_groups:
        for position in positions:
            for axis in (data_axis, coordinate_axis):
                axis.axvline(
                    x=position,
                    linestyle=linestyle,
                    color=color,
                    alpha=GUIDE_ALPHA,
                    lw=linewidth,
                    zorder=3,
                )


def validate_plot_dimensions(
    width: float | None,
    aspect_ratio: float | None,
    interval_offset: float,
    min_coord_width: float,
    dpi: int,
) -> None:
    if width is not None and (not math.isfinite(width) or width <= 0):
        raise ValueError("plot width must be a finite number greater than zero")
    if aspect_ratio is not None and (
        not math.isfinite(aspect_ratio) or aspect_ratio <= 0
    ):
        raise ValueError(
            "plot aspect ratio must be a finite number greater than zero"
        )
    if not math.isfinite(interval_offset) or not 0 <= interval_offset < 1:
        raise ValueError("plot offset must be a finite fraction in [0, 1)")
    if not math.isfinite(min_coord_width) or not 0 <= min_coord_width <= 1:
        raise ValueError("minimum coordinate width must be in [0, 1]")
    if dpi <= 0:
        raise ValueError("plot DPI must be greater than zero")


def resolve_figure_size(
    *,
    default_width: float,
    default_height: float,
    width: float | None,
    aspect_ratio: float | None,
) -> tuple[float, float]:
    resolved_width = default_width if width is None else width
    resolved_height = (
        default_height
        if aspect_ratio is None
        else resolved_width * aspect_ratio
    )
    return resolved_width, resolved_height


def get_cycle_x_padding(
    align_to_combined: bool,
    connection_extension: float,
) -> float:
    """Reserve enough standalone x space for cycle-closing connectors."""
    if align_to_combined:
        return 0.0
    # Closing connectors can extend two lengths beyond the first segment.
    # Keep another extension between that connector and the axes boundary so
    # its vertical stroke is not clipped when the cycle starts at the first
    # displayed coordinate.
    return max(1.0, 3.0 * connection_extension)


def get_cycle_gene_track_ratio(
    resolved_height: float,
    gene_lane_count: int,
    align_to_combined: bool,
) -> float:
    """Size the standalone gene band from the occupied annotation lanes."""
    if align_to_combined:
        return GENE_TRACK_HEIGHT_RATIO
    desired_gene_height_inches = min(
        CYCLE_GENE_TRACK_MAX_HEIGHT_INCHES,
        CYCLE_GENE_TRACK_BASE_HEIGHT_INCHES
        + CYCLE_GENE_TRACK_HEIGHT_PER_LANE_INCHES * max(0, gene_lane_count),
    )
    desired_gene_height_inches = min(
        desired_gene_height_inches,
        resolved_height * CYCLE_GENE_TRACK_MAX_FIGURE_FRACTION,
    )
    return (
        PLOT_DATA_HEIGHT_RATIO
        * desired_gene_height_inches
        / (resolved_height - desired_gene_height_inches)
    )


def get_cycle_figure_margins(
    resolved_height: float,
    hide_genes: bool,
    align_to_combined: bool,
) -> tuple[float, float]:
    """Use physical margin caps so tall cycle plots do not gain whitespace."""
    if align_to_combined:
        return (0.12 if hide_genes else 0.30), 0.88
    default_bottom = 0.12 if hide_genes else 0.34
    max_bottom_inches = 0.6 if hide_genes else 2.0
    bottom = min(default_bottom, max_bottom_inches / resolved_height)
    top = 1.0 - min(0.12, 0.75 / resolved_height)
    return bottom, top


def get_interval_margin(
    num_intervals: int,
    interval_offset: float,
    genomic_plot_width: float = 100.0,
) -> float:
    """Convert a total gap fraction into an equal per-boundary margin."""
    if num_intervals <= 0 or interval_offset <= 0:
        return 0.0
    total_gap_width = (
        genomic_plot_width * interval_offset / (1.0 - interval_offset)
    )
    return total_gap_width / (num_intervals + 1)


def get_coverage_scale_max(
    coverage_values: list[float],
    coverage_scale: Literal["robust", "full"],
) -> float:
    finite_values = np.asarray(
        [
            value
            for value in coverage_values
            if math.isfinite(value) and value >= 0
        ],
        dtype=float,
    )
    if finite_values.size == 0:
        return 0.0
    if coverage_scale == "full" or finite_values.size < 20:
        return float(np.max(finite_values))
    if coverage_scale != "robust":
        raise ValueError("coverage scale must be 'robust' or 'full'")
    return max(
        float(np.percentile(finite_values, 95)),
        float(np.median(finite_values)),
    )


def parse_cycle_selection(value: str | None) -> list[CycleSelection]:
    if value is None or not value.strip():
        return []
    selections: list[CycleSelection] = []
    seen: set[int] = set()
    for token in re.split(r"[\s,]+", value.strip()):
        match = re.fullmatch(r"([pcPC])(\d+)", token)
        if match is None:
            raise ValueError(
                f"invalid cycle selection '{token}'; expected values such as p1,c2"
            )
        cycle_id = int(match.group(2))
        if cycle_id in seen:
            continue
        seen.add(cycle_id)
        selections.append(
            CycleSelection(
                cycle_id=cycle_id,
                expected_is_cyclic=match.group(1).lower() == "c",
            )
        )
    return selections


def resolve_cycle_ids(
    cycles: dict[int, datatypes.ReconstructedCycle],
    selections: list[CycleSelection],
) -> list[int]:
    if not selections:
        return list(cycles)
    cycle_ids: list[int] = []
    for selection in selections:
        if selection.cycle_id not in cycles:
            raise ValueError(
                f"cycle/path ID {selection.cycle_id} was not found"
            )
        cycle = cycles[selection.cycle_id]
        if cycle.is_cyclic != selection.expected_is_cyclic:
            actual_type = "cycle" if cycle.is_cyclic else "path"
            requested_type = "cycle" if selection.expected_is_cyclic else "path"
            raise ValueError(
                f"ID {selection.cycle_id} is a {actual_type}, not a {requested_type}"
            )
        cycle_ids.append(selection.cycle_id)
    return cycle_ids


def parse_cycle_color_file(
    color_file: pathlib.Path | None,
) -> dict[int, str]:
    if color_file is None:
        return {}
    colors: dict[int, str] = {}
    with color_file.open() as infile:
        for line_number, line in enumerate(infile, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = re.split(r"[\s,\t]+", stripped)
            if fields[0].lower() == "item" and len(fields) > 1:
                continue
            if len(fields) != 2:
                raise ValueError(
                    f"{color_file}:{line_number}: expected ITEM COLOR"
                )
            item, color = fields
            match = re.fullmatch(r"(?:[pcPC])?(\d+)", item)
            if match is None:
                raise ValueError(
                    f"{color_file}:{line_number}: invalid item '{item}'"
                )
            cycle_id = int(match.group(1))
            if cycle_id in colors:
                raise ValueError(
                    f"{color_file}:{line_number}: duplicate color for ID {cycle_id}"
                )
            if not mcolors.is_color_like(color):
                raise ValueError(
                    f"{color_file}:{line_number}: invalid color '{color}'"
                )
            colors[cycle_id] = color
    return colors


def assign_cycle_colors(
    cycle_ids: list[int],
    cycles: dict[int, datatypes.ReconstructedCycle],
) -> dict[int, str]:
    return {
        cycle_id: CYCLE_COLOR if cycles[cycle_id].is_cyclic else PATH_COLOR
        for cycle_id in cycle_ids
    }


def add_cycle_axis_labels(
    ax: Axes,
    tick_positions: list[float],
    labels: list[tuple[str, str]],
    fontsize: float,
) -> list[Text]:
    """Draw a regular cycle/path heading above an enlarged CN annotation."""
    if len(tick_positions) != len(labels):
        raise ValueError(
            "cycle tick positions and labels must have equal length"
        )
    ax.set_yticks(tick_positions)
    ax.set_yticklabels([""] * len(tick_positions))
    label_artists: list[Text] = []
    for tick_position, (heading, cn_label) in zip(
        tick_positions,
        labels,
        strict=True,
    ):
        label_artists.append(
            ax.text(
                -0.015,
                tick_position,
                heading,
                transform=ax.get_yaxis_transform(),
                ha="right",
                va="bottom",
                fontsize=fontsize * 1.1,
                weight="normal",
                clip_on=False,
            )
        )
        label_artists.append(
            ax.text(
                -0.015,
                tick_position,
                cn_label,
                transform=ax.get_yaxis_transform(),
                ha="right",
                va="top",
                fontsize=fontsize,
                weight="normal",
                clip_on=False,
            )
        )
    return label_artists


def set_interval_axis_labels(
    ax: Axes,
    intervals_by_chrom: dict[str, list[datatypes.Interval]],
    interval_starts: dict[str, list[float]],
    total_genomic_length: float,
    fontsize: float,
    min_coord_width: float,
) -> None:
    """Place chromosome labels separately from optional endpoint coordinates."""
    chromosome_positions: list[float] = []
    chromosome_labels: list[str] = []
    coordinate_positions: list[float] = []
    coordinate_labels: list[str] = []

    for chrom in breakpoint_utilities.sort_chrom_names(intervals_by_chrom):
        intervals = intervals_by_chrom[chrom]
        if not intervals:
            continue
        displayed_start = interval_starts[chrom][0]
        last_interval = intervals[-1]
        displayed_end = interval_starts[chrom][-1] + (
            (last_interval.end - last_interval.start)
            * 100.0
            / total_genomic_length
        )
        chromosome_positions.append((displayed_start + displayed_end) * 0.5)
        chromosome_labels.append(chrom)

        for interval_index, interval in enumerate(intervals):
            genomic_fraction = (
                interval.end - interval.start
            ) / total_genomic_length
            if genomic_fraction < min_coord_width:
                continue
            interval_start = interval_starts[chrom][interval_index]
            interval_end = interval_start + (
                (interval.end - interval.start) * 100.0 / total_genomic_length
            )
            coordinate_positions.extend((interval_start, interval_end))
            coordinate_labels.extend(
                (f"{interval.start:,}", f"{interval.end:,}")
            )

    ax.set_xticks(chromosome_positions)
    ax.set_xticklabels(chromosome_labels, fontsize=fontsize)
    ax.set_xticks(coordinate_positions, minor=True)
    ax.set_xticklabels(
        coordinate_labels,
        minor=True,
        rotation=90,
        fontsize=fontsize * GENOMIC_COORDINATE_FONT_SCALE,
    )
    ax.tick_params(axis="x", which="major", length=0, pad=3)
    ax.tick_params(axis="x", which="minor", length=3, pad=12)


def get_gene_font_size(
    font_size_multiplier: float,
    base_font_size: float = DEFAULT_GENE_FONT_SIZE,
) -> float:
    """Return a base point size after applying the plot-wide multiplier."""
    if not math.isfinite(font_size_multiplier) or font_size_multiplier < 0:
        raise ValueError(
            "font size multiplier must be a finite, non-negative number"
        )
    if not math.isfinite(base_font_size) or base_font_size < 0:
        raise ValueError("base font size must be a finite, non-negative number")
    scaled_font_size = base_font_size * font_size_multiplier
    if not math.isfinite(scaled_font_size):
        raise ValueError("scaled font size must be finite")
    return scaled_font_size


def scale_axis_elements(ax: Axes, font_size_multiplier: float) -> None:
    """Scale axis tick marks and spines with the plot-wide multiplier."""
    get_gene_font_size(font_size_multiplier)
    axis_names: tuple[Literal["x"], Literal["y"]] = ("x", "y")
    for axis_name in axis_names:
        for tick_kind in ("major", "minor"):
            ax.tick_params(
                axis=axis_name,
                which=tick_kind,
                length=rcParams[f"{axis_name}tick.{tick_kind}.size"]
                * font_size_multiplier,
                width=rcParams[f"{axis_name}tick.{tick_kind}.width"]
                * font_size_multiplier,
            )
    for spine in ax.spines.values():
        spine.set_linewidth(rcParams["axes.linewidth"] * font_size_multiplier)


def hide_figure_text_if_zero(
    fig: Figure,
    font_size_multiplier: float,
) -> None:
    """Fully suppress text instead of relying on zero-point rasterization."""
    if font_size_multiplier == 0:
        for text_artist in fig.findobj(match=Text):
            text_artist.set_visible(False)


def save_plot_figure(fig: Figure, output_fn: str, dpi: int) -> None:
    """Save exact-size PNG/PDF plots and always release their figure."""
    try:
        fig.savefig(output_fn + ".png", dpi=dpi)
        fig.savefig(output_fn + ".pdf")
    finally:
        plt.close(fig)


def combine_graph_and_cycle_plots(
    output_prefix: str,
    dpi: int,
    axis_layout: GenomicAxisLayout | None = None,
    *,
    hide_genes: bool = False,
) -> None:
    """Stack panels and optionally connect their shared genomic guides.

    The standalone graph retains its gene and coordinate track. In the combined
    figure, that lower annotation band is cropped so the shared annotations are
    shown only beneath the cycles panel and the two panels sit close together.
    """
    graph_image = plt.imread(output_prefix + "_graph.png")
    cycle_image = plt.imread(output_prefix + "_cycles.png")
    if axis_layout is not None:
        graph_axis_bottom = (
            GRAPH_MAIN_AXIS_BOTTOM_WITHOUT_GENES
            if hide_genes
            else GRAPH_MAIN_AXIS_BOTTOM_WITH_GENES
        )
        graph_crop_height = max(
            1,
            min(
                graph_image.shape[0],
                int(round(graph_image.shape[0] * (1.0 - graph_axis_bottom)))
                + 1,
            ),
        )
        graph_image = graph_image[:graph_crop_height]
    graph_height, graph_width = graph_image.shape[:2]
    cycle_height, cycle_width = cycle_image.shape[:2]
    combined_width = max(graph_width, cycle_width)
    combined_height = graph_height + cycle_height
    fig = plt.figure(
        figsize=(combined_width / dpi, combined_height / dpi),
        dpi=dpi,
        frameon=False,
    )
    cycle_left = (combined_width - cycle_width) / (2.0 * combined_width)
    graph_left = (combined_width - graph_width) / (2.0 * combined_width)
    cycle_ax = fig.add_axes(
        [
            cycle_left,
            0,
            cycle_width / combined_width,
            cycle_height / combined_height,
        ]
    )
    graph_ax = fig.add_axes(
        [
            graph_left,
            cycle_height / combined_height,
            graph_width / combined_width,
            graph_height / combined_height,
        ]
    )
    for axis, image_data in (
        (graph_ax, graph_image),
        (cycle_ax, cycle_image),
    ):
        axis.imshow(image_data, interpolation="none")
        axis.set_axis_off()

    if axis_layout is not None:
        cycle_fraction = cycle_height / combined_height
        graph_fraction = graph_height / combined_height
        baseline_label_mask_height = (
            GRAPH_BASELINE_LABEL_MASK_HEIGHT * graph_fraction
        )
        for mask_x, mask_width in (
            (0.0, COMBINED_AXIS_LEFT - 0.002),
            (COMBINED_AXIS_RIGHT + 0.002, 1.0 - COMBINED_AXIS_RIGHT),
        ):
            fig.add_artist(
                Rectangle(
                    (mask_x, cycle_fraction),
                    mask_width,
                    baseline_label_mask_height,
                    transform=fig.transFigure,
                    facecolor="white",
                    edgecolor="none",
                    zorder=5,
                )
            )
        guide_bottom = 0.88 * cycle_fraction
        guide_top = cycle_fraction
        guide_groups = (
            (
                axis_layout.chromosome_guides,
                "--",
                "black",
                CHROMOSOME_GUIDE_LINE_WIDTH,
            ),
            (
                axis_layout.interval_guides,
                ":",
                "0.25",
                INTERVAL_GUIDE_LINE_WIDTH,
            ),
        )
        for positions, linestyle, color, linewidth in guide_groups:
            for position in positions:
                normalized_x = position / axis_layout.axis_max
                figure_x = COMBINED_AXIS_LEFT + normalized_x * (
                    COMBINED_AXIS_RIGHT - COMBINED_AXIS_LEFT
                )
                fig.add_artist(
                    Line2D(
                        [figure_x, figure_x],
                        [guide_bottom, guide_top],
                        transform=fig.transFigure,
                        linestyle=linestyle,
                        color=color,
                        alpha=GUIDE_ALPHA,
                        linewidth=linewidth,
                        zorder=10,
                    )
                )
    save_plot_figure(fig, output_prefix + "_combined", dpi)


def get_text_layout_scale(font_size_multiplier: float) -> float:
    """Return bounded canvas expansion for layouts with dense text."""
    get_gene_font_size(font_size_multiplier)
    return min(
        max(1.0, font_size_multiplier),
        MAX_TEXT_LAYOUT_SCALE,
    )


def get_gene_label_padding(
    total_genomic_length: float,
    gene_font_size: float,
    plot_width: float,
    max_label_characters: int,
) -> float:
    """Estimate genomic padding needed to prevent adjacent gene-label overlap."""
    if total_genomic_length <= 0 or max_label_characters <= 0:
        return 0.0
    usable_axis_width_points = max(plot_width * 72.0 * 0.75, 1.0)
    estimated_label_width_points = max_label_characters * gene_font_size * 0.58
    label_width_fraction = (
        estimated_label_width_points / usable_axis_width_points
    )
    return total_genomic_length * min(max(0.02, label_width_fraction), 0.12)


def assign_gene_label_lanes(
    gene_spans: list[tuple[int, str, float, float]],
    gene_font_size: float,
    plot_width: float,
    displayed_axis_span: float,
    axis_width_fraction: float = 0.75,
) -> dict[int, float]:
    """Assign collision-free gene lanes using rendered italic-text widths."""
    if not gene_spans:
        return {}
    usable_axis_width_points = max(
        plot_width * 72.0 * axis_width_fraction,
        1.0,
    )
    italic_font = FontProperties(style="italic")
    padded_spans: list[tuple[float, float, int]] = []
    for gene_id, gene_name, gene_start, gene_end in gene_spans:
        label_width_points = (
            0.0
            if gene_font_size == 0
            else TextPath(
                (0, 0),
                gene_name,
                size=gene_font_size,
                prop=italic_font,
            )
            .get_extents()
            .width
            + GENE_LABEL_HORIZONTAL_PADDING_POINTS
        )
        label_width = (
            label_width_points / usable_axis_width_points * displayed_axis_span
        )
        center = (gene_start + gene_end) * 0.5
        padded_spans.append(
            (
                min(gene_start, center - label_width * 0.5),
                max(gene_end, center + label_width * 0.5),
                gene_id,
            )
        )

    lane_assignments: dict[int, int] = {}
    lane_ends: list[float] = []
    lane_gap = displayed_axis_span * GENE_LABEL_LANE_GAP_FRACTION
    for span_start, span_end, gene_id in sorted(padded_spans):
        available_lane: int | None = next(
            (
                lane_index
                for lane_index, lane_end in enumerate(lane_ends)
                if span_start >= lane_end + lane_gap
            ),
            None,
        )
        if available_lane is None:
            available_lane = len(lane_ends)
            lane_ends.append(float("-inf"))
        lane_assignments[gene_id] = available_lane
        lane_ends[available_lane] = span_end

    lane_heights = np.linspace(0.15, 0.78, len(lane_ends))
    return {
        gene_id: float(lane_heights[lane_index])
        for gene_id, lane_index in lane_assignments.items()
    }


def parse_gene_subset_file(gene_subset_file: pathlib.Path) -> list[str]:
    """Parse newline-, whitespace-, or comma-delimited gene names."""
    genes: list[str] = []
    with gene_subset_file.open() as infile:
        for line in infile:
            for gene_name in re.split(r"[\s,]+", line.strip()):
                if gene_name:
                    genes.append(gene_name)
    return genes


def merge_gene_subsets(
    gene_subset_list: list[str],
    gene_subset_file: pathlib.Path | None,
) -> list[str]:
    genes = list(gene_subset_list)
    if gene_subset_file is not None:
        file_genes = parse_gene_subset_file(gene_subset_file)
        if not file_genes:
            typer.secho(
                "Warning: gene subset file "
                f"{gene_subset_file} is empty; plotting all genes unless "
                "--gene-subset-list was also provided.",
                fg=typer.colors.YELLOW,
                err=True,
            )
        genes.extend(file_genes)

    deduped_genes = []
    seen_genes = set()
    for gene_name in genes:
        if gene_name in seen_genes:
            continue
        seen_genes.add(gene_name)
        deduped_genes.append(gene_name)
    return deduped_genes


def get_discordant_edge_linewidth(
    edge_read_count: float,
    max_read_count: float,
) -> float:
    if max_read_count <= 0:
        return DISCORDANT_EDGE_MIN_LINE_WIDTH
    aa_max_read_count = max(
        max_read_count, AA_DISCORDANT_EDGE_MIN_MAX_READ_COUNT
    )
    return max(
        DISCORDANT_EDGE_MIN_LINE_WIDTH,
        DISCORDANT_EDGE_MAX_LINE_WIDTH
        * min(1.0, max(edge_read_count, 0.0) / aa_max_read_count),
    )


def get_discordant_edge_arc_height(
    plot_distance: float,
    plot_width: float,
    max_segment_cn: float,
) -> float:
    """Return a distance-aware arc apex with publication-scale occupancy."""
    normalized_distance = min(
        max(plot_distance, 0.0) / max(plot_width, 1.0), 1.0
    )
    apex_scale = (
        MIN_ARC_APEX_CN_SCALE
        + (MAX_ARC_APEX_CN_SCALE - MIN_ARC_APEX_CN_SCALE) * normalized_distance
    )
    return max(max_segment_cn, 1.0) * apex_scale


def get_discordant_edge_arc_base(_max_segment_cn: float) -> float:
    return 0.0


def get_graph_axis_limits(
    *,
    max_coverage: float,
    max_segment_cn: float,
    cn_sum_squares: float,
    cn_coverage_cross_product: float,
    max_arc_apex: float,
    max_coverage_cutoff: float = float("inf"),
    coverage_headroom: float = DEFAULT_COVERAGE_HEADROOM,
    arc_top_padding: float = DEFAULT_ARC_TOP_PADDING,
) -> GraphAxisLimits:
    """Fit the twin axes, then expand both proportionally to contain arcs.

    CoRAL's zero-intercept least-squares fit models coverage as
    ``coverage_per_cn * CN``. The initial axis-limit ratio is set to that
    fitted slope. If an arc needs additional CN headroom, both limits are
    multiplied by the same factor so the visual CN/coverage alignment is
    unchanged.
    """
    natural_coverage_ymax = coverage_headroom * max(max_coverage, 0.0)
    if math.isfinite(max_coverage_cutoff):
        natural_coverage_ymax = min(
            natural_coverage_ymax,
            max(max_coverage_cutoff, 0.0),
        )
    fitted_coverage_ymax = max(natural_coverage_ymax, 1.0)

    if cn_sum_squares > 0 and cn_coverage_cross_product > 0:
        coverage_per_cn = cn_coverage_cross_product / cn_sum_squares
        fitted_cn_ymax = fitted_coverage_ymax / coverage_per_cn
    else:
        fitted_cn_ymax = max(
            coverage_headroom * max(max_segment_cn, 0.0),
            1.0,
        )

    required_cn_ymax = max(
        fitted_cn_ymax,
        max(max_arc_apex, 0.0) * (1.0 + arc_top_padding),
    )
    expansion_factor = required_cn_ymax / fitted_cn_ymax
    return GraphAxisLimits(
        coverage_ymax=fitted_coverage_ymax * expansion_factor,
        cn_ymax=required_cn_ymax,
        expansion_factor=expansion_factor,
    )


def get_graph_coverage_label(bam_path: pathlib.Path | None) -> str:
    if bam_path is None:
        return "Graph average coverage"
    return "BAM coverage"


def get_graph_legend_output_prefix(output_prefix: str) -> pathlib.Path:
    return pathlib.Path(f"{output_prefix}_legend")


def get_graph_legend_handles(coverage_label: str) -> list[Patch | Line2D]:
    legend_handles: list[Patch | Line2D] = [
        Patch(facecolor="silver", edgecolor="none", label=coverage_label),
        Line2D(
            [0],
            [0],
            color="black",
            lw=6,
            label="Predicted segment CN",
        ),
    ]
    for orientation, color in DISCORDANT_EDGE_COLORS.items():
        legend_handles.append(
            Line2D(
                [0],
                [0],
                color=color,
                lw=2,
                label=f"Discordant edge {orientation}",
            )
        )
    legend_handles.extend(
        (
            Patch(
                facecolor=CYCLE_SEGMENT_FACE_COLOR,
                edgecolor=CYCLE_SEGMENT_EDGE_COLOR,
                linewidth=CYCLE_SEGMENT_LINE_WIDTH,
                label="Genomic sequence edge",
            ),
            Line2D(
                [0],
                [0],
                color=CYCLE_COLOR,
                lw=CYCLE_CONNECTION_LINE_WIDTH,
                label="Cycle SV edge",
            ),
            Line2D(
                [0],
                [0],
                color=PATH_COLOR,
                lw=CYCLE_CONNECTION_LINE_WIDTH,
                label="Path SV edge",
            ),
        )
    )
    return legend_handles


def _draw_endpoint_stub(
    ax: Axes,
    x: float,
    y: float,
    strand: str,
    *,
    length: float = 0.13,
    color: str = "black",
    linewidth: float = 4.0,
) -> None:
    """Draw the genomic segment incident on a breakpoint endpoint."""
    segment_end = x - length if strand == "+" else x + length
    ax.plot(
        [x, segment_end],
        [y, y],
        color=color,
        lw=linewidth,
        solid_capstyle="butt",
        clip_on=False,
    )
    ax.plot(
        [x],
        [y],
        marker="o",
        markersize=4.5,
        color=color,
        clip_on=False,
    )


def _draw_orientation_example(
    ax: Axes,
    orientation: str,
    color: str | tuple[float, float, float],
    font_size_multiplier: float,
) -> None:
    """Draw one miniature discordant-edge junction."""
    left_strand, right_strand = orientation
    left_x, right_x, y = 0.28, 0.72, 0.34
    _draw_endpoint_stub(ax, left_x, y, left_strand)
    _draw_endpoint_stub(ax, right_x, y, right_strand)
    arc = Arc(
        ((left_x + right_x) / 2, y),
        right_x - left_x,
        0.48,
        theta1=0,
        theta2=180,
        color=color,
        lw=3,
    )
    ax.add_patch(arc)
    label_font_size = get_gene_font_size(
        font_size_multiplier,
        DEFAULT_LEGEND_FONT_SIZE,
    )
    ax.text(
        left_x,
        0.42,
        left_strand,
        ha="center",
        va="bottom",
        weight="bold",
        fontsize=label_font_size,
    )
    ax.text(
        right_x,
        0.42,
        right_strand,
        ha="center",
        va="bottom",
        weight="bold",
        fontsize=label_font_size,
    )
    ax.text(
        0.5,
        0.82,
        orientation,
        ha="center",
        va="center",
        weight="bold",
        fontsize=label_font_size,
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


def write_graph_legend(
    output_prefix: str,
    coverage_label: str,
    *,
    font_size_multiplier: float = 1.0,
    dpi: int = 300,
) -> None:
    fontsize = get_gene_font_size(
        font_size_multiplier,
        DEFAULT_LEGEND_FONT_SIZE,
    )
    title_fontsize = get_gene_font_size(font_size_multiplier, 14.0)
    small_fontsize = get_gene_font_size(font_size_multiplier, 9.0)
    legend_output_prefix = get_graph_legend_output_prefix(output_prefix)
    legend_output_prefix.parent.mkdir(parents=True, exist_ok=True)
    png_path = legend_output_prefix.with_suffix(".png")
    pdf_path = legend_output_prefix.with_suffix(".pdf")
    layout_scale = get_text_layout_scale(font_size_multiplier)
    fig = plt.figure(
        figsize=(8.8 * layout_scale, 4.6 * layout_scale),
        facecolor="white",
    )
    layout = fig.add_gridspec(
        4,
        2,
        height_ratios=(1.15, 0.22, 0.9, 0.9),
        hspace=0.05,
        wspace=0.18,
    )

    summary_ax = fig.add_subplot(layout[0, :])
    summary_ax.axis("off")
    summary_ax.set_title(
        "Legend",
        fontsize=title_fontsize,
        weight="bold",
        pad=4,
    )
    summary_ax.legend(
        handles=get_graph_legend_handles(coverage_label),
        loc="center",
        frameon=False,
        fontsize=fontsize,
        title="Discordant edge width ∝ read count",
        title_fontsize=fontsize,
        ncol=2,
        columnspacing=2.2,
        handlelength=3.0,
    )

    coordinate_ax = fig.add_subplot(layout[1, :])
    coordinate_ax.set_xlim(0, 1)
    coordinate_ax.set_ylim(0, 1)
    coordinate_ax.axis("off")
    coordinate_ax.add_patch(
        FancyArrowPatch(
            (0.27, 0.64),
            (0.73, 0.64),
            arrowstyle="-|>",
            mutation_scale=11,
            lw=1.2,
            color="dimgray",
        )
    )
    coordinate_ax.text(
        0.5,
        0.74,
        "reference coordinate increases",
        ha="center",
        va="bottom",
        fontsize=small_fontsize,
        color="dimgray",
    )
    for grid_cell, orientation in zip(
        (layout[2, 0], layout[2, 1], layout[3, 0], layout[3, 1]),
        ("+-", "++", "-+", "--"),
    ):
        example_ax = fig.add_subplot(grid_cell)
        _draw_orientation_example(
            example_ax,
            orientation,
            DISCORDANT_EDGE_COLORS[orientation],
            font_size_multiplier,
        )

    fig.subplots_adjust(left=0.04, right=0.96, top=0.96, bottom=0.04)
    hide_figure_text_if_zero(fig, font_size_multiplier)
    try:
        fig.savefig(png_path, dpi=dpi, bbox_inches="tight")
        fig.savefig(pdf_path, bbox_inches="tight")
    finally:
        plt.close(fig)


# makes a gene object from parsed refGene data
# this stores global properties for the gene
class Gene:
    def __init__(self, gchrom, gstart, gend, gdata):
        self.gchrom = gchrom
        self.gstart = gstart
        self.gend = gend
        self.gname = gdata[-4]
        self.strand = gdata[3]
        self.height = 0.5
        # self.highlight_name = highlight_name
        estarts = [int(x) for x in gdata[9].rsplit(",") if x]
        eends = [int(x) for x in gdata[10].rsplit(",") if x]
        self.eposns = list(zip(estarts, eends))

    def __str__(self):
        return f"Gene Name: {self.gname}, Chromosome: {self.gchrom}, Start: {self.gstart}, End: {self.gend}, Strand: {self.strand}"


@dataclass
class GraphViz:
    """ """

    lr_bamfh: Optional[pysam.AlignmentFile] = None
    bam: Optional[bam_types.BAMWrapper] = None
    graph: datatypes.BreakpointGraph | None = None

    graph_amplified_intervals: dict[
        core_types.ChrTag, list[datatypes.Interval]
    ] = field(default_factory=lambda: defaultdict(list))
    num_amplified_intervals: int = 0
    cycle_amplified_intervals: dict[
        core_types.ChrTag, list[datatypes.Interval]
    ] = field(default_factory=lambda: defaultdict(list))
    discordant_edges: list[datatypes.BreakpointEdge] = field(
        default_factory=list
    )

    cycles: dict[int, datatypes.ReconstructedCycle] = field(
        default_factory=dict
    )
    genes: DefaultDict[str, intervaltree.IntervalTree] = field(
        default_factory=lambda: defaultdict(intervaltree.IntervalTree)
    )
    plot_bounds: tuple[str, int, int] | None = None

    def open_bam(self, bam_fn: str) -> None:
        self.lr_bamfh = pysam.AlignmentFile(bam_fn, "rb")
        self.bam = bam_types.BAMWrapper(bam_fn, "rb")

    @property
    def sequence_edges_by_chr(self) -> dict[str, list[datatypes.SequenceEdge]]:
        seq_edges_by_chr: dict[str, list[datatypes.SequenceEdge]] = defaultdict(
            list
        )
        for seq_edge in self.graph.sequence_edges:  # type: ignore[union-attr]
            seq_edges_by_chr[seq_edge.chr].append(seq_edge)
        return seq_edges_by_chr

    def parse_genes(
        self,
        ref_genome: core_types.ReferenceGenome,
        gene_subset_list=None,
        restrict_to_bushman=False,
        refgene_file: pathlib.Path | None = None,
    ) -> None:
        if refgene_file is None:
            if ref_genome == core_types.ReferenceGenome.other:
                logger.warning(
                    "No --refgene-file provided for --ref other; "
                    "skipping gene track."
                )
                return
            ref_gene_filepath: Any = (
                importlib.resources.files(supplemental_data)
                / f"refGene_{ref_genome}.txt"
            )
            use_custom = False
        else:
            ref_gene_filepath = refgene_file
            use_custom = True

        bushman_set = set()
        if restrict_to_bushman:
            bushman_filepath = (
                importlib.resources.files(supplemental_data)
                / "Bushman_group_allOnco_May2018.tsv"
            )
            with bushman_filepath.open("r") as infile:
                _ = next(infile)
                for line in infile:
                    if not (bushman_fields := line.rstrip().rsplit()):
                        continue
                    bushman_set.add(bushman_fields[-1].strip('"'))

        seen_names = set()
        with open(ref_gene_filepath) as infile:
            for line in infile:
                if not (fields := line.rstrip().rsplit()):
                    continue
                curr_chrom: str = fields[2]
                if (
                    not use_custom
                    and ref_genome
                    in {
                        core_types.ReferenceGenome.hg19,
                        core_types.ReferenceGenome.hg38,
                    }
                    and not curr_chrom.startswith("chr")
                ):
                    curr_chrom = "chr" + curr_chrom

                tstart = int(fields[4])
                tend = int(fields[5])
                gname = fields[-4]
                is_other_feature = gname.startswith(("LOC", "LINC", "MIR"))
                if (restrict_to_bushman and gname not in bushman_set) or (
                    gene_subset_list and gname not in gene_subset_list
                ):
                    continue

                if gname not in seen_names and not is_other_feature:
                    seen_names.add(gname)
                    curr_gene = Gene(curr_chrom, tstart, tend, fields)
                    self.genes[curr_chrom][tstart:tend] = curr_gene

        if not self.genes:
            logger.warning(
                "No genes loaded from %s. Chromosome names in the gene file "
                "may not match the BAM/graph. Gene track will be empty.",
                ref_gene_filepath,
            )

    def update_graph_intervals(self) -> None:
        for chrom in self.sequence_edges_by_chr:
            lstart, lend = -2, -2
            if chrom not in self.graph_amplified_intervals:
                self.graph_amplified_intervals[chrom] = []
            for seq_edge in self.sequence_edges_by_chr[chrom]:
                start = seq_edge.start
                end = seq_edge.end
                if start != lend + 1:
                    if lstart >= 0:
                        self.graph_amplified_intervals[chrom].append(
                            datatypes.Interval(chrom, lstart, lend)
                        )
                        self.num_amplified_intervals += 1
                    lstart = start
                    lend = end
                else:
                    lend = end
            self.graph_amplified_intervals[chrom].append(
                datatypes.Interval(chrom, lstart, lend)
            )
            self.num_amplified_intervals += 1

    def merge_intervals(
        self, interval_list: list[tuple[int, int]], padding: float = 0.0
    ) -> list[tuple[int, int]]:
        # takes a list of interval tuples (pos1, pos2). Assumes from same chrom
        # return a list of interval tuples that are merged if overlapping or directly adjacent, or within padding distance
        sorted_intervals = sorted(interval_list)
        merged = [sorted_intervals[0]]
        for current in sorted_intervals[1:]:
            prev = merged[-1]
            if current[0] <= prev[1] + padding:
                merged[-1] = (prev[0], max(prev[1], current[1]))
            else:
                merged.append(current)

        return merged

    def update_cycle_amplified_intervals(
        self,
        cycle_ids: list[int] | None = None,
        *,
        cycle_only: bool = False,
        graph_given: bool = False,
        align_to_graph: bool = False,
    ) -> None:
        """Derive amplified intervals from (selected) cycles"""
        self.num_amplified_intervals = 0
        self.cycle_amplified_intervals = defaultdict(list)
        if cycle_ids is None:
            cycle_ids = list(self.cycles.keys())
        if cycle_only:
            cycle_ids = [
                cycle_id
                for cycle_id in self.cycles
                if self.cycles[cycle_id].is_cyclic
            ]
        if align_to_graph:
            if not graph_given:
                raise ValueError("cycle alignment requires a breakpoint graph")
            for chrom, intervals in self.graph_amplified_intervals.items():
                self.cycle_amplified_intervals[chrom] = list(intervals)
        elif graph_given:  # if the graph file is given, use this to set the amplified intervals
            for cycle_id in cycle_ids:
                for cycle_seg in self.cycles[cycle_id].segments:
                    for graph_intv in self.graph_amplified_intervals[
                        cycle_seg.chr
                    ]:
                        if not graph_intv.encompasses(cycle_seg):
                            continue

                        if (
                            graph_intv
                            not in self.cycle_amplified_intervals[cycle_seg.chr]
                        ):
                            self.cycle_amplified_intervals[
                                cycle_seg.chr
                            ].append(graph_intv)
                        break

        else:  # if the graph file is not given extract from the cycles file
            # collect a list of intervals for each chrom
            cycle_ivald = defaultdict(list)
            for cycle_id in cycle_ids:
                for cycle_seg in self.cycles[cycle_id].segments:
                    cycle_ivald[cycle_seg.chr].append(
                        (cycle_seg.start, cycle_seg.end)
                    )

                # merge
                for chrom, ival_list in cycle_ivald.items():
                    merged = self.merge_intervals(ival_list, padding=10000)
                    self.cycle_amplified_intervals[chrom] = [
                        datatypes.Interval(chrom, start, end)
                        for start, end in merged
                    ]

        for chr_tag in self.cycle_amplified_intervals:
            self.cycle_amplified_intervals[chr_tag] = sorted(
                self.cycle_amplified_intervals[chr_tag]
            )
            self.num_amplified_intervals += len(
                self.cycle_amplified_intervals[chr_tag]
            )

    def set_gene_heights(
        self, rel_genes: list[Gene], padding: float = 0.0
    ) -> None:
        if not rel_genes:
            return
        gname_to_gobj = {x.gname: x for x in rel_genes}
        # merge intervals
        intervals = [(x.gstart, x.gend) for x in rel_genes]
        merged = self.merge_intervals(intervals, padding=padding)

        gene_ival_t = intervaltree.IntervalTree[str]()
        for x in rel_genes:
            gene_ival_t.addi(x.gstart, x.gend, x.gname)

        for mi in merged:
            ghits: list[intervaltree.Interval] = gene_ival_t[mi[0] : mi[1]]
            gene_heights = np.linspace(0.15, 0.75, len(ghits))
            for g, h in zip(ghits, gene_heights):
                gname_to_gobj[g.data].height = h

    def set_displayed_gene_heights(
        self,
        intervals_by_chrom: dict[str, list[datatypes.Interval]],
        interval_starts: dict[str, list[float]],
        total_genomic_length: float,
        plot_width: float,
        gene_font_size: float,
        displayed_axis_span: float,
        axis_width_fraction: float = 0.75,
    ) -> int:
        """Assign lanes after genomic intervals have been mapped to the canvas."""
        displayed_genes: dict[int, tuple[Gene, float, float]] = {}
        for chrom, intervals in intervals_by_chrom.items():
            for interval_index, interval in enumerate(intervals):
                if self.plot_bounds and (
                    chrom != self.plot_bounds[0]
                    or interval.end < self.plot_bounds[1]
                    or interval.start > self.plot_bounds[2]
                ):
                    continue
                for gene_interval in self.genes[chrom][
                    interval.start : interval.end
                ]:
                    gene = gene_interval.data
                    cut_start = max(interval.start, gene.gstart)
                    cut_end = min(interval.end, gene.gend)
                    if self.plot_bounds:
                        cut_start = max(cut_start, self.plot_bounds[1])
                        cut_end = min(cut_end, self.plot_bounds[2])
                    if cut_start >= cut_end:
                        continue
                    gene_start = interval_starts[chrom][interval_index] + (
                        (cut_start - interval.start)
                        * 100.0
                        / total_genomic_length
                    )
                    gene_end = interval_starts[chrom][interval_index] + (
                        (cut_end - interval.start)
                        * 100.0
                        / total_genomic_length
                    )
                    gene_id = id(gene)
                    if gene_id in displayed_genes:
                        prior_gene, prior_start, prior_end = displayed_genes[
                            gene_id
                        ]
                        displayed_genes[gene_id] = (
                            prior_gene,
                            min(prior_start, gene_start),
                            max(prior_end, gene_end),
                        )
                    else:
                        displayed_genes[gene_id] = (gene, gene_start, gene_end)

        assignments = assign_gene_label_lanes(
            [
                (gene_id, gene.gname, gene_start, gene_end)
                for gene_id, (
                    gene,
                    gene_start,
                    gene_end,
                ) in displayed_genes.items()
            ],
            gene_font_size,
            plot_width,
            displayed_axis_span,
            axis_width_fraction,
        )
        for gene_id, (gene, _, _) in displayed_genes.items():
            gene.height = assignments[gene_id]
        return len(set(assignments.values()))

    def plot_graph(
        self,
        title: str,
        output_fn: str,
        fontsize: float = DEFAULT_PLOT_FONT_SIZE,
        dpi: int = 300,
        max_cov_cutoff: float = float("inf"),
        quality_threshold: float = 0,
        gene_font_size: float = DEFAULT_GENE_FONT_SIZE,
        *,
        hide_genes: bool = False,
        font_size_multiplier: float = 1.0,
        plot_width: float | None = None,
        aspect_ratio: float | None = None,
        interval_offset: float = DEFAULT_INTERVAL_OFFSET,
        min_coord_width: float = DEFAULT_MIN_COORD_WIDTH,
        coverage_scale: Literal["robust", "full"] = "robust",
        align_to_combined: bool = False,
    ) -> None:
        """Plot discordant edges and coverage on sequence edges in breakpoint
        graph."""
        validate_plot_dimensions(
            plot_width,
            aspect_ratio,
            interval_offset,
            min_coord_width,
            dpi,
        )
        width, height = resolve_figure_size(
            default_width=DEFAULT_GRAPH_WIDTH,
            default_height=DEFAULT_GRAPH_HEIGHT,
            width=plot_width,
            aspect_ratio=aspect_ratio,
        )
        axis_layout = build_genomic_axis_layout(
            self.graph_amplified_intervals,
            interval_offset,
        )
        margin_between_intervals = get_interval_margin(
            self.num_amplified_intervals,
            interval_offset,
        )
        fig = plt.figure(figsize=(width, height))
        if not hide_genes:
            gs = gridspec.GridSpec(
                2,
                1,
                height_ratios=[
                    PLOT_DATA_HEIGHT_RATIO,
                    GENE_TRACK_HEIGHT_RATIO,
                ],
            )
        else:
            gs = gridspec.GridSpec(
                2,
                1,
                height_ratios=[PLOT_DATA_HEIGHT_RATIO, 0.000001],
            )
        ax = fig.add_subplot(gs[0, 0])
        fig.subplots_adjust(
            left=COMBINED_AXIS_LEFT if align_to_combined else 0.13,
            right=COMBINED_AXIS_RIGHT if align_to_combined else 0.88,
            bottom=0.30 if not hide_genes else 0.12,
            top=0.88,
            hspace=0,
        )
        ax.set_title(title, fontsize=fontsize)
        ax2 = ax.twinx()
        ax.set_zorder(2)
        ax2.set_zorder(1)
        ax.patch.set_visible(False)
        ax2.patch.set_visible(False)
        ax3 = fig.add_subplot(gs[1, 0], sharex=ax)
        ax3.set_zorder(3)
        ax3.patch.set_visible(False)
        for plot_axis in (ax, ax2, ax3):
            scale_axis_elements(plot_axis, font_size_multiplier)
        # ax.yaxis.set_label_coords(-0.05, 0.25)
        # ax2.yaxis.set_label_coords(1.05, 0.33)
        ax.xaxis.set_visible(False)
        ax2.xaxis.set_visible(False)
        ax3.yaxis.set_visible(False)
        ax3.spines["left"].set_visible(False)
        ax3.spines["right"].set_visible(False)
        ax3.spines["top"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax2.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax2.spines["left"].set_visible(False)
        ax2.spines["bottom"].set_linewidth(max(1.5, 2.5 * font_size_multiplier))

        # Draw sequence edges
        total_len_amp = axis_layout.total_genomic_length
        # sorted_chrs = sorted(self.intervals_from_graph.keys(), key = lambda chr: CHR_TAG_TO_IDX[chr])
        zoom_factor = 1.0
        if self.plot_bounds:
            zoom_factor = (
                float(self.plot_bounds[2] - self.plot_bounds[1]) / total_len_amp
            )
        sorted_chrs = breakpoint_utilities.sort_chrom_names(
            self.graph_amplified_intervals.keys()
        )
        amplified_intervals_start = axis_layout.interval_starts
        sequence_edge_plot_positions = {}
        cn_bar_specs: list[tuple[float, float, float]] = []
        ymax = 0
        cn_sum_squares = 0.0
        cn_coverage_cross_product = 0.0
        for chrom in sorted_chrs:
            for seq in self.sequence_edges_by_chr[chrom]:
                x1 = genomic_position_to_axis(
                    seq.start,
                    chrom,
                    self.graph_amplified_intervals,
                    axis_layout,
                )
                x2 = genomic_position_to_axis(
                    seq.end,
                    chrom,
                    self.graph_amplified_intervals,
                    axis_layout,
                )
                sequence_edge_plot_positions[id(seq)] = (x1, x2)
                if self.plot_bounds:
                    if chrom != self.plot_bounds[0]:
                        continue  # Skip if chromosome doesn't match plot bounds

                    if not (
                        seq.end >= self.plot_bounds[1]
                        and seq.start <= self.plot_bounds[2]
                    ):
                        continue  # Skip if interval doesn't overlap with plot bounds

                y = seq.cn
                cn_sum_squares += seq.cn**2
                cn_coverage_cross_product += seq.cn * seq.lr_nc
                ymax = max(y, ymax)

                cn_bar_specs.append((x1, x2, y))

        # Draw amplified interval separators
        if not self.plot_bounds:
            draw_interval_guides(ax, ax3, axis_layout)

        # Draw discordant edges
        assert self.graph
        max_bp_read_count = max(
            [bp.lr_count for bp in self.graph.discordant_edges],
            default=0,
        )
        arc_plot_width = axis_layout.axis_max
        if self.plot_bounds:
            arc_plot_width = (
                (self.plot_bounds[2] - self.plot_bounds[1])
                * 100.0
                / total_len_amp
            )
        discordant_edge_plots = []
        for bp in self.graph.discordant_edges:
            chr1 = bp.node1.chr
            pos1 = bp.node1.pos
            chr2 = bp.node2.chr
            pos2 = bp.node2.pos
            ort = f"{bp.node1.strand}{bp.node2.strand}"
            if (
                chr1 in self.graph_amplified_intervals
                and chr2 in self.graph_amplified_intervals
            ):
                bp_x1 = genomic_position_to_axis(
                    pos1,
                    chr1,
                    self.graph_amplified_intervals,
                    axis_layout,
                )
                bp_x2 = genomic_position_to_axis(
                    pos2,
                    chr2,
                    self.graph_amplified_intervals,
                    axis_layout,
                )
                # check if either bp overlaps before plotting
                if self.plot_bounds:
                    # Check if both breakpoints belong to the same chromosome
                    # as in plot bounds
                    hit1 = (
                        chr1 == self.plot_bounds[0]
                        and self.plot_bounds[1] <= pos1 <= self.plot_bounds[2]
                    )
                    hit2 = (
                        chr2 == self.plot_bounds[0]
                        and self.plot_bounds[1] <= pos2 <= self.plot_bounds[2]
                    )
                    if not hit1 and not hit2:
                        continue

                discordant_edge_plots.append((bp, bp_x1, bp_x2, ort))

            else:
                print("Could not place " + str(bp))
                continue

        arc_base_y = get_discordant_edge_arc_base(ymax)
        max_arc_top_y = 0.0
        discordant_edge_arc_specs = []
        for bp, bp_x1, bp_x2, ort in discordant_edge_plots:
            arc_height = get_discordant_edge_arc_height(
                abs(bp_x2 - bp_x1),
                arc_plot_width,
                ymax,
            )
            max_arc_top_y = max(max_arc_top_y, arc_base_y + arc_height)
            discordant_edge_arc_specs.append(
                (bp, bp_x1, bp_x2, ort, arc_height)
            )

        ax2.axhline(
            2.0,
            color="#b2182b",
            alpha=0.65,
            linewidth=0.8,
            zorder=2,
        )
        ax2.set_ylabel("CN", fontsize=fontsize)
        ax2.tick_params(axis="y", labelsize=fontsize)

        # Draw coverage within amplified intervals
        max_cov = 0
        coverage_values: list[float] = []
        coverage_bar_specs: list[tuple[float, float, float]] = []
        if self.bam is not None:
            for chrom in sorted_chrs:
                for inti in range(len(self.graph_amplified_intervals[chrom])):
                    graph_intv = self.graph_amplified_intervals[chrom][inti]
                    if self.plot_bounds:
                        if chrom != self.plot_bounds[0]:
                            continue  # Skip if chromosome doesn't match plot bounds

                        if not (
                            graph_intv.end >= self.plot_bounds[1]
                            and graph_intv.start <= self.plot_bounds[2]
                        ):
                            continue  # Skip if interval doesn't overlap with plot bounds

                    window_size = 150
                    ival_len = graph_intv.end - graph_intv.start
                    if self.plot_bounds:
                        ival_len = self.plot_bounds[2] - self.plot_bounds[1]

                    if ival_len >= 1_000_000:
                        window_size = 10_000
                    elif ival_len >= 100_000:
                        window_size = 1_000

                    for w in range(
                        graph_intv.start, graph_intv.end, window_size
                    ):
                        intv = Interval(chrom, w, w + window_size)
                        cov = (
                            self.bam.count_raw_coverage(
                                intv,
                                quality_threshold=quality_threshold,
                                read_callback_type="nofilter",
                            )
                            / window_size
                        )
                        coverage_values.append(cov)
                        max_cov = max(cov, max_cov)
                        x = (
                            amplified_intervals_start[chrom][inti]
                            + (w - graph_intv.start) * 100.0 / total_len_amp
                        )
                        coverage_bar_specs.append(
                            (
                                x,
                                window_size * 100.0 / total_len_amp,
                                cov,
                            )
                        )
                    w = graph_intv.end - (
                        (graph_intv.end - graph_intv.start + 1) % window_size
                    )
                    if w < graph_intv.end:
                        cov = (
                            self.bam.count_raw_coverage(
                                Interval(chrom, w, w + window_size),
                                quality_threshold=quality_threshold,
                                read_callback_type="nofilter",
                            )
                            * 1.0
                            / window_size
                        )
                        coverage_values.append(cov)
                        max_cov = max(cov, max_cov)
                        x = (
                            amplified_intervals_start[chrom][inti]
                            + (w - graph_intv.start) * 100.0 / total_len_amp
                        )
                        coverage_bar_specs.append(
                            (
                                x,
                                window_size * 100.0 / total_len_amp,
                                cov,
                            )
                        )
        else:
            for chrom in sorted_chrs:
                for seq in self.sequence_edges_by_chr[chrom]:
                    if self.plot_bounds:
                        if chrom != self.plot_bounds[0]:
                            continue  # Skip if chromosome doesn't match plot bounds

                        if not (
                            seq.end >= self.plot_bounds[1]
                            and seq.start <= self.plot_bounds[2]
                        ):
                            continue  # Skip if interval doesn't overlap with plot bounds

                    x1, x2 = sequence_edge_plot_positions[id(seq)]
                    coverage_values.append(seq.lr_nc)
                    max_cov = max(seq.lr_nc, max_cov)
                    coverage_bar_specs.append(
                        (
                            x1,
                            x2 - x1,
                            seq.lr_nc,
                        )
                    )
        max_cov = get_coverage_scale_max(coverage_values, coverage_scale)
        axis_limits = get_graph_axis_limits(
            max_coverage=max_cov,
            max_segment_cn=ymax,
            cn_sum_squares=cn_sum_squares,
            cn_coverage_cross_product=cn_coverage_cross_product,
            max_arc_apex=max_arc_top_y,
            max_coverage_cutoff=max_cov_cutoff,
        )
        if (
            math.isfinite(max_cov_cutoff)
            and axis_limits.coverage_ymax > max_cov_cutoff
        ):
            logger.warning(
                "Expanded the coverage axis above --max-coverage %.3f "
                "to %.3f so discordant-edge arcs remain visible while "
                "preserving the fitted coverage/CN alignment.",
                max_cov_cutoff,
                axis_limits.coverage_ymax,
            )
        ax2.set_ylim(0, axis_limits.cn_ymax)
        ax.set_ylabel("Coverage", fontsize=fontsize)
        ax.set_ylim(0, axis_limits.coverage_ymax)
        ax.tick_params(axis="y", labelsize=fontsize)
        align_cn_axis_ticks(ax, ax2, axis_limits)
        for bar_x, bar_width, coverage in coverage_bar_specs:
            ax.add_patch(
                Rectangle(
                    (bar_x, 0),
                    bar_width,
                    coverage,
                    facecolor="silver",
                    edgecolor="none",
                    zorder=4,
                )
            )
        coverage_per_cn = axis_limits.coverage_ymax / axis_limits.cn_ymax
        for cn_x1, cn_x2, segment_cn in cn_bar_specs:
            ax.hlines(
                segment_cn * coverage_per_cn,
                cn_x1,
                cn_x2,
                color="black",
                lw=5,
                zorder=6,
            )
        for bp, bp_x1, bp_x2, ort, raw_arc_height in discordant_edge_arc_specs:
            ax2.add_patch(
                Arc(
                    ((bp_x1 + bp_x2) * 0.5, arc_base_y),
                    abs(bp_x2 - bp_x1),
                    2 * raw_arc_height,
                    theta1=0,
                    theta2=180,
                    color=DISCORDANT_EDGE_COLORS[ort],
                    lw=get_discordant_edge_linewidth(
                        bp.lr_count,
                        max_bp_read_count,
                    ),
                    alpha=0.8,
                    zorder=2,
                )
            )

        # draw genes below plot
        if not hide_genes:
            displayed_axis_span = (
                (self.plot_bounds[2] - self.plot_bounds[1])
                * 100.0
                / total_len_amp
                if self.plot_bounds
                else 100
                + (self.num_amplified_intervals + 1) * margin_between_intervals
            )
            self.set_displayed_gene_heights(
                self.graph_amplified_intervals,
                amplified_intervals_start,
                total_len_amp,
                width,
                gene_font_size,
                displayed_axis_span,
                (
                    COMBINED_AXIS_RIGHT - COMBINED_AXIS_LEFT
                    if align_to_combined
                    else 0.75
                ),
            )
            for chrom in sorted_chrs:
                for inti in range(len(self.graph_amplified_intervals[chrom])):
                    graph_intv = self.graph_amplified_intervals[chrom][inti]
                    if self.plot_bounds:
                        if chrom != self.plot_bounds[0]:
                            continue  # Skip if chromosome doesn't match plot bounds

                        if not (
                            graph_intv.end >= self.plot_bounds[1]
                            and graph_intv.start <= self.plot_bounds[2]
                        ):
                            continue  # Skip if interval doesn't overlap with plot bounds

                    rel_genes = [
                        x.data
                        for x in self.genes[chrom][
                            graph_intv.start : graph_intv.end
                        ]
                    ]
                    for gene_obj in rel_genes:  # plot line for the gene
                        height = gene_obj.height
                        cut_gs = max(graph_intv.start, gene_obj.gstart)
                        cut_ge = min(graph_intv.end, gene_obj.gend)
                        # if self.plot_bounds:
                        #     cut_gs = max(self.plot_bounds[1], cut_gs)
                        #     cut_ge = min(self.plot_bounds[2], cut_ge)

                        gene_start = (
                            amplified_intervals_start[chrom][inti]
                            + (cut_gs - graph_intv.start)
                            * 100.0
                            / total_len_amp
                        )
                        gene_end = (
                            amplified_intervals_start[chrom][inti]
                            + (cut_ge - graph_intv.start)
                            * 100.0
                            / total_len_amp
                        )
                        ax3.hlines(
                            height,
                            gene_start,
                            gene_end,
                            color="cornflowerblue",
                            lw=4.5,
                        )  # Draw horizontal bars for genes
                        if self.plot_bounds:
                            if (
                                cut_ge < self.plot_bounds[1]
                                or cut_gs > self.plot_bounds[2]
                            ):
                                continue

                            cut_gs = max(self.plot_bounds[1], cut_gs)
                            cut_ge = min(self.plot_bounds[2], cut_ge)
                            gene_start = (
                                amplified_intervals_start[chrom][inti]
                                + (cut_gs - graph_intv.start)
                                * 100.0
                                / total_len_amp
                            )
                            gene_end = (
                                amplified_intervals_start[chrom][inti]
                                + (cut_ge - graph_intv.start)
                                * 100.0
                                / total_len_amp
                            )

                        if gene_font_size > 0:
                            ax3.text(
                                (gene_start + gene_end) / 2,
                                height + 0.05,
                                gene_obj.gname,
                                ha="center",
                                va="bottom",
                                fontsize=gene_font_size,
                                style="italic",
                                clip_on=False,
                                zorder=10,
                            )

                        if gene_obj.strand == "+":
                            ax3.plot(
                                gene_start,
                                height,
                                marker=">",
                                color="black",
                                markersize=7,
                            )
                        elif gene_obj.strand == "-":
                            ax3.plot(
                                gene_end,
                                height,
                                marker="<",
                                color="black",
                                markersize=7,
                            )

                        for (
                            exon_start,
                            exon_end,
                        ) in gene_obj.eposns:  # plot exon bars
                            if (
                                not exon_end > graph_intv.start
                                or not exon_start < graph_intv.end
                            ):
                                continue

                            cut_es = max(graph_intv.start, exon_start)
                            cut_ee = min(graph_intv.end, exon_end)
                            # if self.plot_bounds:
                            #     cut_es = max(self.plot_bounds[1], cut_es)
                            #     cut_ee = min(self.plot_bounds[2], cut_ee)
                            #
                            exon_start_pos = (
                                amplified_intervals_start[chrom][inti]
                                + (cut_es - graph_intv.start)
                                * 100.0
                                / total_len_amp
                            )
                            exon_end_pos = (
                                amplified_intervals_start[chrom][inti]
                                + (cut_ee - graph_intv.start)
                                * 100.0
                                / total_len_amp
                            )

                            exon_min_width = (
                                0.2 * zoom_factor
                            )  # Adjust the minimum width as needed
                            exon_width = exon_end_pos - exon_start_pos
                            if exon_width < exon_min_width:
                                diff = (exon_min_width - exon_width) / 2
                                exon_start_pos -= diff
                                exon_end_pos += diff

                            ax3.hlines(
                                height,
                                exon_start_pos,
                                exon_end_pos,
                                color="black",
                                lw=7.5,
                            )

        # Ticks and labels
        xtickpos = []
        xticklabels = []
        if not self.plot_bounds:
            axis_max = axis_layout.axis_max
            ax.set_xlim(0, axis_max)
            ax2.set_xlim(0, axis_max)
            ax3.set_xlim(0, axis_max)
            set_interval_axis_labels(
                ax3,
                self.graph_amplified_intervals,
                amplified_intervals_start,
                total_len_amp,
                fontsize,
                min_coord_width,
            )

        else:  # self.plot_bounds are given
            # look up the segment
            pchrom, pstart, pend = self.plot_bounds
            nint_chr = len(self.graph_amplified_intervals[pchrom])
            relint = None
            rint_ = None
            for inti in range(nint_chr):
                istart, iend = (
                    self.graph_amplified_intervals[pchrom][inti].start,
                    self.graph_amplified_intervals[pchrom][inti].end,
                )
                if istart <= pstart <= iend:
                    relint = inti
                    rint_ = self.graph_amplified_intervals[pchrom][inti]
                    break

            if relint is None:
                print(
                    f"Could not identify region {pchrom}:{pstart}-{pend} in graph regions. Region should be fully contained in graph.",
                )

            else:
                plot_start = (
                    amplified_intervals_start[pchrom][relint]
                    + (pstart - rint_.start) * 100.0 / total_len_amp
                )
                plot_end = (
                    amplified_intervals_start[pchrom][relint]
                    + (pend - rint_.start) * 100.0 / total_len_amp
                )
                xtickpos.append(plot_start)
                xtickpos.append(plot_end)
                xticklabels.append(pchrom + ":" + str(pstart))
                xticklabels.append(pchrom + ":" + str(pend))
                ax3.set_xticks(xtickpos)
                ax3.set_xticklabels(
                    xticklabels,
                    size=14.0 * font_size_multiplier,
                )

                ax.set_xlim(plot_start, plot_end)
                ax2.set_xlim(plot_start, plot_end)
                ax3.set_xlim(plot_start, plot_end)

        ax3.yaxis.set_major_formatter(ticker.NullFormatter())
        ax3.set_ylim(0, 1)
        fig.subplots_adjust(hspace=0)
        hide_figure_text_if_zero(fig, font_size_multiplier)
        save_plot_figure(fig, output_fn, dpi)

    def close_bam(self):
        try:
            self.lr_bamfh.close()
        except AttributeError:
            pass

    def plot_cycles(
        self,
        title: str,
        output_fn: str,
        num_cycles: int | None = None,
        fontsize: float = DEFAULT_PLOT_FONT_SIZE,
        dpi: int = 300,
        gene_font_size: float = DEFAULT_GENE_FONT_SIZE,
        *,
        cycle_only: bool = False,
        hide_genes: bool = False,
        font_size_multiplier: float = 1.0,
        cycle_ids: list[int] | None = None,
        configured_cycle_colors: dict[int, str] | None = None,
        plot_width: float | None = None,
        aspect_ratio: float | None = None,
        interval_offset: float = DEFAULT_INTERVAL_OFFSET,
        min_coord_width: float = DEFAULT_MIN_COORD_WIDTH,
        align_to_combined: bool = False,
    ) -> None:
        """Plot cycles & paths returned from cycle decomposition"""
        validate_plot_dimensions(
            plot_width,
            aspect_ratio,
            interval_offset,
            min_coord_width,
            dpi,
        )
        cycles_to_plot = (
            list(self.cycles) if cycle_ids is None else list(cycle_ids)
        )
        if cycle_ids is None and num_cycles is not None:
            cycles_to_plot = [
                cycle_id
                for cycle_id in self.cycles
                if int(cycle_id) <= num_cycles
            ]
        if cycle_only:
            cycles_to_plot = [
                cycle_id
                for cycle_id in cycles_to_plot
                if self.cycles[cycle_id].is_cyclic
            ]
        if cycle_ids is None:
            cycles_to_plot = sorted(cycles_to_plot)
        if not cycles_to_plot:
            raise ValueError("no cycles or paths remain after filtering")
        if configured_cycle_colors:
            logger.warning(
                "Ignoring per-item cycle colors; cycle sequence edges and "
                "SV edges now use fixed role colors."
            )
        height = sum(
            [
                3 * len(self.cycles[cycle_id].segments) - 1
                for cycle_id in cycles_to_plot
            ]
        ) + 9 * (len(cycles_to_plot) - 1)
        layout_scale = get_text_layout_scale(font_size_multiplier)
        width, resolved_height = resolve_figure_size(
            default_width=DEFAULT_CYCLE_WIDTH,
            default_height=max(4, height * 0.25) * layout_scale,
            width=plot_width,
            aspect_ratio=aspect_ratio,
        )
        axis_layout = build_genomic_axis_layout(
            self.cycle_amplified_intervals,
            interval_offset,
        )
        margin_between_intervals = get_interval_margin(
            self.num_amplified_intervals,
            interval_offset,
        )
        total_len_amp = axis_layout.total_genomic_length
        amplified_intervals_start = axis_layout.interval_starts
        displayed_axis_span = (
            100
            + (self.num_amplified_intervals + 1) * margin_between_intervals
        )
        gene_lane_count = 0
        if not hide_genes:
            gene_lane_count = self.set_displayed_gene_heights(
                self.cycle_amplified_intervals,
                amplified_intervals_start,
                total_len_amp,
                width,
                gene_font_size,
                displayed_axis_span,
                (
                    COMBINED_AXIS_RIGHT - COMBINED_AXIS_LEFT
                    if align_to_combined
                    else 0.71
                ),
            )
        fig = plt.figure(figsize=(width, resolved_height))
        if not hide_genes:
            gene_track_ratio = get_cycle_gene_track_ratio(
                resolved_height,
                gene_lane_count,
                align_to_combined,
            )
            gs = gridspec.GridSpec(
                2,
                1,
                height_ratios=[PLOT_DATA_HEIGHT_RATIO, gene_track_ratio],
            )
        else:
            gs = gridspec.GridSpec(
                2,
                1,
                height_ratios=[PLOT_DATA_HEIGHT_RATIO, 0.000001],
            )
        ax = fig.add_subplot(gs[0, 0])
        ax.set_title(title, fontsize=fontsize)
        ax.xaxis.set_visible(False)
        ax3 = fig.add_subplot(gs[1, 0], sharex=ax)
        ax3.set_zorder(3)
        ax3.patch.set_visible(False)
        for plot_axis in (ax, ax3):
            scale_axis_elements(plot_axis, font_size_multiplier)
        ax3.yaxis.set_visible(False)
        ax3.spines["left"].set_visible(False)
        ax3.spines["right"].set_visible(False)
        ax3.spines["top"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Compute the x coordinates for each amplified interval
        # sorted_chrs = sorted(self.intervals_from_cycle.keys(), key = lambda chr: CHR_TAG_TO_IDX[chr])
        sorted_chrs = breakpoint_utilities.sort_chrom_names(
            self.cycle_amplified_intervals.keys()
        )

        # Draw amplified interval separators
        draw_interval_guides(ax, ax3, axis_layout)

        # Draw cycles
        y_cur = -2
        extension = 1.5
        cycleticks = []
        cycleticklabels: list[tuple[str, str]] = []
        for cycle_id in cycles_to_plot:
            sv_connection_color = (
                CYCLE_COLOR if self.cycles[cycle_id].is_cyclic else PATH_COLOR
            )
            patch_start = len(ax.patches)
            collection_start = len(ax.collections)
            ystart_cycle_id = y_cur
            cycle_min_x = float("inf")
            cycle_max_x = 0.0
            for i, segment in enumerate(self.cycles[cycle_id].segments):
                # Segment i
                interval_idx = 0
                while (
                    segment.start
                    > self.cycle_amplified_intervals[segment.chr][
                        interval_idx
                    ].end
                ):
                    interval_idx += 1
                x1 = (
                    amplified_intervals_start[segment.chr][interval_idx]
                    + (
                        segment.start
                        - self.cycle_amplified_intervals[segment.chr][
                            interval_idx
                        ].start
                    )
                    * 100.0
                    / total_len_amp
                )
                cycle_min_x = min(x1, cycle_min_x)
                xlen = (segment.end - segment.start) * 100.0 / total_len_amp
                cycle_max_x = max(x1 + xlen, cycle_max_x)
                rect = Rectangle(
                    (x1, y_cur),
                    xlen,
                    1,
                    facecolor=CYCLE_SEGMENT_FACE_COLOR,
                    linewidth=CYCLE_SEGMENT_LINE_WIDTH,
                    edgecolor=CYCLE_SEGMENT_EDGE_COLOR,
                )
                ax.add_patch(rect)

                # Connections between segment i and i + 1
                if i < len(self.cycles[cycle_id].segments) - 1:
                    nseg = self.cycles[cycle_id].segments[i + 1]
                    interval_idx_n = 0
                    while (
                        nseg.start
                        > self.cycle_amplified_intervals[nseg.chr][
                            interval_idx_n
                        ].end
                    ):
                        interval_idx_n += 1
                    if segment.strand == "+" and nseg.strand == "-":
                        x2 = x1 + xlen
                        x2n: float = amplified_intervals_start[nseg.chr][
                            interval_idx_n
                        ]
                        x2n += (
                            (
                                nseg.end
                                - self.cycle_amplified_intervals[nseg.chr][
                                    interval_idx_n
                                ].start
                            )
                            * 100.0
                            / total_len_amp
                        )
                        ax.vlines(
                            x=max(x2, x2n) + extension,
                            ymin=y_cur + 0.5,
                            ymax=y_cur - 1.5,
                            colors="b",
                            lw=2,
                        )
                        ax.hlines(
                            y=y_cur + 0.5,
                            xmin=x2,
                            xmax=max(x2, x2n) + extension,
                            colors="b",
                            lw=2,
                        )
                        ax.hlines(
                            y=y_cur - 1.5,
                            xmin=x2n,
                            xmax=max(x2, x2n) + extension,
                            colors="b",
                            lw=2,
                        )
                        y_cur -= 2
                    elif segment.strand == "-" and nseg.strand == "+":
                        x1n: float = amplified_intervals_start[nseg.chr][
                            interval_idx_n
                        ]
                        x1n += (
                            (
                                nseg.start
                                - self.cycle_amplified_intervals[nseg.chr][
                                    interval_idx_n
                                ].start
                            )
                            * 100.0
                            / total_len_amp
                        )
                        ax.vlines(
                            x=min(x1, x1n) - extension,
                            ymin=y_cur + 0.5,
                            ymax=y_cur - 1.5,
                            colors="b",
                            lw=2,
                        )
                        ax.hlines(
                            y=y_cur + 0.5,
                            xmin=min(x1, x1n) - extension,
                            xmax=x1,
                            colors="b",
                            lw=2,
                        )
                        ax.hlines(
                            y=y_cur - 1.5,
                            xmin=min(x1, x1n) - extension,
                            xmax=x1n,
                            colors="b",
                            lw=2,
                        )
                        y_cur -= 2
                    elif segment.strand == "+" and nseg.strand == "+":
                        x2 = x1 + xlen
                        x1n: float = amplified_intervals_start[nseg.chr][
                            interval_idx_n
                        ]
                        x1n += (
                            (
                                nseg.start
                                - self.cycle_amplified_intervals[nseg.chr][
                                    interval_idx_n
                                ].start
                            )
                            * 100.0
                            / total_len_amp
                        )
                        if x2 <= x1n:
                            ax.hlines(
                                y=y_cur + 0.5,
                                xmin=x2,
                                xmax=x1n,
                                colors="b",
                                lw=2,
                            )
                        else:
                            ax.vlines(
                                x=x2 + extension,
                                ymin=y_cur - 0.5,
                                ymax=y_cur + 0.5,
                                colors="b",
                                lw=2,
                            )
                            ax.vlines(
                                x=x1n - extension,
                                ymin=y_cur - 1.5,
                                ymax=y_cur - 0.5,
                                colors="b",
                                lw=2,
                            )
                            ax.hlines(
                                y=y_cur + 0.5,
                                xmin=x2,
                                xmax=x2 + extension,
                                colors="b",
                                lw=2,
                            )
                            ax.hlines(
                                y=y_cur - 0.5,
                                xmin=x1n - extension,
                                xmax=x2 + extension,
                                colors="b",
                                lw=2,
                            )
                            ax.hlines(
                                y=y_cur - 1.5,
                                xmin=x1n - extension,
                                xmax=x1n,
                                colors="b",
                                lw=2,
                            )
                            y_cur -= 2
                    else:  # seg[3] == '-' and nseg[3] == '-'
                        x2n: float = amplified_intervals_start[nseg.chr][
                            interval_idx_n
                        ]
                        x2n += (
                            (
                                nseg.end
                                - self.cycle_amplified_intervals[nseg.chr][
                                    interval_idx_n
                                ].start
                            )
                            * 100.0
                            / total_len_amp
                        )
                        if x1 >= x2n:
                            ax.hlines(
                                y=y_cur + 0.5,
                                xmin=x2n,
                                xmax=x1,
                                colors="b",
                                lw=2,
                            )
                        else:
                            ax.vlines(
                                x=x1 - extension,
                                ymin=y_cur - 0.5,
                                ymax=y_cur + 0.5,
                                colors="b",
                                lw=2,
                            )
                            ax.vlines(
                                x=x2n + extension,
                                ymin=y_cur - 1.5,
                                ymax=y_cur - 0.5,
                                colors="b",
                                lw=2,
                            )
                            ax.hlines(
                                y=y_cur + 0.5,
                                xmin=x1 - extension,
                                xmax=x1,
                                colors="b",
                                lw=2,
                            )
                            ax.hlines(
                                y=y_cur - 0.5,
                                xmin=x1 - extension,
                                xmax=x2n + extension,
                                colors="b",
                                lw=2,
                            )
                            ax.hlines(
                                y=y_cur - 1.5,
                                xmin=x2n,
                                xmax=x2n + extension,
                                colors="b",
                                lw=2,
                            )
                            y_cur -= 2

            # First and last segments
            if not self.cycles[cycle_id].is_cyclic:  # Paths
                seg = self.cycles[cycle_id].segments[0]
                interval_idx = 0
                while (
                    seg.start
                    > self.cycle_amplified_intervals[seg.chr][interval_idx].end
                ):
                    interval_idx += 1
                if seg.strand == "+":
                    x1: float = (
                        amplified_intervals_start[seg.chr][interval_idx]
                        + (
                            seg.start
                            - self.cycle_amplified_intervals[seg.chr][
                                interval_idx
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    ax.hlines(
                        y=ystart_cycle_id + 0.5,
                        xmin=x1 - 2 * extension,
                        xmax=x1,
                        colors="b",
                        lw=2,
                    )
                else:
                    x2: float = (
                        amplified_intervals_start[seg.chr][interval_idx]
                        + (
                            seg.end
                            - self.cycle_amplified_intervals[seg.chr][
                                interval_idx
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    ax.hlines(
                        y=ystart_cycle_id + 0.5,
                        xmin=x2,
                        xmax=x2 + 2 * extension,
                        colors="b",
                        lw=2,
                    )
                seg = self.cycles[cycle_id].segments[-1]
                interval_idx = 0
                while (
                    seg.start
                    > self.cycle_amplified_intervals[seg.chr][interval_idx].end
                ):
                    interval_idx += 1
                if seg.strand == "+":
                    x2: float = amplified_intervals_start[seg.chr][interval_idx]
                    x2 += (
                        (
                            seg.end
                            - self.cycle_amplified_intervals[seg.chr][
                                interval_idx
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    ax.hlines(
                        y=y_cur + 0.5,
                        xmin=x2,
                        xmax=x2 + 2 * extension,
                        colors="b",
                        lw=2,
                    )
                else:
                    x1: float = amplified_intervals_start[seg.chr][interval_idx]
                    x1 += (
                        (
                            seg.start
                            - self.cycle_amplified_intervals[seg.chr][
                                interval_idx
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    ax.hlines(
                        y=y_cur + 0.5,
                        xmin=x1 - 2 * extension,
                        xmax=x1,
                        colors="b",
                        lw=2,
                    )
            else:  # Cycles
                # xmin_ = 0.5
                xmin_ = cycle_min_x - extension
                xmax_ = cycle_max_x + extension

                if len(self.cycles[cycle_id].segments) > 1:
                    xmin_ -= extension
                    xmax_ += extension

                # xmax_ = 99.5 + (self.num_amplified_intervals + 1) * margin_between_intervals
                seg1 = self.cycles[cycle_id].segments[0]
                interval_idx1 = 0
                while (
                    seg1.start
                    > self.cycle_amplified_intervals[seg1.chr][
                        interval_idx1
                    ].end
                ):
                    interval_idx1 += 1
                seg2 = self.cycles[cycle_id].segments[-1]
                interval_idx2 = 0
                while (
                    seg2.start
                    > self.cycle_amplified_intervals[seg2.chr][
                        interval_idx2
                    ].end
                ):
                    interval_idx2 += 1
                if seg1.strand == "-" and seg2.strand == "+":
                    x2: float = amplified_intervals_start[seg1.chr][
                        interval_idx1
                    ]
                    x2 += (
                        (
                            seg1.end
                            - self.cycle_amplified_intervals[seg1.chr][
                                interval_idx1
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    x2n: float = amplified_intervals_start[seg2.chr][
                        interval_idx2
                    ]
                    x2n += (
                        (
                            seg2.end
                            - self.cycle_amplified_intervals[seg2.chr][
                                interval_idx2
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    ax.vlines(
                        x=xmax_,
                        ymin=y_cur + 0.5,
                        ymax=ystart_cycle_id + 0.5,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=ystart_cycle_id + 0.5,
                        xmin=x2,
                        xmax=xmax_,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=y_cur + 0.5, xmin=x2n, xmax=xmax_, colors="b", lw=2
                    )
                elif seg1.strand == "+" and seg2.strand == "-":
                    x1: float = amplified_intervals_start[seg1.chr][
                        interval_idx1
                    ]
                    x1 += (
                        (
                            seg1.start
                            - self.cycle_amplified_intervals[seg1.chr][
                                interval_idx1
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    x1n: float = amplified_intervals_start[seg2.chr][
                        interval_idx2
                    ]
                    x1n += (
                        (
                            seg2.start
                            - self.cycle_amplified_intervals[seg2.chr][
                                interval_idx2
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    ax.vlines(
                        x=xmin_,
                        ymin=y_cur + 0.5,
                        ymax=ystart_cycle_id + 0.5,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=ystart_cycle_id + 0.5,
                        xmin=xmin_,
                        xmax=x1,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=y_cur + 0.5, xmin=xmin_, xmax=x1n, colors="b", lw=2
                    )
                elif seg1.strand == "-" and seg2.strand == "-":
                    x2: float = amplified_intervals_start[seg1.chr][
                        interval_idx1
                    ]
                    x2 += (
                        (
                            seg1.end
                            - self.cycle_amplified_intervals[seg1.chr][
                                interval_idx1
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    x1n: float = amplified_intervals_start[seg2.chr][
                        interval_idx2
                    ]
                    x1n += (
                        (
                            seg2.end
                            - self.cycle_amplified_intervals[seg2.chr][
                                interval_idx2
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    ax.vlines(
                        x=xmax_,
                        ymin=y_cur - 0.5,
                        ymax=ystart_cycle_id + 0.5,
                        colors="b",
                        lw=2,
                    )
                    ax.vlines(
                        x=x1n - extension,
                        ymin=y_cur - 0.5,
                        ymax=y_cur + 0.5,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=ystart_cycle_id + 0.5,
                        xmin=x2,
                        xmax=xmax_,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=y_cur + 0.5,
                        xmin=x1n - extension,
                        xmax=x1n,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=y_cur - 0.5,
                        xmin=x1n - extension,
                        xmax=xmax_,
                        colors="b",
                        lw=2,
                    )
                else:
                    x1: float = amplified_intervals_start[seg1.chr][
                        interval_idx1
                    ]
                    x1 += (
                        (
                            seg1.start
                            - self.cycle_amplified_intervals[seg1.chr][
                                interval_idx1
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    x2n: float = amplified_intervals_start[seg2.chr][
                        interval_idx2
                    ]
                    x2n += (
                        (
                            seg2.end
                            - self.cycle_amplified_intervals[seg2.chr][
                                interval_idx2
                            ].start
                        )
                        * 100.0
                        / total_len_amp
                    )
                    ax.vlines(
                        x=xmin_,
                        ymin=y_cur - 0.5,
                        ymax=ystart_cycle_id + 0.5,
                        colors="b",
                        lw=2,
                    )
                    ax.vlines(
                        x=x2n + extension,
                        ymin=y_cur - 0.5,
                        ymax=y_cur + 0.5,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=ystart_cycle_id + 0.5,
                        xmin=xmin_,
                        xmax=x1,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=y_cur + 0.5,
                        xmin=x2n,
                        xmax=x2n + extension,
                        colors="b",
                        lw=2,
                    )
                    ax.hlines(
                        y=y_cur - 0.5,
                        xmin=xmin_,
                        xmax=x2n + extension,
                        colors="b",
                        lw=2,
                    )

            for collection in ax.collections[collection_start:]:
                collection.set_color(sv_connection_color)
                collection.set_linewidth(CYCLE_CONNECTION_LINE_WIDTH)
            for patch in ax.patches[patch_start:]:
                patch.set_facecolor(CYCLE_SEGMENT_FACE_COLOR)
                patch.set_edgecolor(CYCLE_SEGMENT_EDGE_COLOR)
                patch.set_linewidth(CYCLE_SEGMENT_LINE_WIDTH)

            # Separators between cycles; ticks
            ax.hlines(
                y=y_cur - 2,
                xmin=-1,
                xmax=101
                + (self.num_amplified_intervals + 1) * margin_between_intervals,
                colors="0.45",
                linewidth=0.8,
            )
            cycleticks.append((y_cur + ystart_cycle_id) * 0.5)
            if self.cycles[cycle_id].is_cyclic:
                cycleticklabels.append(
                    (
                        f"Cycle {cycle_id}",
                        f"CN = {self.cycles[cycle_id].overall_cn:.2f}",
                    )
                )
            else:
                cycleticklabels.append(
                    (
                        f"Path {cycle_id}",
                        f"CN = {self.cycles[cycle_id].overall_cn:.2f}",
                    )
                )
            y_cur -= 4

        if not hide_genes:
            for chrom in sorted_chrs:
                for int_idx in range(
                    len(self.cycle_amplified_intervals[chrom])
                ):
                    intv = self.cycle_amplified_intervals[chrom][int_idx]
                    rel_genes = [
                        x.data for x in self.genes[chrom][intv.start : intv.end]
                    ]
                    for gene_obj in rel_genes:  # plot gene lines
                        height = gene_obj.height
                        cut_gs = max(intv.start, gene_obj.gstart)
                        cut_ge = min(intv.end, gene_obj.gend)
                        gene_start = (
                            amplified_intervals_start[chrom][int_idx]
                            + (cut_gs - intv.start) * 100.0 / total_len_amp
                        )
                        gene_end = (
                            amplified_intervals_start[chrom][int_idx]
                            + (cut_ge - intv.start) * 100.0 / total_len_amp
                        )
                        ax3.hlines(
                            height,
                            gene_start,
                            gene_end,
                            color="cornflowerblue",
                            lw=4.5,
                        )  # Draw horizontal bars for genes
                        if gene_font_size > 0:
                            ax3.text(
                                (gene_start + gene_end) / 2,
                                height + 0.05,
                                gene_obj.gname,
                                ha="center",
                                va="bottom",
                                fontsize=gene_font_size,
                                style="italic",
                                clip_on=False,
                                zorder=10,
                            )

                        if gene_obj.strand == "+":
                            ax3.plot(
                                gene_start,
                                height,
                                marker=">",
                                color="black",
                                markersize=7,
                            )
                        elif gene_obj.strand == "-":
                            ax3.plot(
                                gene_end,
                                height,
                                marker="<",
                                color="black",
                                markersize=7,
                            )

                        for (
                            exon_start,
                            exon_end,
                        ) in gene_obj.eposns:  # plot exon bars
                            if (
                                not exon_end > intv.start
                                or not exon_start < intv.end
                            ):
                                continue

                            cut_es = max(intv.start, exon_start)
                            cut_ee = min(intv.end, exon_end)
                            exon_start_pos = (
                                amplified_intervals_start[chrom][int_idx]
                                + (cut_es - intv.start) * 100.0 / total_len_amp
                            )
                            exon_end_pos = (
                                amplified_intervals_start[chrom][int_idx]
                                + (cut_ee - intv.start) * 100.0 / total_len_amp
                            )

                            exon_min_width = (
                                0.2  # Adjust the minimum width as needed
                            )
                            exon_width = exon_end_pos - exon_start_pos
                            if exon_width < exon_min_width:
                                diff = (exon_min_width - exon_width) / 2
                                exon_start_pos -= diff
                                exon_end_pos += diff

                            ax3.hlines(
                                height,
                                exon_start_pos,
                                exon_end_pos,
                                color="black",
                                lw=7.5,
                            )

        # Ticks and labels
        axis_max = axis_layout.axis_max
        x_padding = get_cycle_x_padding(align_to_combined, extension)
        ax.set_xlim(-x_padding, axis_max + x_padding)
        ax.set_ylim(y_cur + 2, 0)
        ax3.set_xlim(-x_padding, axis_max + x_padding)
        set_interval_axis_labels(
            ax3,
            self.cycle_amplified_intervals,
            amplified_intervals_start,
            total_len_amp,
            fontsize,
            min_coord_width,
        )

        add_cycle_axis_labels(ax, cycleticks, cycleticklabels, fontsize)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)

        ax3.yaxis.set_major_formatter(ticker.NullFormatter())
        ax3.set_ylim(0, 1)
        bottom_margin, top_margin = get_cycle_figure_margins(
            resolved_height,
            hide_genes,
            align_to_combined,
        )
        fig.subplots_adjust(
            left=COMBINED_AXIS_LEFT,
            right=COMBINED_AXIS_RIGHT if align_to_combined else 0.94,
            bottom=bottom_margin,
            top=top_margin,
            hspace=0,
        )
        hide_figure_text_if_zero(fig, font_size_multiplier)
        save_plot_figure(fig, output_fn, dpi)


@core_utils.profile_fn_with_call_counter
def plot_amplicon(
    ref: core_types.ReferenceGenome,
    bam_path: pathlib.Path | None,
    graph_file: io.TextIOWrapper | None,
    cycle_file: io.TextIOWrapper | None,
    output_prefix: str,
    num_cycles: int | None,
    max_coverage: float,
    min_mapq: float,
    gene_subset_list: list[str],
    gene_fontsize: float,
    region: str | None,
    *,
    should_plot_graph: bool,
    should_plot_cycles: bool,
    should_hide_genes: bool,
    should_restrict_to_bushman_genes: bool,
    should_plot_only_cyclic_walks: bool,
    font_size_multiplier: float = 1.0,
    refgene_file: pathlib.Path | None = None,
    gene_subset_file: pathlib.Path | None = None,
    legend_output_prefix: str | None = None,
    plot_width: float | None = None,
    aspect_ratio: float | None = None,
    interval_offset: float = DEFAULT_INTERVAL_OFFSET,
    min_coord_width: float = DEFAULT_MIN_COORD_WIDTH,
    coverage_scale: Literal["robust", "full"] = "robust",
    cycle_list: str | None = None,
    cycle_color_file: pathlib.Path | None = None,
    combine_graph_cycles: bool = False,
    dpi: int = 300,
) -> None:
    validate_plot_dimensions(
        plot_width,
        aspect_ratio,
        interval_offset,
        min_coord_width,
        dpi,
    )
    cycle_selections = parse_cycle_selection(cycle_list)
    if cycle_selections and num_cycles is not None:
        raise ValueError("--cycle-list and --num-cycles are mutually exclusive")
    if should_plot_only_cyclic_walks and any(
        not selection.expected_is_cyclic for selection in cycle_selections
    ):
        raise ValueError(
            "path selections cannot be combined with --only-cyclic-paths"
        )
    if combine_graph_cycles and not (should_plot_graph and should_plot_cycles):
        raise ValueError(
            "combined output requires both a graph and a cycles file"
        )
    if should_plot_graph:
        if not graph_file:
            print("Please specify the breakpoint graph file to plot.")
            sys.exit(1)

    if should_plot_cycles and not cycle_file:
        print("Please specify the cycle file, in *.bed format, to plot.")
        sys.exit(1)

    g = GraphViz()
    effective_gene_font_size = get_gene_font_size(
        font_size_multiplier,
        gene_fontsize,
    )
    effective_plot_font_size = get_gene_font_size(
        font_size_multiplier,
        DEFAULT_PLOT_FONT_SIZE,
    )
    merged_gene_subset = merge_gene_subsets(gene_subset_list, gene_subset_file)
    g.parse_genes(
        ref,
        set(merged_gene_subset),
        should_restrict_to_bushman_genes,
        refgene_file,
    )
    if should_plot_graph:
        bp_graph = parse_breakpoint_graph(graph_file)  # type: ignore[arg-type]
        if bam_path:
            g.open_bam(bam_path)
        g.graph = bp_graph
        if region:
            pchrom = region.split(":")[0]
            pb1, pb2 = region.split(":")[1].rsplit("-")
            g.plot_bounds = (pchrom, int(pb1), int(pb2))
        g.update_graph_intervals()
        gtitle = output_prefix
        if "/" in output_prefix:
            gtitle = output_prefix.split("/")[-1]
        graph_title = (
            f"{gtitle}-Amplicon Graph" if combine_graph_cycles else gtitle
        )
        g.plot_graph(
            graph_title,
            output_prefix + "_graph",
            max_cov_cutoff=max_coverage,
            quality_threshold=min_mapq,
            hide_genes=should_hide_genes,
            gene_font_size=effective_gene_font_size,
            fontsize=effective_plot_font_size,
            font_size_multiplier=font_size_multiplier,
            plot_width=plot_width,
            aspect_ratio=aspect_ratio,
            interval_offset=interval_offset,
            min_coord_width=min_coord_width,
            coverage_scale=coverage_scale,
            align_to_combined=combine_graph_cycles,
            dpi=dpi,
        )
        write_graph_legend(
            output_prefix
            if legend_output_prefix is None
            else legend_output_prefix,
            get_graph_coverage_label(bam_path),
            font_size_multiplier=font_size_multiplier,
            dpi=dpi,
        )

    if should_plot_cycles:
        recon_cycles = parse_cycle_file(
            cycle_file,
            output_prefix,
            None if cycle_selections else num_cycles,
        )
        for cycle_id, cycle in recon_cycles.items():
            g.cycles[cycle_id] = cycle
        cycle_ids_ = resolve_cycle_ids(g.cycles, cycle_selections)
        if should_plot_only_cyclic_walks:
            cycle_ids_ = [
                cycle_id
                for cycle_id in cycle_ids_
                if g.cycles[cycle_id].is_cyclic
            ]
        configured_cycle_colors: dict[int, str] = {}
        if cycle_color_file is not None:
            logger.warning(
                "Ignoring --cycle-color-file; cycle sequence edges and SV "
                "edges now use fixed role colors."
            )

        graph_given_ = graph_file is not None
        if graph_given_:
            g.graph = bp_graph
            if not g.graph_amplified_intervals:
                g.update_graph_intervals()
        g.update_cycle_amplified_intervals(
            cycle_ids=cycle_ids_,
            cycle_only=should_plot_only_cyclic_walks,
            graph_given=graph_given_,
            align_to_graph=combine_graph_cycles,
        )
        gtitle = output_prefix
        if "/" in output_prefix:
            gtitle = output_prefix.split("/")[-1]
        cycles_title = f"{gtitle}-Cycles" if combine_graph_cycles else gtitle
        g.plot_cycles(
            cycles_title,
            output_prefix + "_cycles",
            num_cycles=num_cycles,
            cycle_only=should_plot_only_cyclic_walks,
            hide_genes=should_hide_genes,
            gene_font_size=effective_gene_font_size,
            fontsize=effective_plot_font_size,
            font_size_multiplier=font_size_multiplier,
            cycle_ids=cycle_ids_,
            configured_cycle_colors=configured_cycle_colors,
            plot_width=plot_width,
            aspect_ratio=aspect_ratio,
            interval_offset=interval_offset,
            min_coord_width=min_coord_width,
            align_to_combined=combine_graph_cycles,
            dpi=dpi,
        )
    if combine_graph_cycles:
        combine_graph_and_cycle_plots(
            output_prefix,
            dpi,
            build_genomic_axis_layout(
                g.graph_amplified_intervals,
                interval_offset,
            ),
            hide_genes=should_hide_genes,
        )
    g.close_bam()
    if graph_file:
        print(
            f"Visualization completed for {colorama.Fore.LIGHTCYAN_EX}"
            f"{graph_file.name} ({cycle_file.name if cycle_file else 'no cycles'}"  # type: ignore[union-attr]
            f"{colorama.Style.RESET_ALL})"
        )
    elif cycle_file:
        print(
            f"Visualization completed for {colorama.Fore.LIGHTCYAN_EX}{cycle_file.name}"
        )


def draw_cycle(
    g: BreakpointGraph,
    cycle: list[int],
    cycle_id: int,
    cycle_only: bool,
    graph_given: bool,
) -> None:
    pass
