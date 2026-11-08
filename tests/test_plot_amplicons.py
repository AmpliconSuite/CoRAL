from __future__ import annotations

import math
import pathlib
from types import SimpleNamespace

from matplotlib import pyplot as plt
from typer.testing import CliRunner

from coral import core_types, plot_amplicons


def test_parse_gene_subset_file_accepts_text_and_csv(
    tmp_path: pathlib.Path,
) -> None:
    gene_file = tmp_path / "genes.csv"
    gene_file.write_text("EGFR\nMYC,MDM2\n  CDK4  CDKN2A\n")

    assert plot_amplicons.parse_gene_subset_file(gene_file) == [
        "EGFR",
        "MYC",
        "MDM2",
        "CDK4",
        "CDKN2A",
    ]


def test_merge_gene_subsets_dedupes_cli_and_file(
    tmp_path: pathlib.Path,
) -> None:
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("MYC\nEGFR\n")

    assert plot_amplicons.merge_gene_subsets(["TP53", "MYC"], gene_file) == [
        "TP53",
        "MYC",
        "EGFR",
    ]


def test_empty_gene_subset_file_warns_and_defaults_to_all(
    tmp_path: pathlib.Path,
    capsys: object,
) -> None:
    gene_file = tmp_path / "empty_genes.txt"
    gene_file.write_text("\n,   \n")

    assert plot_amplicons.merge_gene_subsets([], gene_file) == []
    assert "empty; plotting all genes" in capsys.readouterr().err


def test_discordant_edge_linewidth_scales_with_read_count() -> None:
    assert plot_amplicons.get_discordant_edge_linewidth(1.0, 4.0) == 1.0
    assert plot_amplicons.get_discordant_edge_linewidth(2.0, 4.0) == 2.0
    assert plot_amplicons.get_discordant_edge_linewidth(4.0, 4.0) == 4.0
    assert plot_amplicons.get_discordant_edge_linewidth(10.0, 10.0) == 4.0
    assert plot_amplicons.get_discordant_edge_linewidth(4.0, 0.0) == 0.30


def test_discordant_edge_arcs_follow_aa_plotted_distance_convention() -> None:
    max_segment_cn = 4.0

    assert plot_amplicons.get_discordant_edge_arc_base(max_segment_cn) == 0.0
    short_arc = plot_amplicons.get_discordant_edge_arc_height(
        10.0,
        1000.0,
        max_segment_cn,
    )
    long_arc = plot_amplicons.get_discordant_edge_arc_height(
        900.0,
        1000.0,
        max_segment_cn,
    )
    full_width_arc = plot_amplicons.get_discordant_edge_arc_height(
        1000.0,
        1000.0,
        max_segment_cn,
    )

    assert math.isclose(
        short_arc,
        max_segment_cn
        * (1.25 + (plot_amplicons.MAX_ARC_APEX_CN_SCALE - 1.25) * 0.01),
    )
    assert short_arc < long_arc
    assert math.isclose(
        long_arc,
        max_segment_cn
        * (1.25 + (plot_amplicons.MAX_ARC_APEX_CN_SCALE - 1.25) * 0.9),
    )
    assert math.isclose(
        full_width_arc,
        plot_amplicons.MAX_ARC_APEX_CN_SCALE * max_segment_cn,
    )


def test_interval_offset_is_a_total_gap_budget() -> None:
    for interval_count in (1, 5, 40):
        margin = plot_amplicons.get_interval_margin(interval_count, 0.1)
        total_gap = margin * (interval_count + 1)
        assert math.isclose(total_gap / (100.0 + total_gap), 0.1)


def test_figure_size_supports_width_and_height_to_width_ratio() -> None:
    assert plot_amplicons.resolve_figure_size(
        default_width=12.0,
        default_height=5.0,
        width=None,
        aspect_ratio=None,
    ) == (12.0, 5.0)
    assert plot_amplicons.resolve_figure_size(
        default_width=12.0,
        default_height=5.0,
        width=8.0,
        aspect_ratio=0.5,
    ) == (8.0, 4.0)


def test_cycle_layout_keeps_closing_edges_and_tall_figures_compact() -> None:
    assert plot_amplicons.get_cycle_x_padding(False, 1.5) == 4.5
    assert plot_amplicons.get_cycle_x_padding(True, 1.5) == 0.0
    assert plot_amplicons.CYCLE_GENE_TRACK_BASE_HEIGHT_INCHES == 0.8
    assert plot_amplicons.CYCLE_GENE_TRACK_HEIGHT_PER_LANE_INCHES == 0.4
    assert plot_amplicons.CYCLE_GENE_TRACK_MAX_HEIGHT_INCHES == 3.0
    assert plot_amplicons.CYCLE_GENE_TRACK_MAX_FIGURE_FRACTION == 0.45

    sparse_ratio = plot_amplicons.get_cycle_gene_track_ratio(24.0, 1, False)
    crowded_ratio = plot_amplicons.get_cycle_gene_track_ratio(
        24.0,
        5,
        False,
    )
    sparse_height = 24.0 * sparse_ratio / (
        plot_amplicons.PLOT_DATA_HEIGHT_RATIO + sparse_ratio
    )
    crowded_height = 24.0 * crowded_ratio / (
        plot_amplicons.PLOT_DATA_HEIGHT_RATIO + crowded_ratio
    )
    assert math.isclose(sparse_height, 1.2)
    assert math.isclose(crowded_height, 2.8)
    assert sparse_ratio < crowded_ratio

    gene_ratio = plot_amplicons.get_cycle_gene_track_ratio(53.0, 6, False)
    gene_height = 53.0 * gene_ratio / (
        plot_amplicons.PLOT_DATA_HEIGHT_RATIO + gene_ratio
    )
    assert math.isclose(gene_height, 3.0)

    bottom, top = plot_amplicons.get_cycle_figure_margins(
        53.0,
        False,
        False,
    )
    assert math.isclose(bottom * 53.0, 2.0)
    assert math.isclose((1.0 - top) * 53.0, 0.75)


def test_interval_coordinates_use_smaller_font_than_chromosomes() -> None:
    fig, ax = plt.subplots()
    try:
        plot_amplicons.set_interval_axis_labels(
            ax,
            {"chr1": [plot_amplicons.datatypes.Interval("chr1", 10, 110)]},
            {"chr1": [0.0]},
            100.0,
            18.0,
            0.0,
        )
        assert {label.get_fontsize() for label in ax.get_xticklabels()} == {
            18.0
        }
        assert {
            label.get_fontsize()
            for label in ax.get_xticklabels(minor=True)
        } == {18.0 * plot_amplicons.GENOMIC_COORDINATE_FONT_SCALE}
    finally:
        plt.close(fig)


def test_plot_dimension_validation_rejects_invalid_values() -> None:
    invalid_arguments = (
        (0.0, None, 0.1, 0.0, 300),
        (None, 0.0, 0.1, 0.0, 300),
        (None, None, 1.0, 0.0, 300),
        (None, None, 0.1, 1.1, 300),
        (None, None, 0.1, 0.0, 0),
    )
    for arguments in invalid_arguments:
        try:
            plot_amplicons.validate_plot_dimensions(*arguments)
        except ValueError:
            pass
        else:
            raise AssertionError(f"accepted invalid plot options: {arguments}")


def test_robust_coverage_scale_ignores_isolated_spike() -> None:
    coverage = [10.0] * 99 + [1000.0]
    assert plot_amplicons.get_coverage_scale_max(coverage, "robust") == 10.0
    assert plot_amplicons.get_coverage_scale_max(coverage, "full") == 1000.0


def test_cycle_selection_parses_order_and_validates_type() -> None:
    selections = plot_amplicons.parse_cycle_selection("p1, c2 p1")
    assert selections == [
        plot_amplicons.CycleSelection(1, False),
        plot_amplicons.CycleSelection(2, True),
    ]
    cycles = {
        1: SimpleNamespace(is_cyclic=False),
        2: SimpleNamespace(is_cyclic=True),
    }
    assert plot_amplicons.resolve_cycle_ids(cycles, selections) == [1, 2]

    try:
        plot_amplicons.resolve_cycle_ids(
            cycles,
            [plot_amplicons.CycleSelection(1, True)],
        )
    except ValueError as error:
        assert "is a path" in str(error)
    else:
        raise AssertionError("accepted a path selected as a cycle")


def test_cycle_colors_are_assigned_by_walk_type(
    tmp_path: pathlib.Path,
) -> None:
    color_file = tmp_path / "cycle_colors.tsv"
    color_file.write_text("ITEM COLOR\np1 #0072B2\nc2 darkorange\n")

    assert plot_amplicons.parse_cycle_color_file(color_file) == {
        1: "#0072B2",
        2: "darkorange",
    }
    cycles = {
        1: SimpleNamespace(is_cyclic=False),
        2: SimpleNamespace(is_cyclic=True),
        3: SimpleNamespace(is_cyclic=True),
    }
    assert plot_amplicons.assign_cycle_colors([2, 1, 3], cycles) == {
        2: plot_amplicons.CYCLE_COLOR,
        1: plot_amplicons.PATH_COLOR,
        3: plot_amplicons.CYCLE_COLOR,
    }
    assert plot_amplicons.CYCLE_COLOR == "#5F78B5"


def test_default_cycle_plot_style_matches_publication_figures() -> None:
    style = plot_amplicons.DEFAULT_CYCLE_PLOT_STYLE

    assert style == plot_amplicons.CyclePlotStyle(
        segment_line_width=1.8,
        connection_line_width=1.5,
        cycle_color="#5F78B5",
        path_color="#D55E00",
        segment_face_color="#F2D2A2",
        segment_edge_color="#3F3F3F",
    )
    assert plot_amplicons.CYCLE_SEGMENT_LINE_WIDTH == style.segment_line_width
    assert (
        plot_amplicons.CYCLE_CONNECTION_LINE_WIDTH
        == style.connection_line_width
    )
    assert plot_amplicons.CYCLE_COLOR == style.cycle_color
    assert plot_amplicons.PATH_COLOR == style.path_color
    assert plot_amplicons.CYCLE_SEGMENT_FACE_COLOR == style.segment_face_color
    assert plot_amplicons.CYCLE_SEGMENT_EDGE_COLOR == style.segment_edge_color


def test_cycle_axis_labels_use_regular_heading_and_larger_cn() -> None:
    fig, ax = plot_amplicons.plt.subplots()
    try:
        labels = plot_amplicons.add_cycle_axis_labels(
            ax,
            [1.0],
            [("Cycle 1", "CN = 15.14")],
            fontsize=10.0,
        )

        heading, cn_label = labels
        assert heading.get_text() == "Cycle 1"
        assert heading.get_fontweight() == "normal"
        assert cn_label.get_text() == "CN = 15.14"
        assert cn_label.get_fontweight() == "normal"
        assert cn_label.get_fontsize() == 10.0
        assert heading.get_fontsize() > cn_label.get_fontsize()
    finally:
        plot_amplicons.plt.close(fig)


def test_genomic_axis_layout_uses_one_coordinate_transform() -> None:
    intervals = {
        "chr1": [
            plot_amplicons.datatypes.Interval("chr1", 100, 200),
            plot_amplicons.datatypes.Interval("chr1", 300, 500),
        ],
        "chr2": [plot_amplicons.datatypes.Interval("chr2", 50, 150)],
    }
    layout = plot_amplicons.build_genomic_axis_layout(intervals, 0.08)

    assert len(layout.interval_guides) == 1
    assert len(layout.chromosome_guides) == 1
    assert math.isclose(
        plot_amplicons.genomic_position_to_axis(
            425,
            "chr1",
            intervals,
            layout,
        ),
        layout.interval_starts["chr1"][1]
        + 125 * 100.0 / layout.total_genomic_length,
    )


def test_interval_guides_reach_the_coordinate_axis() -> None:
    intervals = {
        "chr1": [
            plot_amplicons.datatypes.Interval("chr1", 100, 200),
            plot_amplicons.datatypes.Interval("chr1", 300, 400),
        ]
    }
    layout = plot_amplicons.build_genomic_axis_layout(intervals, 0.08)
    fig, (data_axis, coordinate_axis) = plot_amplicons.plt.subplots(2, 1)
    try:
        plot_amplicons.draw_interval_guides(
            data_axis,
            coordinate_axis,
            layout,
        )

        assert len(data_axis.lines) == 1
        assert len(coordinate_axis.lines) == 1
        assert coordinate_axis.lines[0].get_linewidth() == (
            plot_amplicons.INTERVAL_GUIDE_LINE_WIDTH
        )
    finally:
        plot_amplicons.plt.close(fig)


def test_minimum_coordinate_width_filters_only_endpoint_labels() -> None:
    fig, ax = plot_amplicons.plt.subplots()
    try:
        intervals = {
            "chr1": [
                plot_amplicons.datatypes.Interval("chr1", 100, 110),
                plot_amplicons.datatypes.Interval("chr1", 200, 290),
            ]
        }
        plot_amplicons.set_interval_axis_labels(
            ax,
            intervals,
            {"chr1": [0.0, 20.0]},
            total_genomic_length=100.0,
            fontsize=8.0,
            min_coord_width=0.2,
        )

        assert [label.get_text() for label in ax.get_xticklabels()] == ["chr1"]
        assert [
            label.get_text() for label in ax.get_xticklabels(minor=True)
        ] == ["200", "290"]
    finally:
        plot_amplicons.plt.close(fig)


def test_graph_axis_limits_preserve_coral_fit_when_arcs_expand_plot() -> None:
    limits = plot_amplicons.get_graph_axis_limits(
        max_coverage=40.0,
        max_segment_cn=4.0,
        cn_sum_squares=20.0,
        cn_coverage_cross_product=200.0,
        max_arc_apex=8.0,
    )

    assert math.isclose(limits.cn_ymax, 8.4)
    assert math.isclose(limits.coverage_ymax, 84.0)
    assert math.isclose(limits.expansion_factor, 1.68)
    assert math.isclose(
        limits.coverage_ymax / limits.cn_ymax,
        10.0,
    )
    assert math.isclose(40.0 / limits.coverage_ymax, 4.0 / limits.cn_ymax)


def test_graph_axis_limits_keep_original_fit_without_arcs() -> None:
    limits = plot_amplicons.get_graph_axis_limits(
        max_coverage=40.0,
        max_segment_cn=4.0,
        cn_sum_squares=20.0,
        cn_coverage_cross_product=200.0,
        max_arc_apex=0.0,
    )

    assert math.isclose(limits.coverage_ymax, 50.0)
    assert math.isclose(limits.cn_ymax, 5.0)
    assert math.isclose(limits.expansion_factor, 1.0)


def test_cn_ticks_align_with_coverage_ticks() -> None:
    fig, coverage_axis = plot_amplicons.plt.subplots()
    cn_axis = coverage_axis.twinx()
    limits = plot_amplicons.GraphAxisLimits(700.0, 45.0, 1.0)
    try:
        coverage_axis.set_ylim(0, limits.coverage_ymax)
        coverage_axis.set_yticks([0.0, 200.0, 400.0, 600.0])
        cn_axis.set_ylim(0, limits.cn_ymax)

        plot_amplicons.align_cn_axis_ticks(
            coverage_axis,
            cn_axis,
            limits,
        )

        expected_ticks = [
            0.0,
            200.0 * 45.0 / 700.0,
            400.0 * 45.0 / 700.0,
            600.0 * 45.0 / 700.0,
        ]
        assert all(
            math.isclose(observed, expected)
            for observed, expected in zip(
                cn_axis.get_yticks(), expected_ticks, strict=True
            )
        )
    finally:
        plot_amplicons.plt.close(fig)


def test_graph_axis_limits_handle_missing_coverage_fit() -> None:
    limits = plot_amplicons.get_graph_axis_limits(
        max_coverage=0.0,
        max_segment_cn=4.0,
        cn_sum_squares=16.0,
        cn_coverage_cross_product=0.0,
        max_arc_apex=8.0,
    )

    assert limits.cn_ymax > 8.0
    assert limits.coverage_ymax > 1.0
    assert math.isclose(
        limits.coverage_ymax / limits.expansion_factor,
        1.0,
    )


def test_graph_axis_limits_prioritize_arc_fit_over_coverage_cutoff() -> None:
    limits = plot_amplicons.get_graph_axis_limits(
        max_coverage=40.0,
        max_segment_cn=4.0,
        cn_sum_squares=20.0,
        cn_coverage_cross_product=200.0,
        max_arc_apex=8.0,
        max_coverage_cutoff=30.0,
    )

    assert math.isclose(limits.cn_ymax, 8.4)
    assert math.isclose(limits.coverage_ymax, 84.0)
    assert limits.coverage_ymax > 30.0
    assert math.isclose(
        limits.coverage_ymax / limits.cn_ymax,
        10.0,
    )


def test_font_size_multiplier_scales_base_sizes() -> None:
    assert plot_amplicons.get_gene_font_size(2.0) == 24.0
    assert plot_amplicons.get_gene_font_size(0.5) == 6.0
    assert plot_amplicons.get_gene_font_size(0.0) == 0.0
    assert plot_amplicons.get_gene_font_size(2.0, 10.0) == 20.0


def test_font_size_multiplier_rejects_invalid_values() -> None:
    for invalid_value in (-0.1, math.inf, -math.inf, math.nan, 1e308):
        try:
            plot_amplicons.get_gene_font_size(invalid_value)
        except ValueError:
            pass
        else:
            raise AssertionError(
                f"accepted invalid multiplier: {invalid_value}"
            )


def test_font_size_multiplier_scales_axis_marks_and_spines() -> None:
    fig, ax = plot_amplicons.plt.subplots()
    try:
        plot_amplicons.scale_axis_elements(ax, 2.0)
        x_tick = ax.xaxis.get_major_ticks()[0]
        assert x_tick.tick1line.get_markersize() == (
            plot_amplicons.rcParams["xtick.major.size"] * 2.0
        )
        assert x_tick.tick1line.get_markeredgewidth() == (
            plot_amplicons.rcParams["xtick.major.width"] * 2.0
        )
        assert ax.spines["bottom"].get_linewidth() == (
            plot_amplicons.rcParams["axes.linewidth"] * 2.0
        )
        plot_amplicons.scale_axis_elements(ax, 0.5)
        assert x_tick.tick1line.get_markersize() == (
            plot_amplicons.rcParams["xtick.major.size"] * 0.5
        )
        assert ax.spines["bottom"].get_linewidth() == (
            plot_amplicons.rcParams["axes.linewidth"] * 0.5
        )
    finally:
        plot_amplicons.plt.close(fig)


def test_zero_font_size_multiplier_hides_text_and_axes() -> None:
    fig, ax = plot_amplicons.plt.subplots()
    try:
        ax.set_title("title")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.text(0.5, 0.5, "annotation")
        plot_amplicons.scale_axis_elements(ax, 0.0)
        plot_amplicons.hide_figure_text_if_zero(fig, 0.0)
        assert all(
            not text.get_visible()
            for text in fig.findobj(match=plot_amplicons.Text)
        )
        x_tick = ax.xaxis.get_major_ticks()[0]
        assert x_tick.tick1line.get_markersize() == 0.0
        assert x_tick.tick1line.get_markeredgewidth() == 0.0
        assert all(spine.get_linewidth() == 0.0 for spine in ax.spines.values())
    finally:
        plot_amplicons.plt.close(fig)


def test_save_plot_figure_closes_figure(tmp_path: pathlib.Path) -> None:
    fig = plot_amplicons.plt.figure()
    figure_number = fig.number

    plot_amplicons.save_plot_figure(fig, str(tmp_path / "closed"), dpi=20)

    assert figure_number not in plot_amplicons.plt.get_fignums()
    assert (tmp_path / "closed.png").exists()
    assert (tmp_path / "closed.pdf").exists()


def test_combine_graph_and_cycle_plots_stacks_exact_width(
    tmp_path: pathlib.Path,
) -> None:
    output_prefix = str(tmp_path / "amplicon")
    graph_fig = plot_amplicons.plt.figure(figsize=(2, 1))
    cycle_fig = plot_amplicons.plt.figure(figsize=(2, 2))
    plot_amplicons.save_plot_figure(graph_fig, output_prefix + "_graph", 20)
    plot_amplicons.save_plot_figure(cycle_fig, output_prefix + "_cycles", 20)

    plot_amplicons.combine_graph_and_cycle_plots(output_prefix, 20)

    combined = plot_amplicons.plt.imread(output_prefix + "_combined.png")
    assert combined.shape[:2] == (60, 40)
    assert (tmp_path / "amplicon_combined.pdf").exists()


def test_aligned_combined_plot_crops_graph_annotation_band(
    tmp_path: pathlib.Path,
) -> None:
    output_prefix = str(tmp_path / "amplicon")
    graph_fig = plot_amplicons.plt.figure(figsize=(2, 1))
    cycle_fig = plot_amplicons.plt.figure(figsize=(2, 2))
    plot_amplicons.save_plot_figure(graph_fig, output_prefix + "_graph", 20)
    plot_amplicons.save_plot_figure(cycle_fig, output_prefix + "_cycles", 20)
    axis_layout = plot_amplicons.GenomicAxisLayout(
        total_genomic_length=1.0,
        interval_starts={},
        axis_max=1.0,
        chromosome_guides=(),
        interval_guides=(),
    )

    plot_amplicons.combine_graph_and_cycle_plots(
        output_prefix,
        20,
        axis_layout,
    )

    combined = plot_amplicons.plt.imread(output_prefix + "_combined.png")
    expected_graph_height = (
        int(
            round(20 * (1.0 - plot_amplicons.GRAPH_MAIN_AXIS_BOTTOM_WITH_GENES))
        )
        + 1
    )
    assert combined.shape[:2] == (40 + expected_graph_height, 40)


def test_combined_cycle_axis_keeps_all_graph_intervals() -> None:
    graph_viz = plot_amplicons.GraphViz()
    graph_viz.graph_amplified_intervals = {
        "chr1": [plot_amplicons.datatypes.Interval("chr1", 100, 200)],
        "chr2": [plot_amplicons.datatypes.Interval("chr2", 300, 450)],
    }

    graph_viz.update_cycle_amplified_intervals(
        cycle_ids=[],
        graph_given=True,
        align_to_graph=True,
    )

    assert graph_viz.cycle_amplified_intervals == (
        graph_viz.graph_amplified_intervals
    )
    assert graph_viz.num_amplified_intervals == 2


def test_dense_text_layout_expansion_is_bounded() -> None:
    assert plot_amplicons.get_text_layout_scale(0.0) == 1.0
    assert plot_amplicons.get_text_layout_scale(0.5) == 1.0
    assert plot_amplicons.get_text_layout_scale(1.0) == 1.0
    assert plot_amplicons.get_text_layout_scale(2.0) == 2.0
    assert plot_amplicons.get_text_layout_scale(20.0) == 2.0


def test_gene_label_padding_scales_with_font_and_canvas() -> None:
    compact = plot_amplicons.get_gene_label_padding(1000.0, 9.0, 12.0, 8)
    large_font = plot_amplicons.get_gene_label_padding(1000.0, 18.0, 12.0, 8)
    wide_canvas = plot_amplicons.get_gene_label_padding(1000.0, 9.0, 24.0, 8)

    assert compact > 20.0
    assert large_font > compact
    assert wide_canvas < compact


def test_gene_label_lanes_separate_close_labels_across_intervals() -> None:
    lanes = plot_amplicons.assign_gene_label_lanes(
        [
            (1, "ALDH2", 90.0, 90.2),
            (2, "PTPN11", 90.3, 90.5),
            (3, "BCL7A", 99.0, 99.2),
        ],
        gene_font_size=12.0,
        plot_width=12.0,
        displayed_axis_span=100.0,
    )

    assert lanes[1] != lanes[2]
    assert lanes[3] in {lanes[1], lanes[2]}


def test_gene_label_lanes_expand_instead_of_forcing_overlap() -> None:
    lanes = plot_amplicons.assign_gene_label_lanes(
        [(gene_id, f"GENE{gene_id}", 50.0, 50.1) for gene_id in range(1, 7)],
        gene_font_size=9.0,
        plot_width=12.0,
        displayed_axis_span=100.0,
        axis_width_fraction=0.65,
    )

    assert len(set(lanes.values())) == 6


def test_gene_heights_use_current_coral_lane_positions() -> None:
    genes = [
        SimpleNamespace(gname="A", gstart=0, gend=100, height=0.0),
        SimpleNamespace(gname="B", gstart=10, gend=90, height=0.0),
    ]

    plot_amplicons.GraphViz().set_gene_heights(genes)

    heights = sorted(gene.height for gene in genes)
    assert math.isclose(heights[0], 0.15)
    assert math.isclose(heights[1], 0.75)


def test_graph_plot_does_not_require_bam(
    monkeypatch: object,
    tmp_path: pathlib.Path,
) -> None:
    records = {}

    def fail_open_bam(_self: object, bam_path: object) -> None:
        raise AssertionError(f"unexpected BAM open: {bam_path}")

    def noop_parse_genes(
        _self: object,
        *_args: object,
        **_kwargs: object,
    ) -> None:
        return None

    def record_plot_graph(
        self: plot_amplicons.GraphViz,
        *_args: object,
        **kwargs: object,
    ) -> None:
        records["bam"] = self.bam
        records["sequence_edges"] = len(self.graph.sequence_edges)
        records["discordant_edges"] = len(self.graph.discordant_edges)
        records["fontsize"] = kwargs["fontsize"]
        records["gene_font_size"] = kwargs["gene_font_size"]
        records["font_size_multiplier"] = kwargs["font_size_multiplier"]

    monkeypatch.setattr(plot_amplicons.GraphViz, "open_bam", fail_open_bam)
    monkeypatch.setattr(
        plot_amplicons.GraphViz, "parse_genes", noop_parse_genes
    )
    monkeypatch.setattr(
        plot_amplicons.GraphViz, "plot_graph", record_plot_graph
    )

    graph_path = pathlib.Path("tests/data/amplicon1_graph.txt")
    with graph_path.open() as graph_file:
        plot_amplicons.plot_amplicon(
            core_types.ReferenceGenome.hg38,
            None,
            graph_file,
            None,
            str(tmp_path / "graph_only"),
            None,
            math.inf,
            0.0,
            [],
            12.0,
            None,
            should_plot_graph=True,
            should_plot_cycles=False,
            should_hide_genes=True,
            should_restrict_to_bushman_genes=False,
            should_plot_only_cyclic_walks=False,
            font_size_multiplier=0.5,
        )

    assert records == {
        "bam": None,
        "sequence_edges": 4,
        "discordant_edges": 2,
        "fontsize": 9.0,
        "gene_font_size": 6.0,
        "font_size_multiplier": 0.5,
    }


def test_graph_plot_uses_bam_when_provided(
    monkeypatch: object,
    tmp_path: pathlib.Path,
) -> None:
    records = {}
    bam_path = tmp_path / "reads.bam"

    def record_open_bam(
        _self: object,
        observed_bam_path: pathlib.Path,
    ) -> None:
        records["opened_bam_path"] = observed_bam_path

    def noop_parse_genes(
        _self: object,
        *_args: object,
        **_kwargs: object,
    ) -> None:
        return None

    def record_plot_graph(
        _self: object,
        *_args: object,
        **_kwargs: object,
    ) -> None:
        records["plot_called"] = True

    monkeypatch.setattr(plot_amplicons.GraphViz, "open_bam", record_open_bam)
    monkeypatch.setattr(
        plot_amplicons.GraphViz, "parse_genes", noop_parse_genes
    )
    monkeypatch.setattr(
        plot_amplicons.GraphViz, "plot_graph", record_plot_graph
    )

    graph_path = pathlib.Path("tests/data/amplicon1_graph.txt")
    with graph_path.open() as graph_file:
        plot_amplicons.plot_amplicon(
            core_types.ReferenceGenome.hg38,
            bam_path,
            graph_file,
            None,
            str(tmp_path / "with_bam"),
            None,
            math.inf,
            0.0,
            [],
            12.0,
            None,
            should_plot_graph=True,
            should_plot_cycles=False,
            should_hide_genes=True,
            should_restrict_to_bushman_genes=False,
            should_plot_only_cyclic_walks=False,
        )

    assert records == {
        "opened_bam_path": bam_path,
        "plot_called": True,
    }


def test_cycle_plot_receives_scaled_font_sizes(
    monkeypatch: object,
    tmp_path: pathlib.Path,
) -> None:
    records = {}

    def noop_parse_genes(
        _self: object,
        *_args: object,
        **_kwargs: object,
    ) -> None:
        return None

    def record_plot_cycles(
        _self: object,
        *_args: object,
        **kwargs: object,
    ) -> None:
        records["fontsize"] = kwargs["fontsize"]
        records["gene_font_size"] = kwargs["gene_font_size"]
        records["font_size_multiplier"] = kwargs["font_size_multiplier"]

    monkeypatch.setattr(
        plot_amplicons.GraphViz,
        "parse_genes",
        noop_parse_genes,
    )
    monkeypatch.setattr(
        plot_amplicons.GraphViz,
        "plot_cycles",
        record_plot_cycles,
    )

    cycle_path = pathlib.Path("tests/data/amplicon1_cycles.txt")
    with cycle_path.open() as cycle_file:
        plot_amplicons.plot_amplicon(
            core_types.ReferenceGenome.hg38,
            None,
            None,
            cycle_file,
            str(tmp_path / "cycles_only"),
            None,
            math.inf,
            0.0,
            [],
            12.0,
            None,
            should_plot_graph=False,
            should_plot_cycles=True,
            should_hide_genes=True,
            should_restrict_to_bushman_genes=False,
            should_plot_only_cyclic_walks=False,
            font_size_multiplier=2.0,
        )

    assert records == {
        "fontsize": 36.0,
        "gene_font_size": 24.0,
        "font_size_multiplier": 2.0,
    }


def test_graph_legend_output_prefix_uses_sample_prefix() -> None:
    legend_prefix = plot_amplicons.get_graph_legend_output_prefix("skbr3/out")
    assert legend_prefix == pathlib.Path("skbr3/out_legend")


def test_write_graph_legend_creates_visual_explainer(
    tmp_path: pathlib.Path,
) -> None:
    output_prefix = str(tmp_path / "sample")

    plot_amplicons.write_graph_legend(
        output_prefix,
        "Graph average coverage",
        dpi=72,
    )

    legend_prefix = plot_amplicons.get_graph_legend_output_prefix(output_prefix)
    assert legend_prefix.with_suffix(".png").is_file()
    assert legend_prefix.with_suffix(".pdf").is_file()


def test_graph_legend_distinguishes_sequence_and_sv_edge_roles() -> None:
    handles = plot_amplicons.get_graph_legend_handles("Coverage")
    labels = [handle.get_label() for handle in handles]

    assert "Genomic sequence edge" in labels
    assert "Cycle SV edge" in labels
    assert "Path SV edge" in labels
    cycle_sv_handle = next(
        handle for handle in handles if handle.get_label() == "Cycle SV edge"
    )
    path_sv_handle = next(
        handle for handle in handles if handle.get_label() == "Path SV edge"
    )
    assert cycle_sv_handle.get_color() == plot_amplicons.CYCLE_COLOR
    assert path_sv_handle.get_color() == plot_amplicons.PATH_COLOR
    assert cycle_sv_handle.get_linewidth() == 1.5
    assert path_sv_handle.get_linewidth() == 1.5


def test_graph_legend_is_refreshed_when_scale_changes(
    tmp_path: pathlib.Path,
) -> None:
    output_prefix = str(tmp_path / "legend")
    plot_amplicons.write_graph_legend(
        output_prefix,
        "Graph average coverage",
        font_size_multiplier=0.5,
        dpi=20,
    )
    first_png = (tmp_path / "legend_legend.png").read_bytes()

    plot_amplicons.write_graph_legend(
        output_prefix,
        "Graph average coverage",
        font_size_multiplier=2.0,
        dpi=20,
    )
    second_png = (tmp_path / "legend_legend.png").read_bytes()

    assert first_png != second_png


def test_plot_cli_passes_gene_subset_file(
    monkeypatch: object,
    tmp_path: pathlib.Path,
) -> None:
    from coral import cli

    gene_file = tmp_path / "genes.csv"
    gene_file.write_text("EGFR,MYC\n")
    records = {}

    def record_plot_amplicon(*_args: object, **kwargs: object) -> None:
        records["gene_subset_file"] = kwargs["gene_subset_file"]

    monkeypatch.setattr(
        cli.plot_amplicons, "plot_amplicon", record_plot_amplicon
    )

    result = CliRunner().invoke(
        cli.coral_app,
        [
            "plot",
            "--ref",
            "hg38",
            "--graph",
            "tests/data/amplicon1_graph.txt",
            "--output-prefix",
            str(tmp_path / "plot" / "out"),
            "--gene-subset-file",
            str(gene_file),
        ],
    )

    assert result.exit_code == 0, result.output
    assert records["gene_subset_file"] == gene_file


def test_plot_cli_passes_font_size_multiplier(
    monkeypatch: object,
    tmp_path: pathlib.Path,
) -> None:
    from coral import cli

    records = {}

    def record_plot_amplicon(*_args: object, **kwargs: object) -> None:
        records["font_size_multiplier"] = kwargs["font_size_multiplier"]

    monkeypatch.setattr(
        cli.plot_amplicons,
        "plot_amplicon",
        record_plot_amplicon,
    )

    result = CliRunner().invoke(
        cli.coral_app,
        [
            "plot",
            "--ref",
            "hg38",
            "--graph",
            "tests/data/amplicon1_graph.txt",
            "--output-prefix",
            str(tmp_path / "plot" / "out"),
            "--font-size",
            "0.5",
        ],
    )

    assert result.exit_code == 0, result.output
    assert records["font_size_multiplier"] == 0.5


def test_plot_cli_passes_publication_layout_and_cycle_options(
    monkeypatch: object,
    tmp_path: pathlib.Path,
) -> None:
    from coral import cli

    records = {}
    color_file = tmp_path / "colors.tsv"
    color_file.write_text("c2 #0072B2\n")

    def record_plot_amplicon(*_args: object, **kwargs: object) -> None:
        records.update(kwargs)

    monkeypatch.setattr(
        cli.plot_amplicons,
        "plot_amplicon",
        record_plot_amplicon,
    )

    result = CliRunner().invoke(
        cli.coral_app,
        [
            "plot",
            "--ref",
            "hg38",
            "--graph",
            "tests/data/amplicon1_graph.txt",
            "--cycles",
            "tests/data/amplicon1_cycles.txt",
            "--output-prefix",
            str(tmp_path / "plot" / "out"),
            "--width",
            "8",
            "--aspect-ratio",
            "0.5",
            "--offset",
            "0.08",
            "--min-coord-width",
            "0.2",
            "--coverage-scale",
            "full",
            "--cycle-list",
            "c2",
            "--cycle-color-file",
            str(color_file),
            "--combined",
            "--dpi",
            "600",
        ],
    )

    assert result.exit_code == 0, result.output
    assert records["plot_width"] == 8.0
    assert records["aspect_ratio"] == 0.5
    assert records["interval_offset"] == 0.08
    assert records["min_coord_width"] == 0.2
    assert records["coverage_scale"] == "full"
    assert records["cycle_list"] == "c2"
    assert records["cycle_color_file"] == color_file
    assert records["combine_graph_cycles"] is True
    assert records["dpi"] == 600


def test_plot_cli_reports_plot_option_errors_without_traceback(
    monkeypatch: object,
    tmp_path: pathlib.Path,
) -> None:
    from coral import cli

    def reject_plot(*_args: object, **_kwargs: object) -> None:
        raise ValueError("ID 13 is a path, not a cycle")

    monkeypatch.setattr(cli.plot_amplicons, "plot_amplicon", reject_plot)
    result = CliRunner().invoke(
        cli.coral_app,
        [
            "plot",
            "--ref",
            "hg38",
            "--graph",
            "tests/data/amplicon1_graph.txt",
            "--output-prefix",
            str(tmp_path / "invalid"),
        ],
        standalone_mode=False,
    )

    assert isinstance(result.exception, cli.typer.BadParameter)
    assert "ID 13 is a path, not a cycle" in str(result.exception)


def test_plot_cli_preserves_explicit_gene_fontsize_by_default(
    monkeypatch: object,
    tmp_path: pathlib.Path,
) -> None:
    from coral import cli

    records = {}

    def record_plot_amplicon(*args: object, **_kwargs: object) -> None:
        records["gene_fontsize"] = args[9]

    monkeypatch.setattr(
        cli.plot_amplicons,
        "plot_amplicon",
        record_plot_amplicon,
    )

    result = CliRunner().invoke(
        cli.coral_app,
        [
            "plot",
            "--ref",
            "hg38",
            "--graph",
            "tests/data/amplicon1_graph.txt",
            "--output-prefix",
            str(tmp_path / "plot" / "out"),
            "--gene-fontsize",
            "20",
        ],
    )

    assert result.exit_code == 0, result.output
    assert records["gene_fontsize"] == 20.0
    assert "--gene-fontsize will be ignored" not in result.output


def test_plot_cli_explicit_global_fontsize_overrides_gene_fontsize(
    monkeypatch: object,
    tmp_path: pathlib.Path,
) -> None:
    from coral import cli

    records = {}

    def record_plot_amplicon(*args: object, **kwargs: object) -> None:
        records["gene_fontsize"] = args[9]
        records["font_size_multiplier"] = kwargs["font_size_multiplier"]

    monkeypatch.setattr(
        cli.plot_amplicons,
        "plot_amplicon",
        record_plot_amplicon,
    )

    result = CliRunner().invoke(
        cli.coral_app,
        [
            "plot",
            "--ref",
            "hg38",
            "--graph",
            "tests/data/amplicon1_graph.txt",
            "--output-prefix",
            str(tmp_path / "plot" / "out"),
            "--gene-fontsize",
            "20",
            "--font-size",
            "1",
        ],
    )

    assert result.exit_code == 0, result.output
    assert records["gene_fontsize"] == plot_amplicons.DEFAULT_GENE_FONT_SIZE
    assert records["font_size_multiplier"] == 1.0
    assert "--gene-fontsize will be ignored" in result.output


def test_plot_cli_rejects_invalid_font_size_multipliers() -> None:
    from coral import cli

    for invalid_value in ("-0.1", "nan", "inf", "-inf", "1e308"):
        result = CliRunner().invoke(
            cli.coral_app,
            [
                "plot",
                "--ref",
                "hg38",
                "--graph",
                "tests/data/amplicon1_graph.txt",
                "--output-prefix",
                "unused",
                "--font-size",
                invalid_value,
            ],
            standalone_mode=False,
        )
        assert isinstance(result.exception, cli.typer.BadParameter)


def test_plot_all_cli_passes_gene_subset_file(
    monkeypatch: object,
    tmp_path: pathlib.Path,
) -> None:
    from coral import cli

    gene_file = tmp_path / "genes.csv"
    gene_file.write_text("EGFR,MYC\n")
    records = []

    def record_plot_amplicon(*_args: object, **kwargs: object) -> None:
        records.append(kwargs["gene_subset_file"])

    monkeypatch.setattr(
        cli.plot_amplicons, "plot_amplicon", record_plot_amplicon
    )

    result = CliRunner().invoke(
        cli.coral_app,
        [
            "plot_all",
            "--ref",
            "hg38",
            "--graph-dir",
            "tests/data",
            "--output-prefix",
            str(tmp_path / "plot_all" / "out"),
            "--gene-subset-file",
            str(gene_file),
        ],
    )

    assert result.exit_code == 0, result.output
    assert records
    assert all(record == gene_file for record in records)


def test_plot_all_cli_passes_font_size_multiplier(
    monkeypatch: object,
    tmp_path: pathlib.Path,
) -> None:
    from coral import cli

    records = []

    def record_plot_amplicon(*_args: object, **kwargs: object) -> None:
        records.append(kwargs["font_size_multiplier"])

    monkeypatch.setattr(
        cli.plot_amplicons,
        "plot_amplicon",
        record_plot_amplicon,
    )

    result = CliRunner().invoke(
        cli.coral_app,
        [
            "plot_all",
            "--ref",
            "hg38",
            "--graph-dir",
            "tests/data",
            "--output-prefix",
            str(tmp_path / "plot_all" / "out"),
            "--font-size",
            "2",
        ],
    )

    assert result.exit_code == 0, result.output
    assert records
    assert all(record == 2.0 for record in records)


def test_plot_all_cli_warns_once_when_global_fontsize_overrides_gene_fontsize(
    monkeypatch: object,
    tmp_path: pathlib.Path,
) -> None:
    from coral import cli

    records = []

    def record_plot_amplicon(*_args: object, **kwargs: object) -> None:
        records.append(kwargs["gene_fontsize"])

    monkeypatch.setattr(
        cli.plot_amplicons,
        "plot_amplicon",
        record_plot_amplicon,
    )

    result = CliRunner().invoke(
        cli.coral_app,
        [
            "plot_all",
            "--ref",
            "hg38",
            "--graph-dir",
            "tests/data",
            "--output-prefix",
            str(tmp_path / "plot_all" / "out"),
            "--gene-fontsize",
            "20",
            "--font-size",
            "0.5",
        ],
    )

    assert result.exit_code == 0, result.output
    assert records
    assert all(
        record == plot_amplicons.DEFAULT_GENE_FONT_SIZE for record in records
    )
    assert result.output.count("--gene-fontsize will be ignored") == 1


def test_plot_all_cli_passes_shared_legend_prefix(
    monkeypatch: object,
    tmp_path: pathlib.Path,
) -> None:
    from coral import cli

    output_prefix = str(tmp_path / "plot_all" / "out")
    records = []

    def record_plot_amplicon(*_args: object, **kwargs: object) -> None:
        records.append(kwargs["legend_output_prefix"])

    monkeypatch.setattr(
        cli.plot_amplicons, "plot_amplicon", record_plot_amplicon
    )

    result = CliRunner().invoke(
        cli.coral_app,
        [
            "plot_all",
            "--ref",
            "hg38",
            "--graph-dir",
            "tests/data",
            "--output-prefix",
            output_prefix,
        ],
    )

    assert result.exit_code == 0, result.output
    assert records
    assert all(record == output_prefix for record in records)
