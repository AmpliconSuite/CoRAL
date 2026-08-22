from __future__ import annotations

import math
import pathlib
from types import SimpleNamespace

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
    assert plot_amplicons.get_discordant_edge_linewidth(4.0, 0.0) == 1.0


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
        max_segment_cn * (1.25 + 0.75 * 0.01),
    )
    assert short_arc < long_arc
    assert math.isclose(
        long_arc,
        max_segment_cn * (1.25 + 0.75 * 0.9),
    )
    assert math.isclose(full_width_arc, 2.0 * max_segment_cn)


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
            raise AssertionError(f"accepted invalid multiplier: {invalid_value}")


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


def test_dense_text_layout_expansion_is_bounded() -> None:
    assert plot_amplicons.get_text_layout_scale(0.0) == 1.0
    assert plot_amplicons.get_text_layout_scale(0.5) == 1.0
    assert plot_amplicons.get_text_layout_scale(1.0) == 1.0
    assert plot_amplicons.get_text_layout_scale(2.0) == 2.0
    assert plot_amplicons.get_text_layout_scale(20.0) == 2.0


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
        record == plot_amplicons.DEFAULT_GENE_FONT_SIZE
        for record in records
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
