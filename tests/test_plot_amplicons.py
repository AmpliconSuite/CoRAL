from __future__ import annotations

import math
import pathlib

from typer.testing import CliRunner

from coral import cli, core_types, plot_amplicons


def test_parse_gene_subset_file_accepts_text_and_csv(tmp_path: pathlib.Path) -> None:
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


def test_discordant_edge_linewidth_scales_with_cn() -> None:
    assert plot_amplicons.get_discordant_edge_linewidth(1.0, 2.0) == 1.5
    assert plot_amplicons.get_discordant_edge_linewidth(2.0, 2.0) == 3.0
    assert plot_amplicons.get_discordant_edge_linewidth(4.0, 2.0) == 3.0
    assert plot_amplicons.get_discordant_edge_linewidth(4.0, 0.0) == 1.0


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
        **_kwargs: object,
    ) -> None:
        records["bam"] = self.bam
        records["sequence_edges"] = len(self.graph.sequence_edges)
        records["discordant_edges"] = len(self.graph.discordant_edges)

    monkeypatch.setattr(plot_amplicons.GraphViz, "open_bam", fail_open_bam)
    monkeypatch.setattr(plot_amplicons.GraphViz, "parse_genes", noop_parse_genes)
    monkeypatch.setattr(plot_amplicons.GraphViz, "plot_graph", record_plot_graph)

    graph_path = pathlib.Path("sample_data/test4/amplicon1_graph.txt")
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
        )

    assert records == {
        "bam": None,
        "sequence_edges": 4,
        "discordant_edges": 2,
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
    monkeypatch.setattr(plot_amplicons.GraphViz, "parse_genes", noop_parse_genes)
    monkeypatch.setattr(plot_amplicons.GraphViz, "plot_graph", record_plot_graph)

    graph_path = pathlib.Path("sample_data/test4/amplicon1_graph.txt")
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


def test_plot_cli_passes_gene_subset_file(
    monkeypatch: object,
    tmp_path: pathlib.Path,
) -> None:
    gene_file = tmp_path / "genes.csv"
    gene_file.write_text("EGFR,MYC\n")
    records = {}

    def record_plot_amplicon(*_args: object, **kwargs: object) -> None:
        records["gene_subset_file"] = kwargs["gene_subset_file"]

    monkeypatch.setattr(cli.plot_amplicons, "plot_amplicon", record_plot_amplicon)

    result = CliRunner().invoke(
        cli.coral_app,
        [
            "plot",
            "--ref",
            "hg38",
            "--graph",
            "sample_data/test4/amplicon1_graph.txt",
            "--output-prefix",
            str(tmp_path / "plot" / "out"),
            "--gene-subset-file",
            str(gene_file),
        ],
    )

    assert result.exit_code == 0, result.output
    assert records["gene_subset_file"] == gene_file


def test_plot_all_cli_passes_gene_subset_file(
    monkeypatch: object,
    tmp_path: pathlib.Path,
) -> None:
    gene_file = tmp_path / "genes.csv"
    gene_file.write_text("EGFR,MYC\n")
    records = []

    def record_plot_amplicon(*_args: object, **kwargs: object) -> None:
        records.append(kwargs["gene_subset_file"])

    monkeypatch.setattr(cli.plot_amplicons, "plot_amplicon", record_plot_amplicon)

    result = CliRunner().invoke(
        cli.coral_app,
        [
            "plot_all",
            "--ref",
            "hg38",
            "--graph-dir",
            "sample_data/test4",
            "--output-prefix",
            str(tmp_path / "plot_all" / "out"),
            "--gene-subset-file",
            str(gene_file),
        ],
    )

    assert result.exit_code == 0, result.output
    assert records
    assert all(record == gene_file for record in records)
