#!/usr/bin/env python3
from __future__ import annotations

import enum
import logging
import pickle

import click

from coral import cycle2bed, cycle_decomposition, hsr, plot_amplicons
from coral.breakpoint import infer_breakpoint_graph
from coral.cnv_seed import run_seeding
from coral.constants import CNGAP_MAX, CNSIZE_MIN, GAIN


class Solver(enum.Enum):
    GUROBI = "gurobi"
    SCIP = "scip"
    # Coin-OR solvers
    BONMIN = "bonmin"
    COUENNE = "couenne"


@click.group(help="Long-read amplicon reconstruction pipeline and associated utilities.")
def cli_mode() -> None:
    pass


def validate_cns_file(ctx: click.Context, param_name: str, cns_file: click.File) -> click.File:
    cns_filename = cns_file.name
    if cns_filename.endswith(".bed") or cns_filename.endswith(".cns"):
        return cns_file
    raise click.BadParameter("Invalid cn-seg file format! (Only .bed and .cns formats are supported.)")


@cli_mode.command(help="Filter and merge amplified intervals.")
@click.option(
    "--cn-seg",
    type=str,
    required=True,
    help="Long read segmented whole genome CN calls (.bed or CNVkit .cns file).",
)
@click.option(
    "--out",
    type=str,
    default="<INPUT_CNS_BASENAME>_CNV_SEEDS.bed",
    help="OPTIONAL: Prefix filename for output bed file. Default: <INPUT_CNS_BASENAME>_CNV_SEEDS.bed",
)
@click.option(
    "--gain",
    type=float,
    default=GAIN,
    help="OPTIONAL: CN gain threshold for interval to be considered as a seed. Default: 6.0",
)
@click.option(
    "--min-seed-size",
    type=float,
    default=CNSIZE_MIN,
    help="OPTIONAL: Minimum size (in bp) for interval to be considered as a seed. Default: 100000",
)
@click.option(
    "--max-seg-gap",
    type=float,
    default=CNGAP_MAX,
    help="OPTIONAL: Maximum gap size (in bp) to merge two proximal intervals. Default: 300000",
)
@click.pass_context  # TODO: do we need to actually log arguments or can we remove this for housekeeping?
def seed(ctx: click.Context, cn_seg: str, out: str, gain: float, min_seed_size: float, max_seg_gap: float) -> None:
    print(f"Performing seeding mode with options: {ctx.params}")
    run_seeding(cn_seg, out, gain, min_seed_size, max_seg_gap)


@cli_mode.command(help="Reconstruct focal amplifications")
# TODO; unify shared params with hsr via click subgroups
@click.option("--lr-bam", type=str, required=True, help="Sorted indexed (long read) bam file.")
@click.option("--cnv-seed", type=click.File("r"), required=True, help="Bed file of CNV seed intervals.")
@click.option(
    "--cn-seg",
    type=click.File("r"),
    required=True,
    callback=validate_cns_file,
    help="Long read segmented whole genome CN calls (.bed or CNVkit .cns file).",
)
@click.option("--output-prefix", type=str, required=True, help="Prefix of output files.")
@click.option("--output-bp", is_flag=True, default=False, help="If specified, only output the list of breakpoints.")
@click.option(
    "--skip-cycle-decomp",
    is_flag=True,
    default=False,
    help="If specified, only reconstruct and output the breakpoint graph for all amplicons.",
)
@click.option(
    "--output_all_path_constraints",
    is_flag=True,
    default=False,
    help="If specified, output all path constraints in *.cycles file.",
)
@click.option(
    "--min-bp-support",
    type=float,
    default=1.0,
    help="Ignore breakpoints with less than (min_bp_support * normal coverage) long read support.",
)
@click.option(
    "--cycle-decomp-alpha",
    type=float,
    default=0.01,
    help="Parameter used to balance CN weight and path constraints in greedy cycle extraction.",
)
@click.option(
    "--cycle-decomp-time-limit",
    type=int,
    default=7200,
    help="Maximum running time (in seconds) reserved for integer program solvers.",
)
@click.option("--cycle-decomp-threads", type=int, help="Number of threads reserved for integer program solvers.")
@click.option(
    "--postprocess-greedy-sol",
    is_flag=True,
    default=False,
    help="Postprocess the cycles/paths returned in greedy cycle extraction.",
)
@click.option("--log-file", type=str, default="", help="Name of log file.")
@click.option("--solver", type=click.Choice([solver.value for solver in Solver]), default=Solver.GUROBI.value)
@click.pass_context
def reconstruct(
    ctx: click.Context,
    lr_bam: str,
    cnv_seed: click.File,
    cn_seg: click.File,
    output_prefix: str,
    output_bp: bool,
    skip_cycle_decomp: bool,
    output_all_path_constraints: bool,
    min_bp_support: float,
    cycle_decomp_alpha: float,
    cycle_decomp_time_limit: int,
    cycle_decomp_threads: int,
    postprocess_greedy_sol: bool,
    log_file: str,
    solver: str,
) -> None:
    print(f"Performing reconstruction with options: {ctx.params}")
    logging.basicConfig(
        filename=f"{output_prefix}/infer_breakpoint_graph.log" or "infer_breakpoint_graph.log",
        filemode="w",
        level=logging.DEBUG,
        format="%(asctime)s:%(levelname)-4s [%(filename)s:%(lineno)d] %(message)s",
    )
    b2bn = infer_breakpoint_graph.reconstruct_graph(
        lr_bam,
        cnv_seed,
        cn_seg,
        output_prefix,
        output_bp,
        skip_cycle_decomp,
        output_all_path_constraints,
        min_bp_support,
        cycle_decomp_alpha,
        cycle_decomp_time_limit,
        cycle_decomp_threads,
        postprocess_greedy_sol,
        log_file,
    )
    if not (output_bp or skip_cycle_decomp):
        cycle_decomposition.reconstruct_cycles(
            output_prefix,
            output_all_path_constraints,
            cycle_decomp_alpha,
            cycle_decomp_time_limit,
            cycle_decomp_threads,
            solver,
            postprocess_greedy_sol,
            b2bn,
        )
    b2bn.closebam()
    print("\nCompleted reconstruction.")


@cli_mode.command(name="decompose", help="Pass existing breakpoint file directly to LP solver.")
@click.option("--output-prefix", type=str, required=True, help="Prefix of output files.")
@click.option(
    "--output_all_path_constraints",
    is_flag=True,
    default=False,
    help="If specified, output all path constraints in *.cycles file.",
)
@click.option(
    "--cycle-decomp-alpha",
    type=float,
    default=0.01,
    help="Parameter used to balance CN weight and path constraints in greedy cycle extraction.",
)
@click.option(
    "--cycle-decomp-time-limit",
    type=int,
    default=7200,
    help="Maximum running time (in seconds) reserved for integer program solvers.",
)
@click.option("--cycle-decomp-threads", type=int, help="Number of threads reserved for integer program solvers.")
@click.option("--bp-graph", type=click.File("rb"), help="Existing BP graph file.")
@click.option(
    "--postprocess-greedy-sol",
    is_flag=True,
    default=False,
    help="Postprocess the cycles/paths returned in greedy cycle extraction.",
)
@click.option("--solver", type=click.Choice([solver.value for solver in Solver]), default=Solver.GUROBI.value)
def cycle_decomposition_mode(
    output_prefix: str,
    output_all_path_constraints: bool,
    cycle_decomp_alpha: float,
    cycle_decomp_time_limit: int,
    cycle_decomp_threads: int,
    bp_graph: click.File,
    postprocess_greedy_sol: bool,
    solver: str,
):
    bb = infer_breakpoint_graph.BamToBreakpointNanopore(None, [pickle.load(bp_graph)])
    cycle_decomposition.reconstruct_cycles(
        output_prefix,
        output_all_path_constraints,
        cycle_decomp_alpha,
        cycle_decomp_time_limit,
        cycle_decomp_threads,
        solver,
        postprocess_greedy_sol=False,
        bb=bb,
    )


@cli_mode.command(name="hsr", help="Detect possible integration points of ecDNA HSR amplifications.")
@click.option("--lr-bam", type=str, required=True, help="Sorted indexed (long read) bam file.")
@click.option(
    "--cycles", type=str, required=True, help="AmpliconSuite-formatted cycles file"
)  # TODO: update type to file
@click.option(
    "--cn-seg",
    type=str,
    required=True,
    help="Long read segmented whole genome CN calls (.bed or CNVkit .cns file).",
)
@click.option("--output-prefix", type=str, required=True, help="Prefix of output files.")
@click.option("--normal-cov", type=float, required=True, help="Estimated diploid coverage.")
@click.option("--bp-match-cutoff", type=int, default=100, help="Breakpoint matching cutoff.")
@click.option(
    "--bp-match-cutoff-clustering", type=int, default=2000, help="Crude breakpoint matching cutoff for clustering."
)
@click.pass_context
def hsr_mode(
    ctx: click.Context,
    lr_bam: str,
    cycles: str,
    cn_seg: str,
    output_prefix: str,
    normal_cov: float,
    bp_match_cutoff: int,
    bp_match_cutoff_clustering: int,
) -> None:
    print(f"Performing HSR mode with options: {ctx.params}")
    hsr.locate_hsrs(lr_bam, cycles, cn_seg, output_prefix, normal_cov, bp_match_cutoff, bp_match_cutoff_clustering)


@cli_mode.command(name="plot", help="Generate plots of amplicon cycles and/or graph from AA-formatted output files")
@click.option("--ref", type=click.Choice(["hg19", "hg38", "GRCh38", "mm10", "GRCh37"]), required=True)
@click.option("--bam", type=str, help="Sorted indexed (long read) bam file.")
@click.option("--graph", type=str, help="AmpliconSuite-formatted *.graph file")  # TODO: update type to file
@click.option("--cycles", type=str, help="AmpliconSuite-formatted cycles file")  # TODO: update type to file
@click.option("--output-prefix", "-o", type=str, required=True, help="Prefix of output files.")
@click.option("--plot-graph", is_flag=True, default=False, help="Visualize breakpoint graph.")
@click.option("--plot-cycles", is_flag=True, default=False, help="Visualize (selected) cycles.")
@click.option("--only-cyclic-paths", is_flag=True, default=False, help="Only plot cyclic paths from cycles file.")
@click.option("--num-cycles", type=int, help="Only plot the first NUM_CYCLES cycles.")
@click.option(
    "--max-coverage", type=float, default=float("inf"), help="Limit the maximum visualized coverage in the graph."
)
@click.option("--min-mapq", type=float, default=0, help="Minimum mapping quality to count read in coverage plotting.")
@click.option(
    "--gene-subset-list",
    type=list[str],
    default=[],
    help="List of genes to visualize (will show all by default).",
)
@click.option("--hide-genes", is_flag=True, default=False, help="Do not show gene track.")
@click.option("--gene-fontsize", type=float, default=12.0, help="Change size of gene font.")
@click.option(
    "--bushman-genes", is_flag=True, default=False, help="Reduce gene set to the Bushman cancer-related gene set"
)
@click.option(
    "--region",
    type=str,
    help="(Graph plotting) Specifically visualize only this region, argument formatted as 'chr1:pos1-pos2'.",
)
@click.pass_context
def plot_mode(
    ctx: click.Context,
    ref: str,
    bam: str,
    graph: str,
    cycles: str,
    output_prefix: str,
    plot_graph: bool,
    plot_cycles: bool,
    only_cyclic_paths: bool,
    num_cycles: bool,
    max_coverage: float,
    min_mapq: float,
    gene_subset_list: list[str],
    hide_genes: bool,
    gene_fontsize: float,
    bushman_genes: bool,
    region: str,
) -> None:
    print(f"Performing plot mode with options: {ctx.params}")
    plot_amplicons.plot_amplicons(
        ref,
        bam,
        graph,
        cycles,
        output_prefix,
        plot_graph,
        plot_cycles,
        only_cyclic_paths,
        num_cycles,
        max_coverage,
        min_mapq,
        gene_subset_list,
        hide_genes,
        gene_fontsize,
        bushman_genes,
        region,
    )


@cli_mode.command(name="cycle2bed", help="Convert cycle files in AA format to bed format.")
@click.option("--cycle-file", type=str, required=True, help="Input AA-formatted cycle file.")
@click.option("--output-file", type=str, required=True, help="Output file name.")
@click.option("--num-cycles", type=int, help="Only plot the first NUM_CYCLES cycles.")
@click.option(
    "--rotate-to-min",
    is_flag=True,
    default=False,
    help="Output cycles starting from the canonically smallest segment with positive strand.",
)
@click.pass_context
def cycle2bed_mode(ctx: click.Context, cycle_file: str, output_file: str, num_cycles: int, rotate_to_min: bool) -> None:
    print(f"Performing cycle to bed mode with options: {ctx.params}")
    if rotate_to_min:
        cycle2bed.convert_cycles_to_bed(cycle_file, output_file, True, num_cycles)
    else:
        cycle2bed.convert_cycles_to_bed(cycle_file, output_file, False, num_cycles)
