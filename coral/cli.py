from __future__ import annotations

import logging
import os
import pathlib
import pickle
from typing import Annotated

import pysam
import typer

from coral import (
    bam_types,
    cycle2bed,
    cycle_decomposition,
    datatypes,
    hsr,
    plot_amplicons,
    plot_cn,
)
from coral.breakpoint import infer_breakpoint_graph
from coral.breakpoint.parse_graph import parse_breakpoint_graph
from coral.cnv_seed import run_seeding
from coral.datatypes import Solver

coral_app = typer.Typer(
    help="Long-read amplicon reconstruction pipeline and associated utilities."
)
logger = logging.getLogger(__name__)


def validate_cns_file(cns_file: typer.FileText) -> typer.FileText:
    cns_filename = cns_file.name
    if cns_filename.endswith(".bed") or cns_filename.endswith(".cns"):
        return cns_file
    raise typer.BadParameter(
        "Invalid cn-seg file format! (Only .bed and .cns formats are supported.)"
    )


# Note: typer.Arguments are required, typer.Options are optional
BamArg = Annotated[
    pathlib.Path, typer.Option(help="Sorted indexed (long read) bam file.")
]
CnvSeedArg = Annotated[
    typer.FileText, typer.Option(help="Bed file of CNV seed intervals.")
]
CnSegArg = Annotated[
    typer.FileText,
    typer.Option(
        help="Long read segmented whole genome CN calls (.bed or CNVkit .cns file).",
        callback=validate_cns_file,
    ),
]
OutputDirArg = Annotated[str, typer.Option(help="Directory of output files.")]
OutputPrefixArg = Annotated[str, typer.Option(help="Prefix of output files.")]
OutputPCFlag = Annotated[
    bool,
    typer.Option(
        help="If specified, output all path constraints in *.cycles file."
    ),
]
PostProcessFlag = Annotated[
    bool,
    typer.Option(
        help="Postprocess the cycles/paths returned in greedy cycle extraction."
    ),
]
SolverArg = Annotated[Solver, typer.Option(help="LP solver to use.")]
ThreadsArg = Annotated[
    int,
    typer.Option(
        help="Number of threads reserved for integer program solvers."
    ),
]
TimeLimitArg = Annotated[
    int,
    typer.Option(
        help="Maximum running time (in seconds) reserved for integer program solvers."
    ),
]
AlphaArg = Annotated[
    float,
    typer.Option(
        help="Parameter used to balance CN weight and path constraints in greedy cycle extraction."
    ),
]


@coral_app.command(help="Filter and merge amplified intervals.")
def seed(
    ctx: typer.Context,
    cn_seg: CnSegArg,
    output_prefix: OutputPrefixArg,
    gain: Annotated[
        float,
        typer.Option(
            help="CN gain threshold for interval to be considered as a seed."
        ),
    ] = 6.0,
    min_seed_size: Annotated[
        float,
        typer.Option(
            help="Minimum size (in bp) for interval to be considered as a seed."
        ),
    ] = 100000,
    max_seg_gap: Annotated[
        float,
        typer.Option(
            help="Maximum gap size (in bp) to merge two proximal intervals."
        ),
    ] = 300000,
) -> None:
    print(f"Performing seeding mode with options: {ctx.params}")
    run_seeding(cn_seg, output_prefix, gain, min_seed_size, max_seg_gap)


@coral_app.command(help="Reconstruct focal amplifications")
def reconstruct(
    ctx: typer.Context,
    output_dir: OutputDirArg,
    lr_bam: BamArg,
    cnv_seed: CnvSeedArg,
    cn_seg: CnSegArg,
    log_file: Annotated[str, typer.Option(help="Name of log file.")] = "",
    cycle_decomp_alpha: AlphaArg = 0.01,
    cycle_decomp_time_limit: TimeLimitArg = 7200,
    cycle_decomp_threads: ThreadsArg = -1,
    solver: SolverArg = Solver.GUROBI,
    output_all_path_constraints: OutputPCFlag = False,
    postprocess_greedy_sol: PostProcessFlag = False,
    output_bp: Annotated[
        bool,
        typer.Option(help="If specified, only output the list of breakpoints."),
    ] = False,
    skip_cycle_decomp: Annotated[
        bool,
        typer.Option(
            help="If specified, only reconstruct and output the breakpoint graph for all amplicons.",
        ),
    ] = False,
    min_bp_support: Annotated[
        float,
        typer.Option(
            help="Ignore breakpoints with less than (min_bp_support * normal coverage) long read support."
        ),
    ] = 1.0,
) -> None:
    print(f"Performing reconstruction with options: {ctx.params}")

    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=f"{output_dir}/infer_breakpoint_graph.log",
        filemode="w+",
        level=logging.DEBUG,
        format="%(asctime)s:%(levelname)-4s [%(filename)s:%(lineno)d] %(message)s",
    )
    logging.getLogger("pyomo").setLevel(logging.INFO)
    b2bn = infer_breakpoint_graph.reconstruct_graphs(
        lr_bam, cnv_seed, cn_seg, output_dir, output_bp, min_bp_support
    )
    solver_options = datatypes.SolverOptions(
        num_threads=cycle_decomp_threads,
        time_limit_s=cycle_decomp_time_limit,
        output_dir=output_dir,
        model_prefix="pyomo",
        solver=solver,
    )
    if not (output_bp or skip_cycle_decomp):
        cycle_decomposition.reconstruct_cycles(
            b2bn,
            solver_options,
            cycle_decomp_alpha,
            should_postprocess_greedy_sol=postprocess_greedy_sol,
            output_all_path_constraints=output_all_path_constraints,
        )
    b2bn.closebam()
    print("\nCompleted reconstruction.")


@coral_app.command(
    name="cycle", help="Pass existing breakpoint file directly to LP solver."
)
def cycle_decomposition_mode(
    bp_graph: Annotated[
        typer.FileText, typer.Option(help="Existing BP graph file.")
    ],
    output_dir: OutputDirArg,
    alpha: AlphaArg = 0.01,
    time_limit_s: TimeLimitArg = 7200,
    threads: ThreadsArg = -1,
    solver: SolverArg = Solver.GUROBI,
    output_all_path_constraints: OutputPCFlag = False,
    postprocess_greedy_sol: PostProcessFlag = False,
) -> None:
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        filename=f"{output_dir}/infer_breakpoint_graph.log",
        filemode="w+",
        level=logging.DEBUG,
        format="%(asctime)s:%(levelname)-4s [%(filename)s:%(lineno)d] %(message)s",
    )
    logging.getLogger("pyomo").setLevel(logging.ERROR)

    parsed_bp_graph = parse_breakpoint_graph(bp_graph)
    amplicon_idx = int(bp_graph.name.split("_")[0].split("amplicon")[1])
    parsed_bp_graph.amplicon_idx = amplicon_idx

    solver_options = datatypes.SolverOptions(
        num_threads=threads,
        time_limit_s=time_limit_s,
        output_dir=output_dir,
        model_prefix="pyomo",
        solver=solver,
    )
    cycle_decomposition.cycle_decomposition_single_graph(
        parsed_bp_graph,
        solver_options,
        alpha,
        should_postprocess=postprocess_greedy_sol,
        output_all_path_constraints=output_all_path_constraints,
    )
    # cycle_decomposition.reconstruct_cycles(
    #     output_dir,
    #     output_all_path_constraints,
    #     alpha,
    #     time_limit_s,
    #     threads,
    #     solver,
    #     should_postprocess_greedy_sol=postprocess_greedy_sol,
    #     bb=bb,
    # )


@coral_app.command(
    name="hsr",
    help="Detect possible integration points of ecDNA HSR amplifications.",
)
def hsr_mode(
    ctx: typer.Context,
    lr_bam: BamArg,
    cn_seg: CnSegArg,
    output_prefix: OutputPrefixArg,
    cycles: Annotated[
        typer.FileText,
        typer.Option(help="AmpliconSuite-formatted cycles file."),
    ],
    normal_cov: Annotated[
        float, typer.Option(help="Estimated diploid coverage.")
    ],
    bp_match_cutoff: Annotated[
        int, typer.Option(help="Breakpoint matching cutoff.")
    ] = 100,
    bp_match_cutoff_clustering: Annotated[
        int,
        typer.Option(help="Crude breakpoint matching cutoff for clustering."),
    ] = 2000,
) -> None:
    print(f"Performing HSR mode with options: {ctx.params}")
    hsr.locate_hsrs(
        lr_bam,
        cycles,
        cn_seg,
        output_prefix,
        normal_cov,
        bp_match_cutoff,
        bp_match_cutoff_clustering,
    )


@coral_app.command(
    name="plot",
    help="Generate plots of amplicon cycles and/or graph from AA-formatted output files",
)
def plot_mode(
    ctx: typer.Context,
    ref: Annotated[str, typer.Option(help="Reference genome.")],
    graph: Annotated[
        typer.FileText | None,
        typer.Option(help="AmpliconSuite-formatted graph file (*_graph.txt)."),
    ],
    bam: BamArg,
    output_prefix: OutputPrefixArg,
    cycle_file: Annotated[
        typer.FileText | None,
        typer.Option(
            help="AmpliconSuite-formatted cycles file (*_cycles.txt)."
        ),
    ] = None,
    num_cycles: Annotated[
        int | None, typer.Option(help="Only plot the first NUM_CYCLES cycles.")
    ] = None,
    region: Annotated[
        str | None,
        typer.Option(
            help="Specifically visualize only this region, argument formatted as 'chr1:pos1-pos2'."
        ),
    ] = None,
    plot_graph: Annotated[
        bool, typer.Option(help="Visualize breakpoint graph.")
    ] = True,
    plot_cycles: Annotated[
        bool, typer.Option(help="Visualize (selected) cycles.")
    ] = False,
    only_cyclic_paths: Annotated[
        bool, typer.Option(help="Only plot cyclic paths from cycles file.")
    ] = False,
    max_coverage: Annotated[
        float,
        typer.Option(
            help="Limit the maximum visualized coverage in the graph."
        ),
    ] = float("inf"),
    min_mapq: Annotated[
        float,
        typer.Option(
            help="Minimum mapping quality to count read in coverage plotting."
        ),
    ] = 0.0,
    gene_subset_list: Annotated[
        list[str],
        typer.Option(
            help="List of genes to visualize (will show all by default)."
        ),
    ] = [],
    hide_genes: Annotated[
        bool, typer.Option(help="Do not show gene track.")
    ] = False,
    gene_fontsize: Annotated[
        float, typer.Option(help="Change size of gene font.")
    ] = 12.0,
    bushman_genes: Annotated[
        bool,
        typer.Option(
            help="Reduce gene set to the Bushman cancer-related gene set."
        ),
    ] = False,
) -> None:
    print(f"Performing plot mode with options: {ctx.params}")
    if "/" in output_prefix:
        os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
    plot_amplicons.plot_amplicons(
        ref,
        bam,
        graph,
        cycle_file,
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


@coral_app.command(
    name="cycle2bed", help="Convert cycle files in AA format to bed format."
)
def cycle2bed_mode(
    ctx: typer.Context,
    cycle_file: Annotated[
        typer.FileText,
        typer.Option(
            help="AmpliconSuite-formatted cycles file (*_cycles.txt)."
        ),
    ],
    output_file: Annotated[str, typer.Option(help="Output file name.")],
    num_cycles: Annotated[
        int | None, typer.Option(help="Only plot the first NUM_CYCLES cycles.")
    ] = None,
    rotate_to_min: Annotated[
        bool,
        typer.Option(
            help="Output cycles starting from the canonically smallest segment with positive strand."
        ),
    ] = False,
) -> None:
    print(f"Performing cycle to bed mode with options: {ctx.params}")
    cycle2bed.convert_cycles_to_bed(
        cycle_file, output_file, rotate_to_min, num_cycles
    )


@coral_app.command(
    name="plot_cn",
    help="Generate CN plots from .cnr (Copy Number Ratio) files generated by CNVkit.",
)
def plot_cn_mode(
    ctx: typer.Context,
    cnr: Annotated[
        typer.FileText,
        typer.Option(help="CNVkit-generated .cnr (Copy Number Ratio) file."),
    ],
    output_dir: OutputDirArg,
    name: Annotated[str, typer.Option(help="Name of sample.")],
) -> None:
    print(f"Performing plot mode with options: {ctx.params}")
    plot_cn.plot_cnr(cnr, output_dir, name)
