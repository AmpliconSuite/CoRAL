from __future__ import annotations

import logging
import os
import pathlib
from typing import Annotated

import colorama
import typer

from coral import (
    core_types,
    cycle2bed,
    cycle_decomposition,
    datatypes,
    global_state,
    hsr,
    output,
    plot_amplicons,
    plot_cn,
    summary,
)
from coral.breakpoint import infer_breakpoint_graph
from coral.breakpoint.parse_graph import (
    get_all_graphs_from_dir,
    parse_breakpoint_graph,
)
from coral.cnv_seed import run_seeding
from coral.core_utils import (
    get_reconstruction_paths_from_shared_dir,
)
from coral.datatypes import CycleDecompOptions, OutputPCOptions, Solver
from coral.output import cycle_output
from coral.scoring import score_simulation

colorama.init()
coral_app = typer.Typer(
    help="Long-read amplicon reconstruction pipeline and associated utilities.",
    pretty_exceptions_show_locals=False,  # Prints all local variables in the
    # error traceback, which is typically kind of insane with WGS data
)
logger = logging.getLogger(__name__)


def validate_cns_file(cns_file: typer.FileText) -> typer.FileText:
    cns_filename = cns_file.name
    if cns_filename.endswith((".bed", ".cns")):
        return cns_file
    raise typer.BadParameter(
        "Invalid cn-seg file format! (Only .bed and .cns formats are supported.)"
    )


# Note: typer.Arguments are required, typer.Options are optional
BamArg = Annotated[
    pathlib.Path | None,
    typer.Option(help="Sorted indexed (long read) bam file."),
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
OutputPrefixArg = Annotated[str, typer.Option(help="Prefix of output files.")]
ReconstructLogArg = Annotated[
    pathlib.Path | None,
    typer.Option(
        help="Name of the main *.log file, which can be used to trace the status "
        "of reconstruct run(s)."
    ),
]
OutputPCArg = Annotated[
    OutputPCOptions,
    typer.Option(
        help="Options for outputting path constraints in *.graph and *.cycles files."
    ),
]
PostProcessFlag = Annotated[
    bool,
    typer.Option(
        help="Postprocess the cycles/paths returned in greedy cycle extraction."
    ),
]
CycleDecompArg = Annotated[
    CycleDecompOptions,
    typer.Option(
        help="Cycle decomposition strategies, either directly perform greedy cycle extraction, "
        "or first try cycle minimization and if not feasible, switch to greedy."
    ),
]
IgnorePathConstraintsFlag = Annotated[
    bool,
    typer.Option(
        help="If specified, ignore path constraints in cycle decomposition."
    ),
]
SolverArg = Annotated[Solver, typer.Option(help="LP solver to use.")]
ThreadsArg = Annotated[
    int,
    typer.Option(
        help="Number of threads reserved for integer program solvers."
    ),
]
SolverTimeLimitArg = Annotated[
    int,
    typer.Option(
        help="Maximum running time (in seconds) reserved for integer program solvers."
    ),
]
GlobalTimeLimitArg = Annotated[
    int,
    typer.Option(
        help="Maximum running time (in seconds) reserved for full sample analysis."
    ),
]
AlphaArg = Annotated[
    float,
    typer.Option(
        help="Parameter used to balance CN weight and path constraints in greedy cycle extraction."
    ),
]
ReferenceGenomeArg = Annotated[
    core_types.ReferenceGenome,
    typer.Option(help="Reference genome."),
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
    print(
        f"{colorama.Style.DIM}{colorama.Fore.LIGHTYELLOW_EX}"
        f"Performing seeding mode with options: {ctx.params}"
        f"{colorama.Style.RESET_ALL}"
    )
    if "/" in output_prefix:
        os.makedirs(os.path.dirname(output_prefix), exist_ok=True)
    run_seeding(cn_seg, output_prefix, gain, min_seed_size, max_seg_gap)


@coral_app.command(help="Reconstruct focal amplifications")
def reconstruct(
    ctx: typer.Context,
    output_prefix: OutputPrefixArg,
    lr_bam: BamArg,
    cnv_seed: CnvSeedArg,
    cn_seg: CnSegArg,
    global_time_limit: GlobalTimeLimitArg = 21600,  # 6 hrs in seconds
    cycle_decomp_mode: CycleDecompArg = CycleDecompOptions.MAX_WEIGHT,
    cycle_decomp_alpha: AlphaArg = 0.01,
    solver: SolverArg = Solver.GUROBI,
    solver_time_limit: SolverTimeLimitArg = 7200,  # 2 hrs in seconds
    solver_threads: ThreadsArg = -1,
    output_path_constraints: OutputPCArg = OutputPCOptions.LONGEST,
    postprocess_greedy_sol: PostProcessFlag = False,
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
    ignore_path_constraints: IgnorePathConstraintsFlag = False,
    profile: Annotated[
        bool, typer.Option(help="Profile resource usage.")
    ] = False,
    log_file: ReconstructLogArg = None,
) -> None:
    print(
        f"{colorama.Style.DIM}{colorama.Fore.LIGHTYELLOW_EX}"
        f"Performing reconstruction with options: {ctx.params}"
        f"{colorama.Style.RESET_ALL}"
    )

    solver_output_dir = pathlib.Path.cwd()
    solver_output_prefix = output_prefix
    if output_prefix.rfind("/") > 0:
        solver_output_dir = output_prefix[: output_prefix.rfind("/")]
        solver_output_prefix = output_prefix[output_prefix.rfind("/") + 1 :]
    pathlib.Path(f"{solver_output_dir}/models").mkdir(
        parents=True, exist_ok=True
    )
    # validate_cycle_flags(ctx)

    log_fn = f"{output_prefix}_reconstruct.log"
    if log_file:
        log_fn = log_file
    logging.basicConfig(
        filename=log_fn,
        filemode="w+",
        level=logging.DEBUG,
        format="%(asctime)s:%(levelname)-4s [%(filename)s:%(lineno)d] %(message)s",
    )
    logging.getLogger("pyomo").setLevel(logging.INFO)
    global_state.STATE_PROVIDER.should_profile = profile
    global_state.STATE_PROVIDER.output_prefix = output_prefix
    global_state.STATE_PROVIDER.time_limit_s = global_time_limit

    b2bn = infer_breakpoint_graph.reconstruct_graphs(
        lr_bam,
        cnv_seed,
        cn_seg,
        output_prefix,
        output_path_constraints,
        min_bp_support,
    )
    solver_options = datatypes.SolverOptions(
        num_threads=solver_threads,
        time_limit_s=solver_time_limit,
        output_dir=solver_output_dir,
        output_prefix=solver_output_prefix,
        model_prefix="pyomo",
        solver=solver,
    )
    if not skip_cycle_decomp:
        cycle_decomposition.reconstruct_cycles(
            b2bn.lr_graph,
            solver_options,
            cycle_decomp_mode,
            cycle_decomp_alpha,
            should_postprocess_greedy_sol=postprocess_greedy_sol,
            pc_output_option=output_path_constraints,
            ignore_path_constraints=ignore_path_constraints,
        )

    b2bn.closebam()
    if profile:
        summary.output.add_resource_usage_summary(solver_options)
    print("\nCompleted reconstruction.")


@coral_app.command(
    name="cycle", help="Pass existing breakpoint file directly to LP solver."
)
def cycle_decomposition_mode(
    ctx: typer.Context,
    bp_graph: Annotated[
        typer.FileText, typer.Option(help="Existing BP graph file.")
    ],
    output_prefix: OutputPrefixArg,
    alpha: AlphaArg = 0.01,
    solver_time_limit: SolverTimeLimitArg = 7200,
    threads: ThreadsArg = -1,
    solver: SolverArg = Solver.GUROBI,
    global_time_limit: GlobalTimeLimitArg = 21600,
    cycle_decomp_mode: CycleDecompArg = CycleDecompOptions.MAX_WEIGHT,
    output_path_constraints: OutputPCArg = OutputPCOptions.LONGEST,
    postprocess_greedy_sol: PostProcessFlag = False,
    ignore_path_constraints: IgnorePathConstraintsFlag = False,
    profile: Annotated[
        bool, typer.Option(help="Profile resource usage.")
    ] = False,
    log_file: ReconstructLogArg = None,
) -> None:
    print(
        f"{colorama.Style.DIM}{colorama.Fore.LIGHTYELLOW_EX}"
        f"Performing reconstruction with options: {ctx.params}"
        f"{colorama.Style.RESET_ALL}"
    )

    solver_output_dir = pathlib.Path.cwd()
    solver_output_prefix = output_prefix
    if output_prefix.rfind("/") > 0:
        solver_output_dir = pathlib.Path(
            output_prefix[: output_prefix.rfind("/")]
        )
        solver_output_prefix = output_prefix[output_prefix.rfind("/") + 1 :]
    pathlib.Path(f"{solver_output_dir}/models").mkdir(
        parents=True, exist_ok=True
    )

    log_fn = f"{output_prefix}_cycle_decomposition.log"
    if log_file:
        log_fn = log_file
    logging.basicConfig(
        filename=log_fn,
        filemode="w+",
        level=logging.DEBUG,
        format="%(asctime)s:%(levelname)-4s [%(filename)s:%(lineno)d] %(message)s",
    )
    logging.getLogger("pyomo").setLevel(logging.ERROR)
    global_state.STATE_PROVIDER.should_profile = profile
    global_state.STATE_PROVIDER.output_prefix = output_prefix
    global_state.STATE_PROVIDER.time_limit_s = global_time_limit

    parsed_bp_graph = parse_breakpoint_graph(bp_graph)
    amplicon_idx = int(bp_graph.name.split("_")[-2].split("amplicon")[1])
    parsed_bp_graph.amplicon_idx = amplicon_idx - 1
    parsed_bp_graph.infer_amplicon_intervals(amplicon_idx - 1)

    solver_options = datatypes.SolverOptions(
        num_threads=threads,
        time_limit_s=solver_time_limit,
        output_dir=solver_output_dir,
        output_prefix=solver_output_prefix,
        model_prefix="pyomo",
        solver=solver,
    )
    if cycle_decomp_mode == datatypes.CycleDecompOptions.MAX_WEIGHT:
        cycle_decomposition.cycle_decomposition_single_graph(
            parsed_bp_graph,
            solver_options,
            alpha,
            should_postprocess=postprocess_greedy_sol,
            should_force_greedy=True,
            pc_output_option=output_path_constraints,
            ignore_path_constraints=ignore_path_constraints,
        )
    else:
        cycle_decomposition.cycle_decomposition_single_graph(
            parsed_bp_graph,
            solver_options,
            alpha,
            should_postprocess=postprocess_greedy_sol,
            should_force_greedy=False,
            pc_output_option=output_path_constraints,
            ignore_path_constraints=ignore_path_constraints,
        )
    if profile:
        summary.output.add_resource_usage_summary(solver_options)
    print("\nCompleted cycle reconstruction.")


@coral_app.command(
    name="cycle_all",
    help="Pass all breakpoint files in directory to LP solver.",
)
def cycle_decomposition_all_mode(
    ctx: typer.Context,
    bp_dir: Annotated[
        pathlib.Path,
        typer.Option(help="Directory containing breakpoint graph files."),
    ],
    output_prefix: OutputPrefixArg,
    alpha: AlphaArg = 0.01,
    solver_time_limit: SolverTimeLimitArg = 7200,
    threads: ThreadsArg = -1,
    solver: SolverArg = Solver.GUROBI,
    global_time_limit: GlobalTimeLimitArg = 21600,
    cycle_decomp_mode: CycleDecompArg = CycleDecompOptions.MAX_WEIGHT,
    output_path_constraints: OutputPCArg = OutputPCOptions.LONGEST,
    postprocess_greedy_sol: PostProcessFlag = False,
    profile: Annotated[
        bool, typer.Option(help="Profile resource usage.")
    ] = False,
    ignore_path_constraints: IgnorePathConstraintsFlag = False,
    log_file: ReconstructLogArg = None,
) -> None:
    print(
        f"{colorama.Style.DIM}{colorama.Fore.LIGHTYELLOW_EX}"
        f"Performing reconstruction with options: {ctx.params}"
        f"{colorama.Style.RESET_ALL}"
    )

    solver_output_dir = pathlib.Path.cwd()
    solver_output_prefix = output_prefix
    if output_prefix.rfind("/") > 0:
        solver_output_dir = output_prefix[: output_prefix.rfind("/")]
        solver_output_prefix = output_prefix[output_prefix.rfind("/") + 1 :]
    pathlib.Path(f"{solver_output_dir}/models").mkdir(
        parents=True, exist_ok=True
    )

    log_fn = f"{output_prefix}_cycle_decomposition.log"
    if log_file:
        log_fn = log_file
    logging.basicConfig(
        filename=log_fn,
        filemode="w+",
        level=logging.DEBUG,
        format="%(asctime)s:%(levelname)-4s [%(filename)s:%(lineno)d] %(message)s",
    )
    logging.getLogger("pyomo").setLevel(logging.ERROR)
    logging.getLogger("gurobipy").setLevel(logging.ERROR)

    global_state.STATE_PROVIDER.should_profile = profile
    global_state.STATE_PROVIDER.output_prefix = output_prefix
    global_state.STATE_PROVIDER.time_limit_s = global_time_limit

    bp_graphs = get_all_graphs_from_dir(bp_dir)
    for graph in bp_graphs:
        graph.infer_amplicon_intervals(graph.amplicon_idx - 1)

    solver_options = datatypes.SolverOptions(
        num_threads=threads,
        time_limit_s=solver_time_limit,
        output_dir=solver_output_dir,
        output_prefix=solver_output_prefix,
        model_prefix="pyomo",
        solver=solver,
    )
    if cycle_decomp_mode == datatypes.CycleDecompOptions.MAX_WEIGHT:
        cycle_decomposition.cycle_decomposition_all_graphs(
            bp_graphs,
            solver_options,
            alpha,
            should_postprocess=postprocess_greedy_sol,
            should_force_greedy=True,
            pc_output_option=output_path_constraints,
            ignore_path_constraints=ignore_path_constraints,
        )
    else:
        cycle_decomposition.cycle_decomposition_all_graphs(
            bp_graphs,
            solver_options,
            alpha,
            should_postprocess=postprocess_greedy_sol,
            should_force_greedy=False,
            pc_output_option=output_path_constraints,
            ignore_path_constraints=ignore_path_constraints,
        )
    if profile:
        summary.output.add_resource_usage_summary(solver_options)
    print("\nCompleted cycle reconstruction.")


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
    print(
        f"{colorama.Style.DIM}{colorama.Fore.LIGHTYELLOW_EX}"
        f"Performing HSR mode with options: {ctx.params}"
        f"{colorama.Style.RESET_ALL}"
    )
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
    ref: ReferenceGenomeArg,
    output_prefix: OutputPrefixArg,
    graph: Annotated[
        typer.FileText | None,
        typer.Option(help="AmpliconSuite-formatted graph file (*_graph.txt)."),
    ] = None,
    bam: BamArg = None,
    cycles: Annotated[
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
    print(
        f"{colorama.Style.DIM}{colorama.Fore.LIGHTYELLOW_EX}"
        f"Performing plot mode with options: {ctx.params}"
        f"{colorama.Style.RESET_ALL}"
    )
    if "/" in output_prefix:
        os.makedirs(os.path.dirname(output_prefix), exist_ok=True)

    plot_graph = False
    plot_cycles = False
    if graph:
        plot_graph = True
    if cycles:
        plot_cycles = True
    if not graph and not cycles:
        raise typer.BadParameter("Must input graph or cycles file to plot.")
    plot_amplicons.plot_amplicon(
        ref,
        bam,
        graph,
        cycles,
        output_prefix,
        num_cycles,
        max_coverage,
        min_mapq,
        gene_subset_list,
        gene_fontsize,
        region,
        should_plot_graph=plot_graph,
        should_plot_cycles=plot_cycles,
        should_hide_genes=hide_genes,
        should_restrict_to_bushman_genes=bushman_genes,
        should_plot_only_cyclic_walks=only_cyclic_paths,
    )


@coral_app.command(
    name="plot_all",
    help="Generate plots for all amplicons in a given directory.",
)
def plot_all_mode(
    ctx: typer.Context,
    ref: ReferenceGenomeArg,
    output_prefix: OutputPrefixArg,
    reconstruction_dir: Annotated[
        pathlib.Path | None,
        typer.Option(help="Reconstruction directory."),
    ] = None,
    cycles_dir: Annotated[
        pathlib.Path | None,
        typer.Option(help="Cycle directory."),
    ] = None,
    graph_dir: Annotated[
        pathlib.Path | None,
        typer.Option(help="Graph directory."),
    ] = None,
    bam: BamArg = None,
    only_cyclic_paths: Annotated[
        bool, typer.Option(help="Only plot cyclic paths from cycles file.")
    ] = False,
    num_cycles: Annotated[
        int | None, typer.Option(help="Only plot the first NUM_CYCLES cycles.")
    ] = None,
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
    region: Annotated[
        str | None,
        typer.Option(
            help="Specifically visualize only this region, argument formatted as 'chr1:pos1-pos2'."
        ),
    ] = None,
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
    profile: Annotated[
        bool, typer.Option(help="Profile resource usage.")
    ] = False,
) -> None:
    print(
        f"{colorama.Style.DIM}{colorama.Fore.LIGHTYELLOW_EX}"
        f"Performing plot_all mode with options: {ctx.params}"
        f"{colorama.Style.RESET_ALL}"
    )
    if "/" in output_prefix:
        pathlib.Path(output_prefix).mkdir(parents=True, exist_ok=True)
    global_state.STATE_PROVIDER.should_profile = profile

    # TODO: make this into a typer validation function, re-use in score mode
    shared_dir_set = reconstruction_dir is not None
    separate_dirs_set = (cycles_dir is not None) or (graph_dir is not None)
    if shared_dir_set == separate_dirs_set:
        raise typer.BadParameter(
            "Must specify either a shared reconstruction directory or "
            "separate graph or cycles directories."
        )

    if shared_dir_set:
        reconstruction_paths = get_reconstruction_paths_from_shared_dir(
            reconstruction_dir  # type: ignore[arg-type]
        )
        if not reconstruction_paths:
            raise typer.BadParameter(
                f"No reconstruction files found in the given directory: {reconstruction_dir}"
            )
        for graph_path, cycle_path in reconstruction_paths:
            with graph_path.open("r") as graph_file:
                cycle_file = None if cycle_path is None else cycle_path.open("r")
                amplicon_idx = int(
                    graph_path.name.split("_")[-2].split("amplicon")[1]
                )
                plot_amplicons.plot_amplicon(
                    ref,
                    bam,
                    graph_file,
                    cycle_file,
                    output_prefix=f"{output_prefix}_amplicon{amplicon_idx}",
                    num_cycles=num_cycles,
                    max_coverage=max_coverage,
                    min_mapq=min_mapq,
                    gene_subset_list=gene_subset_list,
                    gene_fontsize=gene_fontsize,
                    region=region,
                    should_plot_graph=True,
                    should_plot_cycles=True if cycle_file is not None else False,
                    should_hide_genes=hide_genes,
                    should_restrict_to_bushman_genes=bushman_genes,
                    should_plot_only_cyclic_walks=only_cyclic_paths,
                )
                if cycle_file is not None:
                    cycle_file.close()
    else:
        if graph_dir is not None:
            for graph_path in graph_dir.glob("*_graph.txt"):
                with graph_path.open("r") as graph_file:
                    amplicon_idx = int(
                        graph_path.name.split("_")[-2].split("amplicon")[1]
                    )
                    plot_amplicons.plot_amplicon(
                        ref,
                        bam,
                        graph_file,
                        None,
                        output_prefix=f"{output_prefix}_amplicon{amplicon_idx}",
                        num_cycles=num_cycles,
                        max_coverage=max_coverage,
                        min_mapq=min_mapq,
                        gene_subset_list=gene_subset_list,
                        gene_fontsize=gene_fontsize,
                        region=region,
                        should_plot_graph=True,
                        should_plot_cycles=False,
                        should_hide_genes=hide_genes,
                        should_restrict_to_bushman_genes=bushman_genes,
                        should_plot_only_cyclic_walks=only_cyclic_paths,
                    )
        if cycles_dir is not None:
            for cycles_path in cycles_dir.glob("*_cycles.txt"):
                with cycles_path.open("r") as cycle_file:
                    amplicon_idx = int(
                        cycles_path.name.split("_")[-2].split("amplicon")[1]
                    )
                    plot_amplicons.plot_amplicon(
                        ref,
                        bam,
                        None,
                        cycle_file,
                        output_prefix=f"{output_prefix}_amplicon{amplicon_idx}",
                        num_cycles=num_cycles,
                        max_coverage=max_coverage,
                        min_mapq=min_mapq,
                        gene_subset_list=gene_subset_list,
                        gene_fontsize=gene_fontsize,
                        region=region,
                        should_plot_graph=False,
                        should_plot_cycles=True,
                        should_hide_genes=hide_genes,
                        should_restrict_to_bushman_genes=bushman_genes,
                        should_plot_only_cyclic_walks=only_cyclic_paths,
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
    print(
        f"{colorama.Style.DIM}{colorama.Fore.LIGHTYELLOW_EX}"
        f"Performing cycle to bed mode with options: {ctx.params}"
        f"{colorama.Style.RESET_ALL}"
    )
    cycle2bed.convert_cycles_to_bed(
        cycle_file, output_file, rotate_to_min, num_cycles, print_command = True
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
    output_prefix: OutputPrefixArg,
    name: Annotated[str, typer.Option(help="Name of sample.")],
) -> None:
    print(
        f"{colorama.Style.DIM}{colorama.Fore.LIGHTYELLOW_EX}"
        f"Performing plot mode with options: {ctx.params}"
        f"{colorama.Style.RESET_ALL}"
    )
    plot_cn.plot_cnr(cnr, output_prefix, name)


@coral_app.command(
    name="score", help="Score cycle reconstructions vs. ground-truth."
)
def score_mode(
    ctx: typer.Context,
    ground_truth: Annotated[
        pathlib.Path, typer.Option(help="Ground-truth directory.")
    ],
    output_prefix: OutputPrefixArg,
    reconstruction_dir: Annotated[
        pathlib.Path, typer.Option(help="Reconstruction directory.")
    ],
    cycles_dir: Annotated[
        pathlib.Path | None,
        typer.Option(help="Cycle directory."),
    ] = None,
    tolerance: Annotated[
        int, typer.Option(help="Breakpoint matching tolerance.")
    ] = 100,
    to_skip: Annotated[
        list[str] | None, typer.Option(help="List of datasets to skip.")
    ] = None,
) -> None:
    print(
        f"{colorama.Style.DIM}{colorama.Fore.LIGHTYELLOW_EX}"
        f"Performing score mode with options: {ctx.params}"
        f"{colorama.Style.RESET_ALL}"
    )
    pathlib.Path(f"{output_prefix}").mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        filename=f"{output_prefix}/score_reconstruction.log",
        filemode="w+",
        level=logging.DEBUG,
        format="%(asctime)s:%(levelname)-4s [%(filename)s:%(lineno)d] %(message)s",
    )
    score_simulation.score_simulations(
        ground_truth,
        reconstruction_dir=reconstruction_dir,
        cycles_dir=cycles_dir,  # type: ignore[arg-type]
        output_prefix=output_prefix,
        tolerance=tolerance,
        to_skip=to_skip if to_skip else [],
    )


@coral_app.command(
    name="plot_resources",
    help="Generate plots of resource usage.",
)
def plot_resource_usage(
    ctx: typer.Context,
    reconstruction_dir: Annotated[
        pathlib.Path, typer.Option(help="Reconstruction directory.")
    ],
    output_prefix: OutputPrefixArg,
) -> None:
    print(
        f"{colorama.Style.DIM}{colorama.Fore.LIGHTYELLOW_EX}"
        f"Performing plot resource usage mode with options: {ctx.params}"
        f"{colorama.Style.RESET_ALL}"
    )
    pathlib.Path(f"{output_prefix}").mkdir(parents=True, exist_ok=True)

    summary.parsing.plot_resource_usage(reconstruction_dir, output_prefix)
