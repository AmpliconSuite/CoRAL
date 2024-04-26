#!/usr/bin/env python3
import argparse

from cnv_seed import *
# importing of most modules is done after run mode is selected to minimize chances of colliding function or variable names, existing now and in the future.

def seed_mode(args):
    print("Performing seeding mode with options:", args)
    run_seeding(args)


def reconstruct_mode(args):
    print("Performing reconstruction with options:", args)
    import infer_breakpoint_graph
    infer_breakpoint_graph.reconstruct(args)


def hsr_mode(args):
    print("Performing HSR mode with options:", args)
    import hsr
    hsr.locate_hsrs(args)


def plot_mode(args):
    print("Performing plot mode with options:", args)
    import plot_amplicons
    plot_amplicons.plot_amplicons(args)


def cycle2bed_mode(args):
    print("Performing cycle to bed mode with options:", args)
    import cycle2bed
    cycle2bed.convert_cycles_to_bed(args.cycle_fn, args.output_fn, args.num_cycles)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Long-read amplicon reconstruction pipeline and associated utilities.")
    # Subparsers for different modes (seed, reconstruct, hsr, plot, cycle2bed)
    subparsers = parser.add_subparsers(dest="mode", help="Select mode")

    # SEEDING MODE ARGS
    seed_parser = subparsers.add_parser("seed", help="Filter and merge amplified intervals.")
    seed_parser.add_argument('--cn_seg', help="Input copy number segment file from Cnvkit", type=str,
                             required=True)
    seed_parser.add_argument('--out',
                        help="OPTIONAL: Prefix filename for output bed file. Default: <INPUT_CNS_BASENAME>_CNV_SEEDS.bed",
                        type=str, default='')
    seed_parser.add_argument('--gain',
                        help="OPTIONAL: CN gain threshold for interval to be considered as a seed. Default: 6",
                        type=float, default=GAIN)
    seed_parser.add_argument('--min_seed_size',
                        help="OPTIONAL: Minimum size (in bp) for interval to be considered as a seed. Default: 100000",
                        type=int, default=CNSIZE_MIN)
    seed_parser.add_argument('--max_seg_gap',
                        help="OPTIONAL: Maximum gap size (in bp) to merge two proximal intervals. Default: 300000",
                        type=int, default=CNGAP_MAX)

    # RECONSTRUCTION MODE ARGUMENTS
    reconstruct_parser = subparsers.add_parser("reconstruct", help="Reconstruct focal amplifications")
    reconstruct_parser.add_argument("--lr_bam", help="Sorted indexed (long read) bam file.", required=True)
    reconstruct_parser.add_argument("--seed", help="File including seed intervals.", required=True)
    reconstruct_parser.add_argument("--output_prefix", help="Prefix of output files.", required=True)
    reconstruct_parser.add_argument("--cnseg", help="Long read CNV segmentation file.", required=True)
    reconstruct_parser.add_argument("--output_bp", help="If specified, only output the list of breakpoints.",
                                    action='store_true')
    reconstruct_parser.add_argument("--output_all_path_constraints",
                        help="If specified, output all path constraints in *.cycles file.", action='store_true')
    reconstruct_parser.add_argument("--min_bp_support",
                        help="Ignore breakpoints with less than (min_bp_support * normal coverage) long read support.",
                        type=float, default=1.0)
    reconstruct_parser.add_argument("--cycle_decomp_alpha",
                        help="Parameter used to balance CN weight and path constraints in greedy cycle extraction.",
                        type=float)
    reconstruct_parser.add_argument("--cycle_decomp_time_limit",
                        help="Maximum running time (in seconds) reserved for integer program solvers.", type=int)
    reconstruct_parser.add_argument("--cycle_decomp_threads", help="Number of threads reserved for integer program solvers.",
                        type=int)
    reconstruct_parser.add_argument("--postprocess_greedy_sol",
                        help="Postprocess the cycles/paths returned in greedy cycle extraction.", action='store_true')
    reconstruct_parser.add_argument("--log_fn", help="Name of log file.")

    # HSR MODE ARGS
    hsr_parser = subparsers.add_parser("hsr", help="Detect possible integration points of ecDNA HSR amplifications.")
    hsr_parser.add_argument("--lr_bam", help="Sorted indexed long read bam file.", required=True)
    hsr_parser.add_argument("--cycles", help="Cycle file in *.bed format.", required=True)
    hsr_parser.add_argument("--cn_segments", help="Cnvkit *.cns file.", required=True)
    hsr_parser.add_argument("--out_prefix", help="Prefix of output file name.", required=True)
    hsr_parser.add_argument("--normal_cov", help="Estimated diploid coverage.", required=True)
    hsr_parser.add_argument("--bp_match_cutoff", help="Breakpoint matching cutoff.", type=int, default=100)
    hsr_parser.add_argument("--bp_match_cutoff_clustering",
                        help="Crude breakpoint matching cutoff for clustering.", type=int, default=2000)

    # PLOT MODE ARGS
    plot_parser = subparsers.add_parser("plot", help="Generate plots of amplicon cycles and/or graph from AA-formatted output files")
    plot_parser.add_argument("--ref", help="Name of reference genome used",
                             choices=["hg19", "hg38", 'GRCh38', "mm10"], required=True)
    plot_parser.add_argument("--bam", help="Sorted & indexed bam file.")
    plot_parser.add_argument("--graph", help="AmpliconSuite-formatted *.graph file.")
    plot_parser.add_argument("--cycles", help="AmpliconSuite-formatted cycles file, in *.bed format.")
    plot_parser.add_argument("--output_prefix", "-o", help="Prefix of output files.", required=True)
    plot_parser.add_argument("--plot_graph", help="Visualize breakpoint graph.", action='store_true')
    plot_parser.add_argument("--plot_cycles", help="Visualize (selected) cycles.", action='store_true')
    plot_parser.add_argument("--only_cyclic_paths", help="Only plot cyclic paths from cycles file",
                             action='store_true')
    plot_parser.add_argument("--num_cycles", help="Only plot the first NUM_CYCLES cycles.", type=int)
    plot_parser.add_argument("--max_coverage", help="Limit the maximum visualized coverage in the graph",
                             type=float, default=float('inf'))
    plot_parser.add_argument("--min_mapq", help="Minimum mapping quality to count read in coverage plotting",
                             type=float, default=0)
    plot_parser.add_argument("--gene_subset_list", help="List of genes to visualize (will show all by default)",
                             nargs='+', default=[])
    plot_parser.add_argument("--hide_genes", help="Do not show gene track", action='store_true',
                             default=False)
    plot_parser.add_argument("--gene_fontsize", help="Change size of gene font", type=float, default=12)
    plot_parser.add_argument("--bushman_genes", help="Reduce gene set to the Bushman cancer-related gene set",
                             action='store_true', default=False)
    plot_parser.add_argument("--region",
                             help="(Graph plotting) Specifically visualize only this region, argument formatted as 'chr1:pos1-pos2'.")

    # CYCLE2BED MODE ARGS
    c2b_parser = subparsers.add_parser("cycle2bed", help="Convert cycle files in AA format to bed format.")
    c2b_parser.add_argument("--cycle_fn", help="Input AA-formatted cycle file.", required=True)
    c2b_parser.add_argument("--output_fn", help="Output file name.", required=True)
    c2b_parser.add_argument("--num_cycles", help="If specified, only convert the first NUM_CYCLES cycles.",
                            type=int)


    args = parser.parse_args()

    if args.mode == "seed":
        seed_mode(args)
    elif args.mode == "reconstruct":
        reconstruct_mode(args)
    elif args.mode == "hsr":
        hsr_mode(args)
    elif args.mode == "plot":
        plot_mode(args)
    elif args.mode == "cycle2bed":
        cycle2bed_mode(args)
    else:
        parser.print_help()


