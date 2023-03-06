"""
An algorithm for refining an existing breakpoint graph with long-read data.
"""
from typing import Union

import typer

from long_read_aa.breakpoints import BreakpointGraph, BreakpointEdge, SequenceEdge

app = typer.Typer()


def initialize_breakpoint_graph(amplicon_graph: str, amplicon_cycles: str) -> BreakpointGraph.BreakpointGraph:
    """Initializes a breakpoint graph with existing AmpliconArchitect data.

    Args:
        amplicon_graph: AmpliconArchitect amplicon graph
        amplicon_cycles: AmpliconArchitect amplicon cycles

    Returns:
        A BreakpointGraph.
    """
    return BreakpointGraph.BreakpointGraph()

def fetch_long_read_chimeric_alignments(breakpoint_graph: BreakpointGraph.BreakpointGraph, long_read_bam_file: str) -> BreakpointGraph.BreakpointGraph:
    """Fetches long-read chimeric alignments.

    Returns:
        An updated BreakpointGraph.
    """
    return breakpoint_graph


@app.command()
def refine_breakpoint_graph(
    short_read_bam_file: str = typer.Argument(..., help="Path to short-read BAM file."),
    long_read_bam_file: str = typer.Argument(..., help="Path to long-read BAM file."),
    amplicon_graph: str = typer.Argument(..., help="Path to amplicon graph file."),
    amplicon_cycles: str = typer.Argument(..., help="Path to amplicon cycles file."),
    short_read_copy_number_segments: str = typer.Option(None, help="Path to copy-number segements (.cns) for short-read data)."),
    long_read_copy_number_segments: str = typer.Option(None, help="Path to copy-number segements (.cns) for long-read data)."),
    max_bp_match_cutoff: int = typer.Option(100, help="Maximum distance for matching breakpoints."),
    min_cluster_cutoff: int = typer.Option(3, help="Hard cutoff for considering a long-read breakpoint cluster."),
    max_bp_distance_cutoff: int = typer.Argument(2000, help="Maximum distance for adding a breakpoint to a cluster."),
    small_deletion_cutoff: int = typer.Argument(10000, help="Cutoff for handling small deletions during refinement."),
    min_deletion_length: int = typer.Argument(10000, help="MInimum length of all deletion breakpoints returned by AmpliconArchitect.")
):
    """Refines an existing breakpoint graph with long reads.

    Args:
        short_read_bam_file: Path to short-read BAM file.
        long_read_bam_file: Path to long-read BAM file.
        amplicon_graph: Path to amplicon graph file, as inferred with AmpliconArchitect.
        amplicon_cycles: Path to ampicon cycles file, as inferred with AmpliconArchitect.
        max_bp_match_cutoff: Maximum distance for matching breakpoints.
        min_cluster_cutoff: Minimum cluster size for considering a long read breakpoint cluster.
        max_bp_distance_cutoff: Maximum distance allowed for adding a breakpiont to a cluster.
        small_deletion_cutoff: Cutoff for hadnlign small deletions during refinement.
        min_deletion_length: Minimum length of deletions allowed, as returned by AmpliconArchitect.
    """

    min_bp_match_cutoff = [] # Each breakpoint has a unique matching cutoff

    read_length = dict() # Map read name -> read length
    chimeric_alignments = dict() # Map read name -> chimeric alignments (two or more records for one read)
    large_indel_alignments = dict() # Map read name -> alignments with one record per read but large indels showing in CIGAR string

    amplicon_intervals = [] # Amplicon intervals
    seq_edges = [] # Sequence edges
    discordant_edges = [] # Discordant edges
    concordant_edges = [] # Concordant edges

    discordant_edges_pos = dict() 
    small_del_indices = [] # Indices of +- breakpoints on the same chr and < small_del_cutoff in discordant_edges
    source_edges = [] # Special AA discordant edges connected to a source node 
    new_bp_list_ = [] # Long read only breakpoints (discordant edges)
    new_bp_stats_ = [] # Statistics of long read only breakpoints (discordant edges)

    sr_length = 0.0 # Average short read length
    min_sr_alignment_length = 30.0 # Minimum alignment length of short reads - used for CN assignment
    max_sr_insert = 2000
    normal_cov_sr = 0 # Normal short read coverage - used for CN assignment
    normal_cov_lr = 0 # Normal long read coverage - used for CN assignment

    breakpoint_graph = initialize_breakpoint_graph(amplicon_graph, amplicon_cycles, short_read_copy_number_segments, long_read_copy_number_segments)

    fetch_long_read_chimeric_alignments(breakpoint_graph, long_read_bam_file, copy=False) 
