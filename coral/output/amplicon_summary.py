import io
import logging

from coral.breakpoint import infer_breakpoint_graph

logger = logging.getLogger(__name__)


def output_amplicon_info(
    amplicon_idx: int,
    bb: infer_breakpoint_graph.LongReadBamToBreakpointMetadata,
    output_file: io.TextIOWrapper,
) -> None:
    amplicon_intervals = [
        ai
        for ai in bb.amplicon_intervals
        if bb.ccid2id[ai.amplicon_id] == amplicon_idx + 1
    ]

    output_file.write("AmpliconID = " + str(amplicon_idx))
    output_file.write("#Intervals = " + str(len(amplicon_intervals)))
    output_file.write("AmpliconIntervals = ")
