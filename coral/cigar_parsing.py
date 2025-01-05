"""
Convert cigar strings from SA:Z: field into chimeric alignments
A chimeric alignment is represented by the following 3 lists
qint: 0-based intervals on query (reads) for each local alignment
rint: 0-based intervals on reference genome for each local alignment
qual: mapping qualities for each local alignment

TBD:
implement merge_alignment
"""

from __future__ import annotations

import logging
import re

from coral import datatypes
from coral.datatypes import (
    ChimericAlignment,
    CigarAlignment,
    CigarBounds,
    Strand,
)

logger = logging.getLogger(__name__)


def convert_pbmm2_to_bwa_mem(cigar_str):
    """
    Convert pbmm2 or "= X" cigar string to the format used by minimap2, bwa mem.

    This function is not currently used but may be useful to have in the future.

    Args:
        cigar_str: CIGAR string in format *=*X
    Returns:
        updated cigar string in the *S*M format
    """
    # Regular expression to capture CIGAR elements (number + character)
    cigar_elements = re.findall(r"(\d+)([A-Z=])", cigar_str)

    converted_cigar = []
    current_M_length = 0

    for length, op in cigar_elements:
        length = int(length)
        if op == "=" or op == "X":
            # For = and X, treat them as matches (M in BWA-MEM)
            current_M_length += length
        else:
            # For other operations, push the accumulated M, then the operation
            if current_M_length > 0:
                converted_cigar.append(f"{current_M_length}M")
                current_M_length = 0
            converted_cigar.append(f"{length}{op}")

    # After loop, if there's any accumulated M left, add it
    if current_M_length > 0:
        converted_cigar.append(f"{current_M_length}M")

    # logging.debug("Original cigar: %s" %cigar_str)
    # logging.debug("Converted cigar: %s" %(''.join(converted_cigar)))
    return "".join(converted_cigar)


def query_ends_from_cigar(cigar_str: str, strand: str) -> CigarAlignment:
    """
    Retrieve alignment ends and alignment length from cigar string

    Args:
        cigar_str: CIGAR string
        strand: Strand of the read ("+" or "-")

    Returns:
        The start and end position on the read on positive strand, and the alignment length on
        the reference genome.
    """
    # Consumable ops on the reference genome (M, D, N, =, X)
    ref_consumable_ops = {"M", "D", "N", "=", "X"}
    query_consumable_ops = {"M", "I", "=", "X"}  # Consumable ops on the query

    query_start = 0
    query_consumed = 0
    ref_consumed = 0  # This will track the length on the reference genome

    # Parse the CIGAR string using regex to extract (length, operation) tuples
    cigar = re.findall(r"(\d+)([A-Z=])", cigar_str)

    # If the strand is negative, reverse the CIGAR operations
    if strand == "-":
        cigar = cigar[::-1]

    # Process CIGAR to find query and reference alignment lengths
    for length, cigar_op in cigar:
        length = int(length)  # Convert length to an integer

        # Handle soft/hard clipping which affects query_start but not reference
        if cigar_op == "S" or cigar_op == "H":  # Soft/Hard clip
            if query_consumed == 0:  # Before any alignment operation
                query_start += length

        # Handle query alignment (including matches and insertions)
        if cigar_op in query_consumable_ops:
            query_consumed += length

        # Handle reference alignment (e.g., matches, deletions, skipped regions)
        if cigar_op in ref_consumable_ops:
            ref_consumed += length

    query_end = query_start + query_consumed
    ref_alignment_length = ref_consumed  # Reference length consumed

    return CigarAlignment(query_start, query_end, ref_alignment_length)


def alignment_from_satags(
    sa_list: list[str], read_name: str
) -> list[ChimericAlignment]:
    """
    Convert "SA:Z" a list of strings into a new chimeric alignment.
    Require at least one (soft) clip and one match for each canonical alignment
    record in a chimeric alignment.
        If not, trigger a warning message in logger

    Args:
        sa_list: A list of "SA:Z" tags from bam
        read_name: Read name
    Returns:
        chimeric alignments in the form of qint, rint and qual list
        Alignments sorted according to the starting positions on the read on
        positive strand
    """
    alignments = []
    for sa in sa_list:
        t = sa.split(",")
        if "S" not in t[3] or ("M" not in t[3] and "=" not in t[3]):
            # Require a chimeric alignment record having at least some
            # (soft)clips and matches
            logger.warning(
                "Found chimeric alignment without match or soft clips."
            )
            logger.warning(f"All CIGAR strings: {sa_list}")
            return []

        qs, qe, al = query_ends_from_cigar(t[3], t[2])
        chr_tag, pos, strand = t[0], int(t[1]), t[2]
        # Convert to 0-based coordinates
        if strand == "+":
            read_interval = datatypes.ReadInterval(
                chr_tag, pos - 1, pos + al - 2, Strand.FORWARD, read_name
            )
        else:
            read_interval = datatypes.ReadInterval(
                chr_tag, pos + al - 2, pos - 1, Strand.REVERSE, read_name
            )
        alignments.append(
            ChimericAlignment(
                CigarBounds(qs, qe),
                read_interval,
                int(t[4]),
                float(t[-1]) / (qe - qs),
            )
        )
    return alignments
