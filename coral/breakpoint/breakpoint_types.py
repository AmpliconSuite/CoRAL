"""Type aliases and data containers for breakpoint graph elements."""

from __future__ import annotations

import io
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List

import numpy as np

from coral.datatypes import CNSInterval, CNSIntervalTree


@dataclass
class CNSSegData:
    """Container class for mapping CN (copy number) per chromosome segment, as parsed from a .cns or .bed file."""

    tree: Dict[
        str, CNSIntervalTree
    ]  # IntervalTree mapping chromosome segments to their (per-chromosome) index in input file
    # Raw list of all chromosome intervals given, of the form
    # (chromosome, start_idx, end_idx)
    intervals: List[CNSInterval]
    # List of chromosome intervals with their respective CNs,
    # grouped by chromosome number
    intervals_by_chr: Dict[str, List[CNSInterval]]
    log2_cn: List[float]  # Log2 of CN for each chromosome segment

    @classmethod
    def from_file(cls, file: io.TextIOWrapper) -> CNSSegData:
        """Parse CN interval data from a .cns or .bed file."""
        tree, log2_cns = {}, []
        intervals = []
        intervals_by_chr: Dict[str, List[CNSInterval]] = defaultdict(list)
        is_cns = file.name.endswith(".cns")

        idx = 0
        for line in file:
            fields = line.strip().split()
            if (
                not fields or fields[0] == "chromosome" or line.startswith("#")
            ):  # Skip potential header row
                continue
            chr_tag, start, end = fields[0], int(fields[1]), int(fields[2])

            # Calc log2 CN for .cns
            log2_cn = (
                float(fields[4]) if is_cns else np.log2(float(fields[3]) / 2.0)
            )
            raw_cn = 2 * (2**log2_cn) if is_cns else float(fields[3])
            log2_cns.append(log2_cn)
            cns_interval = CNSInterval(chr_tag, start, end - 1, raw_cn)

            intervals.append(cns_interval)
            intervals_by_chr[chr_tag].append(cns_interval)
            if chr_tag not in tree:
                tree[chr_tag] = CNSIntervalTree()
                idx = 0
            tree[chr_tag][start:end] = idx
            idx += 1

        return cls(tree, intervals, intervals_by_chr, log2_cns)
