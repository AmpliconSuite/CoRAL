import io
from dataclasses import dataclass
from typing import Any, Dict, List

import intervaltree
import numpy as np

from coral.types import CnsInterval


@dataclass
class CNSSegData:
    tree: intervaltree.IntervalTree
    intervals: List[CnsInterval]
    intervals_by_chr: Dict[str, List[List[int]]]
    log2_cn: List[float]

    @classmethod
    def from_file(cls, file: io.TextIOWrapper) -> CNSSegData:
        tree, log2_cns = {}, []
        intervals = []
        intervals_by_chr: Dict[str, List[List[Any]]] = {}
        is_cns = file.name.endswith(".cns")

        next(file)  # Skip header row
        idx = 0
        for line in file:
            fields = line.strip().split()
            chr_tag, start, end = fields[0], int(fields[1]), int(fields[2])
            intervals.append([chr_tag, start, end - 1])
            if chr_tag not in tree:
                tree[chr_tag] = intervaltree.IntervalTree()
                intervals_by_chr[chr_tag] = []
                idx = 0
            tree[chr_tag][start:end] = idx
            idx += 1

            # Calc log2 CN for .cns
            log2_cn = float(fields[4]) if is_cns else np.log2(float(fields[3]) / 2.0)
            raw_cn = 2 * (2**log2_cn) if is_cns else float(fields[3])

            log2_cns.append(log2_cn)
            intervals_by_chr[chr_tag].append([chr_tag, start, end - 1, raw_cn])

        return cls(tree, intervals, intervals_by_chr, log2_cns)
