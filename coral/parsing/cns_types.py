from __future__ import annotations

import io
import logging
from dataclasses import dataclass
from typing import Any, Dict, List, NamedTuple, Protocol

import intervaltree
import numpy as np

logger = logging.getLogger(__name__)


class RawCNInterval(NamedTuple):
    chr_tag: str
    start: int
    end: int


class CnsIntervalProtocol(Protocol):
    @property
    def begin(self) -> int: ...  # Enforce typing on intervaltree.Interval
    @property
    def end(self) -> int: ...  # Enforce typing on intervaltree.Interval
    @property
    def data(self) -> int: ...  # CN seg file index


class CnsInterval(CnsIntervalProtocol): ...


@dataclass
class CNSegIntervals:
    tree: intervaltree.IntervalTree  # Interval tree structure for each chromosome
    intervals: List[RawCNInterval]  # CN segments
    intervals_by_chr: Dict[str, List[List[int]]]
    log2_cn: List[float]  # log2 CN for each CN segment

    def get_cn_intervals(self, chr_tag: str, pos: int) -> set[CnsInterval]:
        """Get the CN segment intervals for a given position"""
        return self.tree[chr_tag][pos]

    def get_cn_interval(self, chr_tag: str, pos: int) -> CnsInterval:
        """Get (assumed) single CN segment interval for a given position"""
        return next(iter(self.tree[chr_tag][pos]))

    def get_cn_idx(self, chr_tag: str, pos: int) -> int:
        """Get the CN segment index for a given position"""
        return self.get_cn_interval(chr_tag, pos).data

    @classmethod
    def from_file(cls, file: io.TextIOWrapper) -> CNSegIntervals:
        tree, log2_cns = {}, []
        intervals = []
        intervals_by_chr: Dict[str, List[List[Any]]] = {}
        is_cns = file.name.endswith(".cns")

        next(file)  # Skip header row
        idx = 0
        for line in file:
            fields = line.strip().split()
            chr_tag, start, end = fields[0], int(fields[1]), int(fields[2])
            intervals.append(RawCNInterval(chr_tag, start, end - 1))
            if chr_tag not in tree:
                tree[chr_tag] = intervaltree.IntervalTree()
                intervals_by_chr[chr_tag] = []
                idx = 0
                logger.info(f"resetting idx")
            tree[chr_tag][start:end] = idx
            idx += 1

            # Calc log2 CN for .cns
            log2_cn = float(fields[4]) if is_cns else np.log2(float(fields[3]) / 2.0)
            raw_cn = 2 * (2**log2_cn) if is_cns else float(fields[3])

            log2_cns.append(log2_cn)
            intervals_by_chr[chr_tag].append([chr_tag, start, end - 1, raw_cn])

        return cls(tree, intervals, intervals_by_chr, log2_cns)
