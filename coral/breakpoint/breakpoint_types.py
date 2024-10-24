from __future__ import annotations

import logging
from typing import Literal, NamedTuple

from coral import types
from coral.parsing.cns_types import CnsInterval

logger = logging.getLogger(__name__)


class SourceEdge(NamedTuple):
    kind: Literal["source"]  # Unnecessary but left alone to maintain original tuple indexing in code # TODO: remove + update idx
    cn_idx: int  # Always -1 in code, # TODO: remove + update idx
    strand1: types.Strand
    chr1: str  # Chromosome
    pos1: int  # Position
    o1: int  # Offset
    sr_count: int  # Short Read Count
    sr_flag: str  # Short Read Flag
    sr_cn: float  # Short Read CN
    lr_cn: float  # Long Read CN
    cn: float  # CN
