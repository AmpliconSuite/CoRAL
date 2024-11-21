"""Primitive types, associated containers, and useful aliases within Coral."""

from typing import Any, List, Literal

AmpliconInterval = List[Any]  # tuple[str, int, int, int]
CnsInterval = Any  # tuple[str, int, int]

EdgeType = Literal["e", "c", "d", "s", "t", "ns", "nt"]
EdgeIdx = int
EdgeCount = int
AmpliconWalk = dict[tuple[EdgeType, EdgeIdx], EdgeCount]
