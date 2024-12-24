"""Primitive types, associated containers, and useful aliases within Coral."""

from typing import Any, List, Literal

EdgeType = Literal["e", "c", "d", "s", "t", "ns", "nt"]
EdgeIdx = int
EdgeCount = int
AmpliconWalk = dict[tuple[EdgeType, EdgeIdx], EdgeCount]
