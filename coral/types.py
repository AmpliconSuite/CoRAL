from typing import Literal

CigarString = str
Strand = Literal["+", "-"]

AmpliconInterval = tuple[str, int, int, int]
CnsInterval = tuple[str, int, int]
Edge = tuple[int, int]
