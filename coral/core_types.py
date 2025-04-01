"""Primitive types, associated containers, and useful aliases within Coral."""

import enum

ChrTag = str  # Chromosome identifier, in the form `chr<name>`
CNSIdx = int
ReadName = str
BPIdx = int


class ReferenceGenome(enum.StrEnum):
    """Reference genome."""

    hg19 = "hg19"
    hg38 = "hg38"
    mm10 = "mm10"
