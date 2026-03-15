"""Primitive types, associated containers, and useful aliases within Coral."""

import enum


def chr_sort_key(name: str) -> tuple[int, int, str]:
    """Natural sort key for chromosome names, compatible with any reference genome.

    Sorts autosomes numerically (chr1 < chr2 < ... < chr22), then sex
    chromosomes (X, Y), then mitochondrial (M/MT), then any other contig
    lexicographically. Handles both 'chr'-prefixed and bare names.
    """
    stem = name.removeprefix("chr")
    if stem in ("M", "MT"):
        return (2, 0, "")
    if stem == "X":
        return (1, 0, "")
    if stem == "Y":
        return (1, 1, "")
    try:
        return (0, int(stem), "")
    except ValueError:
        return (1, 2, stem)

ChrTag = str  # Chromosome identifier, in the form `chr<name>`
CNSIdx = int
ReadName = str
BPIdx = int


class ReferenceGenome(enum.StrEnum):
    """Reference genome."""

    hg19 = "hg19"
    hg38 = "hg38"
    mm10 = "mm10"
