"""Primitive types, associated containers, and useful aliases within Coral."""

import enum
import sys

if sys.version_info >= (3, 11):
    from enum import StrEnum
else:
    class StrEnum(str, enum.Enum):  # type: ignore[no-redef]
        pass


def chr_sort_key(name: str) -> tuple[int, int, str]:
    """Natural sort key for chromosome names, compatible with any reference genome.

    Sorts autosomes numerically (chr1 < chr2 < ... < chr22), then sex
    chromosomes (X, Y), then mitochondrial (M/MT), then any other contig
    lexicographically. Handles both 'chr'-prefixed and bare names.
    """
    stem = name[3:] if name[:3].lower() == "chr" else name
    if stem.upper() in ("M", "MT"):
        return (2, 0, "")
    if stem.upper() == "X":
        return (1, 0, "")
    if stem.upper() == "Y":
        return (1, 1, "")
    try:
        return (0, int(stem), "")
    except ValueError:
        return (1, 2, stem)


def is_canonical_chr(name: str) -> bool:
    """Whether `name` is a primary chromosome (autosome, X, Y, or MT).

    Genome-agnostic replacement for the former hardcoded hg38 ``CHR_TAG_TO_IDX``
    membership test: returns True for chr1..chrN, chrX, chrY, chrM/chrMT (with
    or without the ``chr`` prefix), and False for decoy/alt/unplaced contigs
    (e.g. ``chrUn_*``, ``*_random``, ``*_alt``, ``HLA-*``, ``chrEBV``).
    """
    return chr_sort_key(name)[:2] != (1, 2)


ChrTag = str  # Chromosome identifier, in the form `chr<name>`
CNSIdx = int
ReadName = str
BPIdx = int


class ReferenceGenome(StrEnum):
    """Reference genome."""

    hg19 = "hg19"
    hg38 = "hg38"
    mm10 = "mm10"
    t2t = "t2t"
    other = "other"
