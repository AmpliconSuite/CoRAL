"""Tests for the x-axis labels of the HSR integration-site plot.

``coral hsr`` takes chromosome sizes from the BAM header, so the number of
chromosomes plotted varies (e.g., chrM is included when present). The labels
were hard-coded to 1-22, X, Y, which made matplotlib raise when the counts
differed ("The number of FixedLocator locations (25) ... does not match the
number of labels (24)").
"""

from __future__ import annotations

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import pytest

from coral import hsr
from coral.constants import CHR_SIZES

HG38_WITH_CHRM = {**CHR_SIZES, "chrM": 16569}
MM10_LIKE = {f"chr{i}": 1000 for i in [*range(1, 20), "X", "Y", "M"]}
NO_PREFIX = {"1": 1000, "2": 1000, "X": 1000, "MT": 16569}


@pytest.mark.parametrize(
    "chr_sizes",
    [CHR_SIZES, HG38_WITH_CHRM, MM10_LIKE, NO_PREFIX],
    ids=["hg38", "hg38_chrM", "mm10", "no_prefix"],
)
def test_one_label_per_plotted_chromosome(chr_sizes: dict[str, int]) -> None:
    labels = hsr.chr_tick_labels(chr_sizes)
    assert len(labels) == len(chr_sizes)

    # The same call `locate_hsrs` makes; raised ValueError before the fix
    # whenever the BAM header did not have exactly 24 standard chromosomes.
    fig = plt.figure()
    try:
        plt.xticks(range(len(chr_sizes)), labels)
    finally:
        plt.close(fig)


def test_labels_strip_chr_prefix_and_keep_order() -> None:
    assert hsr.chr_tick_labels(HG38_WITH_CHRM)[-3:] == ["X", "Y", "M"]
    assert hsr.chr_tick_labels(NO_PREFIX) == ["1", "2", "X", "MT"]
