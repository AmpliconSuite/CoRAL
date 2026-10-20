"""Exercise seed regrouping after chromosome-arm CN filtering."""

import pathlib

import pytest

from coral.cnv_seed import run_seeding


@pytest.mark.parametrize(
    ("segments", "expected"),
    [
        pytest.param(
            [(50_000, 10), (400_000, 6), (200_000, 10)],
            [(1_450_000, 1_649_999)],
            id="short-first-block-does-not-discard-next-peak",
        ),
        pytest.param(
            [(200_000, 10), (400_000, 6), (200_000, 10)],
            [(1_000_000, 1_199_999), (1_600_000, 1_799_999)],
            id="count-first-segment-after-writing-a-block",
        ),
        pytest.param(
            [(200_000, 10), (400_000, 6), (60_000, 10), (60_000, 12)],
            [(1_000_000, 1_199_999), (1_600_000, 1_719_999)],
            id="count-all-segments-in-subsequent-block",
        ),
        pytest.param(
            [(50_000, 10), (400_000, 10), (200_000, 10)],
            [(1_000_000, 1_649_999)],
            id="surviving-bridge-keeps-block-merged",
        ),
        pytest.param(
            [(50_000, 10), (400_000, 6), (200_000, 7)],
            [],
            id="below-arm-cutoff-peak-remains-filtered",
        ),
        pytest.param(
            [(200_000, 7), (400_000, 6), (200_000, 7)],
            [],
            id="all-candidate-segments-filtered",
        ),
        pytest.param(
            [(200_000, 10), (400_000, 6), (50_000, 10)],
            [(1_000_000, 1_199_999)],
            id="short-final-block-remains-filtered",
        ),
    ],
)
def test_regroup_after_arm_cn_filtering(
    tmp_path: pathlib.Path,
    segments: list[tuple[int, int]],
    expected: list[tuple[int, int]],
) -> None:
    # A fully covered 10 Mb p arm with background CN=4 gives a cutoff of
    # gain + arm CN - 2 = 8. CN=6 bridges enter the initial candidate but
    # are removed by arm filtering, leaving gaps exceeding 300 kb.
    centromere = tmp_path / "centromere.bed"
    centromere.write_text("chr1\t10000000\t11000000\n")
    rows = ["chr1\t0\t1000000\t4\n"]
    start = 1_000_000
    for length, cn in segments:
        end = start + length
        rows.append(f"chr1\t{start}\t{end}\t{cn}\n")
        start = end
    rows.append(f"chr1\t{start}\t10000000\t4\n")
    calls = tmp_path / "input.bed"
    calls.write_text("".join(rows))

    with calls.open() as cn_seg_file:
        run_seeding(
            cn_seg_file,
            str(tmp_path / "sample"),
            gain=6,
            min_seed_size=100_000,
            max_seg_gap=300_000,
            centromere_file=centromere,
            chr_sizes={"chr1": 20_000_000},
        )

    assert (tmp_path / "sample_CNV_SEEDS.bed").read_text() == "".join(
        f"chr1\t{start}\t{end}\n" for start, end in expected
    )
