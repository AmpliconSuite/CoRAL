from __future__ import annotations

import copy
import importlib.resources
import logging
import pathlib

import typer

from coral import supplemental_data
from coral.constants import CHR_SIZES, CNSIZE_MAX
from coral.datatypes import ChrArmInfo, CNInterval, Interval, SingleArmInfo

logger = logging.getLogger(__name__)


def parse_centromere_arms(
    centromere_file: pathlib.Path | None = None,
    chr_sizes: dict[str, int] | None = None,
) -> dict[str, ChrArmInfo]:
    """Parse a centromere BED file into per-chromosome arm info.

    Each contig may have any number of BED entries. All entries for a contig
    must form a single contiguous block (overlapping or directly abutting); if
    two or more distinct non-adjacent regions are found for the same contig, a
    ValueError is raised. Contigs absent from the file are simply not included
    in the returned dict — callers should use .get() and handle the missing
    case (e.g. viral contigs or unplaced scaffolds with no centromere).

    Args:
        centromere_file: Path to a BED file (chr, start, end; extra columns
            are ignored). Defaults to the bundled GRCh38 file.
        chr_sizes: Mapping of chromosome name to length, used to compute
            q-arm sizes. Defaults to the hardcoded hg38 constants.
    """
    if chr_sizes is None:
        logger.warning(
            "No BAM file provided to seed mode -- falling back to hardcoded "
            "hg38 chromosome sizes. If your data is not aligned to hg38, "
            "provide --lr-bam to read chromosome sizes from the BAM header."
        )
        chr_sizes = CHR_SIZES

    if centromere_file is not None:
        ctx = open(centromere_file)
    else:
        ctx = (  # type: ignore[assignment]
            importlib.resources.files(supplemental_data) / "GRCh38_centromere.bed"
        ).open("r")

    regions_by_chr: dict[str, list[Interval]] = {}
    with ctx as fp:
        for line in fp:
            line = line.strip()
            if not line:
                continue
            pieces = line.split()
            chrom = pieces[0]
            regions_by_chr.setdefault(chrom, []).append(
                Interval(chrom, int(pieces[1]), int(pieces[2]))
            )

    chr_arms: dict[str, ChrArmInfo] = {}
    for chrom, regions in regions_by_chr.items():
        regions.sort(key=lambda r: r.start)
        for i in range(len(regions) - 1):
            if regions[i].end < regions[i + 1].start:
                raise ValueError(
                    f"Centromere file contains multiple distinct non-adjacent "
                    f"regions for '{chrom}': "
                    f"{chrom}:{regions[i].start}-{regions[i].end} and "
                    f"{chrom}:{regions[i + 1].start}-{regions[i + 1].end}. "
                    f"All centromere regions for a contig must overlap or abut."
                )
        centro = Interval(chrom, regions[0].start, max(r.end for r in regions))
        chr_size = chr_sizes.get(chrom, centro.end)
        chr_arms[chrom] = ChrArmInfo(
            interval=centro,
            p_arm=SingleArmInfo(
                Interval(chrom, 0, centro.start),
                size=centro.start,
            ),
            q_arm=SingleArmInfo(
                Interval(chrom, centro.end, chr_size),
                size=max(0, chr_size - centro.end),
            ),
        )

    return chr_arms


def aggregate_arm_cn(arm: SingleArmInfo) -> float:
    ccn = 2.0
    if arm.total_length < 0.5 * arm.size:
        return ccn

    sum_cns_len = 0
    for cns in sorted(arm.segs, key=lambda cns: cns.cn):
        ccn = cns.cn
        sum_cns_len += len(cns)
        # TODO: shouldn't we use 0.5 for both thresholds?
        if sum_cns_len >= 0.49 * arm.total_length:
            break
    return ccn


def run_seeding(
    cn_seg_file: typer.FileText,
    output_prefix: str,
    gain: float,
    min_seed_size: float,
    max_seg_gap: float,
    centromere_file: pathlib.Path | None = None,
    chr_sizes: dict[str, int] | None = None,
) -> None:
    """Generate seed intervals from file containing WGS CN calls.

    Breakpoint graph reconstruction initially requires a set of focally
    amplified seed intervals, from which breakpoint edges are explored. This
    method produces these intervals using the given parameters as heuristic
    cutoffs.

    Args:
        cn_seg_file: File containing long-read segmented whole genome CN calls.
        output_prefix: Prefix for output file.
        gain: Minimum CN threshold for an interval to be considered as a seed.
        min_seed_size: Minimum size (in base pairs) of a seed interval.
        max_seg_gap: Maximum gap size (in base pairs) between two adjacent
            potential seed intervals for them to be merged into a single seed.
            If merged, min_seed_size is enforced on the combined interval.
        centromere_file: Path to a centromere BED file (paired p/q arm format).
            Defaults to the bundled GRCh38 file.
        chr_sizes: Chromosome name-to-length mapping. Defaults to the hardcoded
            hg38 constants. Pass the result of core_utils.build_chr_sizes_from_bam()
            when working with a non-hg38 genome.

    """

    chr_arms = parse_centromere_arms(centromere_file, chr_sizes)
    cnv_seeds: list[list[CNInterval]] = []
    cur_seed: list[CNInterval] = []
    for line in cn_seg_file:
        s = line.strip().split()
        if s[0] != "chromosome":
            chr_tag, start, end = s[0], int(s[1]), int(s[2])
            arm_info = chr_arms.get(chr_tag)
            arm_intv = arm_info.interval if arm_info is not None else None
            if cn_seg_file.name.endswith(".cns"):
                cn = 2 * (2 ** float(s[4]))
            elif cn_seg_file.name.endswith(".bed"):
                cn = float(s[-1])  # CN is always the last column in BED input
            else:
                logger.error(cn_seg_file.name + "\n")
                raise SystemExit("Invalid cn_seg file format!\n")
            # Require absolute CN >= max(gain, cn_cutoff_chrarm)
            cn_intv = CNInterval(chr_tag, start, end, cn)
            # Exclude segments overlapping the centromere; if no centromere is
            # defined for this contig, all segments are eligible.
            if cn >= gain and (
                arm_intv is None
                or end <= arm_intv.start
                or start >= arm_intv.end
            ):
                if (
                    len(cur_seed) > 0
                    and chr_tag == cur_seed[-1].chr
                    and start - cur_seed[-1].end <= max_seg_gap
                ):
                    cur_seed.append(cn_intv)
                elif len(cur_seed) == 0:
                    cur_seed = [cn_intv]
                else:
                    cnv_seeds.append(cur_seed)
                    cur_seed = [cn_intv]
            if arm_info is not None and arm_intv is not None:
                if end <= arm_intv.start:
                    arm_info.p_arm.segs.append(cn_intv)
                elif start >= arm_intv.end:
                    arm_info.q_arm.segs.append(cn_intv)

    # Add final seed if non-empty
    if cur_seed:
        cnv_seeds.append(cur_seed)

    for arm_info in chr_arms.values():
        arm_info.p_arm.ccn = aggregate_arm_cn(arm_info.p_arm)
        arm_info.q_arm.ccn = aggregate_arm_cn(arm_info.q_arm)

    if output_prefix:
        output_filename = f"{output_prefix}_CNV_SEEDS.bed"
    else:
        output_filename = cn_seg_file.name.replace(".cns", "CNV_SEEDS.bed")

    with open(output_filename, "w") as fp:
        for seed_intvs in cnv_seeds:
            sum_seed_len = sum([len(cns) for cns in seed_intvs])
            cn_cutoff_chrarm = gain
            seed_chr_tag = seed_intvs[-1].chr
            arm_info = chr_arms.get(seed_chr_tag)
            if sum_seed_len > CNSIZE_MAX:
                cn_cutoff_chrarm = 1.2 * gain
            if arm_info is not None:
                if seed_intvs[-1].end <= arm_info.interval.start:  # p arm
                    cn_cutoff_chrarm += arm_info.p_arm.ccn - 2.0
                elif seed_intvs[0].start >= arm_info.interval.end:  # q arm
                    cn_cutoff_chrarm += arm_info.q_arm.ccn - 2.0
                else:
                    logger.warning(
                        f"Seed on {seed_chr_tag} overlaps centromere region; "
                        f"skipping arm CN correction."
                    )
            for ci in range(len(seed_intvs))[::-1]:
                if seed_intvs[ci].cn < cn_cutoff_chrarm:
                    del seed_intvs[ci]
            if len(seed_intvs) > 0:
                lastseg: Interval | None = None
                sum_seed_len = 0
                for cns in seed_intvs:
                    if lastseg and cns.start - lastseg.end <= max_seg_gap:
                        sum_seed_len += cns.end - cns.start
                        lastseg.end = cns.end
                    elif lastseg is None:
                        lastseg = copy.deepcopy(cns)
                        sum_seed_len += len(cns)
                    elif sum_seed_len >= min_seed_size:
                        fp.write(
                            f"{lastseg.chr}\t{lastseg.start}\t{lastseg.end-1}\n"
                        )
                        sum_seed_len = 0
                        lastseg = cns
                if sum_seed_len >= min_seed_size:
                    fp.write(
                        f"{lastseg.chr}\t{lastseg.start}\t{lastseg.end-1}\n"  # type: ignore
                    )

    if not cnv_seeds:
        print(f"No seed intervals found with CN>={gain}.")
    print("Created " + output_filename)
