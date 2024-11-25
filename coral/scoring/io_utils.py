"""
Utilities for reading and processing breakpoint graphs and cycles file from
AmpliconArchitect.
"""

import re

import numpy as np
import pandas as pd
import pyranges


def read_cycles_intervals_to_bed(file):
    intervals = []
    segments = []
    cycles = []

    with open(file, "r") as f:
        for line in f:
            if line.startswith("Interval"):
                _, _, chrom, start, end = line.split("\t")
                intervals.append((chrom, start, end))

            if line.startswith("Segment"):
                _, _id, chrom, start, end = line.split("\t")
                segments.append((_id, chrom, start, end))

    intervals_pr = pyranges.PyRanges(
        chromosomes=[interval[0] for interval in intervals],
        starts=[interval[1] for interval in intervals],
        ends=[interval[2] for interval in intervals],
    )
    segments_pr = pyranges.PyRanges(
        chromosomes=[segment[1] for segment in segments],
        starts=[segment[2] for segment in segments],
        ends=[segment[3] for segment in segments],
    )
    return intervals_pr, segments_pr


def read_cycles_file_to_bed(input_file, num_cycles=None):
    all_segs = dict()
    cycles = dict()
    with open(input_file) as fp:
        for line in fp:
            t = line.strip().split()
            if t[0] == "Segment":
                all_segs[t[1]] = [t[2], int(t[3]), int(t[4])]
            if t[0][:5] == "Cycle":
                st = t[0].split(";")
                cycle_id = 1
                cycle_weight = 1.0
                cycle_segs = ["0+", "0-"]
                for s in st:
                    s = s.split("=")
                    if s[0] == "Cycle":
                        cycle_id = s[1]
                    if s[0] == "Copy_count":
                        cycle_weight = float(s[1])
                    if s[0] == "Segments":
                        cycle_segs = s[1].split(",")

                if len(cycle_segs) == 1:
                    segs = re.split("(\\+|\\-)", cycle_segs[0])
                    new_segs = []
                    for i in range(len(segs) // 2):
                        new_segs.append(segs[2 * i] + segs[2 * i + 1])
                    cycle_segs = new_segs
                iscyclic = cycle_segs[0] != "0+" or cycle_segs[-1] != "0-"
                cycle = []
                for seg in cycle_segs:
                    segi = seg[:-1]
                    segdir = seg[-1]
                    if int(segi) > 0:
                        if cycle == []:
                            cycle.append(all_segs[segi] + [segdir])
                        else:
                            if (
                                cycle[-1][-1] == "+"
                                and segdir == "+"
                                and cycle[-1][0] == all_segs[segi][0]
                                and cycle[-1][2] + 1 == all_segs[segi][1]
                            ):
                                cycle[-1][2] = all_segs[segi][2]
                            elif (
                                cycle[-1][-1] == "-"
                                and segdir == "-"
                                and cycle[-1][0] == all_segs[segi][0]
                                and cycle[-1][1] - 1 == all_segs[segi][2]
                            ):
                                cycle[-1][1] = all_segs[segi][1]
                            else:
                                cycle.append(all_segs[segi] + [segdir])
                if (
                    cycle[-1][-1] == "+"
                    and cycle[0][-1] == "+"
                    and cycle[-1][0] == cycle[0][0]
                    and cycle[-1][2] + 1 == cycle[0][1]
                ):
                    cycle[0][1] = cycle[-1][1]
                    del cycle[-1]
                if (
                    cycle[-1][-1] == "-"
                    and cycle[0][-1] == "+"
                    and cycle[-1][0] == cycle[0][0]
                    and cycle[-1][1] - 1 == cycle[0][2]
                ):
                    cycle[0][2] = cycle[-1][2]
                    del cycle[-1]
                cycles[int(cycle_id)] = [iscyclic, cycle_weight, cycle]

    bed_df = pd.DataFrame(
        columns=[
            "Chromosome",
            "Start",
            "End",
            "orientation",
            "cycle_id",
            "iscyclic",
            "weight",
        ]
    )
    _num_cycles = len(cycles) if num_cycles is None else num_cycles
    for i in range(1, _num_cycles + 1):
        for seg in cycles[i][2]:
            new_row = pd.DataFrame(
                [
                    [
                        seg[0],
                        seg[1],
                        seg[2],
                        seg[3],
                        i,
                        cycles[i][0],
                        cycles[i][1],
                    ]
                ],
                columns=bed_df.columns,
            )

            bed_df = pd.concat([bed_df, new_row])

    return pyranges.PyRanges(bed_df)


def read_breakpoint_graph(file):
    sequence_edges = []
    breakpoint_edges = []

    with open(file, "r") as f:
        for line in f:
            if line.startswith("sequence"):
                _, left, right, copycount, coverage, size = line.split("\t")[:6]

                left_chrom, left_pos = (
                    left.split(":")[0],
                    int(left.split(":")[1].split("-")[0]),
                )
                right_chrom, right_pos = (
                    right.split(":")[0],
                    int(right.split(":")[1].split("+")[0]),
                )

                sequence_edges.append(
                    [
                        left_chrom,
                        left_pos,
                        right_chrom,
                        right_pos,
                        copycount,
                        coverage,
                        size.strip(),
                    ]
                )

            if line.startswith("discordant"):
                values = line.split("\t")[:5]
                values = [val for val in values if val != ""]

                edge_class, breakpoint_edge, copycount, read_pairs = values[:4]

                left_breakpoint, right_breakpoint = breakpoint_edge.split("->")

                left_chrom, left_pos, left_orient = (
                    left_breakpoint.split(":")[0],
                    int(re.split(r"[+-]", left_breakpoint.split(":")[1])[0]),
                    re.split(r":\d+", left_breakpoint)[1],
                )
                right_chrom, right_pos, right_orient = (
                    right_breakpoint.split(":")[0],
                    int(re.split(r"[+-]", right_breakpoint.split(":")[1])[0]),
                    re.split(r":\d+", right_breakpoint)[1],
                )

                if left_chrom > right_chrom:
                    breakpoint_edges.append(
                        [
                            edge_class,
                            right_chrom,
                            right_pos,
                            right_orient,
                            left_chrom,
                            left_pos,
                            left_orient,
                            copycount,
                            read_pairs.strip(),
                        ]
                    )
                elif left_chrom == right_chrom and left_pos > right_pos:
                    breakpoint_edges.append(
                        [
                            edge_class,
                            right_chrom,
                            right_pos,
                            right_orient,
                            left_chrom,
                            left_pos,
                            left_orient,
                            copycount,
                            read_pairs.strip(),
                        ]
                    )
                else:
                    breakpoint_edges.append(
                        [
                            edge_class,
                            left_chrom,
                            left_pos,
                            left_orient,
                            right_chrom,
                            right_pos,
                            right_orient,
                            copycount,
                            read_pairs.strip(),
                        ]
                    )

    sequence_edges = pd.DataFrame(
        sequence_edges,
        columns=[
            "Chrom1",
            "Start",
            "Chrom2",
            "End",
            "CopyCount",
            "Coverage",
            "Size",
        ],
    )
    breakpoint_edges = pd.DataFrame(
        breakpoint_edges,
        columns=[
            "Class",
            "Chrom1",
            "Pos1",
            "Orient1",
            "Chrom2",
            "Pos2",
            "Orient2",
            "CopyCount",
            "ReadPairs",
        ],
    )

    return sequence_edges, breakpoint_edges


def bin_genome(true_segments, reconstructed_segments, margin_size=10000):
    """
    Bin the intervals by the breakpoints union.
    Warning: Setting a margin size can effect the output of different distances
    """

    df_bins = pd.DataFrame(
        np.concatenate(
            (
                true_segments[["Chromosome", "Start"]].values,
                true_segments[["Chromosome", "End"]].values,
                reconstructed_segments[["Chromosome", "Start"]].values,
                reconstructed_segments[["Chromosome", "End"]].values,
            ),
            axis=0,
        )
    )

    df_bins.columns = ["Chromosome", "Start"]
    # df_bins['Start'] = df_bins['Start']
    df_bins = df_bins.drop_duplicates().sort_values(by=["Chromosome", "Start"])

    # rotate with 1 up the start column
    df_bins_suffix = df_bins.tail(-1)
    df_bins_suffix = pd.concat([df_bins_suffix, df_bins.head(1)])

    df_bins.reset_index(drop=True, inplace=True)
    df_bins_suffix.reset_index(drop=True, inplace=True)
    df_bins = pd.concat([df_bins, df_bins_suffix], axis=1, ignore_index=True)
    df_bins.columns = ["Chromosome", "Start", "Chromosome2", "End"]

    # keep only rows with same chr and non-negative distance
    df_bins = df_bins[
        (df_bins["Chromosome"] == df_bins["Chromosome2"])
        & (abs(df_bins["Start"] - df_bins["End"]) >= 0)
    ]

    df_bins = df_bins[(df_bins["End"] - df_bins["Start"]) > 0]
    df_bins["length"] = abs(df_bins["End"] - df_bins["Start"])
    chrlist = df_bins["Chromosome"].drop_duplicates().tolist()

    return df_bins[["Chromosome", "Start", "End", "length"]], chrlist


def get_seq_edges_from_decoil(results, base_coverage=13.0):
    sequence_edges = []
    for _, entry in results.iterrows():
        sequence_edges.append(
            [
                entry.Chromosome,
                entry.Start,
                entry.Chromosome,
                entry.End,
                (entry.End - entry.Start),
                entry.weight,
            ]
        )
        # sequence_edges.append([entry.Chromosome, entry.Start, entry.Chromosome, entry.End, (entry.End-entry.Start), entry.coverage/base_coverage])

    breakpoint_edges = []
    flip_orient = lambda x: "-" if x == "+" else "+"

    for _, cycle_df in results.groupby("cycle_id"):
        for i in range(len(cycle_df)):
            entry_i = cycle_df.iloc[i, :]

            if i == (len(cycle_df) - 1):
                j = 0
            else:
                j = i + 1

            entry_j = cycle_df.iloc[j, :]

            left_chrom, left_start, left_end, left_orient = entry_i[
                ["Chromosome", "Start", "End", "orientation"]
            ]
            right_chrom, right_start, right_end, right_orient = entry_j[
                ["Chromosome", "Start", "End", "orientation"]
            ]

            right_pos = right_start if right_orient == "+" else right_end
            left_pos = left_end if left_orient == "+" else left_start

            if left_chrom > right_chrom:
                breakpoint_edges.append(
                    [
                        "discordant",
                        right_chrom,
                        right_pos,
                        flip_orient(right_orient),
                        left_chrom,
                        left_pos,
                        left_orient,
                    ]
                )
            elif left_chrom == right_chrom and left_pos > right_pos:
                breakpoint_edges.append(
                    [
                        "discordant",
                        right_chrom,
                        right_pos,
                        flip_orient(right_orient),
                        left_chrom,
                        left_pos,
                        left_orient,
                    ]
                )
            else:
                breakpoint_edges.append(
                    [
                        "discordant",
                        left_chrom,
                        left_pos,
                        left_orient,
                        right_chrom,
                        right_pos,
                        flip_orient(right_orient),
                    ]
                )

    breakpoint_edges = pd.DataFrame(
        breakpoint_edges,
        columns=[
            "Class",
            "Chrom1",
            "Pos1",
            "Orient1",
            "Chrom2",
            "Pos2",
            "Orient2",
        ],
    )

    sequence_edges = pd.DataFrame(
        sequence_edges,
        columns=[
            "Chrom1",
            "Start",
            "Chrom2",
            "End",
            "Size",
            "CopyCount",
        ],
    )

    return breakpoint_edges, sequence_edges
