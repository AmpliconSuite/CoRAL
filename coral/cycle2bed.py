#!/usr/bin/env python3
from __future__ import annotations

from typing import Any

import typer

from coral.constants import CHR_TAG_TO_IDX, INVERT_STRAND_DIRECTION


def convert_cycles_to_bed(
    cycle_file: typer.FileText,
    output_fn: str,
    rotate_to_min: bool = False,
    num_cycles: int | None = None,
):
    """Convert an AA-formatted .txt file into equivalent .bed representation."""
    all_segs: dict[str, list[str | int]] = dict()
    cycles: dict[int, list[Any]] = dict()
    for line in cycle_file:
        t = line.strip().split()
        if t[0] == "Segment":
            all_segs[t[1]] = [t[2], int(t[3]), int(t[4])]
        if t[0][:5] == "Cycle":
            st = t[0].split(";")
            cycle_id = 1
            cycle_weight = 1.0
            cycle_segs = ["0+", "0-"]
            for s in st:
                s_name, s_value = s.split("=")
                if s_name == "Cycle":
                    cycle_id = s_value  # type: ignore[assignment]
                if s_name == "Copy_count":
                    cycle_weight = float(s_value)
                if s_name == "Segments":
                    cycle_segs = s_value.split(",")
            iscyclic = cycle_segs[0] != "0+" or cycle_segs[-1] != "0-"
            cycle: list[list[str | int]] = []
            for seg in cycle_segs:
                segi = seg[:-1]
                segdir = seg[-1]
                if int(segi) > 0:
                    if cycle == []:
                        cycle.append(all_segs[segi] + [segdir])
                    elif (
                        cycle[-1][-1] == "+"
                        and segdir == "+"
                        and cycle[-1][0] == all_segs[segi][0]
                        and cycle[-1][2] + 1 == all_segs[segi][1]  # type: ignore[operator]
                    ):
                        cycle[-1][2] = all_segs[segi][2]
                    elif (
                        cycle[-1][-1] == "-"
                        and segdir == "-"
                        and cycle[-1][0] == all_segs[segi][0]
                        and cycle[-1][1] - 1 == all_segs[segi][2]  # type: ignore[operator]
                    ):
                        cycle[-1][1] = all_segs[segi][1]
                    else:
                        cycle.append(all_segs[segi] + [segdir])
            if (
                cycle[-1][-1] == "+"
                and cycle[0][-1] == "+"
                and cycle[-1][0] == cycle[0][0]
                and cycle[-1][2] + 1 == cycle[0][1]  # type: ignore[operator]
            ):
                cycle[0][1] = cycle[-1][1]
                del cycle[-1]
            if (
                cycle[-1][-1] == "-"
                and cycle[0][-1] == "+"
                and cycle[-1][0] == cycle[0][0]
                and cycle[-1][1] - 1 == cycle[0][2]  # type: ignore[operator]
            ):
                cycle[0][2] = cycle[-1][2]
                del cycle[-1]
            if rotate_to_min and len(cycle) > 1:
                if iscyclic:
                    argmin_idx = cycle.index(
                        min(
                            cycle,
                            key=lambda seg: (CHR_TAG_TO_IDX[seg[0]], seg[1]),  # type: ignore[index]
                        ),
                    )
                    if cycle[argmin_idx][-1] == "+":
                        cycle = cycle[argmin_idx:] + cycle[:argmin_idx]
                    else:
                        cycle = (
                            cycle[: argmin_idx + 1][::-1]
                            + cycle[argmin_idx + 1 :][::-1]
                        )
                        for idx in range(len(cycle)):
                            cycle[idx][-1] = INVERT_STRAND_DIRECTION(
                                cycle[idx][-1]
                            )  # type: ignore[operator]
                elif CHR_TAG_TO_IDX[cycle[-1][0]] < CHR_TAG_TO_IDX[  # type: ignore[index]
                    cycle[0][0]  # type: ignore[index]
                ] or (
                    CHR_TAG_TO_IDX[cycle[-1][0]] == CHR_TAG_TO_IDX[cycle[0][0]]  # type: ignore[index]
                    and cycle[-1][1] < cycle[-1][1]  # type: ignore[index, operator]
                ):
                    cycle = cycle[::-1]
                    if cycle[0][-1] == "-":
                        for idx in range(len(cycle)):
                            cycle[idx][-1] = INVERT_STRAND_DIRECTION(
                                cycle[idx][-1]
                            )  # type: ignore[operator]
            cycles[int(cycle_id)] = [iscyclic, cycle_weight, cycle]

    print("Creating bed-converted cycles file: " + output_fn)
    with open(output_fn, "w") as fp:
        fp.write("#chr\tstart\tend\torientation\tcycle_id\tiscyclic\tweight\n")
        full_num_cycles = len(cycles)
        if num_cycles:
            num_cycles = min(full_num_cycles, num_cycles)
        else:
            num_cycles = full_num_cycles

        for i in range(1, num_cycles + 1):
            for seg in cycles[i][2]:
                fp.write(
                    f"{seg[0]}\t{seg[1]}\t{seg[2]}\t{seg[3]}\t{i}\t{cycles[i][0]}\t{cycles[i][1]}\n"
                )
