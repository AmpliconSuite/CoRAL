from __future__ import annotations

from coral import datatypes


def merge_clusters(
    interval_list: list[datatypes.Interval],
    extend: int = 0,
    margin: float = 0.0,
) -> list[tuple[datatypes.Interval, list[datatypes.Interval]]]:
    ml: list[tuple[datatypes.Interval, list[datatypes.Interval]]] = []
    ci = None
    cl = []
    ai = 0
    cend = len(interval_list)
    for a in interval_list[::-1]:
        ai += 1
        if ci is None or not a.intersects(ci, extend, margin):
            cstart = len(interval_list) - ai + 1
            cl = interval_list[cstart:cend]
            if ci is not None:
                ml.append((ci, cl))
            ci = a
            cl = []
        ci = ci.merge(a, extend)
    cstart = 0
    cl = interval_list[cstart:cend]
    if ci is not None:
        ml.append((ci, cl))
    return ml[::-1]
