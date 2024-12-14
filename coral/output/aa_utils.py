from coral import datatypes


def merge_clusters(
    interval_list: list[datatypes.Interval], extend=0, margin=0.0
):
    ml = []
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
    cl = self[cstart:cend]
    if ci is not None:
        ml.append((ci, cl))
    return ml[::-1]
