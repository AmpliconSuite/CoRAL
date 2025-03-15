from __future__ import annotations

from typing import Protocol, TypedDict

import pyomo.environ as pyo
import pyomo.opt.results.results_ as pyo_results


# Add type-hinting since Pyomo doesn't have type-hinting
class PyomoSolverResults(Protocol):
    status: pyo.SolverStatus
    termination_condition: pyo.TerminationCondition
    solution: pyo_results.SolverResults
    gap: float | None


class PyomoResults(Protocol):
    solver: PyomoSolverResults
    solution: list[pyo_results.SolverResults]


class PyomoDirectGurobiOptions(TypedDict):
    threads: int
    NonConvex: int
    timelimit: int
