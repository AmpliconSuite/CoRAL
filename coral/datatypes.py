import enum
import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, NamedTuple, Optional, Set, cast

import pyomo.core
import pyomo.environ as pyo
import pyomo.opt
import pyomo.util.infeasible
from coral.breakpoint.breakpoint_graph import BreakpointGraph


@dataclass
class EdgeToCN:
    sequence: Dict[int, float] = field(default_factory=dict)
    concordant: Dict[int, float] = field(default_factory=dict)
    discordant: Dict[int, float] = field(default_factory=dict)
    source: Dict[int, float] = field(default_factory=dict)

    @staticmethod
    def from_graph(bp_graph: BreakpointGraph):
        return EdgeToCN(
            sequence={i: edge[-1] for i, edge in enumerate(bp_graph.sequence_edges)},
            concordant={i: edge[-1] for i, edge in enumerate(bp_graph.concordant_edges)},
            discordant={i: edge[-1] for i, edge in enumerate(bp_graph.discordant_edges)},
            source={i: edge[-1] for i, edge in enumerate(bp_graph.source_edges)},
        )


@dataclass
class ParsedLPSolution:
    solver_status: str
    termination_condition: pyo.TerminationCondition
    total_weights_included: float = 0.0
    cycles: List[List[Any]] = field(default_factory=lambda: [[], []])  # cycles, paths
    cycle_weights: List[List[Any]] = field(default_factory=lambda: [[], []])  # cycles, paths
    path_constraints_satisfied: List[List[Any]] = field(default_factory=lambda: [[], []])  # cycles, paths
    path_constraints_satisfied_set: Set[int] = field(default_factory=set)


class InitialSolution(NamedTuple):
    cycles: List
    cycle_weights: List
    satisfied_path_constraints: List


class Solver(str, enum.Enum):
    GUROBI = "gurobi"
    SCIP = "scip"
    # Coin-OR solvers
    BONMIN = "bonmin"
    COUENNE = "couenne"


class ReferenceGenome(str, enum.Enum):
    hg19 = "hg19"
    hg38 = "hg38"
    GRCh38 = "GRCh38"
    mm10 = "mm10"
    GRCh37 = "GRCh37"
