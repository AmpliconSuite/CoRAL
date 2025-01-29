"""Implements a class for BreakpointGraph. Will serve as a container for the
sequence edges, corcordant edges, discordant edges, and breakpoints inferred
from a given sequencing file.
"""

from __future__ import annotations

import logging
import warnings
from collections import defaultdict
from dataclasses import dataclass, field

import cvxopt  # type: ignore[import-untyped]
import cvxopt.modeling  # type: ignore[import-untyped]
import numpy as np

from coral import datatypes, types
from coral.breakpoint.breakpoint_utilities import (
    check_valid_discordant_rc_partition,
    enumerate_partitions,
)
from coral.datatypes import (
    AdjacencyMatrix,
    ConcordantEdge,
    Node,
    SequenceEdge,
    Strand,
    WalkData,
)

logger = logging.getLogger(__name__)


@dataclass
class BreakpointGraph:
    """A container object for the breakpoint graphs."""

    amplicon_intervals: list[datatypes.AmpliconInterval] = field(
        default_factory=list
    )
    sequence_edges: list[datatypes.SequenceEdge] = field(default_factory=list)
    concordant_edges: list[datatypes.ConcordantEdge] = field(
        default_factory=list
    )
    discordant_edges: list[datatypes.DiscordantEdge] = field(
        default_factory=list
    )
    source_edges: list[datatypes.SourceEdge] = field(default_factory=list)

    """
	nodes: adjacent list - keys with format Node(chr, pos, orientation);
	vals = [[sequence edges], [concordant edges], [discordant edges], [source edges]]  
	"""
    node_adjacencies: dict[datatypes.Node, AdjacencyMatrix] = field(
        default_factory=lambda: defaultdict(AdjacencyMatrix)
    )
    endnode_adjacencies: defaultdict[datatypes.Node, list[int]] = field(
        default_factory=lambda: defaultdict(list)
    )

    # TODO: add docstring
    path_constraints: list[datatypes.PathConstraint] = field(
        default_factory=list
    )
    longest_path_constraints: list[datatypes.FinalizedPathConstraint] = field(
        default_factory=list
    )

    max_cn: float = 0.0
    amplicon_idx: int = 0

    # Only filled in after valid LP solution
    # TODO: add docstring
    walks: WalkData[datatypes.OptimizationWalk] = field(
        default_factory=lambda: WalkData([], [])
    )
    walk_weights: WalkData[float] = field(
        default_factory=lambda: WalkData([], [])
    )
    path_constraints_satisfied: WalkData[list[int]] = field(
        default_factory=lambda: WalkData([], [])
    )

    @property
    def num_nodes(self) -> int:
        """Number of nodes in the breakpoint graph (excluding s and t)."""
        return len(self.node_adjacencies)

    @property
    def num_edges(self) -> int:
        return (
            self.num_seq_edges
            + self.num_conc_edges
            + self.num_disc_edges
            + 2 * self.num_src_edges
            + 2 * len(self.endnode_adjacencies)
        )

    @property
    def num_seq_edges(self) -> int:
        return len(self.sequence_edges)

    @property
    def num_conc_edges(self) -> int:
        return len(self.concordant_edges)

    @property
    def num_disc_edges(self) -> int:
        return len(self.discordant_edges)

    @property
    def num_src_edges(self) -> int:
        return len(self.source_edges)

    @property
    def num_nonsrc_edges(self) -> int:
        return self.num_seq_edges + self.num_conc_edges + self.num_disc_edges

    def add_endnode(self, node_: datatypes.Node) -> None:
        """Add a new node to the list corresponding to interval ends.

        Args:
                node_: A breakpoint node.

        """
        if node_ not in self.endnode_adjacencies:
            self.endnode_adjacencies[node_] = []
        else:
            warnings.warn("Node corresponding to interval end already exists.")

    def add_sequence_edge(
        self,
        chr: types.ChrTag,
        l: int,
        r: int,
        lr_count: int = -1,
        lr_nc: float = 0.0,
        cn: float = 0.0,
    ) -> None:
        """Add a sequence edge to the graph."""
        node1, node2 = (
            Node(chr, l, Strand.REVERSE),
            Node(chr, r, Strand.FORWARD),
        )
        lseq = len(self.sequence_edges)
        self.node_adjacencies[node1].sequence.append(lseq)
        self.node_adjacencies[node2].sequence.append(lseq)
        self.sequence_edges.append(SequenceEdge(chr, l, r, lr_nc, lr_count, cn))

    def add_concordant_edge(
        self,
        node1: Node,
        node2: Node,
        lr_count: int = -1,
        reads: set[str] | None = None,
        cn: float = 0.0,
    ) -> None:
        """Add a concordant edge to the graph."""
        if (
            node1.chr != node2.chr
            or node2.pos != node1.pos + 1
            or node1.strand != Strand.FORWARD
            or node2.strand != Strand.REVERSE
        ):
            raise Exception("Invalid concordant edge.")
        if reads is None:
            reads = set()
        lc = len(self.concordant_edges)
        self.node_adjacencies[node1].concordant.append(lc)
        self.node_adjacencies[node2].concordant.append(lc)
        self.concordant_edges.append(
            ConcordantEdge(node1, node2, lr_count, reads, cn)
        )

    def add_discordant_edge(
        self,
        node1: Node,
        node2: Node,
        read_support: int | None = None,
        cn: float = 0.0,
        alignments: set[datatypes.BPAlignments] | None = None,
    ) -> None:
        """Add a discordant edge to the breakpoint graph."""

        if not alignments:
            alignments = set()

        ld = len(self.discordant_edges)
        self.node_adjacencies[node1].discordant.append(ld)
        self.node_adjacencies[node2].discordant.append(ld)
        if node1 in self.endnode_adjacencies:
            self.endnode_adjacencies[node1].append(ld)
        if node2 in self.endnode_adjacencies:
            self.endnode_adjacencies[node2].append(ld)
        self.discordant_edges.append(
            datatypes.DiscordantEdge(
                node1,
                node2,
                lr_count=read_support or 0,
                alignments=alignments,
                cn=cn,
            )
        )

    def sort_edges(self) -> None:
        """Sort sequence and concordant edges according to chromosome and
        position, then reset adjacency matrix.
        """
        self.sequence_edges.sort()
        self.concordant_edges.sort()

        for seqi in range(len(self.sequence_edges)):
            sseg = self.sequence_edges[seqi]
            self.node_adjacencies[sseg.start_node].sequence = [seqi]
            self.node_adjacencies[sseg.end_node].sequence = [seqi]
        for ci in range(len(self.concordant_edges)):
            ce = self.concordant_edges[ci]
            self.node_adjacencies[ce.node1].concordant = [ci]
            self.node_adjacencies[ce.node2].concordant = [ci]

    def compute_cn_lr(self, normal_cov_lr: float) -> None:
        """Estimate CN (copy number) for each edge, using only long reads.

        CN estimation is done via maximum likelihood estimation using the joint
        distribution of
            1) observed number of nucleotides on each sequence edge
                - Uses normal distribution, μ = σ^2 = θ_LR * CN(e) * length(e),
                where length is given in base pairs. (S3.1 in Supplementary Methods)
            2) observed read counts on each concordant/discordant edge
                - Uses Poisson distribution, λ = θ_LR * C_e (S3.2)

        Additionally, each node in the breakpoint graph is required to be
        balanced; that is, the CN assignment for each sequence edge (u,v) = the
        sum of CN values from all breakpoint edges connected to that sequence
        edge (S1.1,S3.4).

        These constraints are enforced using convex optimization via cvxpy.

        Args:
            normal_cov_lr: Normal coverage of long reads. (θ_LR)

        Returns:
            None, but assigns all CN values for all breakpoint graph edges.

        """
        lseq = len(self.sequence_edges)
        lc = len(self.concordant_edges)
        ld = len(self.discordant_edges)
        lsrc = len(self.source_edges)
        logger.debug("Adjacent list for estimating CN:")
        for node in self.node_adjacencies:
            logger.debug(
                f"Node {node}; adjacent list = {self.node_adjacencies[node]}."
            )
        nvariables = lseq + lc + ld + lsrc
        logger.debug(f"There are {nvariables} variables for cvxopt.")
        nconstraints = len(
            [
                node
                for node in self.node_adjacencies
                if node not in self.endnode_adjacencies
            ]
        )
        logger.debug(f"There are {nconstraints} constraints for cvxopt.")

        wcn = []
        wlncn = []
        wlrseg = []
        wcn = [0.5 * normal_cov_lr * se.gap for se in self.sequence_edges]
        wcn += [normal_cov_lr for eci in range(lc)]
        wcn += [normal_cov_lr for edi in range(ld)]
        wcn += [0.5 * normal_cov_lr for srci in range(lsrc)]
        wlncn = [-0.5 for seg in self.sequence_edges]
        wlncn += [ce.lr_count * 1.0 for ce in self.concordant_edges]
        wlncn += [de.lr_count * 1.0 for de in self.discordant_edges]
        wlncn += [-0.5 for srci in range(lsrc)]
        wlrseg = [
            (0.5 * se.lr_nc**2 / (normal_cov_lr * se.gap))
            for se in self.sequence_edges
        ]
        wlrseg += [0.0 for eci in range(lc)]
        wlrseg += [0.0 for edi in range(ld)]
        wlrseg += [
            (0.5 * self.source_edges[srci].cn ** 2 / normal_cov_lr)
            for srci in range(lsrc)
        ]
        wcn = cvxopt.matrix(wcn)
        wlncn = cvxopt.matrix(wlncn)
        wlrseg = cvxopt.matrix(wlrseg)

        cidx = 0
        balance_constraints = np.zeros([nconstraints, nvariables])
        for node in self.node_adjacencies:
            if node not in self.endnode_adjacencies:
                for seqi in self.node_adjacencies[node].sequence:
                    balance_constraints[cidx][seqi] = 1
                for eci in self.node_adjacencies[node].concordant:
                    balance_constraints[cidx][lseq + eci] = -1
                for edi in self.node_adjacencies[node].discordant:
                    balance_constraints[cidx][lseq + lc + edi] = -1
                for srci in self.node_adjacencies[node].source:
                    balance_constraints[cidx][lseq + lc + ld + srci] = -1
                cidx += 1
        balance_constraints = cvxopt.matrix(balance_constraints)

        # Convex optimization function required by cvxopt
        def F_normal(x=None, z=None):
            if x is None:
                return 0, cvxopt.matrix(1.0, (nvariables, 1))
            if min(x) <= 0.0:
                return None
            f = (
                cvxopt.modeling.dot(wlrseg, x**-1)
                + cvxopt.modeling.dot(wcn, x)
                - cvxopt.modeling.dot(wlncn, cvxopt.log(x))
            )
            Df = (wcn - cvxopt.mul(wlncn, x**-1) - cvxopt.mul(wlrseg, x**-2)).T
            if z is None:
                return f, Df
            H = cvxopt.spdiag(
                z[0]
                * (cvxopt.mul(wlncn, x**-2) + cvxopt.mul(2.0 * wlrseg, x**-3))
            )
            return f, Df, H

        options = {"maxiters": 1000, "show_progress": False}
        sol = {}
        if nconstraints > 0:
            sol = cvxopt.solvers.cp(
                F_normal,
                A=balance_constraints,
                b=cvxopt.matrix([0.0 for i in range(nconstraints)]),
                kktsolver="ldl",
                options=options,
            )
            if sol["status"] == "optimal" or sol["status"] == "unknown":
                if sol["status"] == "optimal":
                    logger.debug(
                        "Found optimal solution.",
                    )
                else:
                    logger.debug(
                        "Reached maximum num iterations.",
                    )
                logger.debug(
                    f"\tprimal objective = {sol['primal objective']}",
                )
                logger.debug(
                    f"\tdual objective = {sol['dual objective']}",
                )
                logger.debug(
                    f"\tgap = {sol['gap']}",
                )
                logger.debug(
                    f"\trelative gap = {sol['relative gap']}",
                )
                logger.debug(
                    f"\tprimal infeasibility = {sol['primal infeasibility']}",
                )
                logger.debug(
                    f"dual infeasibility = {sol['dual infeasibility']}",
                )
                for seqi in range(lseq):
                    self.sequence_edges[seqi].cn = sol["x"][seqi] * 2
                    self.max_cn = max(sol["x"][seqi] * 2, self.max_cn)
                for ci in range(lc):
                    self.concordant_edges[ci].cn = sol["x"][lseq + ci] * 2
                    self.max_cn = max(sol["x"][lseq + ci] * 2, self.max_cn)
                for di in range(ld):
                    de = self.discordant_edges[di]
                    if de.is_self_loop:
                        self.discordant_edges[di].cn = sol["x"][lseq + lc + di]
                        self.max_cn = max(sol["x"][lseq + lc + di], self.max_cn)
                    else:
                        self.discordant_edges[di].cn = (
                            sol["x"][lseq + lc + di] * 2
                        )
                        self.max_cn = max(
                            sol["x"][lseq + lc + di] * 2, self.max_cn
                        )
                for srci in range(lsrc):
                    self.source_edges[srci].cn = (
                        sol["x"][lseq + lc + ld + srci] * 2
                    )
                    self.max_cn = max(
                        sol["x"][lseq + lc + ld + srci] * 2, self.max_cn
                    )
        else:
            assert lc == 0 and ld == 0 and lsrc == 0
            logger.debug("Skipped convex optimization.")
            for seqi in range(lseq):
                se = self.sequence_edges[seqi]
                cn_seqi = se.lr_nc * 2.0 / (normal_cov_lr * se.gap)
                self.sequence_edges[seqi].cn = cn_seqi
                self.max_cn = max(cn_seqi, self.max_cn)
        self.max_cn += 1.0

    def infer_discordant_edge_multiplicities(
        self, max_multiplicity: int = 5
    ) -> list[int]:
        """Estimate the upper bound of multiplicities for each discordant edge.

        A single discordant edge may be a part of multiple distinct ecDNA
        species, i.e., amplicon walks. We consider multiple supporting reads
        where the CN multiplicity is within the given `max_multiplicity`
        parameter to be representative of the same species, and attempt to
        normalize multiplicities across reads of a discordant edge accordingly.

        Args:
            max_multiplicity: integer, maximum allowed multiplicities in walks

        Return: a list of integers corresponding to the multiplicities of each
            discordant edge
        """
        rc_list = [de.lr_count for de in self.discordant_edges]  # Read counts
        if len(rc_list) == 0:
            return []
        rc_indices = np.argsort(rc_list)
        rc_list = sorted(rc_list)
        if np.log2(rc_list[-1]) - np.log2(rc_list[0]) < 1.0:
            return [1 for i in rc_indices]
        # Minimize clusters, with maximum sum gap
        #   Sum deviations from multiplicities
        num_clusters = 1
        valid_clustering = False
        best_score_all = -10.0
        best_partitions = []
        distinct_all = []

        # Slowly increase # num clusters until we get a valid clustering,
        # since we want to minimize cycles / clusters.
        while not valid_clustering:
            valid_clustering = False
            for partitions in enumerate_partitions(
                num_clusters - 1, 0, len(rc_list) - 1
            ):
                valid_partition = True
                score_all = 0.0
                distinct = []
                for pi in range(len(partitions)):
                    partition = partitions[pi]
                    if (
                        scored_partition := check_valid_discordant_rc_partition(
                            rc_list,
                            partition,
                            max_multiplicity,
                        )
                    ) is None:
                        valid_partition = False
                        break
                    base_ri, score = scored_partition
                    score_all += score
                    distinct.append([partitions[pi][0], base_ri])
                    if pi > 0:
                        score_all += np.log2(
                            rc_list[partitions[pi][0]]
                        ) - np.log2(
                            rc_list[partitions[pi - 1][1]],
                        )
                if valid_partition:
                    valid_clustering = True
                    if score_all > best_score_all:
                        best_score_all = score_all
                        best_partitions = partitions
                        distinct_all = distinct
            if not valid_clustering:
                num_clusters += 1
        multiplicities_sorted = []
        for pi in range(len(best_partitions)):
            partition = best_partitions[pi]
            base_ = distinct_all[pi]
            for i in range(base_[0], base_[1] + 1):
                multiplicities_sorted.append(1)
            base_ri = base_[1] + 1
            if base_ri > partition[1]:
                continue
            base_avg_rc = np.average(rc_list[base_[0] : base_[1] + 1])
            multiplicity = 2
            while rc_list[base_ri] / base_avg_rc >= multiplicity + 0.5:
                multiplicity += 1
            # Note: sometimes int(round()) will convert 1.5 to 1
            # The following procedure works well
            for i in range(base_ri, partition[1] + 1):
                while rc_list[i] / base_avg_rc >= multiplicity + 0.5:
                    multiplicity += 1
                multiplicities_sorted.append(multiplicity)
        return [
            multiplicities_sorted[list(rc_indices).index(i)]
            for i in range(len(rc_list))
        ]
