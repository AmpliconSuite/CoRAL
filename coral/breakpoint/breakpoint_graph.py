"""Implements a class for BreakpointGraph. Will serve as a container for the
sequence edges, corcordant edges, discordant edges, and breakpoints inferred
from a given sequencing file.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from typing import Any

import cvxopt  # type: ignore[import-untyped]
import cvxopt.modeling  # type: ignore[import-untyped]
import numpy as np

from coral import types
from coral.breakpoint.breakpoint_utilities import (
    check_valid_discordant_rc_partition,
    enumerate_partitions,
)
from coral.constants import CHR_TAG_TO_IDX

logger = logging.getLogger(__name__)


@dataclass
class BreakpointGraph:
    """A container object for the breakpoint graphs."""

    amplicon_intervals: list[types.AmpliconInterval] = field(
        default_factory=list
    )
    sequence_edges: list[list[Any]] = field(default_factory=list)
    concordant_edges: list[list[Any]] = field(default_factory=list)
    discordant_edges: list[list[Any]] = field(default_factory=list)
    source_edges: list[list[Any]] = field(default_factory=list)

    """
	nodes: adjacent list - keys with format (chr, pos, orientation); 
	vals = [[sequence edges], [concordant edges], [discordant edges], [source edges]]  
	"""
    nodes: dict[tuple, list[list]] = field(default_factory=dict)
    endnodes: dict[tuple, list[list]] = field(default_factory=dict)
    max_cn: float = 0.0

    @property
    def num_nodes(self) -> int:
        return len(self.nodes)

    @property
    def num_edges(self) -> int:
        return (
            self.num_seq_edges
            + self.num_conc_edges
            + self.num_disc_edges
            + 2 * self.num_src_edges
            + 2 * len(self.endnodes)
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

    def add_node(self, node_):
        """Add a new node to the breakpoint graph.

        Args:
                node_: A breakpoint node.

        """
        if type(node_) != tuple or len(node_) != 3:
            raise Exception(
                "Breakpoint node must be of form (chr, pos, orientation)."
            )
        if node_ in self.nodes:
            pass
        self.nodes[node_] = [[], [], [], []]

    def add_endnode(self, node_):
        """Add a new node to the list corresponding to interval ends.

        Args:
                node_: A breakpoint node.

        """
        if type(node_) != tuple or len(node_) != 3:
            raise Exception(
                "Breakpoint node must be of form (chr, pos, orientation)."
            )
        if node_ not in self.endnodes:
            self.endnodes[node_] = []
        else:
            warnings.warn("Node corresponding to interval end already exists.")

    def del_endnode(self, node_):
        """Delete a node corresponding to interval ends.

        Args:
                node_: A breakpoint node.

        """
        if node_ in self.endnodes:
            del self.endnodes[node_]
        else:
            warnings.warn("Node corresponding to interval end not exists.")

    def del_discordant_endnodes(self):
        """Delete nodes correspond to interval ends and connect to a discordant edges."""
        del_list = []
        for node in self.endnodes:
            if len(self.endnodes[node]) > 0:
                del_list.append(node)
        for node in del_list:
            del self.endnodes[node]

    def add_sequence_edge(
        self, chr, l, r, sr_count=-1, sr_flag="d", lr_count=-1, lr_nc=0, cn=0.0
    ):
        """Add a sequence edge to the graph."""
        if (chr, l, "-") not in self.nodes or (chr, r, "+") not in self.nodes:
            raise Exception("Breakpoint node must be added first.")
        lseq = len(self.sequence_edges)
        self.nodes[(chr, l, "-")][0].append(lseq)
        self.nodes[(chr, r, "+")][0].append(lseq)
        self.sequence_edges.append(
            [chr, l, r, sr_count, sr_flag, lr_count, lr_nc, r - l + 1, cn]
        )

    def add_concordant_edge(
        self,
        chr1,
        pos1,
        o1,
        chr2,
        pos2,
        o2,
        sr_count=-1,
        sr_flag="d",
        lr_count=-1,
        reads=set([]),
        cn=0.0,
    ):
        """Add a concordant edge to the graph."""
        if chr1 != chr2 or pos2 != pos1 + 1 or o1 != "+" or o2 != "-":
            raise Exception("Invalid concordant edge.")
        if (chr1, pos1, o1) not in self.nodes or (
            chr2,
            pos2,
            o2,
        ) not in self.nodes:
            raise Exception("Breakpoint node must be added first.")
        lc = len(self.concordant_edges)
        self.nodes[(chr1, pos1, o1)][1].append(lc)
        self.nodes[(chr2, pos2, o2)][1].append(lc)
        self.concordant_edges.append(
            [
                chr1,
                pos1,
                o1,
                chr2,
                pos2,
                o2,
                sr_count,
                sr_flag,
                lr_count,
                reads,
                cn,
            ],
        )

    def add_discordant_edge(
        self,
        chr1,
        pos1,
        o1,
        chr2,
        pos2,
        o2,
        sr_count=-1,
        sr_flag="d",
        sr_cn=0.0,
        lr_count=-1,
        reads=set([]),
        cn=0.0,
    ):
        """Add a discordant edge to the breakpoint graph."""
        if (chr1, pos1, o1) not in self.nodes or (
            chr2,
            pos2,
            o2,
        ) not in self.nodes:
            raise Exception("Breakpoint node must be added first.")
        ld = len(self.discordant_edges)
        self.nodes[(chr1, pos1, o1)][2].append(ld)
        self.nodes[(chr2, pos2, o2)][2].append(ld)
        if (chr1, pos1, o1) in self.endnodes:
            self.endnodes[(chr1, pos1, o1)].append(ld)
        if (chr2, pos2, o2) in self.endnodes:
            self.endnodes[(chr2, pos2, o2)].append(ld)
        self.discordant_edges.append(
            [
                chr1,
                pos1,
                o1,
                chr2,
                pos2,
                o2,
                sr_count,
                sr_flag,
                sr_cn,
                lr_count,
                reads,
                cn,
            ],
        )

    def del_discordant_edges(self, del_list, bpi_map):
        """Delete a list discordant edges from the graph."""
        sorted_del_list = sorted(del_list, reverse=True)
        for bpi in sorted_del_list:
            del self.discordant_edges[bpi]
        for node in self.endnodes.keys():
            for i in range(len(self.endnodes[node])):
                if self.endnodes[node][i] in sorted_del_list:
                    del self.endnodes[node][i]
                else:
                    self.endnodes[node][i] = bpi_map[self.endnodes[node][i]]
        for node in self.nodes.keys():
            for i in range(len(self.nodes[node][2])):
                if self.nodes[node][2][i] in sorted_del_list:
                    del self.nodes[node][2][i]
                else:
                    self.nodes[node][2][i] = bpi_map[self.nodes[node][2][i]]

    def add_source_edge(
        self,
        chr1,
        pos1,
        o1,
        sr_count=0,
        sr_flag="d",
        sr_cn=0.0,
        lr_cn=0.0,
        cn=0.0,
    ):
        """Adds a source edge to the graph."""
        if (chr1, pos1, o1) not in self.nodes:
            raise Exception("Breakpoint node must be added first.")
        self.nodes[(chr1, pos1, o1)][3].append(len(self.source_edges))
        self.source_edges.append(
            [
                "source",
                -1,
                "-",
                chr1,
                pos1,
                o1,
                sr_count,
                sr_flag,
                sr_cn,
                lr_cn,
                cn,
            ],
        )

    def del_source_edges(self, del_list, srci_map):
        """Delete a list source edges from the graph."""
        sorted_del_list = sorted(del_list, reverse=True)
        for srci in sorted_del_list:
            del self.source_edges[srci]
        for node in self.nodes.keys():
            for i in range(len(self.nodes[node][3])):
                if self.nodes[node][3][i] in sorted_del_list:
                    del self.nodes[node][3][i]
                else:
                    self.nodes[node][3][i] = srci_map[self.nodes[node][3][i]]

    def del_redundant_sequence_edges(self):
        """Delete redundant sequence edges after merging."""
        if len(self.discordant_edges) == 0:
            return
        del_list = []
        for seqi in range(len(self.sequence_edges)):
            sseg = self.sequence_edges[seqi]
            node1 = (sseg[0], sseg[1], "-")
            node2 = (sseg[0], sseg[2], "+")
            s1 = (
                len(self.nodes[node1][1])
                + len(self.nodes[node1][2])
                + len(self.nodes[node1][3])
            )
            s2 = (
                len(self.nodes[node2][1])
                + len(self.nodes[node2][2])
                + len(self.nodes[node2][3])
            )
            if s1 + s2 == 0:
                del_list.append(seqi)
        for seqi in del_list[::-1]:
            ai = self.sequence_edges[seqi][:3]
            if ai in self.amplicon_intervals:
                del self.amplicon_intervals[self.amplicon_intervals.index(ai)]
            node1 = (ai[0], ai[1], "-")
            node2 = (ai[0], ai[2], "+")
            del self.sequence_edges[seqi]
            del self.nodes[node1]
            del self.nodes[node2]
            self.del_endnode(node1)
            self.del_endnode(node2)
        for seqi in range(len(self.self.sequence_edges)):
            sseg = self.sequence_edges[seqi]
            node1 = (sseg[0], sseg[1], "-")
            node2 = (sseg[0], sseg[2], "+")
            self.nodes[node1][0][0] = seqi
            self.nodes[node2][0][0] = seqi

    def merge_edges(self):
        """Merge sequence edges connected only by concordant edges;
        Delete the nodes and concordant edges accordingly.
        """
        c_del_list = []  # The list of concordant edges to be deleted
        seq_del_list = []  # The list of sequence edges to be deleted
        for ci in range(len(self.concordant_edges)):
            ce = self.concordant_edges[ci]
            node1 = (ce[0], ce[1], ce[2])
            node2 = (ce[3], ce[4], ce[5])
            if (
                len(self.nodes[node1][2]) == 0
                and len(self.nodes[node2][2]) == 0
                and len(self.nodes[node1][3]) == 0
                and len(self.nodes[node2][3]) == 0
            ):
                seqi1 = self.nodes[node1][0][0]
                seqi2 = self.nodes[node2][0][0]
                seq_del_list.append(seqi1)
                del self.nodes[node1]
                del self.nodes[node2]
                c_del_list.append(ci)
        seq_del_list = sorted(seq_del_list)
        si = 0
        li = 0
        for i in range(1, len(seq_del_list)):
            if seq_del_list[i] == seq_del_list[li] + 1:
                li += 1
            else:
                seqi1 = seq_del_list[si]
                seqi2 = seq_del_list[li] + 1
                self.sequence_edges[seqi2][1] = self.sequence_edges[seqi1][1]
                self.sequence_edges[seqi2][3] = -1
                self.sequence_edges[seqi2][4] = "f"
                self.sequence_edges[seqi2][-2] = (
                    self.sequence_edges[seqi2][2]
                    - self.sequence_edges[seqi2][1]
                    + 1
                )
                si = i
                li = i
        seqi1 = seq_del_list[si]
        seqi2 = seq_del_list[li] + 1
        self.sequence_edges[seqi2][1] = self.sequence_edges[seqi1][1]
        self.sequence_edges[seqi2][3] = -1
        self.sequence_edges[seqi2][4] = "f"
        self.sequence_edges[seqi2][-2] = (
            self.sequence_edges[seqi2][2] - self.sequence_edges[seqi2][1] + 1
        )
        for seqi in seq_del_list[::-1]:
            del self.sequence_edges[seqi]
        for ci in sorted(c_del_list, reverse=True):
            del self.concordant_edges[ci]
        for seqi in range(len(self.sequence_edges)):
            sseg = self.sequence_edges[seqi]
            node1 = (sseg[0], sseg[1], "-")
            node2 = (sseg[0], sseg[2], "+")
            self.nodes[node1][0][0] = seqi
            self.nodes[node2][0][0] = seqi
        for ci in range(len(self.concordant_edges)):
            ce = self.concordant_edges[ci]
            node1 = (ce[0], ce[1], ce[2])
            node2 = (ce[3], ce[4], ce[5])
            self.nodes[node1][1][0] = ci
            self.nodes[node2][1][0] = ci

    def sort_edges(self):
        """Sort sequence and concordant edges according to chromosome and position
        Reset adjacent list
        """
        self.sequence_edges.sort(
            key=lambda sseg: (CHR_TAG_TO_IDX[sseg[0]], sseg[1])
        )
        self.concordant_edges.sort(
            key=lambda ce: (CHR_TAG_TO_IDX[ce[0]], ce[1])
        )

        for seqi in range(len(self.sequence_edges)):
            sseg = self.sequence_edges[seqi]
            self.nodes[(sseg[0], sseg[1], "-")][0] = [seqi]
            self.nodes[(sseg[0], sseg[2], "+")][0] = [seqi]
        for ci in range(len(self.concordant_edges)):
            ce = self.concordant_edges[ci]
            self.nodes[(ce[0], ce[1], ce[2])][1] = [ci]
            self.nodes[(ce[3], ce[4], ce[5])][1] = [ci]

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
        for node in self.nodes:
            logger.debug(f"Node {node}; adjacent list = {self.nodes[node]}.")
        nvariables = lseq + lc + ld + lsrc
        logger.debug(f"There are {nvariables} variables for cvxopt.")
        nconstraints = len(
            [node for node in self.nodes if node not in self.endnodes]
        )
        logger.debug(f"There are {nconstraints} constraints for cvxopt.")

        wcn = []
        wlncn = []
        wlrseg = []
        wcn = [0.5 * normal_cov_lr * se[7] for se in self.sequence_edges]
        wcn += [normal_cov_lr for eci in range(lc)]
        wcn += [normal_cov_lr for edi in range(ld)]
        wcn += [0.5 * normal_cov_lr for srci in range(lsrc)]
        wlncn = [-0.5 for seg in self.sequence_edges]
        wlncn += [ce[8] * 1.0 for ce in self.concordant_edges]
        wlncn += [de[9] * 1.0 for de in self.discordant_edges]
        wlncn += [-0.5 for srci in range(lsrc)]
        wlrseg = [
            (0.5 * se[6] ** 2 / (normal_cov_lr * se[7]))
            for se in self.sequence_edges
        ]
        wlrseg += [0.0 for eci in range(lc)]
        wlrseg += [0.0 for edi in range(ld)]
        wlrseg += [
            (0.5 * self.source_edges[srci][-1] ** 2 / normal_cov_lr)
            for srci in range(lsrc)
        ]
        wcn = cvxopt.matrix(wcn)
        wlncn = cvxopt.matrix(wlncn)
        wlrseg = cvxopt.matrix(wlrseg)

        cidx = 0
        balance_constraints = np.zeros([nconstraints, nvariables])
        for node in self.nodes:
            if node not in self.endnodes:
                for seqi in self.nodes[node][0]:
                    balance_constraints[cidx][seqi] = 1
                for eci in self.nodes[node][1]:
                    balance_constraints[cidx][lseq + eci] = -1
                for edi in self.nodes[node][2]:
                    balance_constraints[cidx][lseq + lc + edi] = -1
                for srci in self.nodes[node][3]:
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
                    "\tprimal objective = %f" % (sol["primal objective"]),
                )
                logger.debug(
                    "\tdual objective = %f" % (sol["dual objective"]),
                )
                logger.debug(
                    "\tgap = %f" % (sol["gap"]),
                )
                logger.debug(
                    "\trelative gap = %f" % (sol["relative gap"]),
                )
                logger.debug(
                    "\tprimal infeasibility = %f"
                    % (sol["primal infeasibility"]),
                )
                logger.debug(
                    "dual infeasibility = %f" % (sol["dual infeasibility"]),
                )
                for seqi in range(lseq):
                    self.sequence_edges[seqi][-1] = sol["x"][seqi] * 2
                    self.max_cn = max(sol["x"][seqi] * 2, self.max_cn)
                for ci in range(lc):
                    self.concordant_edges[ci][-1] = sol["x"][lseq + ci] * 2
                    self.max_cn = max(sol["x"][lseq + ci] * 2, self.max_cn)
                for di in range(ld):
                    de = self.discordant_edges[di]
                    if de[0] == de[3] and de[1] == de[4] and de[2] == de[5]:
                        self.discordant_edges[di][-1] = sol["x"][lseq + lc + di]
                        self.max_cn = max(sol["x"][lseq + lc + di], self.max_cn)
                    else:
                        self.discordant_edges[di][-1] = (
                            sol["x"][lseq + lc + di] * 2
                        )
                        self.max_cn = max(
                            sol["x"][lseq + lc + di] * 2, self.max_cn
                        )
                for srci in range(lsrc):
                    self.source_edges[srci][-1] = (
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
                cn_seqi = se[6] * 2.0 / (normal_cov_lr * se[7])
                self.sequence_edges[seqi][-1] = cn_seqi
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
        rc_list = [de[9] for de in self.discordant_edges]  # Read counts
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
