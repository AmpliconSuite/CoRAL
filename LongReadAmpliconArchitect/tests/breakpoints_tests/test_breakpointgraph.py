"""
Tests for BreakpointGraph class.
"""
import unittest

from long_read_aa.breakpoints import (
    BreakpointEdge,
    BreakpointGraph,
    SequenceEdge,
)


class TestBreakpointGraph(unittest.TestCase):
    def setUp(self):
        self.breakpoint1 = BreakpointEdge.BreakpointEdge(
            "chr1", 2000, "+", "chr2", 4000, "-", annotation="discordant"
        )
        self.breakpoint2 = BreakpointEdge.BreakpointEdge(
            "chr10", 1201, "-", "chr10", 20000, "+", annotation="concordant"
        )
        self.sequence_edge = SequenceEdge.SequenceEdge("chr1", 10000, 120000)

    def test_breakpoint_edge_constructor(self):
        """Test BreakpointEdge constructor."""

        self.assertEqual(self.breakpoint1.chromosome1, "chr1")
        self.assertEqual(self.breakpoint1.position1, 2000)
        self.assertEqual(self.breakpoint1.orientation1, "+")
        self.assertEqual(self.breakpoint1.chromosome2, "chr2")
        self.assertEqual(self.breakpoint1.position2, 4000)
        self.assertEqual(self.breakpoint1.orientation2, "-")
        self.assertEqual(self.breakpoint1.support, None)
        self.assertEqual(self.breakpoint1.copy_number, None)
        self.assertEqual(self.breakpoint1.annotation, "discordant")

        self.assertEqual(self.breakpoint1.left_breakpoint, ("chr1", 2000, "+"))
        self.assertEqual(self.breakpoint1.right_breakpoint, ("chr2", 4000, "-"))

        # test error catches
        self.assertRaises(
                Exception,
                BreakpointEdge.BreakpointEdge,
                "chr1", 1, "+", "chr2", 2, "-", None, -1, 0
            )
        self.assertRaises(
                Exception,
                BreakpointEdge.BreakpointEdge,
                "chr1", 1, "+", "chr2", 2, "-", None, 0, -1
            )

        # test print
        self.assertEqual(str(self.breakpoint1), "('chr1', 2000, '+') -> ('chr2', 4000, '-')")
    
    def test_sequence_edge_constructor(self):
        """Test SequenceEdge constructor."""

        self.assertEqual(self.sequence_edge.chromosome, "chr1")
        self.assertEqual(self.sequence_edge.start, 10000)
        self.assertEqual(self.sequence_edge.end, 120000)
        self.assertEqual(self.sequence_edge.length, 120000-10000 + 1)

        # test error statements
        self.assertRaises(
            Exception,
            SequenceEdge.SequenceEdge,
            "chr1", 1000, 500
        )
        self.assertRaises(
            Exception,
            SequenceEdge.SequenceEdge,
            "chr1", 1000, 2000, {'long_read': -1}, {'long_read': 0}, 100
        )
        self.assertRaises(
            Exception,
            SequenceEdge.SequenceEdge,
            "chr1", 1000, 2000, {'long_read': 0}, {'long_read': -1}, 100
        )

        self.assertEqual(str(self.sequence_edge), "chr1:10000-120000")


    def test_breakpoint_graph_constructor(self):
        """Test BreakpointGraph constructor."""
        breakpoint_graph = BreakpointGraph.BreakpointGraph(
            sequence_edges=[self.sequence_edge],
            discordant_edges=[self.breakpoint1],
            concordant_edges=[self.breakpoint2],
        )

        self.assertEqual(breakpoint_graph.discordant_edges, [self.breakpoint1])
        self.assertEqual(breakpoint_graph.concordant_edges, [self.breakpoint2])
        self.assertEqual(len(breakpoint_graph.source_edges), 0)
        self.assertEqual(breakpoint_graph.sequence_edges, [self.sequence_edge])

    def test_breakpoint_graph_adjacency_list(self):
        
        breakpoint_graph = BreakpointGraph.BreakpointGraph(
            sequence_edges=[self.sequence_edge],
            discordant_edges=[self.breakpoint1],
            concordant_edges=[self.breakpoint2],
        )

        breakpoint_graph[("chr1", 12232)] = [self.breakpoint1]

        expected_breakpoint_left = self.breakpoint1.left_breakpoint
        expected_breakpoint_right = self.breakpoint1.right_breakpoint

        self.assertEqual(len(breakpoint_graph[("chr1", 12232)]), 1)

        observed_breakpoint = breakpoint_graph[('chr1', 12232)][0]

        self.assertEqual(observed_breakpoint.left_breakpoint, expected_breakpoint_left)
        self.assertEqual(observed_breakpoint.right_breakpoint, expected_breakpoint_right)

        # test errors
        self.assertRaises(
            Exception,
            breakpoint_graph.__setitem__,
            ("chr1", 1111),
            self.breakpoint2
        )

        self.assertRaises(
            Exception,
            breakpoint_graph.__setitem__,
            'chr1'
        )

        self.assertIsNone(breakpoint_graph[('chr1', 1111)])

        self.assertRaises(
            Exception,
            breakpoint_graph.__getitem__,
            1111
        )

        # test delete
        del breakpoint_graph[('chr1', 12232)]
        self.assertFalse(('chr1', 12232) in breakpoint_graph)
    

if __name__ == "__main__":
    unittest.main()
