"""
A module for performing tests on breakpoint_graph/cigar_parsing.py
"""
import unittest

from long_read_aa.breakpoints import cigar_parsing


class TestCigarParsing(unittest.TestCase):
    def setUp(self):
        # we can define any data or objects that will be used throughout the
        # tests here.
        # If this is  not used, please remove the setUp function.
        pass

    def test_cigar_to_pos_SM(self):
        # test the function cigar2posSM

        # define test input
        test_cigar_string = "S100M"

        query_start, query_end, read_length = cigar_parsing.cigar2posMS(
            test_cigar_string, "+", 100
        )

        # define expected output
        expected_query_start = 0
        expected_query_end = 0
        expected_read_length = 0

        self.assertEqual(query_start, expected_query_start)
        self.assertEqual(query_end, expected_query_end)
        self.assertEqual(read_length, expected_read_length)


if __name__ == "__main__":
    unittest.main()
