#!/usr/bin/env python
"""Unit tests for the pure-Python DNAcopy CBS reconstruction."""

import unittest

import numpy as np

from coral.pycbs import (
    _max_tstat,
    _recursive_breakpoints,
    _sequential_boundary,
    _tail_probability,
)


class PythonCbsTests(unittest.TestCase):
    def test_flat_signal_stays_one_segment(self):
        values = np.zeros(100)
        weights = np.ones(100)
        breaks = _recursive_breakpoints(
            values, weights, 0.01, 100, np.random.default_rng(1)
        )
        self.assertEqual(breaks, [100])

    def test_detects_binary_step(self):
        values = np.concatenate([np.zeros(50), np.ones(50)])
        weights = np.ones(100)
        breaks = _recursive_breakpoints(
            values, weights, 0.01, 100, np.random.default_rng(1)
        )
        self.assertEqual(breaks, [50, 100])

    def test_detects_internal_pulse(self):
        values = np.concatenate([np.zeros(30), np.ones(30), np.zeros(40)])
        weights = np.ones(100)
        breaks = _recursive_breakpoints(
            values, weights, 0.01, 100, np.random.default_rng(1)
        )
        self.assertEqual(breaks, [30, 60, 100])

    def test_weighted_arc_statistic_prefers_real_step(self):
        values = np.concatenate([np.zeros(20), np.ones(20)])
        weights = np.linspace(0.5, 2.0, 40)
        observed = _max_tstat(values, weights)
        self.assertGreater(observed.tstat, 0)
        self.assertEqual(observed.start, 0)
        self.assertEqual(observed.end, 20)

    def test_hybrid_tail_probability_is_finite(self):
        probability = _tail_probability(4.0, 0.05, 300)
        self.assertTrue(np.isfinite(probability))
        self.assertGreaterEqual(probability, 0.0)

    def test_sequential_boundary_has_one_entry_per_rejection_count(self):
        boundary = _sequential_boundary(100, 3)
        self.assertEqual(len(boundary), 4)
        self.assertTrue(all(left < right for left, right in zip(boundary, boundary[1:])))


if __name__ == "__main__":
    unittest.main()
