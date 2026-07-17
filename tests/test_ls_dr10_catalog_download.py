"""Offline tests for LS DR10 query geometry."""

import unittest

from ls_dr10_catalog_download import build_query, compute_box_bounds, ra_degrees


class CoordinateGeometryTests(unittest.TestCase):
    def test_ra_is_normalized(self):
        self.assertEqual(ra_degrees("361.25"), 1.25)
        self.assertEqual(ra_degrees("-1"), 359.0)

    def test_query_wraps_across_zero_ra(self):
        ra_min, ra_max, dec_min, dec_max = compute_box_bounds(359.99, 0.0, 0.02)
        query = build_query(ra_min, ra_max, dec_min, dec_max, ["ra", "dec"])
        self.assertIn("ra >=", query)
        self.assertIn("OR ra <=", query)

    def test_declination_bounds_are_clipped_at_pole(self):
        _, _, dec_min, dec_max = compute_box_bounds(0.0, 89.99, 1.0)
        self.assertEqual(dec_max, 90.0)
        self.assertGreaterEqual(dec_min, -90.0)


if __name__ == "__main__":
    unittest.main()
