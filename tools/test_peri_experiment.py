import math
import unittest

from peri_experiment import duplicate_groups, read_sheet_rows, sample_metadata


def sheet_rows(sheet='Micro-H', label='0TO-3SA-H', points=None):
    units = ('\u03bcm\t', 'mN') if sheet.startswith('Micro') else ('nm', '\u00b5N')
    return [(label, None), ('Displacement', '\tForce'), units] + list(points or [])


class ExperimentalSourceTests(unittest.TestCase):
    def test_unloading_order_repeated_displacements_and_origin_are_preserved(self):
        rows = sheet_rows(points=[(0, 0.2), (10, 2), (10, 3), (4, 1), (None, None), (-1, 0.1)])
        curve = read_sheet_rows('Micro-H', rows)[0]
        self.assertEqual(curve.points, [(4, 0., .2), (5, 10., 2.), (6, 10., 3.), (7, 4., 1.), (9, -1., .1)])

    def test_si_conversion_does_not_mix_nano_and_micro_scales(self):
        micro = read_sheet_rows('Micro-H', sheet_rows(points=[(1, 1)]))[0]
        nano = read_sheet_rows('Nano-H', sheet_rows('Nano-H', points=[(1, 1)]))[0]
        self.assertNotEqual(micro.point_digest(), nano.point_digest())
        self.assertEqual(micro.displacement_unit, 'um')
        self.assertEqual(nano.force_unit, 'uN')

    def test_duplicate_labels_are_retained_but_not_independent_curves(self):
        a = read_sheet_rows('Micro-H', sheet_rows(label='5TO-3SA-H', points=[(0, 0), (2, 3)]))[0]
        b = read_sheet_rows('Micro-H', sheet_rows(label='5TO-5SA-H', points=[(0, 0), (2, 3)]))[0]
        self.assertEqual(duplicate_groups([a, b]), [['Micro-H_5TO-3SA-H', 'Micro-H_5TO-5SA-H']])
        self.assertEqual(len(a.points), 2)
        self.assertEqual(len(b.points), 2)

    def test_equal_points_in_different_orders_are_not_same_loading_history(self):
        a = read_sheet_rows('Micro-H', sheet_rows(points=[(0, 0), (1, 2), (0, 1)]))[0]
        b = read_sheet_rows('Micro-H', sheet_rows(points=[(0, 0), (0, 1), (1, 2)]))[0]
        self.assertEqual(duplicate_groups([a, b]), [])

    def test_unknown_time_geometry_and_contact_zero_stay_unknown(self):
        data = sample_metadata('0TO-3SA-H', 'Micro-H')
        self.assertFalse(data['contains_tung_oil_by_formulation'])
        self.assertIsNone(data['specimen_diameter_m'])
        self.assertIsNone(data['loading_rate_m_s'])
        self.assertFalse(data['physical_time_available'])
        self.assertFalse(data['contact_zero_verified'])

    def test_all_curve_pairs_are_consumed_in_one_pass(self):
        def rows():
            yield ('0TO-3SA-H', None, None, '20TO-3SA-H', None)
            yield ('Displacement', 'Force', None, 'Displacement', 'Force')
            yield ('um', 'mN', None, 'um', 'mN')
            yield (0, 0, None, 0, 1)
            yield (2, 3, None, None, None)
        a, b = read_sheet_rows('Micro-H', rows())
        self.assertEqual(len(a.points), 2)
        self.assertEqual(b.points, [(4, 0., 1.)])

    def test_invalid_or_partial_measurements_fail_instead_of_disappearing(self):
        for point in [(1, None), (None, 2), (1, math.nan), (math.inf, 1), (True, 2), ('=A1', 2)]:
            with self.subTest(point=point), self.assertRaisesRegex(ValueError, 'row 4'):
                read_sheet_rows('Micro-H', sheet_rows(points=[point]))

    def test_wrong_units_or_hydration_labels_fail(self):
        rows = sheet_rows(points=[(1, 2)])
        rows[2] = ('nm', 'mN')
        with self.assertRaisesRegex(ValueError, 'units'):
            read_sheet_rows('Micro-H', rows)
        with self.assertRaisesRegex(ValueError, 'inconsistent sample'):
            read_sheet_rows('Micro-H', sheet_rows(label='0TO-3SA-D', points=[(1, 2)]))

    def test_duplicate_column_labels_fail(self):
        with self.assertRaisesRegex(ValueError, 'duplicate label'):
            read_sheet_rows('Micro-H', [('0TO-3SA-H', None, '0TO-3SA-H', None),
                                      ('Displacement', 'Force', 'Displacement', 'Force'),
                                      ('um', 'mN', 'um', 'mN'), (0, 0, 1, 1)])


if __name__ == '__main__':
    unittest.main()
