from fractions import Fraction as Q
import unittest
from peri_normal_projection import probe, row
from peri_normal_projection_audit import exact_projection, check_native


class ExactNormalContactOracleTests(unittest.TestCase):
    def test_moving_wall_has_known_rational_work_kinetic_and_dissipation(self):
        case = probe([2], [0], [0], .2, [row(new=.1)])
        result = exact_projection(case)
        self.assertEqual(result['velocity'], [Q(1, 2)])
        self.assertEqual(result['normal_impulse'], [Q(1)])
        self.assertEqual(result['wall_work'], Q(1, 2))
        self.assertEqual(result['accepted_kinetic'], Q(1, 4))
        self.assertEqual(result['contact_step_dissipation'], Q(1, 4))

    def test_infeasible_motion_and_invalid_initial_geometry_are_distinct(self):
        case = probe([1], [.5], [0], .1, [row(new=.7), row(normal=-1, old=-1, new=-.3)])
        self.assertIsNone(exact_projection(case))
        with self.assertRaises(ValueError): exact_projection(probe([1], [-.1], [0], .1, [row()]))

    def test_wrong_individual_impulses_cannot_pass_by_zero_net_momentum(self):
        case = probe([1, 1], [0, 1], [-1, 1], .1, [row(), row((1,), normal=-1, old=-1, new=-1)])
        native = dict(converged=True, velocity=[0, 0], position=[0, 1], normal_impulse=[1, 1],
                      free_kinetic=1, accepted_kinetic=0, wall_work=0, contact_step_dissipation=1)
        self.assertEqual(check_native(case, native, 0)['status'], 'exact_rational_primal_energy_and_native_dual_KKT_checked')
        with self.assertRaises(ValueError): check_native(case, dict(native, normal_impulse=[0, 0]), 0)


if __name__ == '__main__': unittest.main()
