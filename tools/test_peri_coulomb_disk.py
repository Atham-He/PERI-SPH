import unittest
from decimal import Decimal as D,localcontext
from peri_coulomb_disk_audit import reference

class IndependentDiskReferenceTests(unittest.TestCase):
    def test_exact_two_component_slip_and_energy(self):
        with localcontext() as ctx:
            ctx.prec=95
            r=reference(list(map(D,['1','-.24','-.18','.4','.3','.3'])))
            self.assertEqual(r['residual'],[D(0),D(0)])
            self.assertEqual(r['friction_loss'],D('.15'))
            self.assertEqual(r['maximum_dissipation_error'],0)

    def test_box_is_not_isotropic_disk(self):
        with localcontext() as ctx:
            ctx.prec=95
            r=reference(list(map(D,['1','-.3','-.3','1','1','.3'])))
            self.assertGreater(r['cone_excess'],D('.12'))
            self.assertGreater(abs(r['residual'][0]),D('.08'))

    def test_relative_velocity_cancelled_before_disk_law(self):
        with localcontext() as ctx:
            ctx.prec=95
            r=reference(list(map(D,['1','.1','.05','0','0','.3'])))
            self.assertEqual(r['residual'],[D(0),D(0)])
            self.assertEqual(r['friction_loss'],0)

    def test_exact_boundary_and_origin_choices(self):
        with localcontext() as ctx:
            ctx.prec=95
            r=reference(list(map(D,['1','.5','0','0','0','.5'])))
            self.assertEqual(r['branch'],'slip_exterior')
            self.assertEqual(r['impulse_derivative'],[[D(1),D(0)],[D(0),D(0)]])
            z=reference([D(0)]*5+[D('.3')]);self.assertEqual(z['branch'],'zero_radius_origin')

if __name__=='__main__':unittest.main()
