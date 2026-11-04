import math
import unittest
from decimal import Decimal as D,localcontext
from peri_pd_reference_audit import expected_reference,expected_surface_points,jacobi_eigenvalues
from peri_pd_material_reference import material_step
from peri_gel_affine_audit import reference_step,green
from peri_bounded_gel_audit import cubic_minimum

class ActualPDIndependentReferenceTests(unittest.TestCase):
    def test_independent_jacobi_spectrum_and_zero_mode(self):
        eigen=jacobi_eigenvalues([[2,-1,0],[-1,2,0],[0,0,0]])
        for a,b in zip(eigen,[0,1,3]):self.assertAlmostEqual(a,b,places=12)

    def test_fixed_dense_and_pole_rules_preserve_area_and_second_moment(self):
        for kind in ['surface','dense_surface']:
            points,area=expected_surface_points(kind,0,16,1/64)
            self.assertEqual(len(points),114 if kind=='surface' else 640)
            self.assertTrue(all(x>0 for x in area));self.assertAlmostEqual(sum(area),4*math.pi,places=11)
            self.assertAlmostEqual(sum(w*p[1]**2 for p,w in zip(points,area)),4*math.pi/3,places=11)

    def test_actual_shape_reference_has_exact_rigidity_rank(self):
        ref=expected_reference(0,tuple([4*math.pi/99]*33),1.25)
        self.assertEqual(ref['rigidity_rank'],93);self.assertEqual(len(ref['shape_gradient']),33*33)
        self.assertGreater(min(map(float,ref['shape_eigenvalue_ratio'])),1e-3)

    def test_cached_material_reference_matches_qualified_antiderivative(self):
        with localcontext() as ctx:
            ctx.prec=95
            f0=list(map(D,['1','.1','0','0','1','0','0','0','1']))
            f1=list(map(D,['.99','.12','.01','0','.98','0','.01','0','1.01']))
            z=[D(0)]*9;p=list(map(D,['5','1','.5','.3']))
            a=material_step(f0,f1,z,D(0),D('.02'),p,True)
            b=reference_step(f0,f1,z,D(0),D('.02'),p,True)
            for field in ['green','viscous','average_pk1']:
                self.assertLess(max(abs(x-y) for x,y in zip(a[field],b[field])),D('1e-60'))
            for column,other in zip(a['endpoint_tangent'],b['tangent']):
                self.assertLess(max(abs(x-y) for x,y in zip(column,other)),D('1e-60'))

    def test_continuous_wall_path_requires_more_than_endpoint_feasibility(self):
        with localcontext() as ctx:
            ctx.prec=95
            result=[]
            for n in [8,16]:
                theta=D(n-1)/n;h=D(1)/n;u0=theta*theta*(3-2*theta)
                coefficients=[D(0),1-u0-(6*theta-6*theta*theta)*h,-(3-6*theta)*h*h,2*h**3]
                self.assertEqual(sum(coefficients),0)
                result.append(-cubic_minimum(coefficients))
            self.assertGreater(result[0],D('.01'));self.assertLess(result[1],D('.003'))

if __name__=='__main__':unittest.main()
