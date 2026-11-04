import math
import unittest
from decimal import Decimal as D,localcontext
from peri_gel_affine_audit import reference_step,green,identity_matrix,scale,check


class AffineMaterialReferenceTests(unittest.TestCase):
    def test_identity_tangent_matches_isotropic_closed_form(self):
        case=dict(name='identity',parameters=[5000,1000,500,.3],initial=[1,0,0,0,1,0,0,0,1],green=[0]*9,
            viscous=[0]*9,old_dissipation=0,steps=[dict(dt=.15,F=[1,0,0,0,1,0,0,0,1])])
        h=.5; alpha=-math.expm1(-h)/h; shear=500+500*(1-alpha)/h
        tangent=[]
        for column in range(9):
            H=[float(i==column) for i in range(9)]; trace=H[0]+H[4]+H[8]
            tangent.append([shear*(H[i]+H[(i%3)*3+i//3]-(2*trace/3 if i//3==i%3 else 0))+
                            (2500*trace if i//3==i%3 else 0) for i in range(9)])
        row=dict(type='step',step=1,dt=.15,green=[0]*9,viscous=[0]*9,average_pk1=[0]*9,tangent=tangent,
            stored_old=0,stored_new=0,physical_dissipation=0,total_dissipation=0,work=0,energy_residual=0,
            minimum_determinant=1,tangent_asymmetry=0,tangent_available=True)
        result=check(case,[row,dict(type='terminal',completed=True)],0)
        self.assertLess(result['maximum_normalized_error'],1e-13)

    def test_held_relaxation_has_closed_positive_loss(self):
        with localcontext() as context:
            context.prec=65
            f=[D('1.1'),D(0),D(0),D(0),D('.95'),D(0),D(0),D(0),D(1)]
            p=list(map(D,[5000,1000,500]))+[D('.3')]
            ref=reference_step(f,f,[D(0)]*9,D(0),D('.15'),p,False)
            x=green(f); expected=p[2]*(1-(-D(1)).exp())*sum(v*v for v in x)
            self.assertLess(abs(ref['physical_dissipation']-expected),D('1e-55'))

    def test_affine_F_is_not_linear_Green(self):
        with localcontext() as context:
            context.prec=65
            f=identity_matrix(); final=list(f); final[1]=D('.2')
            p=[D(5000),D(1000),D(500),D('.3')]
            ref=reference_step(f,final,[D(0)]*9,D(0),D('.15'),p,False)
            alpha=(1-(-D('.5')).exp())/D('.5')
            wrong=scale(green(final),1-alpha)
            self.assertGreater(max(abs(a-b) for a,b in zip(ref['viscous'],wrong)),D('1e-5'))

    def test_positive_endpoints_do_not_hide_singular_interior(self):
        with localcontext() as context:
            context.prec=65
            final=[D(-1),D(0),D(0),D(0),D(-1),D(0),D(0),D(0),D(1)]
            with self.assertRaises(ValueError):
                reference_step(identity_matrix(),final,[D(0)]*9,D(0),D('.1'),[D(5000),D(1000),D(500),D('.3')],False)


if __name__=='__main__': unittest.main()
