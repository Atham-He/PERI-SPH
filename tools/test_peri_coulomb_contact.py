import unittest
from decimal import Decimal as D,localcontext
from fractions import Fraction as F
from peri_coulomb_contact_audit import coefficients,fdot,fraction_disk,original_reference,operator_reference


def plane(weights=None):
    return dict(weights=[1] if weights is None else weights,normal=[1,0,0],tangent=[[0,1,0],[0,0,1]],
                old_offset=0,new_offset=0,mu=.4,wall_tangent_velocity=[0,0])


class CoupledIndependentReferenceTests(unittest.TestCase):
    def test_branch_is_decided_before_finite_decimal_conversion(self):
        with localcontext() as ctx:
            ctx.prec=95
            tiny=F(1,10**100)
            inside=fraction_disk(F(1),[F(1),F(0)],[tiny,F(0)],F(1))
            outside=fraction_disk(F(1),[F(1),F(0)],[-tiny,F(0)],F(1))
            self.assertEqual(inside['branch'],'stick_interior')
            self.assertEqual(outside['branch'],'slip_exterior')
            self.assertEqual(inside['residual'][0],D('1e-100'))
            self.assertEqual(inside['velocity_derivative'],[[D(1),D(0)],[D(0),D(1)]])

    def test_signed_gather_transpose_preserves_exact_work(self):
        c=coefficients([plane([-.25,1.25])],6)
        v=[.1,.2,.3,-.4,.5,-.6];p=[.7,-.8,.9]
        gathered=fdot([fdot(row,v) for row in c],p)
        scattered=fdot(v,[fdot([row[i] for row in c],p) for i in range(6)])
        self.assertEqual(gathered,scattered)
        self.assertNotEqual(gathered,fdot(v,[sum(abs(row[i])*p[k] for k,row in enumerate(c)) for i in range(6)]))

    def test_joint_nonsymmetric_nonlinear_jacobian(self):
        with localcontext() as ctx:
            ctx.prec=95
            header=dict(operator_matrix=[[2,.2,-.1],[.05,1.4,.3],[-.12,.1,1.2]],
                        operator_load=[.2,.3,.4],operator_cubic_beta=[.1,.2,.3])
            v=list(map(F,['.1','.3','-.2']));p=list(map(F,['.5','-.1','.04']))
            I=F('1.7');V=F('.8');dt=F('.4');contacts=[plane()]
            def ref(v,p):
                R,J=operator_reference(header,v)
                return original_reference([0]*3,contacts,dt,I,V,v,p,R,J)
            r=ref(v,p);epsilon=F(1,10**22)
            for k in range(6):
                vp,vm,pp,pm=v.copy(),v.copy(),p.copy(),p.copy()
                if k<3:vp[k]+=epsilon*V;vm[k]-=epsilon*V
                else:pp[k-3]+=epsilon*I;pm[k-3]-=epsilon*I
                a,b=ref(vp,pp),ref(vm,pm)
                for row in range(6):
                    derivative=(a['residual'][row]-b['residual'][row])/D('2e-22')
                    self.assertLess(abs(derivative-r['jacobian'][row][k]),D('1e-38'))
            self.assertNotEqual(r['jacobian'][0][1],r['jacobian'][1][0])
            self.assertNotEqual(r['jacobian'][4][3],0)  # normal radius couples to tangent

    def test_original_gather_cancellation_is_not_float_sum(self):
        c=coefficients([plane([1e16,1,-1e16])],9)
        v=[0,1,0]*3
        self.assertEqual(fdot(c[1],v),1)
        premature=0.
        for a,b in zip(c[1],v):premature+=a*b
        self.assertEqual(premature,0)
        with localcontext() as ctx:
            ctx.prec=95
            r=original_reference([0]*9,[plane([1e16,1,-1e16])],1,1,1,v,[1,0,0],[0]*9,[[0]*9 for _ in range(9)])
            self.assertEqual(r['relative_velocity'][0][0],1)

    def test_mass_impact_energy_uses_physical_mass(self):
        with localcontext() as ctx:
            ctx.prec=95
            q=D(5).sqrt();v=[D(0),2*(1-D('.3')/q),1-D('.3')/q]
            p=[D(1),-D('.6')/q,-D('.3')/q];vf=[D(-1),D(2),D(1)]
            kinetic=sum(x*x for x in v)/2;algorithmic=sum((x-y)**2 for x,y in zip(v,vf))/2
            friction=-sum(p[k]*v[k] for k in [1,2])
            self.assertLess(abs(3-kinetic-algorithmic-friction),D('1e-90'))
            self.assertGreater(friction,0)
            self.assertGreater(algorithmic,0)

    def test_analytic_controls_cover_frames_and_positive_friction(self):
        from peri_coulomb_global_controls import controls
        from peri_coulomb_global_audit import analytic_particle
        with localcontext() as ctx:
            ctx.prec=95
            cases={c['name']:c for c in controls()}
            a=analytic_particle(cases['physical_slide_two_components']['input'])
            b=analytic_particle(cases['normalization_preserves_physical_mass']['input'])
            self.assertEqual(a,b)
            b=analytic_particle(cases['coordinate_frame_covariance']['input'])
            self.assertEqual(b[0],[a[0][1],a[0][2],a[0][0]])
            self.assertGreater(a[2]['friction_loss'],0)
            for label in ['physical_strict_stick','physical_common_wall_velocity','physical_zero_mu','physical_separated']:
                self.assertEqual(analytic_particle(cases[label]['input'])[2]['friction_loss'],0)

    def test_new_input_controls_have_original_sizes_and_explicit_scope(self):
        from peri_coulomb_global_controls import controls
        cases=controls()
        self.assertEqual(len({c['name'] for c in cases}),30)
        for c in cases:
            p=c['input'];n=p['n'];m=p['m']
            self.assertEqual(len(p['linear_rhs']),n+3*m)
            self.assertEqual(len(p['trial_impulses']),3*m)
            expected=10+3*n+n*n+m*(n//3+14)+n+3*m+n+3*m
            self.assertEqual(len(c['tokens']),expected)

if __name__=='__main__':unittest.main()
