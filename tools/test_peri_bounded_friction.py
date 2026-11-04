import unittest
from decimal import Decimal as D,localcontext
from fractions import Fraction as F
from peri_bounded_friction_audit import protocol,neumaier
from peri_pd_reference_audit import check_geometry


def case(kind=1):
    return dict(dt=.01,stroke=.01,tangent_stroke=[.003,-.002],protocol_kind=kind)


class FrictionMotionAndHistoryReferenceTests(unittest.TestCase):
    def test_reference_16_step_clock_and_normal_motion_preserved(self):
        with localcontext() as ctx:
            ctx.prec=95
            time=0.;duration=16*.01
            for step in range(16):
                p=protocol(case(0),step,time)
                self.assertEqual(p['phase'],0)
                self.assertEqual(p['actual_dt'],duration-time if step==15 else .01)
                theta=D.from_float(time)/D.from_float(duration)
                self.assertLess(abs(p['old_offset'][1]-(-1+D.from_float(.01)*theta*theta*(3-2*theta))),D('1e-90'))
                time=p['new_time']
            self.assertEqual(time,duration)

    def test_piecewise_cubic_continuity_hold_release_and_tangent_telescoping(self):
        with localcontext() as ctx:
            ctx.prec=95
            time=0.;previous=None;travel=[D(0)]*2;positive=False;negative=False
            for step in range(40):
                p=protocol(case(),step,time)
                self.assertEqual(p['step'],step)
                if previous:
                    self.assertEqual(p['old_offset'],previous['new_offset'])
                    self.assertEqual(p['old_tangent_displacement'],previous['new_tangent_displacement'])
                cubic=p['normal_cubic'][1]
                # Floating accepted time is round(t+dt), while the cubic uses
                # the actual dt parameter. Bound the endpoint-clock rounding
                # with the known smoothstep Lipschitz constant (<=1.5).
                clock=abs(D.from_float(p['new_time'])-D.from_float(p['old_time'])-D.from_float(p['actual_dt']))
                duration=D.from_float(.16 if p['phase']==0 else .4-.24)
                bound=3*D.from_float(.01)*clock/duration+D('1e-90')
                self.assertLessEqual(abs(sum(cubic)-p['new_offset'][1]),bound)
                for a in range(2):
                    change=p['new_tangent_displacement'][1][a]-p['old_tangent_displacement'][1][a]
                    travel[a]+=change
                if step<16:positive|=p['new_tangent_displacement'][1][0]>p['old_tangent_displacement'][1][0]
                elif step<24:
                    self.assertEqual(cubic[1:],[D(0)]*3)
                    self.assertEqual(p['old_tangent_displacement'],p['new_tangent_displacement'])
                else:negative|=p['new_tangent_displacement'][1][0]<p['old_tangent_displacement'][1][0]
                previous=p;time=p['new_time']
            self.assertEqual(time,40*.01)
            self.assertEqual(previous['new_offset'],[D(-1),D(-1)])
            self.assertLess(max(map(abs,travel)),D('1e-90'))
            self.assertTrue(positive and negative)

    def test_hidden_between_endpoint_penetration_is_detected(self):
        from peri_bounded_gel_audit import cubic_minimum
        with localcontext() as ctx:
            ctx.prec=95
            time=0.
            for step in range(16):p=protocol(case(),step,time);time=p['new_time']
            c=p['normal_cubic'][1]
            # A particle exactly on the plate at both endpoints may lie below
            # the true cubic plate between endpoints. Endpoint-only QA misses it.
            gap=[D(0),sum(c[1:])-c[1],-c[2],-c[3]]
            self.assertLess(cubic_minimum(gap),0)
            self.assertLess(abs(sum(gap)),D('1e-90'))
            self.assertAlmostEqual(float(-cubic_minimum(gap)/D.from_float(.01)),.002746921031773374,places=12)

    def test_compensated_history_is_not_replaced_by_plain_accumulation(self):
        s=c=0.
        for x in [1e16,1.,-1e16]:s,c=neumaier(x,s,c)
        self.assertEqual(s,0.)
        self.assertEqual(c,1.)
        self.assertEqual(F.from_float(s)+F.from_float(c),1)

    def test_normal_only_reader_rejects_friction_scope_by_default(self):
        header=dict(actual_volume_nodes=33,reference=dict(generation=1,horizon=1.25),
            normal_only=False,full_S3_complete=False,normalized_benchmark=True)
        with self.assertRaises(AssertionError):check_geometry(dict(rotation=0),header)

    def test_early_tangent_cycle_preserves_normal_and_accepted_clock(self):
        with localcontext() as ctx:
            ctx.prec=95;time=0.;signs=[]
            for step in range(40):
                old=protocol(case(1),step,time);new=protocol(case(2),step,time)
                for field in ['old_time','new_time','actual_dt','old_offset','new_offset','normal_cubic']:
                    self.assertEqual(old[field],new[field])
                change=new['new_tangent_displacement'][1][0]-new['old_tangent_displacement'][1][0]
                signs.append(1 if change>0 else (-1 if change<0 else 0))
                time=new['new_time']
            self.assertEqual(signs,[1]*2+[-1]*4+[1]*2+[0]*32)

    def test_early_cycle_endpoints_are_c1_and_held_afterward(self):
        from peri_bounded_friction_audit import early_tangent_shape
        with localcontext() as ctx:
            ctx.prec=95
            self.assertEqual([early_tangent_shape(t,.01) for t in [0,.02,.06,.08,.4]],[0,1,-1,0,0])
            for t in [.02,.06,.08]:
                e=1e-8
                for direction in [-1,1]:
                    slope=abs((early_tangent_shape(t+direction*e,.01)-early_tangent_shape(t,.01))/D.from_float(e))
                    self.assertLess(slope,D('1e-4'))


if __name__=='__main__':unittest.main()
