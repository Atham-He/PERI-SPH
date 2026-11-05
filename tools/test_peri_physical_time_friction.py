import copy,json,tempfile,unittest
from decimal import localcontext
from pathlib import Path
from peri_bounded_friction_audit import FrictionInspector,inspect_stream
from peri_physical_time_friction_audit import PhysicalTimeFrictionInspector,protocol,schedule_grid

def case(nominal=5000):
    c=dict(name='new_v2_clock_only_unit',protocol_kind=2,stroke=.01,tangent_stroke=[.0006,-.0004],
           physical_time=dict(denominator=1000000,nominal_ticks=nominal,
                              normal_event_ticks=[160000,240000,400000],tangent_event_ticks=[20000,60000,80000]))
    g=schedule_grid(c);c.update(steps=len(g)-1,dt=nominal/1000000,checkpoint_steps=[len(g)-1]);return c

def header(c):
    p=c['physical_time'];g=schedule_grid(c);ends=p['normal_event_ticks'];counts=[0,0,0]
    for t in g[:-1]:counts[next(k for k,e in enumerate(ends) if t<e)]+=1
    return dict(format=2,protocol_kind=2,physical_clock_profile='explicit_rational_events_v2',
                time_denominator=p['denominator'],nominal_ticks=p['nominal_ticks'],
                normal_event_ticks=ends.copy(),tangent_event_ticks=p['tangent_event_ticks'].copy(),
                steps_requested=len(g)-1,nominal_dt=c['dt'],phase_step_counts=counts,
                phase_end_times=[t/p['denominator'] for t in ends],end_time=.4,
                tangent_cycle_cursor='fixed_rational_C1_early_cycle_then_zero',tangent_cycle_factors=[0,1,-1,0],
                tangent_cycle_end_times=[t/p['denominator'] for t in p['tangent_event_ticks']],checkpoint_steps=c['checkpoint_steps'])

class NewPhysicalTimeReaderControls(unittest.TestCase):
    def test_fixed_physical_end_and_noncommensurate_grid(self):
        for nominal,count in [(10000,40),(5000,80),(2500,160),(7000,59),(13000,35)]:
            c=case(nominal);self.assertEqual(c['steps'],count)
            a=object.__new__(PhysicalTimeFrictionInspector);a.case=c
            a.validate_protocol_header(c,header(c));self.assertEqual(float(a.expected_end_time()),.4)
            g=schedule_grid(c)
            for k,t in enumerate(g[:-1]):
                p=protocol(c,k,t/1000000);self.assertEqual(p['new_tick'],g[k+1])
                self.assertEqual(p['old_time']+p['actual_dt'],p['new_time'])

    def test_old_format_keeps_strict_meaning(self):
        c=case();h=header(c);a=object.__new__(PhysicalTimeFrictionInspector)
        bad=copy.deepcopy(h);bad['format']=1
        with self.assertRaises(AssertionError):a.validate_protocol_header(c,bad)
        v1=object.__new__(FrictionInspector)
        with self.assertRaises(AssertionError):v1.validate_protocol_header(c,h)

    def test_changed_physical_time_not_admitted_by_same_step_count(self):
        c=case();h=header(c);a=object.__new__(PhysicalTimeFrictionInspector)
        h['phase_end_times']=[.08,.12,.2]
        with self.assertRaises(AssertionError):a.validate_protocol_header(c,h)

    def test_in_range_wrong_tick_refused(self):
        c=case();a=object.__new__(PhysicalTimeFrictionInspector);a.case=c;a.maximum=0.
        with localcontext() as ctx:
            ctx.prec=95
            p=protocol(c,1,.005);p['old_tick']+=1
            with self.assertRaises(AssertionError):a.interval(p,dict(step=1,time=.005))

    def test_wrong_accepted_time_refused(self):
        with self.assertRaises(AssertionError):protocol(case(),1,.0050001)

    def test_explicit_reader_class_transport(self):
        class TransportOnly:
            def __init__(self,c,h,full):
                assert c['name']=='transport_only' and h['format']==2 and full is False
                self.geometry={'maximum_normalized_cache_error':0};self.maximum=0;self.trial_points=0
                self.statistics={'positive_input_work':0,'algorithmic_loss':0};self.row_slip_signs={}
            def initial(self):return {'time':0},None,None,None
        c=dict(name='transport_only',steps=0,expected_code=0)
        rows=[dict(type='identity'),dict(type='header',format=2,end_time=0),
              dict(type='end',completed=True,accepted_steps=0,last_accepted_time=0)]
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'zero_transport.jsonl';p.write_text('\n'.join(map(json.dumps,rows))+'\n')
            result=inspect_stream(c,p,0,False,inspector_type=TransportOnly)
            self.assertEqual(result['accepted_steps'],0)

if __name__=='__main__':unittest.main()
