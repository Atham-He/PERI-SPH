"""Synthetic transport/byte-admission controls only; never PD trajectory proof."""
import copy,json,tempfile,unittest
from pathlib import Path
from peri_physical_checkpoint_support import (trace_index,compare_new_branch,schedule_grid,initial_contract)

def case():
    return dict(steps=80,checkpoint_steps=[16,40,64,80],physical_time=dict(denominator=1000000,nominal_ticks=5000,
        normal_event_ticks=[160000,240000,400000],tangent_event_ticks=[20000,60000,80000]))

def synthetic(c,first,last):
    g=schedule_grid(c);den=c['physical_time']['denominator']
    h=dict(type='header',format=2,initial_state=dict(time=0),time_denominator=den) if first==1 else dict(
        type='resume_header',selected_generation=first-1,accepted_step=first-1,new_initial_state=False,accepted_time=g[first-1]/den)
    rows=[dict(type='identity',pid=1),h]
    for k in range(first,last+1):
        a,b=g[k-1]/den,g[k]/den
        rows.append(dict(type='step',committed=True,local_step_gates_passed=True,candidate_ready=True,reason='all_applicable_step_gates_passed',
            state=dict(step=k,time=b,last_dt=b-a,synthetic_endpoint=[1,2,3]),
            interval=dict(step=k-1,old_tick=g[k-1],new_tick=g[k],time_denominator=den,old_time=a,new_time=b,actual_dt=b-a),
            solve=dict(trace=[dict(synthetic_trial_T=[k,1,2])])) )
    rows.append(dict(type='end',completed=True,accepted_steps=last,last_accepted_time=g[last]/den,actual_solver_calls=last-first+1))
    return rows

class NewV2PacketControls(unittest.TestCase):
    def write(self,root,name,rows):
        p=Path(root)/name;p.write_text('\n'.join(map(json.dumps,rows))+'\n');return p
    def test_full80_and_split16_whole_packets(self):
        c=case()
        with tempfile.TemporaryDirectory() as d:
            a=trace_index(self.write(d,'all',synthetic(c,1,80)),c,1,80)
            b=trace_index(self.write(d,'split',synthetic(c,1,16)),c,1,16)
            self.assertEqual(compare_new_branch(a,b,c)['different_steps'],[])
    def test_changed_intermediate_T_with_identical_endpoint_requires_reference(self):
        c=case();rows=synthetic(c,1,16);rows[9]['solve']['trace'][0]['synthetic_trial_T'][1]=99
        with tempfile.TemporaryDirectory() as d:
            a=trace_index(self.write(d,'all',synthetic(c,1,80)),c,1,80)
            b=trace_index(self.write(d,'changed',rows),c,1,16)
            self.assertEqual(compare_new_branch(a,b,c)['different_steps'],[8])
    def test_explicit_resume_clock(self):
        c=case()
        with tempfile.TemporaryDirectory() as d:
            a=trace_index(self.write(d,'all',synthetic(c,1,80)),c,1,80)
            b=trace_index(self.write(d,'resume',synthetic(c,17,80)),c,17,80,resume_generation=16)
            self.assertEqual(compare_new_branch(a,b,c)['different_steps'],[])
    def test_range_valid_wrong_tick_refused(self):
        c=case();rows=synthetic(c,1,16);rows[5]['interval']['old_tick']+=1
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(AssertionError):trace_index(self.write(d,'bad',rows),c,1,16)
    def test_hidden_extra_native_tail_refused(self):
        c=case();rows=synthetic(c,1,16)+[dict(type='hidden_extra')]
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(AssertionError):trace_index(self.write(d,'bad',rows),c,1,16)

if __name__=='__main__':unittest.main()
