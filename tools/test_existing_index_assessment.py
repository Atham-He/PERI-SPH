"""NEW assessment controls,synthetic packets only,not prior component reruns."""
import copy,json,tempfile,unittest
from pathlib import Path
from peri_streamed_raw_index import build_streamed_index
from peri_existing_index_assessment import assess
from peri_buffered_metadata import buffered_record
from peri_physical_time_friction_audit import schedule_grid


def fixture(path):
    clock=dict(denominator=1000,nominal_ticks=3,normal_event_ticks=[4,7,11],tangent_event_ticks=[2,5,8])
    case=dict(name='NEW_EXISTING_ARTIFACT_ASSESSMENT_SYNTHETIC_ONLY',physical_time=clock)
    grid=schedule_grid(case);case['steps']=len(grid)-1
    meta=dict(type='identity',pid=3,ppid=2,start_ticks='1',boot_id='00000000-0000-0000-0000-000000000000',argv=['SYNTHETIC_ONLY'])
    header=dict(type='header',format=2,steps_requested=case['steps'],time_denominator=1000,nominal_ticks=3,
        normal_event_ticks=clock['normal_event_ticks'],tangent_event_ticks=clock['tangent_event_ticks'],primal_dofs=99,contact_count=116,unknowns=447,
        reference=dict(position=[[0.,0.,0.]]*33),initial_state=dict(step=0,time=0.,last_dt=0.),initial_friction_history=dict(sum=[0.,-0.]),end_time=.011)
    rows=[meta,header]
    for step,(a,b) in enumerate(zip(grid,grid[1:]),1):
        rows.append(dict(type='step',committed=True,local_step_gates_passed=True,candidate_ready=True,
            state=dict(step=step,time=b/1000,last_dt=b/1000-a/1000),friction_history=dict(sum=[float(step),-0.]),
            interval=dict(old_tick=a,new_tick=b,time_denominator=1000,old_time=a/1000,new_time=b/1000,actual_dt=b/1000-a/1000),
            solve=dict(converged=True,trace=[dict(iteration=k,complete_T_payload=[[float((i+j+k)%11)/11 for j in range(81)] for i in range(33)]) for k in range(2)])))
    rows.append(dict(type='end',completed=True,accepted_steps=case['steps'],last_accepted_time=.011))
    path.write_bytes(b''.join(json.dumps(r).encode()+b'\n' for r in rows));rec=buffered_record(str(path));birth={k:v for k,v in meta.items() if k!='type'}
    index=build_streamed_index(path,case,rec,birth);frontiers=copy.deepcopy(index['index'])
    return index,case,birth,rec,frontiers,rows


class ExistingArtifactControls(unittest.TestCase):
    def test_complete_existing_index_all_packets_verified_without_building(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'raw';i,c,b,r,f,_=fixture(p);events=[]
            q=assess(i,p,c,b,r,f,lambda s,k:events.append((s,k)))
            self.assertEqual(q['accepted_steps'],c['steps']);self.assertEqual(q['point_metadata'],66*c['steps'])
            self.assertEqual([k for s,k in events if s=='accepted_metadata_verified'],list(range(c['steps']+1)))
            self.assertFalse(q['material_equation_qualified'])

    def test_full_T_payload_mutation_with_correct_wholeSHA_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'raw';i,c,b,r,f,rows=fixture(p)
            rows[2]['solve']['trace'][1]['complete_T_payload'][32][80]+=.1
            p.write_bytes(b''.join(json.dumps(row).encode()+b'\n' for row in rows));new=buffered_record(str(p))
            i['raw']=new;i['whole_raw_sha256']=new['sha256']
            with self.assertRaises(ValueError):assess(i,p,c,b,new,f)

    def test_frontier_case_and_native_birth_semantic_rejections(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'raw';i,c,b,r,f,_=fixture(p)
            wrong=copy.deepcopy(f);wrong[2]['native_row_sha']='f'*64
            with self.assertRaisesRegex(ValueError,'frontier bytes'):assess(i,p,c,b,r,wrong)
            case=copy.deepcopy(c);case['name']='WRONG'
            with self.assertRaisesRegex(ValueError,'case typed'):assess(i,p,case,b,r,f)
            birth=dict(b,start_ticks='999')
            with self.assertRaisesRegex(ValueError,'argv/birth'):assess(i,p,c,birth,r,f)

    def test_runtime_false_flags_even_consistent_reindexed_artifact_reject(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'raw';i,c,b,r,f,rows=fixture(p);rows[2]['local_step_gates_passed']=False
            p.write_bytes(b''.join(json.dumps(row).encode()+b'\n' for row in rows));r=buffered_record(str(p))
            i=build_streamed_index(p,c,r,b);f=copy.deepcopy(i['index'])
            with self.assertRaisesRegex(ValueError,'runtime record flags'):assess(i,p,c,b,r,f)

    def test_partial_index_and_hidden_tail_reject(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'raw';i,c,b,r,f,_=fixture(p);partial=copy.deepcopy(i);partial['index'].pop()
            with self.assertRaisesRegex(ValueError,'all accepted'):assess(partial,p,c,b,r,f)
            with p.open('ab') as stream:stream.write(b'{"hidden":true}\n')
            r=buffered_record(str(p));i['raw']=r;i['whole_raw_sha256']=r['sha256']
            with self.assertRaisesRegex(ValueError,'hidden'):assess(i,p,c,b,r,f)
