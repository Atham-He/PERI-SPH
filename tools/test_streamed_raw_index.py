"""New nondivisor synthetic complete-payload index equivalence controls."""
import copy,json,os,tempfile,unittest
from pathlib import Path
from peri_reference_checkpoint import canonical,digest
from peri_durable_reference_stream import build_index
from peri_streamed_raw_index import build_streamed_index
from peri_physical_time_friction_audit import schedule_grid


def fixture(path,dt):
    clock=dict(denominator=1000,nominal_ticks=dt,normal_event_ticks=[4,7,11],tangent_event_ticks=[2,5,8])
    case=dict(name='NEW_SYNTHETIC_NONDIVISOR_INDEX_NOT_PHYSICAL',physical_time=clock)
    grid=schedule_grid(case);case['steps']=len(grid)-1
    meta=dict(type='identity',pid=1,ppid=2,start_ticks='3',boot_id='00000000-0000-0000-0000-000000000000',argv=['SYNTHETIC_ONLY'])
    header=dict(type='header',format=2,steps_requested=case['steps'],time_denominator=1000,nominal_ticks=dt,
        normal_event_ticks=clock['normal_event_ticks'],tangent_event_ticks=clock['tangent_event_ticks'],initial_state=dict(step=0,time=0.),initial_friction_history=dict(sum=[0.,-0.]))
    rows=[meta,header]
    for step,(a,b) in enumerate(zip(grid,grid[1:]),1):
        rows.append(dict(type='step',state=dict(step=step,time=b/1000),committed=True,friction_history=dict(sum=[float(step),-0.]),
            interval=dict(old_tick=a,new_tick=b,time_denominator=1000,old_time=a/1000,new_time=b/1000,actual_dt=b/1000-a/1000),
            solve=dict(trace=[dict(iteration=k,payload=dict(T=[[float((i+j+k)%13)/13 for j in range(81)] for i in range(33)])) for k in range(2)])))
    rows.append(dict(type='end',completed=True,accepted_steps=case['steps'],last_accepted_time=.011))
    path.write_bytes(b''.join(json.dumps(r,ensure_ascii=False).encode()+b'\n' for r in rows))
    return case,{k:v for k,v in meta.items() if k!='type'},dict(path=str(path),size=path.stat().st_size,sha256=digest(path.read_bytes())),rows


class StreamedIndexControls(unittest.TestCase):
    def test_complete_index_all_bytes_equal_legacy_on_four_NEW_grids(self):
        with tempfile.TemporaryDirectory() as d:
            for dt in [1,3,5,20]:
                path=Path(d)/('case'+str(dt)+'.raw');case,identity,rec,_=fixture(path,dt)
                expected=build_index(path,case,rec,identity);stages=[];frontiers=[]
                actual=build_streamed_index(path,case,rec,identity,lambda s,k:stages.append((s,k)),lambda k,e:frontiers.append((k,e)))
                self.assertEqual(canonical(expected),canonical(actual))
                self.assertEqual(len(frontiers),case['steps']+1)
                self.assertEqual([k for s,k in stages if s=='before_row_release'],list(range(1,case['steps']+1)))

    def test_complete_payload_corruption_changes_every_relevant_binding(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'case.raw';case,identity,rec,rows=fixture(path,3)
            old=build_streamed_index(path,case,rec,identity)
            rows[2]['solve']['trace'][1]['payload']['T'][32][80]+=.0001
            path.write_bytes(b''.join(json.dumps(r).encode()+b'\n' for r in rows));rec.update(size=path.stat().st_size,sha256=digest(path.read_bytes()))
            new=build_streamed_index(path,case,rec,identity)
            self.assertNotEqual(old['index'][1]['native_row_sha'],new['index'][1]['native_row_sha'])
            self.assertNotEqual(old['index'][1]['raw_cursor'],new['index'][1]['raw_cursor'])
            self.assertNotEqual(old['whole_raw_sha256'],new['whole_raw_sha256'])

    def test_correctSHA_wrong_tick_rejected_without_final_index(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'case.raw';case,identity,rec,rows=fixture(path,5);rows[2]['interval']['new_tick']+=1
            path.write_bytes(b''.join(json.dumps(r).encode()+b'\n' for r in rows));rec.update(size=path.stat().st_size,sha256=digest(path.read_bytes()))
            frontiers=[]
            with self.assertRaisesRegex(ValueError,'rational cursor'):
                build_streamed_index(path,case,rec,identity,accept=lambda k,e:frontiers.append(k))
            self.assertEqual(frontiers,[0])

    def test_hidden_tail_and_wrong_native_identity_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'case.raw';case,identity,rec,rows=fixture(path,20)
            wrong=copy.deepcopy(identity);wrong['start_ticks']='99'
            with self.assertRaisesRegex(ValueError,'birth'):
                build_streamed_index(path,case,rec,wrong)
            with path.open('ab') as f:f.write(b'{"type":"hidden"}\n')
            rec.update(size=path.stat().st_size,sha256=digest(path.read_bytes()))
            with self.assertRaisesRegex(ValueError,'hidden'):
                build_streamed_index(path,case,rec,identity)
