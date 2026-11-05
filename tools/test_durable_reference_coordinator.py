"""New stream/coordinator structural controls,never a real equation audit."""
import ast
import copy
import json
import os
from pathlib import Path
import sys
import unittest
from peri_reference_checkpoint import Store,canonical,digest,encode,context_descriptor
from peri_durable_reference_stream import DurableSession,build_index,mathematical_image
from peri_durable_reference_launch import run
from durable_reference_fixture import fixture,ToyInspector


class StreamControls(unittest.TestCase):
    def setUp(self):
        self.root=Path(os.environ['DURABLE_COORDINATOR_WORK'])/self._testMethodName;self.root.mkdir()
        self.raw,self.index,self.contract,self.rows=fixture(self.root)
        ToyInspector.contract=self.contract;ToyInspector.calls=[]
        (self.root/'store').mkdir();self.store=Store(self.root/'store',self.contract)

    def session(self,head=None,leaf=None):
        return DurableSession(self.raw,self.index,self.contract,self.store,head,leaf,inspector_type=ToyInspector)

    def tampered(self,rows):
        path=self.root/'altered.raw';path.write_bytes(b''.join(json.dumps(r).encode()+b'\n' for r in rows))
        rec=dict(path=str(path),size=path.stat().st_size,sha256=digest(path.read_bytes()))
        return path,rec

    def test_binary_index_all_prefixes_offsets_and_counts(self):
        self.assertEqual(self.index['trial_points'],8*33)
        data=self.raw.read_bytes()
        for e in self.index['index']:
            c=e['raw_cursor'];self.assertEqual(digest(data[:c['offset']]),c['prefix_sha256'])
        self.assertFalse(self.index['equation_review_qualified'])

    def test_correctSHA_wrong_step_clock_and_native_birth_rejected(self):
        evidence=[]
        for name,mutate in [('step',lambda r:r[3]['state'].update(step=8)),
                            ('ticks',lambda r:r[3]['interval'].update(new_tick=123)),
                            ('dt',lambda r:r[3]['interval'].update(actual_dt=.001)),
                            ('identity',lambda r:r[0].update(start_ticks='987'))]:
            rows=copy.deepcopy(self.rows);mutate(rows);path,rec=self.tampered(rows)
            with self.assertRaises(ValueError) as error:
                build_index(path,self.index['case'],rec,{k:v for k,v in self.rows[0].items() if k!='type'})
            evidence.append(dict(name=name,sha256=rec['sha256'],rejected=True,reason=str(error.exception)))
        (self.root/'correctSHA_controls.json').write_bytes(canonical(evidence))

    def test_hidden_tail_partial_endpoint_and_wrong_wholeSHA_rejected(self):
        rows=copy.deepcopy(self.rows);rows.append(dict(type='hidden'))
        path,rec=self.tampered(rows)
        with self.assertRaisesRegex(ValueError,'hidden'):
            build_index(path,self.index['case'],rec,{k:v for k,v in self.rows[0].items() if k!='type'})
        rows=self.rows[:-1];path,rec=self.tampered(rows)
        with self.assertRaisesRegex(ValueError,'newline'):
            build_index(path,self.index['case'],rec,{k:v for k,v in self.rows[0].items() if k!='type'})
        rec=copy.deepcopy(self.index['raw']);rec['sha256']=digest(b'wrong')
        with self.assertRaisesRegex(ValueError,'whole raw SHA'):
            build_index(self.raw,self.index['case'],rec,{k:v for k,v in self.rows[0].items() if k!='type'})

    def test_actual_loop_complete_commits_full_independent_history(self):
        s=self.session();s.commit('start',cpu=1.)
        for k in range(1,9):s.advance();last=s.commit('s'+str(k),cpu=1.+k)
        s.finish_raw();image=s.image(cpu=9.);s.close()
        self.assertEqual(ToyInspector.calls,['initial']+list(range(1,9)))
        self.assertEqual(last['accepted_step'],8)
        self.assertNotEqual(image['reference_viscous'][0][0],image['old']['viscous'][0][0])
        self.assertEqual(len(list((self.root/'store').glob('g_*'))),9)

    def test_explicit_split_restores_without_initial_or_prefix_step_replay(self):
        s=self.session();s.commit('start',cpu=1.)
        for k in range(1,4):s.advance();mid=s.commit('first'+str(k),cpu=2.)
        saved=s.image(cpu=2.);s.close();ToyInspector.calls=[]
        resumed=self.session(mid['HEAD_sha256'],mid['leaf'])
        self.assertEqual(ToyInspector.calls,[])
        self.assertEqual(encode(resumed.image(cpu=2.),self.contract),encode(saved,self.contract))
        for k in range(4,9):resumed.advance();resumed.commit('resume'+str(k),cpu=3.)
        resumed.finish_raw();resumed.close()
        self.assertEqual(ToyInspector.calls,list(range(4,9)))

    def test_old_tail_selection_preserves_later_files_and_mathematical_state(self):
        s=self.session();s.commit('start',cpu=1.)
        for k in range(1,9):
            s.advance();last=s.commit('whole'+str(k),cpu=2.)
            if k==3:mid=last
        expected=mathematical_image(s.image(cpu=2.),self.contract);s.close()
        preserved={str(p):p.read_bytes() for p in (self.root/'store').rglob('*') if p.is_file() and p.name!='HEAD'}
        resumed=self.session(last['HEAD_sha256'],mid['leaf'])
        for k in range(4,9):resumed.advance();resumed.commit('rollback'+str(k),cpu=4.)
        result=mathematical_image(resumed.image(cpu=4.),self.contract);resumed.close()
        self.assertEqual(expected['mathematical_state_sha256'],result['mathematical_state_sha256'])
        self.assertNotEqual(expected['original_state_sha256'],result['original_state_sha256'])
        self.assertEqual(result['excluded_fields'],['cpu_completed'])
        for p,data in preserved.items():self.assertEqual(Path(p).read_bytes(),data)

    def test_changed_raw_prefix_and_implicit_selection_reject_before_steps(self):
        s=self.session();commit=s.commit('start',cpu=0.);s.close();ToyInspector.calls=[]
        with self.assertRaisesRegex(ValueError,'explicit checkpoint'):
            self.session(commit['HEAD_sha256'])
        data=self.raw.read_bytes();self.raw.write_bytes(data.replace(b'SYNTHETIC_NOT_NATIVE',b'X'*len(b'SYNTHETIC_NOT_NATIVE'),1))
        with self.assertRaisesRegex(ValueError,'prefix SHA'):
            self.session(commit['HEAD_sha256'],commit['leaf'])
        self.assertEqual(ToyInspector.calls,[])

    def test_partial_finalize_and_changed_next_line_fail(self):
        s=self.session();s.commit('start',cpu=0.)
        with self.assertRaisesRegex(ValueError,'partial'):
            s.finish_raw()
        self.index['index'][1]['native_row_sha']=digest(b'wrong')
        with self.assertRaisesRegex(ValueError,'externally frozen'):
            s.advance()
        self.assertEqual(ToyInspector.calls,['initial']);s.close()

    def test_original_finalization_body_is_AST_identical(self):
        ref=Path(__file__).parent/'reference/peri_bounded_friction_audit.py'
        tree=ast.parse(ref.read_text());fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='inspect_stream')
        body=fn.body[0].body;idx=next(i for i,n in enumerate(body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='stat' for t in n.targets))
        new=ast.parse(Path(__file__).with_name('peri_durable_reference_summary.py').read_text())
        out=next(n for n in new.body if isinstance(n,ast.FunctionDef) and n.name=='original_summary')
        self.assertEqual([ast.dump(n) for n in body[idx:]],[ast.dump(n) for n in out.body])


class LaunchControls(unittest.TestCase):
    def setUp(self):
        self.root=Path(os.environ['DURABLE_COORDINATOR_WORK'])/self._testMethodName;self.root.mkdir()
        (self.root/'ledger').mkdir()
        self.contract=dict(schema=1,plan_sha256=digest(b'NEW_SYNTHETIC_COORDINATOR_PLAN'),source_sha256=digest(b'SOURCE'),
            CPU_cap_ns=30*10**9,wall_cap_ns=90*10**9,prior_charges=[],measurement='wait4_rusage_ceil_nanoseconds_v1')
        (self.root/'contract.json').write_bytes(canonical(self.contract))
        self.limits=dict(RSS_bytes=256*2**20,AS_bytes=512*2**20,threads=1,output_bytes=8*2**20,output_roots=[str(self.root/'attempt'),str(self.root/'ledger')],
            reserves=[dict(path=str(self.root),minimum_bytes=150*2**30)])

    def command(self,mode):
        return [sys.executable,'-B',str(Path(__file__).with_name('durable_reference_launch_child.py')),
            str(self.root/'attempt'),mode,str(self.root/'ledger'),str(self.root/'contract.json'),'fresh_control']

    def test_launch_is_reserved_born_released_acknowledged_then_wait4_settled(self):
        from peri_reference_budget import Ledger
        result=run(self.root/'attempt',self.root/'ledger',self.contract,'fresh_control',self.command('success'),None,5*10**9,20*10**9,self.limits)
        self.assertEqual(result['return_code'],0)
        self.assertTrue((self.root/'attempt/acknowledged.json').exists())
        state=Ledger(self.root/'ledger',self.contract).read(result['resource_HEAD_sha256'])
        self.assertEqual(state['attempts']['fresh_control']['phase'],'settled')
        self.assertEqual(state['CPU_ns'],result['CPU_ns'])
        with self.assertRaisesRegex(ValueError,'unused'):
            run(self.root/'attempt',self.root/'ledger',self.contract,'fresh_control',self.command('success'),None,5*10**9,20*10**9,self.limits)

    def test_nonzero_complete_child_records_exact_exit_and_halts(self):
        from peri_reference_budget import Ledger
        result=run(self.root/'attempt',self.root/'ledger',self.contract,'fresh_control',self.command('controlled_failure'),None,5*10**9,20*10**9,self.limits)
        self.assertEqual(result['return_code'],9)
        self.assertTrue(Ledger(self.root/'ledger',self.contract).read(result['resource_HEAD_sha256'])['halted'])

    def test_exit_before_handshake_retains_reservation_and_known_code(self):
        from peri_reference_budget import Ledger
        with self.assertRaises(RuntimeError):
            run(self.root/'attempt',self.root/'ledger',self.contract,'fresh_control',self.command('early_exit'),None,5*10**9,20*10**9,self.limits)
        failure=json.loads((self.root/'attempt/failure.json').read_bytes())
        self.assertEqual(failure['return_code'],7)
        state=Ledger(self.root/'ledger',self.contract).read(failure['resource_HEAD_sha256'])
        self.assertEqual(state['CPU_ns'],5*10**9);self.assertEqual(state['attempts']['fresh_control']['phase'],'authorized')

    def test_disk_reserve_rejects_without_new_launch_or_ledger_event(self):
        self.limits['reserves'][0]['minimum_bytes']=10**30
        with self.assertRaisesRegex(RuntimeError,'disk reserve'):
            run(self.root/'attempt',self.root/'ledger',self.contract,'fresh_control',self.command('success'),None,5*10**9,20*10**9,self.limits)
        self.assertFalse((self.root/'attempt').exists());self.assertEqual(list((self.root/'ledger').iterdir()),[])

    def test_worker_requires_explicit_full_material_admission_before_imports(self):
        from peri_durable_reference_worker import run as worker
        plan=self.root/'denied_plan.json';plan.write_bytes(canonical(dict(material_audit_admitted=False)))
        request=self.root/'denied_request.json'
        request.write_bytes(canonical(dict(material_audit_admitted=False,plan=dict(path=str(plan),size=plan.stat().st_size,sha256=digest(plan.read_bytes())))))
        with self.assertRaisesRegex(ValueError,'not prospectively admitted'):
            worker(request,digest(request.read_bytes()))
        self.assertFalse((self.root/'attempt').exists())
