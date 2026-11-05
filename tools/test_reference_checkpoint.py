import copy,json,os,subprocess,sys,unittest
from pathlib import Path
import decimal as dec
from fractions import Fraction
from peri_reference_checkpoint import (Store,pack,unpack,canonical,digest,encode,decode,snapshot,install,context_descriptor)
from reference_checkpoint_fixture import initial,advance,inspector

WORK=Path(os.environ['REFERENCE_CP_CONTROL_WORK']);WORK.mkdir(exist_ok=True)
def tree(root):return {str(p.relative_to(root)):digest(p.read_bytes()) for p in root.rglob('*') if p.is_file()}

class DurableReferenceControls(unittest.TestCase):
    def setUp(self):
        self.root=WORK/self._testMethodName;self.root.mkdir();self.contract,self.zero,self.case,self.header=initial()
    def state(self,n):
        value=self.zero
        for _ in range(n):value=advance(value,self.contract,self.case,self.header)
        return value
    def store(self):
        p=self.root/'store';p.mkdir();return Store(p,self.contract)
    def test_lossless_decimal_quantum_negativezero_fraction_and_key_types(self):
        value={('typed',Fraction(1,3)):[dec.Decimal('-0E-195'),-0.,dec.Decimal('1.2300000000000000000000'),b'\x00\xff'],17:({-1,1},frozenset({3,5}),True,None)}
        self.assertEqual(pack(unpack(pack(value))),pack(value))
        s=self.state(3);before=context_descriptor();data=encode(s,self.contract);out=decode(data,self.contract)
        self.assertEqual(encode(out,self.contract),data);self.assertEqual(context_descriptor(),before)
        (self.root/'roundtrip.bin').write_bytes(data)
    def test_complete_memory_inventory_snapshot_noadvance_install(self):
        s=self.state(3);audit,capture=inspector(s,self.case,self.header)
        data=encode(snapshot(audit,s['old'],s['reference_viscous'],s['reference_loss'],s['friction_history'],capture,s['raw_cursor'],self.contract,s['cpu_completed']),self.contract)
        destination,target=inspector(self.zero,self.case,self.header);old,z,loss,history=install(decode(data,self.contract),destination,target,self.contract)
        restored=snapshot(destination,old,z,loss,history,target,s['raw_cursor'],self.contract,s['cpu_completed'])
        self.assertEqual(encode(restored,self.contract),data);self.assertEqual(context_descriptor(),s['decimal_context'])
        destination.unregistered_history=1
        with self.assertRaisesRegex(ValueError,'inventory'):snapshot(destination,old,z,loss,history,target,s['raw_cursor'],self.contract,s['cpu_completed'])
    def test_split_and_uninterrupted_all_bytes_equal(self):
        full=self.state(8);store=self.store();s=self.state(3);p=store.publish(s,'split3',None)
        restored,data=store.read(p['leaf'],p['HEAD_sha256'])
        for _ in range(5):restored=advance(restored,self.contract,self.case,self.header)
        self.assertEqual(encode(full,self.contract),encode(restored,self.contract))
        (self.root/'continuous8.bin').write_bytes(encode(full,self.contract));(self.root/'split8.bin').write_bytes(encode(restored,self.contract))
    def test_older_tail_selection_preserves_later_files_and_exact_result(self):
        store=self.store();three=self.state(3);a=store.publish(three,'old3',None);eight=self.state(8);b=store.publish(eight,'later8',a['HEAD_sha256']);before=tree(store.path)
        selected,_=store.read(a['leaf'],b['HEAD_sha256']);c=store.publish(selected,'rollback3',b['HEAD_sha256'])
        for _ in range(5):selected=advance(selected,self.contract,self.case,self.header)
        self.assertEqual(encode(selected,self.contract),encode(eight,self.contract))
        after=tree(store.path)
        self.assertTrue(all(name=='HEAD' or after[name]==value for name,value in before.items()))
        self.assertEqual(c['sequence'],3)
    def test_precision_traps_flags_restore_losslessly(self):
        s=self.state(2);s['decimal_context']['flags']['Clamped']=True;blob=encode(s,self.contract)
        audit,capture=inspector(self.zero,self.case,self.header);dec.getcontext().clear_flags()
        install(decode(blob,self.contract),audit,capture,self.contract)
        self.assertEqual(context_descriptor(),s['decimal_context']);self.assertTrue(dec.getcontext().flags[dec.Clamped])
    def test_correctSHA_semantic_corruption_rejected_without_live_mutation(self):
        s=self.state(3);audit,capture=inspector(s,self.case,self.header)
        changes=[('raw_cursor',lambda x:x['raw_cursor'].update(offset=x['raw_cursor']['offset']+1)),
          ('precision',lambda x:x['decimal_context'].update(prec=94)),
          ('trap',lambda x:x['decimal_context']['traps'].update(Overflow=False)),
          ('clock',lambda x:x['old'].update(time=.0001)),
          ('nonsymmetric_history',lambda x:x['reference_viscous'][0].__setitem__(1,dec.Decimal('.01'))),
          ('negative_loss',lambda x:x['reference_loss'].__setitem__(0,dec.Decimal('-1'))),
          ('missing_statistics',lambda x:x['audit']['statistics'].pop('normal_loss')),
          ('missing_capture',lambda x:x['capture'].update(selected=None)),
          ('trial_counter',lambda x:x['audit'].update(trial_points=99*33)),
          ('unknown_member',lambda x:x.update(unregistered_context=True)),
          ('row_history',lambda x:x['audit']['row_slip_signs'].update({(116,0):{1}})),
          ('CPU',lambda x:x.update(cpu_completed=-1.))]
        original=encode(s,self.contract);context_before=context_descriptor()
        for name,change in changes:
            bad=copy.deepcopy(s);change(bad);payload=pack(bad)
            data=canonical(dict(schema=1,contract_sha256=digest(canonical(self.contract)),payload=payload,payload_sha256=digest(canonical(payload))))
            (self.root/(name+'.bin')).write_bytes(data)
            with self.subTest(name=name):
                with self.assertRaises(ValueError):decode(data,self.contract)
            self.assertEqual(encode(s,self.contract),original);self.assertEqual(context_descriptor(),context_before)
    def test_contract_binding_rejects_new_self_consistent_candidate(self):
        s=self.state(2);other=copy.deepcopy(self.contract);other['raw_sha256']=digest(b'wrongRaw')
        blob=encode(s,other)
        with self.assertRaisesRegex(ValueError,'binding'):decode(blob,self.contract)
    def test_snapshot_byte_tamper_rejected_at_selected_SHA(self):
        store=self.store();p=store.publish(self.state(2),'good2',None)
        f=store.path/p['leaf']/'state.json';data=f.read_bytes();f.write_bytes(data[:-1]+b' ')
        with self.assertRaisesRegex(ValueError,'SHA'):store.read(p['leaf'],p['HEAD_sha256'])
    def test_stale_HEAD_rejected_before_branch_or_state_mutation(self):
        store=self.store();p=store.publish(self.zero,'zero',None);before=tree(store.path)
        with self.assertRaisesRegex(ValueError,'HEAD'):store.publish(self.state(1),'forbidden','0'*64)
        self.assertEqual(tree(store.path),before)
        with self.assertRaisesRegex(ValueError,'ancestor'):store.read('g_automatic_latest',p['HEAD_sha256'])
    def test_uncommitted_directory_never_adopted_by_selection(self):
        store=self.store();p=store.publish(self.zero,'zero',None)
        def hook(stage):
            if stage=='before_HEAD':raise RuntimeError('declared_precommit_stop')
        with self.assertRaisesRegex(RuntimeError,'declared'):store.publish(self.state(1),'orphan',p['HEAD_sha256'],hook)
        self.assertEqual(store.head(),p['HEAD_sha256'])
        with self.assertRaisesRegex(ValueError,'ancestor'):store.read('g_orphan',p['HEAD_sha256'])
    def killed_child(self,stage):
        store=self.store();p=store.publish(self.zero,'zero',None);candidate=self.state(1)
        (self.root/'contract.json').write_bytes(canonical(self.contract));(self.root/'candidate.bin').write_bytes(encode(candidate,self.contract))
        worker=Path(__file__).resolve().parent/'reference_checkpoint_child.py'
        command=[sys.executable,'-B',str(worker),str(self.root),stage,p['HEAD_sha256']]
        child=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        out,err=child.communicate(timeout=30);(self.root/'child.stdout').write_bytes(out);(self.root/'child.stderr').write_bytes(err)
        meta=json.loads(out);self.assertEqual(meta['pid'],child.pid);self.assertEqual(meta['ppid'],os.getpid());self.assertEqual(meta['argv'],command)
        self.assertEqual(child.returncode,-9);self.assertEqual(err,b'');self.assertEqual((self.root/(stage+'.witness')).read_text(),stage)
        (self.root/'child.attempt.json').write_bytes(canonical(dict(command=command,identity=meta,return_code=child.returncode,expected_signal=9)))
        return store,p,candidate
    def test_actual_child_kill_before_HEAD_keeps_old_frontier(self):
        store,p,candidate=self.killed_child('before_HEAD');self.assertEqual(store.head(),p['HEAD_sha256'])
        with self.assertRaisesRegex(ValueError,'ancestor'):store.read('g_child_before_HEAD',p['HEAD_sha256'])
        restored,_=store.read('g_zero',p['HEAD_sha256']);self.assertEqual(encode(restored,self.contract),encode(self.zero,self.contract))
    def test_actual_child_kill_after_HEAD_exposes_only_committed_frontier(self):
        store,p,candidate=self.killed_child('after_HEAD');self.assertNotEqual(store.head(),p['HEAD_sha256'])
        restored,_=store.read('g_child_after_HEAD',store.head());self.assertEqual(encode(restored,self.contract),encode(candidate,self.contract))
        restored0,_=store.read('g_zero',store.head());self.assertEqual(encode(restored0,self.contract),encode(self.zero,self.contract))
    def test_only_initial_named_infinity_sentinel_allowed(self):
        self.assertEqual(decode(encode(self.zero,self.contract),self.contract)['audit']['statistics']['minimum_observed_gap'],dec.Decimal('Infinity'))
        s=self.state(1);s['audit']['statistics']['minimum_observed_gap']=dec.Decimal('Infinity')
        with self.assertRaisesRegex(ValueError,'nonfinite'):encode(s,self.contract)
        with self.assertRaises(ValueError):pack(dec.Decimal('NaN'))
        with self.assertRaises(ValueError):pack(object())

if __name__=='__main__':unittest.main()
