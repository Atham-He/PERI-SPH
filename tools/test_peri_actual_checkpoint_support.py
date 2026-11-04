"""New parser/import boundaries, not another material/raw review or solver run."""
import copy,hashlib,json,struct,tempfile,unittest
from pathlib import Path
from peri_actual_checkpoint_support import trace_index,compare_new_branch,initial_contract,selected,file_tree,assert_preserved,restored_noadvance_metadata
from peri_bounded_session_archive_audit import encode
from peri_bounded_mechanical_stream_audit import packed


class CompletePacketTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        self.rows=[dict(type='identity',pid=123,ppid=122,start_ticks='42',boot_id='synthetic',argv=['synthetic']),
                   dict(type='header',initial_state=dict(time=0.))]
        for step in range(1,41):
            self.rows.append(dict(type='step',committed=True,local_step_gates_passed=True,candidate_ready=True,
                reason='all_applicable_step_gates_passed',state=dict(step=step,time=step*.01),
                interval=dict(step=step-1,old_time=(step-1)*.01,new_time=step*.01,actual_dt=.01),
                solve=dict(trace=[dict(original_nonsymmetric_J=[[1.,2.],[0.,3.]],all_point_T=[1.,-0.])]),
                accepted_physical_endpoint='same'))
        self.rows.append(dict(type='end',completed=True,accepted_steps=40,last_accepted_time=.4,actual_solver_calls=40))
        self.base=self.write('base',self.rows);self.qualified=trace_index(self.base,1,40)
    def tearDown(self):self.temp.cleanup()
    def write(self,name,rows):
        p=self.root/name;p.write_bytes(b''.join(json.dumps(r,separators=(',',':'),allow_nan=False).encode()+b'\n' for r in rows));return p
    def resumed(self):
        return [dict(self.rows[0],pid=456),dict(type='resume_header',selected_generation=8,accepted_step=8,accepted_time=.08,new_initial_state=False)]+copy.deepcopy(self.rows[10:-1])+[dict(self.rows[-1],actual_solver_calls=32)]
    def test_all_full_newton_packets_match_despite_distinct_process_identity(self):
        branch=trace_index(self.write('resume',self.resumed()),9,40,resume_generation=8)
        result=compare_new_branch(self.qualified,branch)
        self.assertEqual(result['different_steps'],[]);self.assertEqual(result['compared_complete_steps'],32)
    def test_changed_trial_tangent_with_identical_endpoints_is_not_imported(self):
        rows=self.resumed();rows[6]['solve']['trace'][0]['all_point_T'][0]+=1
        branch=trace_index(self.write('changed_T',rows),9,40,resume_generation=8)
        result=compare_new_branch(self.qualified,branch)
        self.assertEqual(result['different_steps'],[13]);self.assertFalse(result['original_full_reference_not_repeated'])
    def test_dropped_newton_packet_is_not_hidden_by_equal_endpoint(self):
        rows=self.resumed();rows[2]['solve']['trace']=[]
        branch=trace_index(self.write('dropped_trial',rows),9,40,resume_generation=8)
        self.assertEqual(compare_new_branch(self.qualified,branch)['different_steps'],[9])
    def test_seed_prefix_reuses_only_its_exact_full_packets(self):
        rows=copy.deepcopy(self.rows[:10])+[dict(self.rows[-1],accepted_steps=8,last_accepted_time=.08,actual_solver_calls=8)]
        seed=trace_index(self.write('seed',rows),1,8)
        self.assertEqual(compare_new_branch(self.qualified,seed)['compared_complete_steps'],8)
    def test_initial_header_changes_are_rejected(self):
        rows=copy.deepcopy(self.rows);rows[1]['undeclared_mass_scaling']=True
        branch=trace_index(self.write('bad_header',rows),1,40)
        with self.assertRaises(AssertionError):compare_new_branch(self.qualified,branch)
    def test_structural_refusals(self):
        variants={}
        rows=copy.deepcopy(self.rows);rows[3]['state']['step']=1;variants['duplicate_step']=rows
        rows=copy.deepcopy(self.rows);rows[3]['interval']['old_time']=0.;variants['wrong_clock']=rows
        rows=copy.deepcopy(self.rows);rows[3]['committed']=False;variants['uncommitted']=rows
        rows=copy.deepcopy(self.rows);rows[3]['type']='checkpoint';variants['unexpected_metadata']=rows
        rows=copy.deepcopy(self.rows);rows[-1]['actual_solver_calls']=39;variants['wrong_solve_count']=rows
        rows=copy.deepcopy(self.rows);rows[-1]['last_accepted_time']=.39;variants['wrong_terminal']=rows
        rows=copy.deepcopy(self.rows);rows.append(dict(type='step'));variants['extra_row']=rows
        for name,rows in variants.items():
            with self.subTest(name=name),self.assertRaises((AssertionError,RuntimeError,KeyError)):
                trace_index(self.write(name,rows),1,40)
    def test_partial_line_and_missing_terminal_are_rejected(self):
        p=self.root/'partial';p.write_bytes(self.base.read_bytes()[:-1])
        with self.assertRaises(AssertionError):trace_index(p,1,40)
        with self.assertRaises(RuntimeError):trace_index(self.write('missing',self.rows[:-1]),1,40)
    def test_noadvance_requires_owned_restore_and_preserves_explicit_generation(self):
        rows=[dict(self.rows[0]),dict(type='resume_header',selected_generation=7,accepted_step=0,accepted_time=0.,new_initial_state=False),
              dict(type='restored_noadvance',selected_generation=7,accepted_step=0,accepted_time=0.,actual_solver_calls=0,new_physical_steps=0,
                   owned_state_installed_and_resaved=True,full_restart_qualified=False,archive_sha256='ab'*32)]
        arguments=dict(generation=7,step=0,time=0.,archive_sha256='ab'*32)
        self.assertEqual(restored_noadvance_metadata(self.write('owned_noadvance',rows),**arguments)['actual_solver_calls'],0)
        for name,change in [('read_only',dict(type='validated_noadvance')),('retagged_generation',dict(selected_generation=1)),
                            ('hidden_advance',dict(actual_solver_calls=1)),('not_installed',dict(owned_state_installed_and_resaved=False))]:
            bad=copy.deepcopy(rows);bad[-1].update(change)
            with self.subTest(name=name),self.assertRaises(AssertionError):restored_noadvance_metadata(self.write(name,bad),**arguments)
    def test_selected_generation_and_complete_physical_bytes(self):
        # Immutable t0 fixture is input for this NEW parser-format test only.
        # No old native, material computation or previous raw review is called.
        fixture=Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material/engineering-actual-checkpoint-driver-preflight-v4')
        if not fixture.exists():self.skipTest('bounded zero-advance fixture not mounted')
        value=initial_contract(fixture);archive=(fixture/'initial_candidate.bin').read_bytes();physical=value['physical'];trusted=value['trusted']
        store=self.root/'store';branch=store/'g_fixture';branch.mkdir(parents=True)
        (branch/'state.bin').write_bytes(archive);(branch/'prefix.phy').write_bytes(physical)
        fields=dict(SelectionSchema=packed(2,(),1),Branch=packed(1,(9,),list(b'g_fixture')),
            ArchiveSHA256=packed(1,(32,),list(hashlib.sha256(archive).digest())),
            PhysicalPrefixSHA256=packed(1,(32,),list(hashlib.sha256(physical).digest())),ContractSHA256=packed(1,(32,),list(trusted)))
        (store/'selected.bin').write_bytes(encode(fields,1));expected=dict(archive=archive,physical=physical)
        self.assertTrue(selected(store,'selected.bin',1,expected,trusted)['full_candidate_byte_equal'])
        before=file_tree(store);self.assertEqual(before,file_tree(store));assert_preserved(before,file_tree(store))
        # Recomputed selection hashes do not make altered accepted bytes equal.
        altered=bytearray(physical);altered[-1]^=1;(branch/'prefix.phy').write_bytes(altered)
        fields['PhysicalPrefixSHA256']=packed(1,(32,),list(hashlib.sha256(altered).digest()));(store/'selected.bin').write_bytes(encode(fields,1))
        with self.assertRaises(AssertionError):selected(store,'selected.bin',1,expected,trusted)
        with self.assertRaises(AssertionError):assert_preserved(before,file_tree(store))


if __name__=='__main__':unittest.main()
