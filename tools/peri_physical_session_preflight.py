"""NEW explicit v2 complete-session ZERO-PD-step engineering qualification."""
import argparse,hashlib,json,shutil,subprocess,sys,traceback
from pathlib import Path
from peri_archive_component_build import prepare,compile_and_freeze
from peri_run_body_archive_engineering import reserves,stamp,identity
from peri_actual_checkpoint_campaign import Execution,verified_execution_component
from peri_bounded_body_archive_audit import record
from peri_canonical_archive_audit import parse
from peri_bounded_session_archive_audit import encode
from peri_bounded_mechanical_stream_audit import packed
from peri_elastic_normal import read,save,verify

def assess(root,legacy,build,native):
    data=(root/'initial.phy').read_bytes();trusted=native['contract_sha256'];archive=(root/'initial_candidate.bin').read_bytes()
    _,fields=parse(archive);_,old=parse((legacy/'initial_candidate.bin').read_bytes())
    body={k:v for k,v in fields.items() if k.startswith('body/')}
    assert body=={k:v for k,v in old.items() if k.startswith('body/')}
    assert data[:48]==b'PERIPHY2'+(2).to_bytes(8,'little')+bytes.fromhex(trusted)
    length=int.from_bytes(data[48:56],'little');assert len(data)==56+length
    _,physical=parse(data[56:]);old_physical=(legacy/'initial.phy').read_bytes();_,old_record=parse(old_physical[56:])
    expected_physical=dict(old_record)
    for key,value in [('RecordSchema',2),('OldTick',0),('AcceptedTick',0),('TimeDenominator',1000000)]:
        expected_physical[key]=packed(2,(),value)
    assert physical==expected_physical and data[56:]==encode(expected_physical,1)
    old_contract={k.removeprefix('contract/'):v for k,v in old.items() if k.startswith('contract/')}
    expected=dict(old_contract)
    for name in ['NormalBaseSteps','ProtocolKind','NominalDt','NormalStroke','TwoTangentStrokes']:del expected[name]
    for name,value in [('Schema',2),('PhysicalStreamSchema',2),('TimeDenominator',1000000),('NominalTicks',5000)]:expected[name]=packed(2,(),value)
    for k,(n,t) in enumerate(zip([160000,240000,400000],[20000,60000,80000])):
        expected['NormalEvent'+str(k)]=packed(2,(),n);expected['TangentEvent'+str(k)]=packed(2,(),t)
    for name,value in [('PhysicalNormalStroke',.01),('PhysicalTangentStroke0',.0006),('PhysicalTangentStroke1',-.0004)]:expected[name]=packed(3,(),value)
    expected['CheckpointSteps/count']=packed(2,(),4);expected['CheckpointSteps']=packed(2,(4,),[16,40,64,80])
    static={k.removeprefix('body/'):v for k,v in body.items() if not k.startswith('body/accepted/')}
    fingerprints=dict(SourceManifest=build['source_manifest']['sha256'],ELF=build['ELF']['sha256'],
                      BuildRuntimeManifest=build['build_runtime_manifest']['sha256'],Input=record(root/'original.input')['sha256'],
                      StaticBody=hashlib.sha256(encode(static,1)).hexdigest())
    for name,value in fingerprints.items():expected[name+'SHA256']=packed(1,(32,),list(bytes.fromhex(value)))
    actual_contract={k.removeprefix('contract/'):v for k,v in fields.items() if k.startswith('contract/')}
    assert actual_contract==expected
    assert (root/'contract.bin').read_bytes()==encode(expected,1)
    assert hashlib.sha256(encode(expected,1)).hexdigest()==trusted
    expected_candidate={k:v for k,v in old.items() if not k.startswith(('contract/','output/'))}
    expected_candidate.update({'contract/'+k:v for k,v in expected.items()})
    expected_candidate['accepted_algorithm/AcceptedTick']=packed(2,(),0)
    for name,value in [('Schema',1),('Records',1),('Bytes',len(data)),('LastStep',0)]:expected_candidate['output/'+name]=packed(2,(),value)
    expected_candidate['output/LastTime']=packed(3,(),0.)
    for key,raw in [('HeaderSHA256',data[:48]),('PrefixSHA256',data)]:
        expected_candidate['output/'+key+'/count']=packed(2,(),32)
        expected_candidate['output/'+key]=packed(1,(32,),list(hashlib.sha256(raw).digest()))
    assert fields==expected_candidate and archive==encode(expected_candidate,1)
    assert len(fields)==149 and len(physical)==142 and native['new_physical_steps']==0
    return dict(status='v2_complete_session_initial149fields_physical142fields_independently_byte_matched',
                candidate_fields=149,physical_record_fields=142,candidate_bytes=len(archive),physical_prefix_bytes=len(data),
                contract_sha256=trusted,new_PD_steps=0,initial_body_identical_to_qualified_legacy_t0=True,
                full_nonzero_restart_qualified=False,S4_complete=False,calibration_allowed=False)

def fixture(root,trusted):
    _,fields=parse((root/'initial_candidate.bin').read_bytes());archive=encode(fields,7);physical=(root/'initial.phy').read_bytes()
    store=root/'new_v2_zero_store';branch=store/'g_fixture';branch.mkdir(parents=True)
    (branch/'state.bin').write_bytes(archive);(branch/'prefix.phy').write_bytes(physical)
    selected=encode(dict(SelectionSchema=packed(2,(),1),Branch=packed(1,(9,),list(b'g_fixture')),
                         ArchiveSHA256=packed(1,(32,),list(hashlib.sha256(archive).digest())),
                         PhysicalPrefixSHA256=packed(1,(32,),list(hashlib.sha256(physical).digest())),
                         ContractSHA256=packed(1,(32,),list(bytes.fromhex(trusted)))),7)
    (store/'HEAD').write_bytes(selected);(store/'selected_input.bin').write_bytes(selected)
    before={str(p.relative_to(store)):record(p) for p in store.rglob('*') if p.is_file()}
    save(root/'new_v2_zero_fixture_receipt.json',dict(accepted_step=0,logical_generation=7,inputs=before,new_PD_steps=0))
    return ['--restore-noadvance',str(root),str(store),trusted,'selected_input.bin','7',hashlib.sha256(selected).hexdigest(),'v2zero'],before,archive,physical

def main(root):
    repo=Path(__file__).resolve().parents[1];legacy=root.parent/'engineering-actual-checkpoint-driver-preflight-v5'
    assert root.parent==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists(),'no redispatch'
    reserves(root);assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    component=verified_execution_component(root.parent);root.mkdir(mode=0o700)
    try:
        stage='tests/3d_examples/extra_3d_PD_framework_tests/'
        files=[stage+n for n in ['bounded_gel_physical_checkpoint_probe.cpp','bounded_gel_physical_time_observer.h',
            'coulomb_contact_probe_io.h','bounded_gel_probe_io.h']]
        files += ['tools/'+n for n in ['peri_physical_session_preflight.py','peri_archive_component_build.py','peri_run_body_archive_engineering.py',
            'peri_actual_checkpoint_campaign.py','peri_bounded_body_archive_audit.py','peri_canonical_archive_audit.py',
            'peri_bounded_session_archive_audit.py','peri_bounded_mechanical_stream_audit.py','peri_physical_time_friction_audit.py',
            'peri_bounded_friction_audit.py','peri_bounded_gel_audit.py','test_peri_physical_time_friction.py']]
        prepared=prepare(root,repo,files,stage+'bounded_gel_physical_checkpoint_probe.cpp','physical_time_CP')
        production=repo/'build/peri-latest-20261103/extra_3d_PD_gel_compression/bin/extra_3d_PD_gel_compression'
        assert record(production)['sha256']=='bf398b1415cbded594c4164760989dc91c99f7935ec3c4f047eb7c8bbee88b69'
        baseline=record(legacy/'initial_candidate.bin')
        assert baseline['sha256']=='b920b706b7e522d93e5b84c4edb73d3a8390eb141370d4f0296b3dfd64b3a525'
        anchors=[baseline,record(legacy/'initial.phy'),record(legacy/'contract.bin'),record(legacy/'independent_preflight.json'),record(production)]
        mutations=[('legal_tick_wrong_history','accepted_algorithm/AcceptedTick',1,2,'accepted rational tick/history mismatch'),
                   ('denominator_changed','contract/TimeDenominator',2000000,2,'external trusted contract mismatch'),
                   ('normal_amplitude_changed','contract/PhysicalNormalStroke',.0101,3,'external trusted contract mismatch'),
                   ('tangent_amplitude_changed','contract/PhysicalTangentStroke0',.0007,3,'external trusted contract mismatch'),
                   ('wrong_next_phase','accepted_algorithm/NextPhase',1,2,'exact accepted protocol clock'),
                   ('wrong_output_time','output/LastTime',.001,3,'structural output cursor')]
        plan=dict(prepared_at=stamp(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
                  source_snapshots=prepared['source_snapshots'],compile_command=prepared['compile_command'],link_command=prepared['link_command'],
                  compile_budget=dict(CPU_seconds=180,wall_seconds=240,AS_GiB=2),native_budget=dict(CPU_seconds=30,wall_seconds=45,AS_GiB=1,threads=1),
                  anchors=anchors,execution_component=component,
                  declared_positive_modes=['capabilities','initialize','inspect_detached','legacy_t0_readonly_roundtrip','owned_zero_restore'],
                  mutations=mutations,expected_candidate_fields=149,expected_physical_fields=142,
                  scope='NEW explicit v2 complete schema/actual owned zero-step only;not nonzero physical CP or time convergence',
                  physical_dispatch_allowed=False,new_PD_steps=0,file_bound_MiB=256,total_bound_GiB=1,
                  host_F_reserve_GiB=120,output_reserve_GiB=150,
                  failure_policy='Preserve precise raw/returncode/failure and stop dependent controls;never redispatch this state or relax gates.',S4_complete=False,calibration_allowed=False)
        save(root/'plan.json',plan);save(root/'state.json',dict(status='v2_zero_step_preflight_attempted',identity=identity(__import__('os').getpid()),started_at=stamp()))
        build=compile_and_freeze(prepared)
        with (root/'original.input').open('x') as f:
            f.write('PD_PHYSICAL_TIME_V2 0 8 1000000 5000 160000 240000 400000 20000 60000 80000 0.01 0.0006 -0.0004 0.22 0.22 1 5 1 0.5 0.3 0.1 1.25 0.015625 4 16 40 64 80\n')
        for a,b in [('initial_candidate.bin','legacy_candidate.bin'),('initial.phy','legacy.phy')]:shutil.copy2(legacy/a,root/b)
        trusted_legacy=hashlib.sha256((legacy/'contract.bin').read_bytes()).hexdigest()
        inputs=[record(root/n) for n in ['original.input','legacy_candidate.bin','legacy.phy']]
        save(root/'frozen_execution.json',dict(plan=record(root/'plan.json'),build=build,inputs=inputs,anchors=anchors))
        executor=Execution(root,plan);results=[]
        def execute(label,args,code=0):
            verify(build['support']+build['source_snapshots']+inputs+anchors+[build['ELF'],build['source_manifest'],build['build_runtime_manifest']])
            r=executor.run(label,[build['ELF']['path']]+args,wall=45,cpu=30,gib=1,expected_code=code);results.append(record(root/(label+'.attempt.json')))
            rows=[json.loads(line) for line in (root/(label+'.stdout')).read_text().splitlines()]
            assert rows[0]['type']=='identity' or label=='capabilities'
            if label!='capabilities':assert all(rows[0][k]==v for k,v in r['birth_identity'].items()) and rows[0]['argv']==r['command']
            if code==0:assert r['stderr']['size']==0
            verify(build['support']+build['source_snapshots']+inputs+anchors)
            return rows[-1],r
        c,_=execute('capabilities',['--capabilities']);assert c['physical_time_schema']==2 and not c['full_restart_qualified']
        n,_=execute('initialize',['--initialize',str(root)]);assessment=assess(root,legacy,build,n);trusted=assessment['contract_sha256']
        inputs += [record(root/n) for n in ['initial_candidate.bin','initial.phy','contract.bin']]
        n,_=execute('inspect_detached',['--inspect-detached',str(root),'initial_candidate.bin','1',trusted]);assert n['session_schema']==2 and n['accepted_tick']==0 and n['physical_solver_calls']==0
        n,_=execute('legacy_t0_readonly_roundtrip',['--roundtrip-legacy-t0',str(root),'legacy_candidate.bin','1',trusted_legacy]);assert n['session_schema']==1 and not n['owned_state_installed'] and n['physical_solver_calls']==0
        args,before,archive,physical=fixture(root,trusted)
        n,_=execute('owned_zero_restore',args)
        assert n['owned_state_installed_and_resaved'] is True and n['actual_solver_calls']==n['new_physical_steps']==n['accepted_step']==0
        store=root/'new_v2_zero_store'
        for branch in ['g_fixture','g_v2zero_resume','g_v2zero_s0']:
            assert (store/branch/'state.bin').read_bytes()==archive and (store/branch/'prefix.phy').read_bytes()==physical
        for name,r in before.items():
            if name!='HEAD':assert record(store/name)==r
        _,fields=parse((root/'initial_candidate.bin').read_bytes())
        rejects=[]
        for name,key,value,tag,cause in mutations:
            altered=dict(fields);altered[key]=packed(tag,(),value);data=encode(altered,1)
            target=root/(name+'.bin');target.write_bytes(data);inputs.append(record(target))
            # Canonical SHA has been re-created by the independent encoder.
            _,parsed=parse(data);assert parsed==altered
            unchanged={str(p.relative_to(store)):record(p) for p in store.rglob('*') if p.is_file()}
            n,r=execute(name,['--inspect-detached',str(root),target.name,'1',trusted],1)
            assert r['native_failure_witnesses'] and cause in r['native_failure_witnesses'][0]['message']
            w=r['native_failure_witnesses'][0];assert w['actual_solver_calls']==0 and w['append_witness'] is None
            assert unchanged=={str(p.relative_to(store)):record(p) for p in store.rglob('*') if p.is_file()}
            rejects.append(dict(name=name,expected_cause=cause,input=record(target),no_live_or_output_mutation=True))
        assessment.update(recorded_at=stamp(),source_commit=plan['source_commit'],plan=record(root/'plan.json'),
                          frozen_execution=record(root/'frozen_execution.json'),phases=results,clock_semantic_rejections=rejects,
                          zero_step_owned_install_resave_bytes_identical=True,legacy_t0_readonly_full_candidate_and_prefix_matched=True,
                          actual_compiler_dependency_count=build['actual_compiler_dependency_count'],new_PD_steps=0,
                          full_nonzero_restart_qualified=False,time_convergence_qualified=False)
        save(root/'independent_preflight.json',assessment)
        save(root/'terminal.json',dict(status=assessment['status'],return_code=0,finished_at=stamp(),receipt=record(root/'independent_preflight.json'),new_PD_steps=0,S4_complete=False))
        print(json.dumps(dict(status=assessment['status'],candidate_fields=149,physical_fields=142,correct_SHA_clock_rejections=len(rejects),new_PD_steps=0)))
    except BaseException as error:
        save(root/'failure.json',dict(recorded_at=stamp(),failure=repr(error),traceback=traceback.format_exc(),new_PD_steps=0,full_nonzero_restart_qualified=False));raise

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();main(a.root)
