"""Prospectively frozen NEW v2 history campaign; prepare/run each one-shot.

Qualification is fixed-clock bounded history and original physical-equation
consistency only. The 1% trajectory loss and 2% refinement S4 gates remain
separate, explicitly unqualified; raw metrics are mandatory and never relabeled.
"""
import argparse,hashlib,json,math,os,resource,shutil,subprocess,sys,traceback
from pathlib import Path
# Packaged executable entry must bind its frozen/reference before imports.
if (Path(__file__).resolve().parent/'reference').is_dir():
    sys.path.insert(0,str(Path(__file__).resolve().parent/'reference'))
from peri_actual_checkpoint_campaign import Execution,verified_execution_component

def campaign_case(header_case):
    import copy
    from peri_physical_time_friction_audit import schedule_grid
    case=copy.deepcopy(header_case)
    case.update(name='physical_time_v2_complete_history_80_fixed_end04',expected_code=0,
                require_statistics=['nonzero_tangent0_rows','nonzero_tangent1_rows'])
    assert case['steps']==80 and case['dt']==.005 and case['checkpoint_steps']==[16,40,64,80]
    assert schedule_grid(case)[-1]==400000 and len(schedule_grid(case))==81
    return case

def require_ready(header,initial,semantic,binding):
    assert header['status']=='existing_actual_v2_header_independently_qualified_by_NEW_path_admitted_reader_ZERO_steps'
    assert header['physical_solver_calls']==0 and header['old_failure_preserved'] is True
    assert initial['maximum_normalized_error']<=1e-10 and initial['new_PD_steps']==0
    assert initial['planned_steps']==80 and initial['fixed_end_time']==.4 and initial['phase_step_counts']==[32,16,32]
    assert initial['accepted_initial_body_all_bytes_match'] and len(initial['clock_header_negative_controls'])==6
    assert semantic['status']=='physical_time_v2_readonly_CP_correct_SHA_semantic_refusals_and_valid_binding_passed'
    assert semantic['semantic_rejections']==22 and semantic['positive_binding_controls']==1 and semantic['native_solver_calls']==0
    assert semantic['qualified_physical_nonzero_seed'] is False and semantic['no_private_or_live_or_output_changes']
    assert binding['status']=='new_v2_full_reference_wrapper_fixedpure95_binding_and_explicit_v2_module_admission_qualified_ZERO_material'
    q=binding['observation'];assert q['precision']==95 and q['C_decimal_loaded'] is False and q['material_calculations']==0

def prepare(root):
    from peri_elastic_normal import read,record,verify,save,reserves,stamp
    from peri_physical_checkpoint_support import initial_contract
    from peri_physical_checkpoint_negative import specification
    repo=Path(__file__).resolve().parents[1];base=root.parent
    native=base/'engineering-s4-physical-header-v1';header=base/'analysis-s4-header-v3-admission-path'
    semantic=base/'engineering-s4-v2-complete-semantic-v1';binding=base/'engineering-s4-v2-reference-binding-v1'
    assert base==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists(),'unused new campaign root required'
    h=read(header/'independent_review.json');hi=read(header/'header_independent_review.json')
    sem=read(semantic/'semantic_rejection_review.json');bound=read(binding/'independent_review.json');require_ready(h,hi,sem,bound)
    anchors=[]
    for p in [header,semantic,binding]:
        terminal=read(p/'terminal.json');assert terminal['return_code']==0;verify([terminal['receipt']])
        anchors += [record(p/'terminal.json'),terminal['receipt'],record(p/'plan.json')]
    anchors += [record(header/'header_independent_review.json'),record(semantic/'semantic_rejection_review.json'),
                record(native/'failure.json'),record(base/'analysis-s4-header-v2-existing-raw/failure.json')]
    for receipt in sem['results']:
        a=read(receipt['path']);anchors += [receipt,a['stdout'],a['stderr']]
    component=verified_execution_component(base);build=read(native/'frozen_execution.json')['build']
    # Bind the EXACT modules that were newly preflighted, not just a status.
    for name in ['peri_physical_checkpoint_negative.py','peri_physical_checkpoint_support.py']:
        assert record(repo/'tools'/name)['sha256']==record(semantic/'frozen/reference'/name)['sha256']
    for name in ['peri_physical_checkpoint_reference.py','peri_physical_friction_negative.py','peri_physical_time_friction_audit.py']:
        assert record(repo/'tools'/name)['sha256']==record(binding/'frozen/reference'/name)['sha256']
    resources=reserves(base)
    mem=int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))*1024
    assert mem>=8*2**30,'8GiB memory admission'
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    case=campaign_case(read(native/'plan.json')['case']);trust=h['initial_complete_bytes']['contract_sha256']
    root.mkdir(mode=0o700);frozen=root/'frozen';reference=frozen/'reference';reference.mkdir(parents=True)
    artifacts=anchors+component['immutable_evidence']+build['support']+build['source_snapshots']+[build['ELF'],record(native/'frozen_execution.json')]
    artifacts += read(native/'plan.json')['software']
    for source in sorted((repo/'tools').glob('*.py')):
        target=reference/source.name;shutil.copy2(source,target);artifacts += [record(source),record(target)]
    for name in ['peri_physical_checkpoint_campaign_v2.py','peri_physical_checkpoint_reference.py','peri_fixed_python_raw_audit.py']:
        source=repo/'tools'/name;target=frozen/name;shutil.copy2(source,target);artifacts += [record(source),record(target)]
    source=binding/'frozen/_pydecimal.py';target=frozen/'_pydecimal.py';shutil.copy2(source,target);artifacts += [record(source),record(target)]
    for name in ['source_manifest.json','build_runtime_manifest.json','original.input','initial_candidate.bin','initial.phy','contract.bin']:
        source=native/name;target=root/name;shutil.copy2(source,target);artifacts += [record(source),record(target)]
    exe=frozen/'physical_PD_checkpoint';shutil.copy2(build['ELF']['path'],exe);artifacts.append(record(exe))
    initial_contract(root,case,bytes.fromhex(trust));negative=specification();assert negative==sem['specification']
    plan=dict(prepared_at=stamp(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
              native_producer_source=read(native/'plan.json')['source_commit'],native_producer_ELF=build['ELF'],executable=str(exe),
              case=case,trusted_contract=trust,artifacts=artifacts,resources=resources,negative_specification=negative,execution_component=component,
              native_CPU_seconds=1800,native_wall_seconds=1800,native_AS_GiB=1,threads=1,file_bound_GiB=4,total_output_GiB=4,
              host_reserve_GiB=120,output_reserve_GiB=150,reference_AS_GiB=8,reference_max_CPU_seconds=36000,
              reference_cost_formula='max(7200,ceil(2*actual_native_trial_material_points+2000))CPU;wall=CPU+600;freeze exact admission BEFORE new reference;never extend after failure.',
              cost_basis=dict(old_native40_seconds=212.721386824,old_raw40_bytes=259209107,old_reference_points=6534,
                              old_reference_CPU_seconds=7528.674968096,prospective80_reference_CPU_seconds=28136,
                              native_1800_seconds_is_new_batch_budget_not_retry=True,reference160_estimate54272_NOT_admitted=True),
              frozen_entry_preflight=dict(label='entry_binding',CPU_seconds=30,wall_seconds=45,AS_GiB=1,physical_solver_calls=0),
              segments=dict(continuous=[1,80],split_first=[1,16],split_resume=[17,80],rollback40=[41,80]),
              owned_noadvance_steps=[16,40,64],checkpoint_steps=[16,40,64,80],planned_total_accepted_solves=200,
              normalized_physical_gate=1e-10,original_precision_gate=1e-12,native_Newton_tolerance=1e-12,dense_penetration_stroke_gate=.01,
              full_candidate_and_physical_prefix_byte_equality=True,
              convergence_gates=dict(force_refinement_relative=.02,whole_interval_penetration_stroke=.01,algorithmic_loss_over_positive_input=.01),
              qualification_scope='NEW fixed-clock protocol/full physical-equation reference and complete nonzero history restart only;not time/space/horizon/inertia/rate convergence. Emit actual loss/input,do not call >1% S4 passed.',
              failure_policy='Stop dependent phases and preserve exact raw/errno/append/debug/durable cursor and returncode;no replay/reset/fallback/reference refresh or changes to physical mass,force,speed,gates.',
              S3_prior_qualification_preserved=True,full_restart_qualified=False,S4_complete=False,calibration_allowed=False)
    plan['critical_runtime']=[record(exe)]+[record(root/n) for n in ['source_manifest.json','build_runtime_manifest.json','original.input']]
    plan['critical_runtime'] += [x for x in build['support'] if '.so' in Path(x['path']).name]
    verify(artifacts);save(root/'plan.json',plan)
    print(json.dumps(dict(status='NEW_physical_time_v2_complete_history_campaign_prepared_NOT_dispatched',root=str(root),plan=record(root/'plan.json'))))

def run(root):
    sys.path.insert(0,str(root/'frozen/reference'))
    from peri_elastic_normal import read,record,verify,save,identity,stamp
    from peri_physical_checkpoint_support import derive_continuous,selected,file_tree,assert_preserved,restored_noadvance_metadata,trace_packets
    plan=read(root/'plan.json');verify_entry_binding(root);plan_record=record(root/'plan.json');assert not (root/'state.json').exists(),'never redispatch campaign'
    resource.setrlimit(resource.RLIMIT_AS,(8*2**30,8*2**30))
    save(root/'state.json',dict(status='physical_time_v2_complete_history_campaign_attempted',started_at=stamp(),identity=identity(),plan=record(root/'plan.json')))
    execution=Execution(root,plan);results=[];failure=None;trusted=bytes.fromhex(plan['trusted_contract']);exe=plan['executable']
    def verify_all():verify([plan_record]+plan['artifacts'])
    def native(label,args):
        verify_all();result=execution.run(label,[exe]+args,wall=plan['native_wall_seconds'],cpu=plan['native_CPU_seconds'],gib=plan['native_AS_GiB'])
        with Path(result['stdout']['path']).open() as stream:meta=json.loads(next(stream))
        assert meta['type']=='identity' and meta['argv']==result['command']
        assert all(meta[k]==v for k,v in result['birth_identity'].items()) and result['stderr']['size']==0
        if result['child_identity']:assert all(meta[k]==v for k,v in result['child_identity'].items())
        verify_all();return result
    def head(store):return hashlib.sha256((store/'HEAD').read_bytes()).hexdigest()
    def reference(label,first,last):
        attempt=read(root/(label+'.attempt.json'));points=0
        for row,_ in trace_packets(attempt['stdout']['path'],plan['case'],first,last,resume_generation=None if first==1 else first-1):
            if row['type']=='step':points+=33*len(row['solve']['trace'])
        cpu=max(7200,math.ceil(2*points+2000));assert cpu<=plan['reference_max_CPU_seconds'],'new reference cost exceeds prospective cap;do not attempt'
        save(root/(label+'.reference_admission.json'),dict(recorded_at=stamp(),attempt=record(root/(label+'.attempt.json')),
            actual_native_trial_material_points=points,CPU_seconds=cpu,wall_seconds=cpu+600,AS_GiB=8,
            prior_measured_cost=plan['cost_basis'],formula=plan['reference_cost_formula']))
        execution.run(label+'.reference',['/usr/bin/python3',str(root/'frozen/peri_physical_checkpoint_reference.py'),'--root',str(root),'--label',label],wall=cpu+600,cpu=cpu,gib=8)
        receipt=read(root/(label+'.independent_review.json'))
        assert receipt['status'] in ['new_physical_time_CP_continuous_full_reference_passed_NOT_full_restart','new_physical_time_CP_branch_full_trial_packets_identical_math_evidence_imported','new_physical_time_CP_branch_distinct_trial_context_independently_passed']
        assert 'failure' not in receipt
        if label=='continuous':
            from peri_physical_checkpoint_reference import require_full_reference
            require_full_reference(receipt,plan['case'])
        verify_all();return record(root/(label+'.independent_review.json'))
    def bytes_assessment(label,store,steps,derived):
        checks=[selected(store,'selected_'+label+'_s'+str(step)+'.bin',step,derived['checkpoints'][step],trusted) for step in steps]
        result=dict(status='new_physical_time_v2_all149_candidate_and142_physical_output_bytes_matched',recorded_at=stamp(),native_attempt=record(root/(label+'.attempt.json')),
            checkpoints=checks,full_candidate_byte_equal=True,full_physical_prefix_byte_equal=True,full_restart_qualified=False)
        save(root/(label+'.bytes_review.json'),result);return record(root/(label+'.bytes_review.json'))
    try:
        verify_all();continuous=root/'continuous.store';split=root/'split.store';continuous.mkdir();split.mkdir()
        native('continuous',['--fresh',str(root),str(continuous),plan['trusted_contract'],'80','continuous'])
        derived=derive_continuous(root,root/'continuous.stdout',plan['case'],trusted);results.append(bytes_assessment('continuous',continuous,[16,40,64,80],derived))
        results.append(reference('continuous',1,80))
        for step in plan['owned_noadvance_steps']:
            label='na'+str(step);before=file_tree(continuous)
            save(root/(label+'.inputs.json'),dict(recorded_at=stamp(),selected_input=record(continuous/('selected_continuous_s'+str(step)+'.bin')),expected_HEAD=record(continuous/'HEAD'),prior_files=before))
            attempt=native(label,['--restore-noadvance',str(root),str(continuous),plan['trusted_contract'],'selected_continuous_s'+str(step)+'.bin',str(step),head(continuous),label])
            expected=derived['checkpoints'][step];proof=selected(continuous,'selected_'+label+'_s'+str(step)+'.bin',step,expected,trusted)
            metadata=restored_noadvance_metadata(attempt['stdout']['path'],generation=step,step=step,time=derived['index']['steps'][step-1]['time'],archive_sha256=hashlib.sha256(expected['archive']).hexdigest())
            assert_preserved(before,file_tree(continuous),allowed_changes={'HEAD'})
            save(root/(label+'.independent_noadvance.json'),dict(recorded_at=stamp(),proof=proof,metadata=metadata,prior_files_preserved=True,new_physical_steps=0))
            results.append(record(root/(label+'.independent_noadvance.json')))
        native('split_first',['--fresh',str(root),str(split),plan['trusted_contract'],'16','split_first'])
        results.append(bytes_assessment('split_first',split,[16],derived));results.append(reference('split_first',1,16))
        for label,step,selection in [('split_resume',16,'selected_split_first_s16.bin'),('rollback40',40,'selected_split_resume_s40.bin')]:
            before=file_tree(split);save(root/(label+'.inputs.json'),dict(recorded_at=stamp(),selected_input=record(split/selection),expected_HEAD=record(split/'HEAD'),prior_files=before))
            native(label,['--resume',str(root),str(split),plan['trusted_contract'],selection,str(step),head(split),'80',label])
            assert_preserved(before,file_tree(split),allowed_changes={'HEAD'})
            results.append(bytes_assessment(label,split,[40,64,80] if step==16 else [64,80],derived))
            results.append(reference(label,step+1,80))
        from peri_physical_checkpoint_negative import run_controls
        negative_review=run_controls(root,plan,execution,derived);negative=read(negative_review['path'])
        assert negative['qualified_physical_nonzero_seed'] and negative['accepted_step']==16 and negative['semantic_rejections']==22
        results.append(negative_review)
        verify_all();status='new_physical_time_v2_bounded_PD_complete_history_restart_independently_qualified'
    except Exception as error:
        failure=dict(error=repr(error),traceback=traceback.format_exc());status='new_physical_time_v2_complete_history_campaign_stopped_preserve_failure'
    save(root/'terminal.json',dict(status=status,finished_at=stamp(),results=results,failure=failure,plan=record(root/'plan.json'),
        full_restart_qualified=failure is None,S3_prior_qualification_preserved=True,S4_complete=False,calibration_allowed=False))
    if failure is None:
        save(root/'independent_review.json',dict(status=status,recorded_at=stamp(),results=results,terminal=record(root/'terminal.json'),
            scope='New v2 rational-clock/33node material/reference/neighbor/cache/contact/duals/algorithm/compensation/149state/142physicalprefix;nonzero owned-noadvance16/40/64,split16+64,old40rollback40. Not S4 convergence,production or hardware-powerloss qualification.',
            full_restart_qualified=True,S3_prior_qualification_preserved=True,S4_complete=False,calibration_allowed=False))
    print(json.dumps(dict(status=status,failure=failure)),flush=True);return int(failure is not None)



def selftest_entry(root):
    from peri_elastic_normal import identity
    from peri_physical_checkpoint_reference import require_full_reference
    from peri_physical_checkpoint_negative import specification
    from peri_physical_checkpoint_support import derive_continuous
    assert len(specification()['controls'])==22
    frozen=Path(__file__).resolve().parent
    assert frozen.name=='frozen' and Path(sys.modules['peri_actual_checkpoint_campaign'].__file__).resolve().parent==frozen/'reference'
    assert Path(sys.modules['peri_physical_checkpoint_support'].__file__).resolve().parent==frozen/'reference'
    print(json.dumps(dict(status='physical_time_v2_frozen_entry_import_binding_passed_ZERO_solver',identity=identity(),
          native_solver_calls=0,material_calculations=0,frozen_script=str(Path(__file__).resolve()),reference_path=str(frozen/'reference'))))

def verify_entry_binding(root):
    from peri_elastic_normal import read,verify,record
    proof=read(root/'entry_binding.independent_review.json');a=read(root/'entry_binding.attempt.json')
    assert proof['status']=='physical_time_v2_frozen_entry_import_binding_passed_ZERO_solver'
    verify([proof['plan'],proof['attempt'],proof['stdout'],proof['stderr']])
    assert proof['plan']==record(root/'plan.json') and proof['native_solver_calls']==0
    assert a['return_code']==0 and a['stop_reason'] is None and a['stderr']['size']==0
    q=read(root/'entry_binding.stdout')
    assert q['status']==proof['status'] and q['native_solver_calls']==q['material_calculations']==0
    assert q['identity']['argv']==a['command'] and all(q['identity'][k]==v for k,v in a['birth_identity'].items())
    assert q['frozen_script']==str(root/'frozen/peri_physical_checkpoint_campaign_v2.py')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run','selftest-entry']);p.add_argument('--root',required=True,type=Path);a=p.parse_args()
    if a.action=='prepare':
        fresh=not a.root.exists()
        try:prepare(a.root)
        except BaseException as error:
            if fresh and a.root.exists():
                with (a.root/'preparation_error.json').open('x') as out:json.dump(dict(error=repr(error),traceback=traceback.format_exc(),native_dispatched=False),out,indent=2)
            raise
    elif a.action=='selftest-entry':selftest_entry(a.root)
    else:sys.exit(run(a.root))
