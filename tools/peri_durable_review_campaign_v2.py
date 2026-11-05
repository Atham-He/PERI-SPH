"""NEW persistent full95 readonly audit experiment of unqualified native80.

Never launches native physics,old index builders,old controllers or old reviews.
Accepted equation steps are assessed once in four declared segments. Explicit
owned noadvance checkpoints at1/8/40 perform no new accepted equation steps.
"""
from pathlib import Path
import sys
_ENTRY_DIR=Path(__file__).resolve().parent
if (_ENTRY_DIR/"reference").is_dir():sys.path.insert(0,str(_ENTRY_DIR/"reference"))
import argparse,copy,hashlib,json,math,os,resource,shutil,subprocess,sys,traceback
from peri_elastic_normal import read,record,save,reserves,stamp,identity
from peri_actual_checkpoint_campaign import Execution,verified_execution_component
from peri_buffered_metadata import verify_unique


INDEX_STATUS='new_existing_actual80_index_artifact_independently_qualified_producer_failure_preserved_NOT_equation_review'


def require_index(q):
    assert q['status']==INDEX_STATUS and q['producer_return_code_preserved']==-9 and q['producer_status_changed'] is False
    assert q['index_builder_called'] is False and q['physical_solver_calls']==q['material_calculations']==0
    assert q['precision']==95 and q['C_decimal_loaded'] is False and q['actual_full_material_audit_admitted'] is False
    a=q['assessment'];assert a['accepted_steps']==80 and a['frontiers']==81 and a['point_metadata']==11847 and a['final_time']==.4
    assert a['all_complete_row_payloads_checked'] is True and a['all81prefix_line_offset_state_friction_row_typedSHA_checked'] is True
    assert a['material_equation_qualified'] is False


def budgets(index,limit=36000):
    endpoints=[1,8,40,80];previous=0;out=[]
    for end in endpoints:
        delta=index['index'][end]['trial_points']-index['index'][previous]['trial_points']
        assert type(delta) is int and delta>0
        cpu=2*delta+500
        out.append(dict(label='equations_'+str(end),first=previous+1,last=end,point_metadata=delta,CPU_seconds=cpu,wall_seconds=cpu+300))
        if end<80:out.append(dict(label='owned_noadvance_'+str(end),first=end+1,last=end,point_metadata=0,CPU_seconds=300,wall_seconds=600))
        previous=end
    assert sum(x['point_metadata'] for x in out)==11847
    assert sum(x['CPU_seconds'] for x in out)==26594 and sum(x['CPU_seconds'] for x in out)<=limit
    assert all(x['CPU_seconds']<=limit for x in out)
    return out


def guard_controls():
    good=dict(status=INDEX_STATUS,producer_return_code_preserved=-9,producer_status_changed=False,index_builder_called=False,
        physical_solver_calls=0,material_calculations=0,precision=95,C_decimal_loaded=False,actual_full_material_audit_admitted=False,
        assessment=dict(accepted_steps=80,frontiers=81,point_metadata=11847,final_time=.4,all_complete_row_payloads_checked=True,
            all81prefix_line_offset_state_friction_row_typedSHA_checked=True,material_equation_qualified=False))
    require_index(good);controls=[]
    for label,fn in [('partial_steps',lambda q:q['assessment'].update(accepted_steps=79)),
                     ('missing_fullrow',lambda q:q['assessment'].update(all_complete_row_payloads_checked=False)),
                     ('changed_failed_producer',lambda q:q.update(producer_status_changed=True)),
                     ('wrong_precision',lambda q:q.update(precision=94)),
                     ('equations_claimed_by_index',lambda q:q['assessment'].update(material_equation_qualified=True))]:
        bad=copy.deepcopy(good);fn(bad)
        try:require_index(bad)
        except AssertionError:controls.append(dict(name=label,rejected=True))
        else:raise AssertionError('admission negative accepted:'+label)
    fake=dict(index=[dict(trial_points=round(11847*k/80)) for k in range(81)])
    phases=budgets(fake)
    try:budgets(fake,26593)
    except AssertionError:controls.append(dict(name='cumulative_new_reference_cost_limit',rejected=True))
    else:raise AssertionError('cost limit negative accepted')
    assert len(controls)==6
    return dict(valid_binding=True,negatives=controls,all_equation_points_assessed_once=11847,prospective_NEW_CPU_seconds=26594)


def prepare(root):
    repo=Path(__file__).resolve().parents[1];base=root.parent;qualified=base/'analysis-s4-existing-index-independent-v3'
    assert base==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists()
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    q=read(qualified/'independent_review.json');require_index(q);end=read(qualified/'terminal.json')
    assert end['return_code']==0 and end['receipt']==record(qualified/'independent_review.json')
    reserves(base);execution=verified_execution_component(base);root.mkdir(mode=0o700)
    try:
        frozen=root/'frozen';shutil.copytree(qualified/'frozen',frozen)
        artifacts=[record(p) for p in sorted(frozen.rglob('*')) if p.is_file()]
        for name in ['peri_durable_audit_entry_v2.py','peri_durable_review_campaign_v2.py']:
            source=repo/'tools'/name;target=frozen/name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
        prior=read(qualified/'plan.json');artifacts+=prior['artifacts']+prior['input_records']
        artifacts += [record(qualified/n) for n in ['plan.json','independent_review.json','terminal.json','review.attempt.json']]
        index=read(q['index_artifact']['path']);oldplan=read(prior['original_plan']['path'])
        phases=budgets(index,oldplan['reference_max_CPU_seconds'])
        lineage=[]
        for name,record_name in [('study-v23-physical-time-complete-history','continuous.reference_admission.json'),
                                  ('analysis-s4-actual-durable-index-v1','plan.json'),('analysis-s4-actual-streamed-index-v2','plan.json'),
                                  ('analysis-s4-existing-index-independent-v3','plan.json')]:
            p=base/name/record_name;v=read(p)
            lineage.append(dict(id=name,CPU_ns=int(v['CPU_seconds'])*10**9,wall_ns=int(v['wall_seconds'])*10**9,
                evidence_sha256=record(p)['sha256'],kind='unknown_full_reservation'))
            artifacts.append(record(p))
        # Prior self-reported CPU diagnostics cannot settle unknown final usage.
        # A fresh cold-entry bound is also conservatively retained in the ledger.
        lineage.append(dict(id='new_cold_entry_preflight',CPU_ns=30*10**9,wall_ns=60*10**9,evidence_sha256=record(frozen/'peri_durable_audit_entry_v2.py')['sha256'],kind='unknown_full_reservation'))
        previous=base/'study-v24-existing-native80-durable-full95'
        old_guard=record(previous/'entry_preflight.stdout');old_source=record(previous/'frozen/peri_durable_review_campaign.py')
        import ast
        def fn_ast(path,name):
            tree=ast.parse(Path(path).read_text())
            return ast.dump(next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name))
        for name in ['require_index','budgets','guard_controls']:
            assert fn_ast(old_source['path'],name)==fn_ast(__file__,name),'unchanged guard transfer AST'
        artifacts += [old_guard,old_source,record(previous/'preparation.json')]
        # The failed earlier parent never dispatched a material child.
        assert not (previous/'state.json').exists() and not list((previous/'store').iterdir())
        lineage[-1]['id']='new_controller_and_child_preflights'
        lineage[-1]['CPU_ns']=60*10**9;lineage[-1]['wall_ns']=120*10**9
        lineage.append(dict(id='previous_v24_child_preflight',CPU_ns=30*10**9,wall_ns=60*10**9,evidence_sha256=old_guard['sha256'],kind='unknown_full_reservation'))
        plan=dict(prepared_at=stamp(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
            artifacts=artifacts,reused_guard_qualification=old_guard,guard_source_AST_preserved=True,execution_component=execution,case=index['case'],index_artifact=q['index_artifact'],qualified_index_review=record(qualified/'independent_review.json'),
            raw=index['raw'],original_native_attempt=prior['native_attempt'],phases=phases,prior_resource_reservations=lineage,
            NEW_reference_CPU_seconds=sum(x['CPU_seconds'] for x in phases),NEW_reference_wall_seconds=sum(x['wall_seconds'] for x in phases),
            inherited_max_reference_CPU_seconds=oldplan['reference_max_CPU_seconds'],reference_AS_GiB=8,RSS_GiB=8,threads=1,total_output_GiB=4,
            scope='NEW persistent fixedpure95 complete equation assessment of existing UNQUALIFIED native80. Four disjoint accepted-step segments1/2..8/9..40/41..80,owned noadvance1/8/40,zero native solves.No old controller/state/review restarted.',
            budget_basis='2*actual_native_point_metadata+2000 across four segments(500each),plus3*300CPU owned-noadvance.Next costs frozen before any phase,OLD admissions unchanged;new batch26594<inherited36000. Prior unmeasured final costs retained at full reserved caps.',
            all_original_Inspector_and_finalization_methods_unchanged=True,all10_original_negative_controls_required=True,
            independent_z_loss_context_capture_not_native_reset=True,actual_material_audit_admitted=True,
            S3_prior_qualification_preserved=True,new_physical_CP_restart_qualified=False,full_restart_qualified=False,S4_complete=False,calibration_allowed=False,
            failure_policy='Stop all dependent phases on any applicable failure,preserve exact code/stage/ledger/CP/debugcursor/raw.SizeSHA.No same-root retry,old state resume,fallback/reset,gates/source or resource extension.')
        verify_unique(artifacts);save(root/'plan.json',plan)
        parent_entry=Execution(root,plan).run('controller_preflight',[sys.executable,'-B',str(frozen/'peri_durable_review_campaign_v2.py'),'--controller-preflight'],wall=60,cpu=30,gib=1)
        parent_binding=read(root/'controller_preflight.stdout')
        assert parent_binding['status']=='actual_frozen_controller_reference_bootstrap_qualified' and parent_binding['material_steps']==0
        entry=Execution(root,plan).run('entry_preflight',[sys.executable,'-B',str(frozen/'peri_durable_audit_entry_v2.py'),'--preflight'],wall=60,cpu=30,gib=1)
        binding=read(root/'entry_preflight.stdout');assert binding['status']=='new_frozen_durable_audit_entry_and_admission_guards_qualified'
        assert binding['guard_controls']['valid_binding'] and len(binding['guard_controls']['negatives'])==6 and binding['precision']==95 and not binding['C_decimal_loaded']
        assert binding['physical_solver_calls']==binding['material_calculations']==0
        sha=lambda b:hashlib.sha256(b).hexdigest()
        canonical=lambda v:json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
        software=sha(canonical(artifacts));contract=dict(schema=1,nodes=33,contacts=116,total_steps=80,index=index['index'],decimal_profile=binding['decimal_profile'],
            CPU_cap=float(plan['NEW_reference_CPU_seconds']),case_sha=index['source_case_sha256'],header_sha=index['header_sha256'],raw_sha256=index['raw']['sha256'],software_sha256=software)
        save(root/'state_contract.json',contract)
        resource_contract=dict(schema=1,plan_sha256=record(root/'plan.json')['sha256'],source_sha256=software,
            CPU_cap_ns=sum(e['CPU_ns'] for e in lineage)+plan['NEW_reference_CPU_seconds']*10**9,
            wall_cap_ns=sum(e['wall_ns'] for e in lineage)+plan['NEW_reference_wall_seconds']*10**9,
            prior_charges=lineage,measurement='wait4_rusage_ceil_nanoseconds_v1')
        save(root/'resource_contract.json',resource_contract);(root/'resource_ledger').mkdir();(root/'store').mkdir()
        verify_unique(artifacts)
        save(root/'preparation.json',dict(status='new_durable_full95_readonly_campaign_prepared_NOT_mathematics_passed',prepared_at=stamp(),plan=record(root/'plan.json'),
            state_contract=record(root/'state_contract.json'),resource_contract=record(root/'resource_contract.json'),entry_preflight=record(root/'entry_preflight.attempt.json'),
            cold_entry_qualified=True,actual_controller_entry_qualified=True,controller_preflight=record(root/'controller_preflight.attempt.json'),actual_solver_calls=0,actual_material_steps=0,S4_complete=False))
        print(json.dumps(dict(status='new_durable_full95_readonly_campaign_prepared',root=str(root),NEW_CPU_seconds=plan['NEW_reference_CPU_seconds'],prior_full_CPU_reservations_seconds=sum(e['CPU_ns'] for e in lineage)/10**9)))
    except BaseException as error:
        save(root/'preparation_failure.json',dict(recorded_at=stamp(),failure=repr(error),traceback=traceback.format_exc()));raise


def run(root):
    sys.path.insert(0,str(root/'frozen/reference'))
    from peri_durable_reference_launch import run as launch
    from peri_reference_budget import Ledger
    prep=read(root/'preparation.json');plan=read(root/'plan.json')
    assert prep['cold_entry_qualified'] and prep['actual_controller_entry_qualified'] and not (root/'state.json').exists(),'one fresh campaign only,never resume old controller'
    pins=[prep[k] for k in ['plan','state_contract','resource_contract','entry_preflight','controller_preflight']]
    verify_unique(plan['artifacts']+pins);require_index(read(plan['qualified_index_review']['path']))
    save(root/'state.json',dict(status='new_persistent_full95_readonly_campaign_attempted',started_at=stamp(),identity=identity(),plan=prep['plan']))
    resource.setrlimit(resource.RLIMIT_AS,(8*2**30,8*2**30));contract=read(root/'resource_contract.json')
    ledger=Ledger(root/'resource_ledger',contract);resource_head=None;state_head=None;selected=None;results=[];last_step=0
    try:
        for phase in plan['phases']:
            reserves(root);verify_unique(plan['artifacts']+pins)
            label=phase['label'];directory=root/('phase_'+label)
            request=dict(plan=prep['plan'],index=plan['index_artifact'],state_contract=prep['state_contract'],store=str(root/'store'),
                selected_leaf=selected,expected_checkpoint_HEAD_sha256=state_head,end_step=phase['last'],
                resource_ledger=str(root/'resource_ledger'),resource_contract=prep['resource_contract'],attempt_id=label,
                launch_directory=str(directory),material_audit_admitted=True)
            request_path=root/(label+'.request.json');save(request_path,request);request_record=record(request_path)
            command=[sys.executable,'-B',str(root/'frozen/peri_durable_audit_entry_v2.py'),'--request',str(request_path),'--request-sha',request_record['sha256']]
            limits=dict(RSS_bytes=8*2**30,AS_bytes=8*2**30,threads=1,output_bytes=4*2**30,output_roots=[str(root)],
                reserves=[dict(path='/mnt/f',minimum_bytes=150*2**30)])
            result=launch(directory,root/'resource_ledger',contract,label,command,resource_head,phase['CPU_seconds']*10**9,phase['wall_seconds']*10**9,limits)
            resource_head=result['resource_HEAD_sha256'];assert result['return_code']==0 and result['status']=='complete_one_shot_resource_witness',result
            review=read(directory/'reference_review.json');assert 'failure' not in review and review['accepted_step']==phase['last']
            commits=review['commits'];assert commits and commits[-1]['accepted_step']==phase['last']
            if phase['point_metadata']==0:
                assert len(commits)==1 and last_step==phase['last']
                oldstate=read(root/'store'/selected/'commit.json');newstate=read(root/'store'/commits[-1]['leaf']/'commit.json')
                assert oldstate['snapshot_sha256']==newstate['snapshot_sha256'],'actual owned noadvance entire independent state bytes'
            else:
                assert phase['first']==last_step+1 and len(commits)==(phase['last']-last_step)+(1 if last_step==0 else 0)
                last_step=phase['last']
            state_head=commits[-1]['HEAD_sha256'];selected=commits[-1]['leaf']
            results.append(dict(phase=phase,request=request_record,attempt=record(directory/'terminal.json'),review=record(directory/'reference_review.json'),
                checkpoint_HEAD_sha256=state_head,resource_HEAD_sha256=resource_head))
            save(root/(label+'.qualified_segment.json'),results[-1])
        final=read(root/'phase_equations_80/reference_review.json')
        from peri_physical_checkpoint_reference import require_full_reference
        require_full_reference(final,plan['case']);assert final['assessment']['independent_trial_material_points']==11847
        verify_unique(plan['artifacts']+pins)
        outcome=dict(status='new_existing_native80_complete_original95_equations_and_owned_auditor_noadvance_qualified_NOT_physical_restart',finished_at=stamp(),
            plan=prep['plan'],qualified_index_review=plan['qualified_index_review'],raw=plan['raw'],original_native_attempt=plan['original_native_attempt'],
            segments=results,final_reference=record(root/'phase_equations_80/reference_review.json'),resource_state=ledger.read(resource_head),
            accepted_equation_steps_assessed_once=80,owned_noadvance_steps=[1,8,40],original_negative_controls=10,
            physical_solver_calls=0,old_native_runs_repeated=False,old_fully_qualified_math_repeated=False,
            actual_physical_CP_restart_qualified=False,trajectory_convergence_qualified=False,S4_complete=False,calibration_allowed=False)
        save(root/'independent_review.json',outcome);save(root/'terminal.json',dict(status=outcome['status'],return_code=0,finished_at=stamp(),receipt=record(root/'independent_review.json')))
        print(json.dumps(dict(status=outcome['status'],accepted_steps=80)))
    except BaseException as error:
        save(root/'failure.json',dict(recorded_at=stamp(),failure=repr(error),traceback=traceback.format_exc(),last_complete_segment=results[-1] if results else None,
            resource_HEAD_sha256=ledger.head(),checkpoint_HEAD_sha256=state_head,selected_leaf=selected))
        save(root/'terminal.json',dict(status='new_persistent_full95_readonly_campaign_failed_dependents_stopped',return_code=1,finished_at=stamp(),failure=record(root/'failure.json')));raise


if __name__=='__main__':
    if '--controller-preflight' in sys.argv:
        assert (_ENTRY_DIR/'reference').is_dir()
        assert Path(sys.modules['peri_elastic_normal'].__file__).resolve()==(_ENTRY_DIR/'reference/peri_elastic_normal.py').resolve()
        print(json.dumps(dict(status='actual_frozen_controller_reference_bootstrap_qualified',material_steps=0,physical_solver_calls=0,
            helper_module=str(Path(sys.modules['peri_elastic_normal'].__file__).resolve()))))
        sys.exit(0)
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);p.add_argument('--root',required=True,type=Path);a=p.parse_args()
    (prepare if a.action=='prepare' else run)(a.root)
