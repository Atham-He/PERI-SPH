"""New actual-PD complete-history campaign; prepare/run are each one-shot.

No old native or qualified raw reference is redispatched. A planned fresh8 /
resume32 / old20 rollback is a NEW preregistered physical restart experiment.
"""
import argparse,hashlib,json,math,os,resource,shutil,signal,subprocess,sys,time,traceback
from pathlib import Path


def require_tail_qualification(name,review,normal_anchor=None):
    """Validate the named physical scope, not only a success label.

    This imports a completed review without solving/reviewing it again. In
    particular, zero-mu needs the explicit same-discrete-model normal limit.
    The future1% trajectory-loss target is deliberately NOT a new S3 gate.
    """
    expected={'actual_PD_strict_stick_hold_release':(40,116),'same114sample_zero_mu':(16,228)}
    assert name in expected,'unknown prerequisite scope'
    assert review['status']=='never_attempted_original_tail_independently_passed'
    assert review['backend_selected_before_reference_imports'] is True and review['physical_solver_calls']==0
    assert review['original_failed_study_unchanged'] is True
    q=review['assessment'];steps,rows=expected[name]
    assert q['name']==name and q['accepted_steps']==steps and q['all_trial_point_tangents_checked'] is True
    assert 0<=q['maximum_normalized_error']<=1e-10
    assert q['independent_trial_material_points']>0
    geometry=q['geometry'];assert geometry['nodes']==33 and geometry['physical_DOF']==99 and geometry['normal_rows']==rows
    assert geometry['exact_rigidity_rank']==93 and q['independent_stability']['rigid_modes']==6
    if name=='actual_PD_strict_stick_hold_release':
        assert q['observations']['strict_stick_rows']>0 and q['observations']['previously_active_separation_rows']>0
    else:
        comparison=review['normal_limit_comparison']
        assert comparison['status']=='new_zero_mu_same114sample_normal_limit_passed' and comparison['steps']==16
        assert 0<=comparison['maximum_normalized_error']<=1e-10
        assert comparison['old_solver_dispatched'] is False and comparison['old_review_repeated'] is False
        assert normal_anchor is not None and comparison['old_qualified_review_imported']==normal_anchor['review']


def prepare(root):
    repo=Path(__file__).resolve().parents[1];sys.path.insert(0,str(repo/'tools'))
    from peri_elastic_normal import read,record,verify,save,reserves,stamp
    from peri_actual_checkpoint_support import initial_contract
    base=root.parent;tail=base/'study-v20-unattempted-friction-tails';producer=base/'engineering-actual-checkpoint-driver-preflight-v5';cross=base/'cross-runtime-v19-raw'
    assert base==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists(),'new unused campaign required'
    assert (repo/'tools/peri_actual_checkpoint_negative.py').is_file(),'semantic CP rejection implementation required before admission'
    from peri_actual_checkpoint_negative import specification
    negative_specification=specification();assert negative_specification['schema']==1 and negative_specification['positive_binding_control']
    execution_component=verified_execution_component(base)
    semantic=base/'engineering-CP-semantic-validator';semantic_receipt=read(semantic/'engineering_receipt.json');semantic_terminal=read(semantic/'terminal.json')
    assert semantic_terminal['return_code']==0 and semantic_receipt['semantic_rejections']==16 and semantic_receipt['positive_binding_controls']==1
    semantic_evidence=verified_semantic_component(semantic,semantic_receipt,semantic_terminal)
    verify(semantic_evidence)
    for name in ['peri_actual_checkpoint_negative.py','peri_actual_checkpoint_support.py','peri_canonical_archive_audit.py','peri_bounded_session_archive_audit.py','peri_bounded_mechanical_stream_audit.py']:
        assert record(repo/'tools'/name)['sha256']==record(semantic/'frozen'/name)['sha256'],'semantic component changed;new preflight required'
    terminal=read(tail/'terminal.json');assert terminal['status']=='two_original_never_attempted_friction_tails_independently_qualified' and terminal['failure'] is None
    tail_plan=read(tail/'plan.json');normal_anchor=tail_plan['normal_anchor']
    normal_records=[normal_anchor[k] for k in ['review','stdout','stderr']];verify(normal_records)
    prerequisites=[record(tail/'terminal.json'),record(tail/'plan.json')]+normal_records
    for name in ['actual_PD_strict_stick_hold_release','same114sample_zero_mu']:
        review=read(tail/(name+'.independent_review.json'));require_tail_qualification(name,review,normal_anchor)
        verify([review['attempt'],review['plan']])
        attempt=read(tail/(name+'.attempt.json'));verify([attempt['stdout'],attempt['stderr']])
        prerequisites += [record(tail/(name+'.independent_review.json')),review['attempt'],attempt['stdout'],attempt['stderr']]
    proof=read(producer/'independent_preflight.json');end=read(producer/'terminal.json');build=read(producer/'frozen_execution.json')['build']
    assert end['return_code']==0 and end['receipt']==record(producer/'independent_preflight.json')
    assert proof['owned_restore_noadvance_control']['original_fixture_and_both_new_owned_branches_full_bytes_equal']
    assert proof['owned_restore_noadvance_control']['actual_solver_calls']==0 and not proof['full_restart_qualified']
    native_inputs=[record(producer/name) for name in ['source_manifest.json','build_runtime_manifest.json','original.input','initial_candidate.bin','initial.phy','contract.bin']]
    verify(build['support']+build['source_snapshots']+native_inputs+[build['ELF']]);resources=reserves(base)
    assert int(dict(x.split(':',1) for x in Path('/proc/meminfo').read_text().splitlines())['MemAvailable'].split()[0])*1024>=8*2**30
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    root.mkdir(mode=0o700);frozen=root/'frozen';reference=frozen/'reference';reference.mkdir(parents=True)
    artifacts=prerequisites+build['support']+build['source_snapshots']+native_inputs+[build['ELF'],record(producer/'independent_preflight.json'),record(producer/'terminal.json')]
    artifacts += execution_component['immutable_evidence']+semantic_evidence
    for source in sorted((cross/'frozen/reference').glob('*.py')):
        target=reference/source.name;shutil.copy2(source,target);artifacts += [record(source),record(target)]
    for name in ['peri_actual_checkpoint_campaign.py','peri_actual_checkpoint_reference.py','peri_actual_checkpoint_support.py',
                 'peri_actual_checkpoint_negative.py','peri_canonical_archive_audit.py','peri_bounded_session_archive_audit.py','peri_bounded_mechanical_stream_audit.py']:
        source=repo/'tools'/name;target=(frozen if name in ['peri_actual_checkpoint_campaign.py','peri_actual_checkpoint_reference.py'] else reference)/name
        assert not target.exists();shutil.copy2(source,target);artifacts += [record(source),record(target)]
    for name in ['peri_fixed_python_raw_audit.py','_pydecimal.py']:
        source=cross/'frozen'/name;target=frozen/name;shutil.copy2(source,target);artifacts += [record(source),record(target)]
    for source_record in native_inputs:
        target=root/Path(source_record['path']).name;shutil.copy2(source_record['path'],target);artifacts.append(record(target))
    executable=frozen/'actual_PD_checkpoint';shutil.copy2(build['ELF']['path'],executable);artifacts.append(record(executable))
    benchmark=read(cross/'independent_review.json');assert benchmark['status']=='existing_v19_raw_independently_qualified_by_new_fixed_pure_Python_reference'
    artifacts += [record(cross/'independent_review.json'),record(cross/'terminal.json')]
    initial=initial_contract(root);assert initial['trusted'].hex()==proof['contract_sha256']
    case=dict(name='actual_PD_complete_history_continuous',mode=1,steps=40,protocol_kind=2,dt=.01,stroke=.01,rotation=0,azimuths=8,
        rho0=1.,parameters=[5.,1.,.5,.3],eta=.1,support=1.25,cap_width=1/64,tangent_stroke=[.0006,-.0004],mu=[.22,.22],expected_code=0,
        require_statistics=['nonzero_tangent0_rows','nonzero_tangent1_rows'])
    plan=dict(prepared_at=stamp(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
        native_producer_source=proof['source_commit'],native_producer_ELF=build['ELF'],executable=str(executable),case=case,
        trusted_contract=initial['trusted'].hex(),artifacts=artifacts,resources=resources,
        negative_specification=negative_specification,
        execution_component=execution_component,
        zero_step_semantic_preflight=record(semantic/'engineering_receipt.json'),
        native_CPU_seconds=900,native_wall_seconds=900,native_AS_GiB=1,threads=1,file_bound_GiB=4,total_output_GiB=4,
        host_reserve_GiB=120,output_reserve_GiB=150,reference_AS_GiB=8,reference_max_CPU_seconds=36000,
        reference_cost_formula='max(7200,ceil(2*actual_native_trial_material_points+2000))CPU;wall=CPU+600;freeze admission BEFORE each new assessment. No failed audit extension.',
        cost_basis=dict(measured_fixed_reference_CPU_seconds=benchmark['CPU_seconds'],points=benchmark['assessment']['independent_trial_material_points'],
            native40_previous_seconds=205.775167119,prospective_native_margin='900s includes durable writer/prefix and seal work;not an old failed run extension'),
        phases=['continuous40/fulloriginal95reference+8negatives/allphysicalbytes','owned_noadvance8','owned_noadvance20','owned_noadvance32',
                'planned_fresh_split8/wholepacket+bytes','explicit_resume8to40/wholepacket+bytes','explicit_rollback20to40/wholepacket+bytes',
                'correctSHA_semantic_controls_before_staging'],
        normalized_physical_gate=1e-10,original_precision_gate=1e-12,native_Newton_tolerance=1e-12,dense_penetration_stroke_gate=.01,
        full_candidate_and_physical_prefix_byte_equality=True,future_force_relative=.02,future_algorithmic_loss_over_positive_work=.01,
        failure_policy='Preserve exact exception/returncode/raw/resource/acceptedvsdebug/appendwitness andstop dependentphases. Neverredispatch state,reset/fallback,force/mass/speed/reference changes or relax gates.',
        new_input_scope='New CP-specific two-axis input,zero-advance preflight only previously;plannedsplit replicate is required restart experiment,notoldraw recovery.',
        S3_complete=False,full_restart_qualified=False,S4_complete=False,calibration_allowed=False)
    plan['critical_runtime']=[record(executable)]+[record(root/name) for name in ['source_manifest.json','build_runtime_manifest.json','original.input']]
    plan['critical_runtime'] += [x for x in build['support'] if '.so' in Path(x['path']).name]
    verify(artifacts);save(root/'plan.json',plan)
    print(json.dumps(dict(status='actual_complete_CP_campaign_prepared_NOT_passed',root=str(root))))


class Execution:
    """Bound a child to full OS identity and persist every distinct attempt."""
    def __init__(self,root,plan):
        self.root=root;self.plan=plan
    @staticmethod
    def process(pid):
        process=Path('/proc')/str(pid);fields=(process/'stat').read_text().rsplit(')',1)[1].split()
        return dict(pid=pid,ppid=int(fields[1]),start_ticks=fields[19],
            boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
            argv=(process/'cmdline').read_bytes().decode().strip(chr(0)).split(chr(0)),state=fields[0])
    def run(self,label,command,*,wall,cpu,gib,expected_code=0):
        from peri_elastic_normal import save,record,reserves,stamp
        root=self.root;assert not (root/(label+'.dispatch.json')).exists(),'never repeat a phase'
        save(root/(label+'.dispatch.json'),dict(recorded_at=stamp(),command=command,wall_seconds=wall,CPU_seconds=cpu,AS_GiB=gib,
            expected_code=expected_code,plan=record(root/'plan.json')))
        def limits():
            resource.setrlimit(resource.RLIMIT_CPU,(cpu,cpu));resource.setrlimit(resource.RLIMIT_AS,(gib*2**30,gib*2**30))
            resource.setrlimit(resource.RLIMIT_FSIZE,(4*2**30,4*2**30));resource.setrlimit(resource.RLIMIT_CORE,(0,0))
        start=time.monotonic();peak=threads=0;seen=None;stop=None
        with (root/(label+'.stdout')).open('xb') as out,(root/(label+'.stderr')).open('xb') as err:
            child=subprocess.Popen(command,cwd=root,stdout=out,stderr=err,start_new_session=True,preexec_fn=limits,
                env=dict(os.environ,PYTHONPATH=str(root/'frozen'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_MAX_THREADS='1'))
            # Capture birth before poll()/wait() can reap a fast noadvance child.
            # Its argv may still be pre-exec; birth fields remain stable at exec.
            born=self.process(child.pid);birth={k:born[k] for k in ['pid','ppid','start_ticks','boot_id']}
            parent=self.process(os.getpid());parent.pop('state')
            while child.poll() is None:
                try:
                    current=self.process(child.pid);state=current.pop('state')
                    assert {k:current[k] for k in birth}==birth,'child birth identity changed'
                    if state=='Z':continue
                    if seen:assert current==seen,'child identity changed'
                    elif current['argv']==command:
                        seen=current;save(root/(label+'.launch.json'),dict(recorded_at=stamp(),child=seen,parent=parent,command=command))
                    status=(Path('/proc')/str(child.pid)/'status').read_text().splitlines()
                    peak=max(peak,next(int(x.split()[1]) for x in status if x.startswith('VmRSS:')))
                    threads=max(threads,next(int(x.split()[1]) for x in status if x.startswith('Threads:')))
                    assert peak<=gib*1024**2 and threads<=1,'RSS/thread budget'
                    reserves(root);assert time.monotonic()-start<=wall,'wall budget'
                    assert sum(p.stat().st_size for p in root.rglob('*') if p.is_file())<=4*2**30,'total output budget'
                except (FileNotFoundError,ProcessLookupError):continue
                except Exception as error:stop=repr(error);os.killpg(child.pid,signal.SIGKILL);break
                time.sleep(.25)
            code=child.wait();out.flush();err.flush();os.fsync(out.fileno());os.fsync(err.fileno())
        result=dict(label=label,command=command,child_identity=seen,birth_identity=birth,child_pid=child.pid,parent_identity=parent,
            return_code=code,stop_reason=stop,wall_seconds=time.monotonic()-start,
            peak_RSS_KiB=peak,peak_threads=threads,stdout=record(root/(label+'.stdout')),stderr=record(root/(label+'.stderr')))
        # Keep exact first native error witness, including physical/debug cursor.
        errors=[]
        with (root/(label+'.stderr')).open('rb') as stream:error_prefix=stream.read(1024*1024)
        for line in error_prefix.decode(errors='replace').splitlines():
            try:
                value=json.loads(line)
                if value.get('type')=='checkpoint_driver_failure':errors.append(value)
            except (ValueError,AttributeError):pass
        result['native_failure_witnesses']=errors
        result['failure_witness_prefix_truncated']=result['stderr']['size']>len(error_prefix)
        save(root/(label+'.attempt.json'),result)
        assert code==expected_code and stop is None,(label,'execution',code,expected_code,stop)
        return result


def verified_execution_component(base):
    """Import the scoped guardian proof only if its class/imports are EXACT.

    Legitimate changes to prepare/assessment outside this class do not require
    rerunning an already passed guardian state. A changed guardian does.
    """
    import ast
    from peri_elastic_normal import read,record,verify
    old=Path(base)/'engineering-CP-execution-guardian';receipt=read(old/'engineering_receipt.json');plan=read(old/'plan.json');terminal=read(old/'terminal.json')
    assert terminal['return_code']==0 and terminal['status']==receipt['status'];verify([terminal['receipt'],receipt['plan']])
    assert receipt['control_count']==5 and receipt['native_solver_calls']==0
    frozen=old/'frozen/peri_actual_checkpoint_campaign.py'
    expected=next(x for x in plan['artifacts'] if x['path']==str(frozen));verify([expected])
    def pieces(path):
        text=Path(path).read_text();tree=ast.parse(text);lines=text.splitlines(keepends=True)
        cls=next(x for x in tree.body if isinstance(x,ast.ClassDef) and x.name=='Execution')
        imports=[ast.dump(x,include_attributes=False) for x in tree.body if isinstance(x,(ast.Import,ast.ImportFrom))]
        return ''.join(lines[cls.lineno-1:cls.end_lineno]),imports
    before=pieces(frozen);current=pieces(__file__);assert before==current,'guardian component changed;new preflight required'
    evidence=[record(old/'terminal.json'),terminal['receipt'],receipt['plan']]
    evidence += [x for x in plan['artifacts'] if not x['path'].startswith('/home/hezhy58/SPHinXsys/')]
    for key in ['fast','reference_binding','controlled_failure','controlled_wall']:
        evidence.append(receipt[key]);attempt=read(receipt[key]['path']);evidence += [attempt['stdout'],attempt['stderr']]
    verify(evidence)
    return dict(status='exact_guardian_class_and_top_imports_match_previously_qualified_component',
        class_SHA256=hashlib.sha256(current[0].encode()).hexdigest(),receipt=record(old/'engineering_receipt.json'),frozen_controller=record(frozen),
        control_count=5,old_preflight_not_repeated=True,immutable_evidence=evidence)


def verified_semantic_component(root,receipt,terminal):
    """Import immutable source/attempt/raw evidence without repeating controls."""
    from peri_elastic_normal import read,record,verify
    evidence=[record(root/'terminal.json'),terminal['receipt'],receipt['plan'],receipt['review']]
    verify(evidence);plan=read(receipt['plan']['path']);review=read(receipt['review']['path'])
    assert review['semantic_rejections']==16 and review['positive_binding_controls']==1
    assert review['accepted_step']==0 and not review['qualified_physical_nonzero_seed']
    evidence += [x for x in plan['artifacts'] if not x['path'].startswith('/home/hezhy58/SPHinXsys/')]
    for attempt_record in review['results']:
        verify([attempt_record]);attempt=read(attempt_record['path'])
        evidence += [attempt_record,attempt['stdout'],attempt['stderr']]
    verify(evidence);return evidence


def run(root):
    sys.path.insert(0,str(root/'frozen/reference'))
    from peri_elastic_normal import read,record,verify,save,identity,stamp
    from peri_actual_checkpoint_support import derive_continuous,selected,file_tree,assert_preserved,restored_noadvance_metadata,trace_packets
    plan=read(root/'plan.json');plan_record=record(root/'plan.json');assert not (root/'state.json').exists(),'never redispatch campaign'
    resource.setrlimit(resource.RLIMIT_AS,(8*2**30,8*2**30))
    save(root/'state.json',dict(status='actual_complete_CP_campaign_attempted',started_at=stamp(),identity=identity(),plan=record(root/'plan.json')))
    execution=Execution(root,plan);results=[];failure=None;trusted=bytes.fromhex(plan['trusted_contract']);exe=plan['executable']
    def verify_all():verify([plan_record]+plan['artifacts'])
    def native(label,args):
        verify_all();result=execution.run(label,[exe]+args,wall=900,cpu=900,gib=1)
        with Path(result['stdout']['path']).open() as stream:meta=json.loads(next(stream))
        assert meta['type']=='identity' and meta['argv']==result['command']
        assert all(meta[k]==v for k,v in result['birth_identity'].items()) and result['stderr']['size']==0
        if result['child_identity']:assert all(meta[k]==v for k,v in result['child_identity'].items())
        verify_all();return result
    def head(store):return hashlib.sha256((store/'HEAD').read_bytes()).hexdigest()
    def reference(label,first,last):
        attempt=read(root/(label+'.attempt.json'));points=0
        for row,_ in trace_packets(attempt['stdout']['path'],first,last,resume_generation=None if first==1 else first-1):
            if row['type']=='step':points+=33*len(row['solve']['trace'])
        cpu=max(7200,math.ceil(2*points+2000));assert cpu<=plan['reference_max_CPU_seconds'],'new reference cost exceeds prospective cap;do not attempt'
        save(root/(label+'.reference_admission.json'),dict(recorded_at=stamp(),attempt=record(root/(label+'.attempt.json')),
            actual_native_trial_material_points=points,CPU_seconds=cpu,wall_seconds=cpu+600,AS_GiB=8,
            prior_measured_cost=plan['cost_basis'],formula=plan['reference_cost_formula']))
        execution.run(label+'.reference',['/usr/bin/python3',str(root/'frozen/peri_actual_checkpoint_reference.py'),'--root',str(root),'--label',label],wall=cpu+600,cpu=cpu,gib=8)
        receipt=read(root/(label+'.independent_review.json'))
        assert receipt['status'] in ['new_actual_CP_continuous_full_reference_passed_NOT_full_restart','new_actual_CP_branch_full_trial_packets_identical_math_evidence_imported','new_actual_CP_branch_distinct_trial_context_independently_passed']
        verify_all();return record(root/(label+'.independent_review.json'))
    def bytes_assessment(label,store,steps,derived):
        checks=[selected(store,'selected_'+label+'_s'+str(step)+'.bin',step,derived['checkpoints'][step],trusted) for step in steps]
        result=dict(status='new_actual_CP_all_candidate_and_physical_output_bytes_matched',recorded_at=stamp(),native_attempt=record(root/(label+'.attempt.json')),
            checkpoints=checks,full_candidate_byte_equal=True,full_physical_prefix_byte_equal=True,full_restart_qualified=False)
        save(root/(label+'.bytes_review.json'),result);return record(root/(label+'.bytes_review.json'))
    try:
        verify_all();continuous=root/'continuous.store';split=root/'split.store';continuous.mkdir();split.mkdir()
        native('continuous',['--fresh',str(root),str(continuous),plan['trusted_contract'],'40','continuous'])
        derived=derive_continuous(root,root/'continuous.stdout');results.append(bytes_assessment('continuous',continuous,[8,20,32,40],derived))
        results.append(reference('continuous',1,40))
        for step in [8,20,32]:
            label='na'+str(step);before=file_tree(continuous)
            save(root/(label+'.inputs.json'),dict(recorded_at=stamp(),selected_input=record(continuous/('selected_continuous_s'+str(step)+'.bin')),expected_HEAD=record(continuous/'HEAD'),prior_files=before))
            attempt=native(label,['--restore-noadvance',str(root),str(continuous),plan['trusted_contract'],'selected_continuous_s'+str(step)+'.bin',str(step),head(continuous),label])
            expected=derived['checkpoints'][step];proof=selected(continuous,'selected_'+label+'_s'+str(step)+'.bin',step,expected,trusted)
            metadata=restored_noadvance_metadata(attempt['stdout']['path'],generation=step,step=step,time=derived['index']['steps'][step-1]['time'],archive_sha256=hashlib.sha256(expected['archive']).hexdigest())
            assert_preserved(before,file_tree(continuous),allowed_changes={'HEAD'})
            save(root/(label+'.independent_noadvance.json'),dict(recorded_at=stamp(),proof=proof,metadata=metadata,prior_files_preserved=True,new_physical_steps=0))
            results.append(record(root/(label+'.independent_noadvance.json')))
        native('split_first',['--fresh',str(root),str(split),plan['trusted_contract'],'8','split_first'])
        results.append(bytes_assessment('split_first',split,[8],derived));results.append(reference('split_first',1,8))
        for label,step,selection in [('split_resume',8,'selected_split_first_s8.bin'),('rollback20',20,'selected_split_resume_s20.bin')]:
            before=file_tree(split);save(root/(label+'.inputs.json'),dict(recorded_at=stamp(),selected_input=record(split/selection),expected_HEAD=record(split/'HEAD'),prior_files=before))
            native(label,['--resume',str(root),str(split),plan['trusted_contract'],selection,str(step),head(split),'40',label])
            assert_preserved(before,file_tree(split),allowed_changes={'HEAD'})
            results.append(bytes_assessment(label,split,[20,32,40] if step==8 else [32,40],derived))
            results.append(reference(label,step+1,40))
        from peri_actual_checkpoint_negative import run_controls
        negative_review=run_controls(root,plan,execution,derived);negative=read(negative_review['path'])
        assert negative['qualified_physical_nonzero_seed'] and negative['accepted_step']==8 and negative['semantic_rejections']==16
        results.append(negative_review)
        verify_all();status='actual_bounded_PD_complete_material_reference_contact_history_restart_independently_qualified'
    except Exception as error:
        failure=dict(error=repr(error),traceback=traceback.format_exc());status='actual_complete_CP_campaign_stopped_preserve_failure'
    save(root/'terminal.json',dict(status=status,finished_at=stamp(),results=results,failure=failure,plan=record(root/'plan.json'),
        full_restart_qualified=failure is None,S3_complete=False,S4_complete=False,calibration_allowed=False))
    if failure is None:
        save(root/'independent_review.json',dict(status=status,recorded_at=stamp(),results=results,terminal=record(root/'terminal.json'),
            scope='Actual bounded33node material/reference/neighbor/cache/signedMLS/contact/normal+two-tangentdual/algorithm/compensation/clock/committedphysicaloutput owned-noadvance/split/tailrollback;not S4 convergence/scalability/calibration or hardware powerloss.',
            full_restart_qualified=True,S3_complete=False,S4_complete=False,calibration_allowed=False))
    print(json.dumps(dict(status=status,failure=failure)),flush=True);return int(failure is not None)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run']);parser.add_argument('--root',required=True,type=Path);args=parser.parse_args()
    if args.action=='prepare':
        fresh=not args.root.exists()
        try:prepare(args.root)
        except Exception as error:
            if fresh and args.root.exists():
                with (args.root/'preparation_error.json').open('x') as stream:json.dump(dict(error=repr(error),traceback=traceback.format_exc(),native_dispatched=False),stream,indent=2)
            raise
    else:sys.exit(run(args.root))
