"""Only two never-attempted original native inputs, each gated by fixed review.

Preparation requires a successful complete new fixed-reference v19 audit.
Original attempted solver case/controller and all prior reviews stay immutable.
"""
import argparse, hashlib, json, os, resource, shutil, signal, subprocess, sys, time, traceback
from pathlib import Path


def review(root,name):
    # Select ONE prospectively frozen backend before reference/Fraction imports.
    reference=root/'frozen/reference';sys.path.insert(0,str(reference))
    from peri_fixed_python_raw_audit import select_backend
    backend=select_backend(reference)
    from peri_elastic_normal import read,record,verify,save,identity,stamp
    from peri_bounded_friction_audit import FrictionInspector,inspect_stream
    from peri_bounded_friction_support import compare_normal_anchor
    plan=read(root/'plan.json');case=next(c for c in plan['controls'] if c['name']==name)
    assert not (root/(name+'.review_state.json')).exists(),'no repeated review'
    verify(plan['artifacts']);attempt=read(root/(name+'.attempt.json'));verify([attempt['stdout'],attempt['stderr']])
    assert attempt['return_code']==0 and attempt['stop_reason'] is None
    save(root/(name+'.review_state.json'),dict(status='new_tail_review_attempted',identity=identity(),started_at=stamp()))
    resource.setrlimit(resource.RLIMIT_AS,(8*2**30,8*2**30));resource.setrlimit(resource.RLIMIT_CPU,(21600,21600))
    resource.setrlimit(resource.RLIMIT_FSIZE,(256*2**20,256*2**20));resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    started=time.monotonic();original_step=FrictionInspector.step;cursor={}
    def observed(self,row,old,z,loss,history):
        cursor.update(old_step=old['step'],recorded_step=row['state']['step'],time=row['state']['time'])
        result=original_step(self,row,old,z,loss,history)
        with (root/(name+'.review_progress.jsonl')).open('a') as stream:
            stream.write(json.dumps(dict(recorded_at=stamp(),cursor=cursor,trial_points=self.trial_points,
                maximum_normalized_error=self.maximum,CPU_seconds=time.process_time()),allow_nan=False)+'\n');stream.flush();os.fsync(stream.fileno())
        return result
    FrictionInspector.step=observed
    receipt=dict(status='new_tail_independent_review_failed',identity=identity(),backend=backend.__file__,
        backend_selected_before_reference_imports=True,physical_solver_calls=0,original_failed_study_unchanged=True,
        S3_complete=False,full_restart_qualified=False,S4_complete=False,calibration_allowed=False)
    try:
        receipt['assessment']=inspect_stream(case,attempt['stdout']['path'],0,True)
        assert receipt['assessment']['all_trial_point_tangents_checked']
        if name=='same114sample_zero_mu':receipt['normal_limit_comparison']=compare_normal_anchor(attempt['stdout']['path'],plan['normal_anchor'])
        verify(plan['artifacts']);verify([attempt['stdout'],attempt['stderr']])
        receipt['status']='never_attempted_original_tail_independently_passed'
    except Exception as error:receipt.update(failure=repr(error),failure_traceback=traceback.format_exc(),last_cursor=cursor)
    finally:FrictionInspector.step=original_step
    receipt.update(recorded_at=stamp(),wall_seconds=time.monotonic()-started,CPU_seconds=time.process_time(),
        peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,attempt=record(root/(name+'.attempt.json')),plan=record(root/'plan.json'))
    save(root/(name+'.independent_review.json'),receipt);return int('failure' in receipt)


def prepare(root):
    repo=Path(__file__).resolve().parents[1];base=root.parent;cross=base/'cross-runtime-v19-raw';old=base/'study-v19-early-friction-cycle'
    sys.path.insert(0,str(repo/'tools'))
    from peri_elastic_normal import read,record,verify,save,reserves,stamp
    from peri_bounded_friction_support import controls,encode,anchored_normal
    assert root.parent==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists(),'no redispatch'
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    terminal=read(cross/'terminal.json');qualified=read(cross/'independent_review.json')
    assert terminal['status']=='new_fixed_reference_independently_qualified' and terminal['return_code']==0 and terminal['stop_reason'] is None
    assert qualified['status']=='existing_v19_raw_independently_qualified_by_new_fixed_pure_Python_reference'
    assert qualified['assessment']['accepted_steps']==40 and qualified['assessment']['all_trial_point_tangents_checked']
    assert len(qualified['negative_controls'])==8 and all(x['rejected'] for x in qualified['negative_controls'])
    verify([qualified['plan'],qualified['audit_state'],qualified['raw'],terminal['independent_review']])
    oldterminal=read(old/'terminal.json');assert oldterminal['status']=='stopped_for_review' and len(oldterminal['attempts'])==1
    cases=controls()[3:];assert len(cases)==2
    oldprep=read(old/'preparation.json');native=Path(oldprep['executable'])
    snapshots=[x for x in oldprep['artifacts'] if x['path'].startswith(str(old/'frozen')+'/')]
    assert snapshots and any(x['path']==str(native) for x in snapshots);verify(snapshots)
    # Reuse actual ORIGINAL frozen executable, not a rebuilt current observer.
    anchor=anchored_normal(Path(read(old/'plan.json')['normal_anchor']['review']['path']).parent)
    artifacts=snapshots+anchor['artifacts']+[record(cross/x) for x in ['terminal.json','independent_review.json','plan.json','audit_state.json']]
    artifacts += [record(old/x) for x in ['terminal.json','state.json','plan.json','preparation.json']]+[qualified['raw'],qualified['original_failed_review']]
    for case in cases:
        name=case['name']
        assert not any((old/(name+suffix)).exists() for suffix in ['.attempt.json','.jsonl','.review_state.json','.independent_review.json'])
        assert (old/(name+'.input')).read_text()==encode(case);artifacts.append(record(old/(name+'.input')))
    budget=reserves(root.parent);root.mkdir(mode=0o700);(root/'frozen').mkdir();ref=root/'frozen/reference';ref.mkdir()
    for path in sorted((cross/'frozen/reference').glob('*.py')):
        target=ref/path.name;shutil.copy2(path,target);artifacts += [record(path),record(target)]
    fixed=cross/'frozen/peri_fixed_python_raw_audit.py';shutil.copy2(fixed,ref/fixed.name);artifacts += [record(fixed),record(ref/fixed.name)]
    script=root/'frozen/peri_friction_unattempted_tail.py';shutil.copy2(__file__,script);artifacts += [record(__file__),record(script)]
    original=root/'frozen/original_b492_friction_probe';shutil.copy2(native,original);artifacts += [record(original)]
    for line in subprocess.check_output(['ldd',str(original)],text=True).splitlines():
        assert 'not found' not in line
        for token in line.split():
            if token.startswith('/'):
                library=Path(token);expected=[x for x in oldprep['artifacts'] if Path(x['path']).resolve()==library.resolve()]
                assert expected and all(record(library)['sha256']==x['sha256'] for x in expected),'original native linked library changed'
                artifacts.append(record(library))
    document=repo/'docs/peri/UNATTEMPTED_FRICTION_TAILS_20261105.md'
    shutil.copy2(document,root/'frozen'/document.name);artifacts += [record(document),record(root/'frozen'/document.name)]
    for case in cases:
        destination=root/(case['name']+'.input');shutil.copy2(old/destination.name,destination);artifacts.append(record(destination))
        assert record(destination)['sha256']==record(old/destination.name)['sha256']
    # Freeze actual reference runtime dependencies, including selected Python
    # interpreter/_pydecimal and original dynamically linked native libraries.
    crossplan=read(cross/'plan.json')
    runtime=[x for x in crossplan['artifacts'] if not x['path'].startswith(str(old)+'/') and not x['path'].startswith('/home/')]
    verify(runtime);artifacts += runtime
    assert qualified['CPU_seconds']<21600
    plan=dict(prepared_at=stamp(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
        controls=cases,executable=str(original),script=str(script),artifacts=artifacts,normal_anchor=anchor,
        imported_complete_fixed_review=record(cross/'independent_review.json'),resources=budget,
        native_AS_GiB=1,native_threads=1,native_output_GiB_per_case=4,review_AS_GiB=8,review_threads=1,
        review_CPU_seconds=21600,review_wall_seconds=22200,total_output_GiB=8,host_reserve_GiB=120,output_reserve_GiB=150,
        cost_basis=dict(complete_fixed_reference_CPU_seconds=qualified['CPU_seconds'],
            complete_fixed_reference_trial_points=qualified['assessment']['independent_trial_material_points'],
            rationale='New two unattempted inputs;21600CPU/22200wall same prospectively bounded fixed reference as complete40step audit,not extension of any failed review. Native original900/1500s budgets unchanged. Freeze before first tail dispatch.'),
        reference_backend='fixed_pure_Python_decimal_before_import',new_negative_controls_required=False,
        inherited_eight_negative_controls=qualified['negative_controls'],normalized_physical_gate=1e-10,roundoff_gate=1e-12,dense_penetration_stroke_gate=.01,
        scope='Only original never-attempted strictstick40 then same114zeroMu16; each new native and fixed independent review before next. Old first native/controller/reviews never repeated.',
        failure_policy='Any applicable unexpected native/review/integrity/resource/accuracy failure freezes exact raw/returncode/cursors and stops later tail. No retry/reset/fallback/gate relaxation.',
        old_v19_terminal_changed=False,old_C_runtime_cause_resolved=False,S3_complete=False,full_restart_qualified=False,S4_complete=False,calibration_allowed=False)
    verify(artifacts);save(root/'plan.json',plan);print(json.dumps(dict(status='two_never_attempted_tails_prepared_NOT_passed',root=str(root))))


def run(root):
    sys.path.insert(0,str(root/'frozen/reference'))
    from peri_elastic_normal import read,record,verify,save,reserves,identity,stamp
    from peri_bounded_friction_support import complete_cursor
    plan=read(root/'plan.json');assert not (root/'state.json').exists(),'no repeated controller'
    save(root/'state.json',dict(status='new_two_tail_state_attempted',identity=identity(),started_at=stamp()))
    attempts=[];failure=None
    def bounded(label,command,wall,cpu,gib,file_gib,output_gib):
        def limits():
            resource.setrlimit(resource.RLIMIT_AS,(gib*2**30,gib*2**30));resource.setrlimit(resource.RLIMIT_CPU,(cpu,cpu))
            resource.setrlimit(resource.RLIMIT_FSIZE,(file_gib*2**30,file_gib*2**30));resource.setrlimit(resource.RLIMIT_CORE,(0,0))
        def live(pid):
            p=Path('/proc')/str(pid);fields=(p/'stat').read_text().rsplit(')',1)[1].split()
            return dict(pid=pid,ppid=int(fields[1]),start_ticks=fields[19],
                boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
                argv=(p/'cmdline').read_bytes().decode().strip('\0').split('\0'))
        started=time.monotonic();stop=None;peak=threads=0;seen=None
        with (root/(label+'.stdout')).open('xb') as out,(root/(label+'.stderr')).open('xb') as err:
            child=subprocess.Popen(command,stdout=out,stderr=err,start_new_session=True,preexec_fn=limits,
                env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_MAX_THREADS='1'))
            while child.poll() is None:
                try:
                    current=live(child.pid)
                    if seen:assert current==seen,'live process identity changed'
                    elif current['argv']==command:seen=current;save(root/(label+'.launch.json'),dict(recorded_at=stamp(),command=command,child=seen,parent=identity()))
                    status=(Path('/proc')/str(child.pid)/'status').read_text().splitlines()
                    peak=max(peak,next(int(x.split()[1]) for x in status if x.startswith('VmRSS:')))
                    threads=max(threads,next(int(x.split()[1]) for x in status if x.startswith('Threads:')))
                    assert peak<=gib*1024**2 and threads<=1,'RSS or threads'
                    reserves(root);assert time.monotonic()-started<=wall,'wall budget'
                    assert sum(p.stat().st_size for p in root.glob(label+'*') if p.is_file())<=output_gib*2**30,'case output budget'
                    assert sum(p.stat().st_size for p in root.rglob('*') if p.is_file())<=8*2**30,'total study output budget'
                except FileNotFoundError:continue
                except Exception as error:stop=repr(error);os.killpg(child.pid,signal.SIGKILL);break
                time.sleep(.25)
            code=child.wait();out.flush();err.flush();os.fsync(out.fileno());os.fsync(err.fileno())
        result=dict(return_code=code,stop_reason=stop,wall_seconds=time.monotonic()-started,peak_RSS_KiB=peak,peak_threads=threads,
            child_identity=seen,command=command,stdout=record(root/(label+'.stdout')),stderr=record(root/(label+'.stderr')))
        save(root/(label+'.execution.json'),result);return result
    try:
        for case in plan['controls']:
            verify(plan['artifacts']);reserves(root);name=case['name']
            native=bounded(name,[plan['executable'],str(root/(name+'.input'))],case['native_wall_seconds'],case['native_wall_seconds'],1,4,4)
            cursor,meta=complete_cursor(native['stdout']['path'],native['stderr']['path'],native['return_code'])
            attempt=dict(native,name=name,native_identity=meta,cursor=cursor);save(root/(name+'.attempt.json'),attempt);attempts.append(record(root/(name+'.attempt.json')))
            assert native['return_code']==0 and native['stop_reason'] is None,(name,'native',native['return_code'],native['stop_reason'])
            assert native['child_identity'] and all(meta[k]==v for k,v in native['child_identity'].items())
            assert cursor['durability_confirmed'] and cursor['last_durable_accepted_steps']==case['steps'] and not cursor['partial_record_present']
            assert native['stderr']['size']==0
            reviewed=bounded(name+'.review',[sys.executable,plan['script'],'review','--root',str(root),'--case',name],22200,21600,8,1,4)
            assert reviewed['return_code']==0 and reviewed['stop_reason'] is None,(name,'review',reviewed['return_code'],reviewed['stop_reason'])
            receipt=read(root/(name+'.independent_review.json'));assert receipt['status']=='never_attempted_original_tail_independently_passed'
            verify(plan['artifacts'])
        status='two_original_never_attempted_friction_tails_independently_qualified'
    except Exception as error:
        failure=dict(error=repr(error),traceback=traceback.format_exc());status='new_friction_tail_stopped_preserve_failure'
    save(root/'terminal.json',dict(status=status,finished_at=stamp(),attempts=attempts,failure=failure,plan=record(root/'plan.json'),
        S3_complete=False,full_restart_qualified=False,S4_complete=False,calibration_allowed=False))
    print(json.dumps(dict(status=status,failure=failure)),flush=True);return int(failure is not None)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run','review']);parser.add_argument('--root',required=True,type=Path);parser.add_argument('--case')
    args=parser.parse_args()
    if args.action=='prepare':
        fresh=not args.root.exists()
        try:prepare(args.root)
        except Exception as error:
            if fresh and args.root.exists():
                with (args.root/'preparation_error.json').open('x') as stream:
                    json.dump(dict(error=repr(error),traceback=traceback.format_exc(),native_dispatched=False),stream,indent=2)
                    stream.flush();os.fsync(stream.fileno())
            raise
    elif args.action=='run':sys.exit(run(args.root))
    else:sys.exit(review(args.root,args.case))
