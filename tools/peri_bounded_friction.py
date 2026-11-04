"""New bounded actualPD friction campaign; immutable per-case independent review.

Run verifies each physical case with a distinct full reference process BEFORE
dispatching the next native case. An attempted state is never dispatched twice.
"""
import argparse,json,math,os,resource,shutil,subprocess,sys,time,traceback
from pathlib import Path
from peri_elastic_normal import identity,read,record,reserves,save,stamp,verify
from peri_bounded_friction_audit import inspect_stream
from peri_bounded_friction_support import encode,complete_cursor,anchored_normal,compare_rotated,compare_normal_anchor
from peri_bounded_friction_negative import OneStepCapture,negatives
from peri_coulomb_global import live_identity
from peri_early_friction_support import qualification_prerequisite


def prepare(root,executable,prospective,normal_root):
    assert not root.exists(),'new unused study required'
    repo=Path(__file__).resolve().parents[1]
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=repo).strip(),'committed source required'
    plan=read(prospective)
    assert plan['normalized_physical_gate']==1e-10 and plan['roundoff_gate']==1e-12
    assert plan['dense_penetration_stroke_gate']==.01 and not plan['S3_complete'] and not plan['calibration_allowed']
    assert plan['campaign']=='early_cycle_unattempted_tail','a new explicit continuation plan is required'
    previous_plan=read(Path(plan['continuation']['prior_root'])/'plan.json')
    trusted_normal=Path(previous_plan['normal_anchor']['review']['path']).parent
    if normal_root is not None:
        assert normal_root.resolve()==trusted_normal.resolve(),'normal anchor must equal the sealed previous plan'
    normal_root=trusted_normal
    continuation=qualification_prerequisite(**plan['continuation'])
    budget=reserves(root.parent)
    mem={line.split(':')[0]:int(line.split()[1])*1024 for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith(('MemTotal:','MemAvailable:'))}
    assert mem['MemAvailable']>=8*2**30,('available_memory',mem)
    budget.update(mem,native_threads=1,reader_AS_GiB=8,affinity_cpus=sorted(os.sched_getaffinity(0)))
    anchor=anchored_normal(normal_root)
    global_review=read(plan['original_global_review'])
    assert global_review['status']=='bounded_original_global_normal_two_tangent_coupling_independently_qualified_with_explicit_precision_exclusions'
    assert global_review['supported_original_control_count']==18 and not global_review['S3_complete']
    disk_review=read(plan['disk_review'])
    assert disk_review['status']=='coulomb_disk_kernel_values_derivatives_enclosures_independently_passed'
    root.mkdir();(root/'frozen').mkdir();artifacts=anchor['artifacts']+continuation['artifacts']+[record(prospective)]
    artifacts.extend([record(plan['original_global_review']),record(plan['disk_review'])])
    verify([global_review['prior_independent_review']]);artifacts.append(global_review['prior_independent_review'])
    files=[
        'src/shared/peridynamics/pd_bounded_gel_friction.h','src/shared/peridynamics/pd_bounded_friction_protocol.h',
        'src/shared/peridynamics/pd_early_tangent_cycle.h',
        'src/shared/peridynamics/pd_bounded_normal_scene.h','src/shared/peridynamics/pd_bounded_gel.h',
        'src/shared/peridynamics/pd_fixed_reference.h','src/shared/peridynamics/pd_surface_quadrature.h',
        'src/shared/peridynamics/pd_gel_affine_path.h','src/shared/peridynamics/pd_gel_material.h','src/shared/peridynamics/pd_gel_material.cpp',
        'src/shared/peridynamics/pd_coulomb_contact.h','src/shared/peridynamics/pd_coulomb_disk.h',
        'src/shared/peridynamics/pd_many_contact.h','src/shared/peridynamics/pd_exact_contact_balance.h',
        'tests/3d_examples/extra_3d_PD_framework_tests/bounded_gel_friction_probe.cpp',
        'tests/3d_examples/extra_3d_PD_framework_tests/coulomb_contact_probe_io.h',
        'tests/3d_examples/extra_3d_PD_framework_tests/bounded_gel_probe_io.h',
        'tests/3d_examples/extra_3d_PD_framework_tests/CMakeLists.txt',
        'tools/peri_bounded_friction.py','tools/peri_bounded_friction_audit.py',
        'tools/peri_bounded_friction_support.py','tools/peri_bounded_friction_negative.py',
        'tools/peri_early_friction_support.py','tools/test_peri_early_friction_support.py',
        'tools/test_peri_bounded_friction.py','tools/test_peri_bounded_friction_support.py','tools/peri_coulomb_contact_audit.py',
        'tools/peri_bounded_gel_audit.py','tools/peri_pd_reference_audit.py','tools/peri_pd_material_reference.py',
        'tools/peri_gel_affine_audit.py','tools/peri_normal_projection_audit.py','tools/peri_dense_contact_audit.py',
        'tools/peri_many_contact_audit.py','tools/peri_elastic_normal.py','tools/peri_coulomb_global.py',
        'tools/peri_coulomb_global_controls.py','tools/peri_coulomb_global_audit.py',
        'docs/peri/BOUNDED_FRICTION_QUALIFICATION_20261105.md',
        'docs/peri/EARLY_FRICTION_CYCLE_20261105.md']
    for file in files:
        source=repo/file;target=root/'frozen'/source.name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
    binary=root/'frozen'/executable.name;shutil.copy2(executable,binary)
    artifacts.extend([record(executable),record(binary)])
    build=repo/'build/peri-latest-20261103'
    obj='extra_3d_PD_framework_tests/CMakeFiles/extra_3d_PD_bounded_gel_friction_probe.dir/bounded_gel_friction_probe.cpp.o'
    command=subprocess.check_output(['ninja','-C',str(build),'-t','commands',obj],text=True)
    assert not any(flag in command for flag in ['-ffast-math','-Ofast','-ffinite-math-only'])
    (root/'frozen/compiler_command.txt').write_text(command);artifacts.append(record(root/'frozen/compiler_command.txt'))
    deps=subprocess.check_output(['ninja','-C',str(build),'-t','deps',obj],text=True);paths=[]
    for line in deps.splitlines():
        if line.startswith('    '):
            p=Path(line.strip());p=p if p.is_absolute() else build/p
            if p.is_file():paths.append(p.resolve())
    paths=sorted(set(paths));assert paths and sum(p.stat().st_size for p in paths)<512*2**20
    manifest=[]
    for index,source in enumerate(paths):
        target=root/'frozen/compiler_headers'/str(index)/source.name;target.parent.mkdir(parents=True);shutil.copy2(source,target)
        entry=dict(original=record(source),frozen=record(target));manifest.append(entry);artifacts.extend(entry.values())
    save(root/'frozen/compiler_dependencies.json',manifest);artifacts.append(record(root/'frozen/compiler_dependencies.json'))
    import _decimal
    artifacts.append(record(Path(sys.executable).resolve()))
    for module in [math,_decimal]:
        if getattr(module,'__file__',None):artifacts.append(record(module.__file__))
        else:assert module.__spec__.origin in ['built-in','frozen']
    for line in subprocess.check_output(['ldd',str(executable)],text=True).splitlines():
        for token in line.split():
            if token.startswith('/') and Path(token).is_file():artifacts.append(record(token))
    cases=continuation['continuation_controls']
    cases[0]['new_negative_controls']=True
    assert all(c['native_wall_seconds']<=3600 and c['output_GiB']<=4 for c in cases)
    plan.update(controls=cases,normal_anchor=anchor,continuation_qualification=continuation,
        prepared_at=stamp(),identity=identity(),resources=budget,
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),python_version=sys.version)
    save(root/'plan.json',plan);artifacts.append(record(root/'plan.json'));inputs=[]
    for case in cases:
        path=root/(case['name']+'.input');path.write_text(encode(case));inputs.append(record(path))
    save(root/'preparation.json',dict(artifacts=artifacts,inputs=inputs,executable=str(binary),compiler_dependencies=len(paths)))
    print(json.dumps(dict(status='prepared_actual_PD_friction_not_passed',cases=len(cases),compiler_dependencies=len(paths),source=plan['source_commit'])))


def event(root,value):
    with (root/'events.jsonl').open('a') as stream:
        stream.write(json.dumps(dict(recorded_at=stamp(),identity=identity(),**value),allow_nan=False)+'\n')
        stream.flush();os.fsync(stream.fileno())


def native_limits(case):
    resource.setrlimit(resource.RLIMIT_CPU,(case['native_wall_seconds'],case['native_wall_seconds']))
    resource.setrlimit(resource.RLIMIT_AS,(2**30,2**30))
    resource.setrlimit(resource.RLIMIT_FSIZE,(4*2**30,4*2**30));resource.setrlimit(resource.RLIMIT_CORE,(0,0))


def review_case(root,name):
    prep=read(root/'preparation.json');plan=read(root/'plan.json')
    case=next(c for c in plan['controls'] if c['name']==name)
    resource.setrlimit(resource.RLIMIT_AS,(8*2**30,8*2**30))
    resource.setrlimit(resource.RLIMIT_CPU,(case['reader_CPU_seconds'],case['reader_CPU_seconds']))
    path=root/(name+'.independent_review.json');assert not path.exists(),'case review already exists'
    save(root/(name+'.review_state.json'),dict(status='started',identity=identity(),started_at=stamp()))
    attempt=read(root/(name+'.attempt.json'));initial=record(root/(name+'.attempt.json'))
    result=dict(status='actual_PD_friction_case_review_failed',name=name,identity=identity(),recorded_at=stamp())
    start=time.monotonic()
    try:
        verify([initial]);verify(prep['artifacts']+prep['inputs']);verify([attempt['stdout'],attempt['stderr']])
        assert attempt['return_code']==0 and not attempt['stop_reason']
        with Path(attempt['stdout']['path']).open() as stream:meta=json.loads(next(stream))
        assert meta==attempt['native_identity'] and meta['argv']==attempt['command']
        assert meta['pid']==attempt['child_pid'] and meta['ppid']==attempt['identity']['pid']
        assert meta['boot_id']==attempt['identity']['boot_id'] and meta['start_ticks'].isdigit()
        capture=OneStepCapture() if name=='actual_PD_two_tangent_press' or case.get('new_negative_controls') else None
        result['assessment']=inspect_stream(case,attempt['stdout']['path'],0,True,capture)
        if capture:
            result['negative_controls']=negatives(capture,root/'single_step_negatives',attempt['stdout'])
        if name=='actual_PD_two_tangent_rotated':
            previous=read(root/'actual_PD_two_tangent_press.independent_review.json')
            assert previous['status']=='actual_PD_friction_case_independently_passed'
            result['frame_comparison']=compare_rotated(root/'actual_PD_two_tangent_press.jsonl',attempt['stdout']['path'])
        if name=='same114sample_zero_mu':
            result['normal_limit_comparison']=compare_normal_anchor(attempt['stdout']['path'],plan['normal_anchor'])
        verify([initial]);verify(prep['artifacts']+prep['inputs']);verify([attempt['stdout'],attempt['stderr']])
        result['status']='actual_PD_friction_case_independently_passed'
    except Exception as error:result.update(failure=repr(error),failure_traceback=traceback.format_exc())
    result.update(wall_seconds=time.monotonic()-start,attempt=initial,preparation=record(root/'preparation.json'),
        S3_complete=False,full_restart_qualified=False,calibration_allowed=False)
    save(path,result);print(json.dumps(dict(status=result['status'],name=name,failure=result.get('failure'))))
    return int(result['status']!='actual_PD_friction_case_independently_passed')


def run(root):
    assert not (root/'state.json').exists(),'never redispatch an attempted state'
    prep=read(root/'preparation.json');plan=read(root/'plan.json');attempts=[];failure=None;failure_traceback=None
    save(root/'state.json',dict(status='attempted',started_at=stamp(),identity=identity(),preparation=record(root/'preparation.json')))
    try:
        verify(prep['artifacts']+prep['inputs'])
        for case in plan['controls']:
            reserves(root);name=case['name'];out=root/(name+'.jsonl');err=root/(name+'.stderr')
            command=[prep['executable'],str(root/(name+'.input'))];start=time.monotonic();stop_reason=None;peak=0;threads=0
            with out.open('xb') as stdout,err.open('xb') as stderr:
                child=subprocess.Popen(command,stdout=stdout,stderr=stderr,preexec_fn=lambda:native_limits(case),
                    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
                observed=live_identity(child.pid);event(root,dict(phase='native_started',case=name,command=command,live=observed,pid=child.pid))
                while child.poll() is None:
                    try:
                        status=Path('/proc')/str(child.pid)/'status'
                        for line in status.read_text().splitlines():
                            if line.startswith('VmRSS:'):peak=max(peak,int(line.split()[1]))
                            elif line.startswith('Threads:'):threads=max(threads,int(line.split()[1]))
                    except FileNotFoundError:pass
                    if time.monotonic()-start>case['native_wall_seconds']:stop_reason='native_wall_budget'
                    elif peak>1024*1024:stop_reason='native_RSS_budget'
                    elif threads>1:stop_reason='native_thread_budget'
                    else:
                        try:reserves(root)
                        except AssertionError:stop_reason='disk_reserve'
                    if stop_reason:child.kill();break
                    time.sleep(.25)
                code=child.wait()
            cursor,meta=complete_cursor(out,err,code)
            attempt=dict(name=name,identity=identity(),command=command,return_code=code,stop_reason=stop_reason,
                child_pid=child.pid,native_identity=meta,
                native_wall_seconds=time.monotonic()-start,peak_RSS_KiB=peak,peak_threads=threads,live_observation=observed,
                stdout=record(out),stderr=record(err),cursor=cursor)
            attempts.append(attempt);save(root/(name+'.attempt.json'),attempt)
            event(root,dict(phase='native_terminal',case=name,return_code=code,stop_reason=stop_reason,cursor=cursor))
            assert not stop_reason and code==0,(name,'execution_failed',code,stop_reason)
            assert meta['argv']==command and meta['pid']==child.pid and meta['ppid']==identity()['pid']
            assert meta['boot_id']==identity()['boot_id'] and meta['start_ticks'].isdigit()
            if observed and observed['argv']==command:assert {k:meta[k] for k in observed}==observed
            assert err.stat().st_size==0 and not cursor['partial_record_present'] and not cursor['native_error_cursor_used']
            assert cursor['last_durable_accepted_steps']==case['steps']
            attempt['runtime_gate']=inspect_stream(case,out,code,False)
            # Distinct reference process, completed and immutable before the
            # next independent initial physical case may start.
            review_command=[sys.executable,str(Path(__file__).resolve()),'review-case','--root',str(root),'--case',name]
            log=root/(name+'.review.stdout');error=root/(name+'.review.stderr')
            with log.open('xb') as stdout,error.open('xb') as stderr:
                reviewer=subprocess.Popen(review_command,stdout=stdout,stderr=stderr,
                    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
                event(root,dict(phase='independent_review_started',case=name,command=review_command,pid=reviewer.pid,live=live_identity(reviewer.pid)))
                try:review_code=reviewer.wait(timeout=case['reader_CPU_seconds']+600)
                except subprocess.TimeoutExpired:reviewer.kill();review_code=reviewer.wait()
            event(root,dict(phase='independent_review_terminal',case=name,return_code=review_code,stdout=record(log),stderr=record(error)))
            assert review_code==0,(name,'independent_review_failed',review_code)
            qualified=read(root/(name+'.independent_review.json'))
            assert qualified['status']=='actual_PD_friction_case_independently_passed'
            attempt['independent_review']=record(root/(name+'.independent_review.json'))
        verify(prep['artifacts']+prep['inputs'])
    except Exception as error:failure=repr(error);failure_traceback=traceback.format_exc()
    terminal=dict(status='all_actual_PD_cases_independently_qualified' if failure is None else 'stopped_for_review',
        finished_at=stamp(),attempts=attempts,failure=failure,failure_traceback=failure_traceback,
        S3_complete=False,full_restart_qualified=False,calibration_allowed=False)
    save(root/'terminal.json',terminal)
    if failure is None:
        save(root/'independent_review.json',dict(status='actual_bounded_PD_friction_trajectories_independently_passed',
            recorded_at=stamp(),identity=identity(),case_reviews=[a['independent_review'] for a in attempts],
            terminal=record(root/'terminal.json'),preparation=record(root/'preparation.json'),
            scope='new early-cycle actual PD trajectory plus previously-unattempted stick/zero-mu; prior press/rotation and failed old protocol retained separately; not convergence/scalability or restart',
            S3_complete=False,S4_complete=False,full_restart_qualified=False,calibration_allowed=False))
    print(json.dumps(dict(status=terminal['status'],attempts=len(attempts),failure=failure)))
    return int(failure is not None)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run','review-case']);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--executable',type=Path);p.add_argument('--prospective',type=Path);p.add_argument('--normal-root',type=Path);p.add_argument('--case')
    args=p.parse_args()
    if args.action=='prepare':
        fresh=not args.root.exists()
        try:prepare(args.root,args.executable,args.prospective,args.normal_root)
        except Exception as error:
            if fresh and args.root.exists():save(args.root/'preparation_error.json',dict(recorded_at=stamp(),failure=repr(error),native_dispatched=False,failure_traceback=traceback.format_exc()))
            raise
    elif args.action=='run':sys.exit(run(args.root))
    else:sys.exit(review_case(args.root,args.case))
