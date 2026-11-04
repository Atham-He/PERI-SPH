"""New bounded disk-kernel controls; never redispatch an existing state."""
import argparse,copy,json,math,os,resource,shutil,subprocess,sys,time,traceback
from pathlib import Path
from peri_elastic_normal import identity,read,record,reserves,save,stamp,verify
from peri_coulomb_disk_audit import inspect

def controls():
    cases=[]
    def case(name,tokens,code=0,root=False,fd=False,message=None):
        c=dict(name=name,tokens=tokens.split(),expected_code=code,disk_root=root,finite_difference=fd)
        if message:c['expected_message']=message
        cases.append(c)
    case('stick','1 .1 .05 0 0 .3',root=True,fd=True)
    case('slip_two_components','1 -.24 -.18 .4 .3 .3',root=True,fd=True)
    case('reverse_slip','1 .24 .18 -.4 -.3 .3',root=True,fd=True)
    case('mixed_nonroot_operator','.7 .2 -.35 .5 .1 .4',fd=True)
    case('zero_mu','1 .2 .1 .4 .3 0',fd=True)
    case('zero_normal_separation','0 0 0 .4 .3 .3',root=True)
    case('origin','0 0 0 0 0 .3',root=True)
    case('exact_disk_boundary','1 .5 0 0 0 .5',root=True)
    case('negative_normal_trial','-.1 .1 .2 .3 .4 .3',fd=True)
    case('slip_basis_rotated','1 .18 -.24 -.3 .4 .3',root=True,fd=True)
    case('large_impulse_small_velocity','1e8 1e8 0 1e-8 0 1',fd=True)
    case('tiny_slip','1 -.5 0 1e-20 0 .5',root=True)
    case('subnormal_relative_velocity','1 -1 0 5e-324 0 1',root=True)
    case('square_box_not_disk','1 -.3 -.3 1 1 .3',fd=True)
    case('precision_range_rejection','1e150 -.5e150 0 1e-150 0 .5',code=2)
    case('unrepresentable_work_rejection','1.7976931348623157e308 1.7976931348623157e308 0 -1.7976931348623157e308 0 1',code=2)
    message='Coulomb disk: nonfinite input or negative coefficient'
    case('invalid_mu','1 0 0 0 0 -.1',code=1,message=message)
    case('nonfinite_normal','inf 0 0 0 0 .3',code=1,message=message)
    case('nonfinite_tangent','1 nan 0 0 0 .3',code=1,message=message)
    case('zero_normal_nonzero_tangent_trial','0 .1 -.2 1 2 .5')
    case('subnormal_normal_impulse','5e-324 0 0 .4 .3 .5')
    return cases

def prepare(root,executable,prospective):
    assert not root.exists(),'new unused directory required'
    repo=Path(__file__).resolve().parents[1]
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=repo).strip()
    plan=read(prospective);assert plan['roundoff_gate']==1e-12 and plan['normalized_accuracy_gate']==1e-10
    budget=reserves(root.parent)
    available=int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:')))*1024
    assert available>=2**30;budget['available_memory_bytes']=available
    root.mkdir();(root/'frozen').mkdir();artifacts=[]
    files=['src/shared/peridynamics/pd_coulomb_disk.h','tests/3d_examples/extra_3d_PD_framework_tests/coulomb_disk_probe.cpp',
        'tests/3d_examples/extra_3d_PD_framework_tests/coulomb_disk_tests.cpp','tests/3d_examples/extra_3d_PD_framework_tests/CMakeLists.txt',
        'tools/peri_coulomb_disk.py','tools/peri_coulomb_disk_audit.py','tools/test_peri_coulomb_disk.py','tools/peri_elastic_normal.py',
        'docs/peri/COULOMB_DISK_KERNEL_20261105.md']
    for file in files:
        source=repo/file;target=root/'frozen'/source.name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
    binary=root/'frozen'/executable.name;shutil.copy2(executable,binary)
    artifacts.extend([record(executable),record(binary),record(prospective)])
    build=repo/'build/peri-latest-20261103';obj='extra_3d_PD_framework_tests/CMakeFiles/extra_3d_PD_coulomb_disk_probe.dir/coulomb_disk_probe.cpp.o'
    command=subprocess.check_output(['ninja','-C',str(build),'-t','commands',obj],text=True)
    assert not any(flag in command for flag in ['-ffast-math','-Ofast','-ffinite-math-only'])
    (root/'frozen/compiler_command.txt').write_text(command);artifacts.append(record(root/'frozen/compiler_command.txt'))
    deps=subprocess.check_output(['ninja','-C',str(build),'-t','deps',obj],text=True);paths=[]
    for line in deps.splitlines():
        if line.startswith('    '):
            source=Path(line.strip());source=source if source.is_absolute() else build/source
            if source.is_file():paths.append(source.resolve())
    paths=sorted(set(paths));assert paths and sum(p.stat().st_size for p in paths)<256*2**20
    manifest=[]
    for index,source in enumerate(paths):
        target=root/'frozen/compiler_headers'/str(index)/source.name;target.parent.mkdir(parents=True);shutil.copy2(source,target)
        entry=dict(original=record(source),frozen=record(target));manifest.append(entry);artifacts.extend(entry.values())
    save(root/'frozen/compiler_dependencies.json',manifest);artifacts.append(record(root/'frozen/compiler_dependencies.json'))
    artifacts.append(record(Path(sys.executable).resolve()))
    import _decimal
    modules={}
    for module in [math,_decimal]:
        modules[module.__name__]=module.__spec__.origin
        if getattr(module,'__file__',None):artifacts.append(record(module.__file__))
        else:assert module.__spec__.origin in ['built-in','frozen']
    save(root/'frozen/interpreter_modules.json',modules);artifacts.append(record(root/'frozen/interpreter_modules.json'))
    for line in subprocess.check_output(['ldd',str(executable)],text=True).splitlines():
        for token in line.split():
            if token.startswith('/') and Path(token).is_file():artifacts.append(record(token))
    plan.update(controls=controls(),prepared_at=stamp(),identity=identity(),resources=budget,python_version=sys.version,
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip())
    save(root/'plan.json',plan);artifacts.append(record(root/'plan.json'));inputs=[]
    for case in plan['controls']:
        path=root/(case['name']+'.input');path.write_text(' '.join(case['tokens'])+'\n');inputs.append(record(path))
    save(root/'preparation.json',dict(artifacts=artifacts,inputs=inputs,executable=str(binary),compiler_dependencies=len(paths)))
    print(json.dumps(dict(status='prepared_kernel_only_not_passed',cases=len(plan['controls']),compiler_dependencies=len(paths),source=plan['source_commit'])))

def limits():
    resource.setrlimit(resource.RLIMIT_CPU,(30,30));resource.setrlimit(resource.RLIMIT_AS,(512*2**20,512*2**20))
    resource.setrlimit(resource.RLIMIT_FSIZE,(4*2**20,4*2**20));resource.setrlimit(resource.RLIMIT_CORE,(0,0))

def records(path):return [json.loads(line) for line in Path(path).open() if line.strip()]

def run(root):
    assert not (root/'state.json').exists(),'old state must not dispatch'
    prep=read(root/'preparation.json');plan=read(root/'plan.json');attempts=[];failure=None;failure_traceback=None
    save(root/'state.json',dict(status='attempted',started_at=stamp(),identity=identity(),preparation=record(root/'preparation.json')))
    try:
        verify(prep['artifacts']+prep['inputs'])
        for case in plan['controls']:
            reserves(root);out=root/(case['name']+'.jsonl');err=root/(case['name']+'.stderr')
            command=[prep['executable'],str(root/(case['name']+'.input'))];start=time.monotonic();timeout=False
            with out.open('xb') as stdout,err.open('xb') as stderr:
                child=subprocess.Popen(command,stdout=stdout,stderr=stderr,preexec_fn=limits,env=dict(os.environ,OMP_NUM_THREADS='1'))
                try:code=child.wait(timeout=30)
                except subprocess.TimeoutExpired:child.kill();code=child.wait();timeout=True
            attempt=dict(name=case['name'],command=command,return_code=code,timed_out=timeout,wall_seconds=time.monotonic()-start,
                stdout=record(out),stderr=record(err),max_child_RSS_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,last_accepted_time=0,accepted_steps=0)
            attempts.append(attempt);assert not timeout and attempt['max_child_RSS_KiB']<=512*1024
            raw=records(out);meta=raw[0]
            assert meta['argv']==command and meta['pid']==child.pid and meta['ppid']==identity()['pid']
            assert meta['boot_id']==identity()['boot_id'] and meta['start_ticks'].isdigit()
            attempt['runtime_gate']=inspect(case,raw,code)
        verify(prep['artifacts']+prep['inputs'])
    except Exception as error:failure=repr(error);failure_traceback=traceback.format_exc()
    save(root/'terminal.json',dict(status='runtime_completed_pending_independent_review' if failure is None else 'stopped_for_review',
        finished_at=stamp(),attempts=attempts,failure=failure,failure_traceback=failure_traceback))
    print(json.dumps(dict(status=read(root/'terminal.json')['status'],attempts=len(attempts),failure=failure)))
    return int(failure is not None)

def review(root):
    assert not (root/'independent_review.json').exists(),'immutable review already exists'
    resource.setrlimit(resource.RLIMIT_AS,(2*2**30,2*2**30));resource.setrlimit(resource.RLIMIT_CPU,(300,300))
    prep=read(root/'preparation.json');plan=read(root/'plan.json');terminal=read(root/'terminal.json')
    inputs=[record(root/name) for name in ['state.json','terminal.json','preparation.json']]
    result=dict(status='coulomb_disk_kernel_review_failed',recorded_at=stamp(),identity=identity());assessed=[]
    try:
        verify(inputs);verify(prep['artifacts']+prep['inputs'])
        assert terminal['status']=='runtime_completed_pending_independent_review' and len(terminal['attempts'])==len(plan['controls'])
        raw={}
        for case,attempt in zip(plan['controls'],terminal['attempts']):
            verify([attempt['stdout'],attempt['stderr']]);data=records(attempt['stdout']['path']);raw[case['name']]=data
            assessed.append(inspect(case,data,attempt['return_code']))
        a=raw['slip_two_components'][2];b=raw['slip_basis_rotated'][2]
        for field in ['residual','normal_derivative']:
            assert max(abs(p-q) for p,q in zip(b[field],[-a[field][1],a[field][0]]))<=1e-12
        assert abs(a['friction_loss']-b['friction_loss'])<=1e-12
        negative=[]
        for label,name in [('binary64_cancellation','large_impulse_small_velocity'),('box_instead_of_disk','square_box_not_disk'),
            ('dropped_normal_coupling','mixed_nonroot_operator'),('friction_loss_wrong_sign','slip_two_components'),
            ('forged_precision_acceptance','precision_range_rejection')]:
            data=copy.deepcopy(raw[name]);r=data[2];case=next(c for c in plan['controls'] if c['name']==name)
            if label=='binary64_cancellation':r['residual'][0]=1e8-min(1e8,1e8-1e-8)
            elif label=='box_instead_of_disk':r['residual']=[0,0]
            elif label=='dropped_normal_coupling':r['normal_derivative']=[0,0]
            elif label=='friction_loss_wrong_sign':r['friction_loss']=-r['friction_loss']
            else:r['within_budget']=True;data[-1]['arithmetic_budget_satisfied']=True
            path=root/(label+'.negative.jsonl')
            with path.open('x') as stream:
                for row in data:stream.write(json.dumps(row,allow_nan=False)+'\n')
            try:inspect(case,data,case['expected_code'])
            except AssertionError:negative.append(dict(name=label,rejected=True,artifact=record(path)))
            else:raise AssertionError(('negative accepted',label))
        verify(inputs);verify(prep['artifacts']+prep['inputs'])
        result.update(status='coulomb_disk_kernel_values_derivatives_enclosures_independently_passed',cases=assessed,
            maximum_normalized_error=max(c.get('maximum_normalized_error',0) for c in assessed),negative_controls=negative,basis_rotation_passed=True,
            scope='dimensionless two-tangent disk natural map only; no global contact solve, actual PD friction trajectory or restart proof',
            accepted_physical_steps=0,full_global_friction_qualified=False,full_restart_qualified=False,S3_complete=False,calibration_allowed=False)
    except Exception as error:result.update(failure=repr(error),failure_traceback=traceback.format_exc(),completed_assessments=assessed)
    result.update(state=record(root/'state.json'),terminal=record(root/'terminal.json'),preparation=record(root/'preparation.json'))
    save(root/'independent_review.json',result);print(json.dumps(dict(status=result['status'],cases=len(assessed),failure=result.get('failure'))))
    return int(not result['status'].endswith('_passed'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run','review']);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--executable',type=Path);p.add_argument('--prospective',type=Path);args=p.parse_args()
    if args.action=='prepare':
        fresh=not args.root.exists()
        try:prepare(args.root,args.executable,args.prospective)
        except Exception as error:
            if fresh and args.root.exists():save(args.root/'preparation_error.json',dict(recorded_at=stamp(),failure=repr(error),failure_traceback=traceback.format_exc(),native_dispatched=False))
            raise
    else:sys.exit(globals()[args.action](args.root))
