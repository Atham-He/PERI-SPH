"""Frozen new global algebra campaign; never redispatch an attempted state."""
import argparse,copy,json,math,os,resource,shutil,subprocess,sys,time,traceback
from pathlib import Path
from peri_elastic_normal import identity,read,record,reserves,save,stamp,verify
from peri_coulomb_global_controls import controls
from peri_coulomb_global_audit import inspect
from peri_coulomb_contact_audit import close_tree,decimal,fdot


def prepare(root,executable,prospective):
    assert not root.exists(),'new unused directory required'
    repo=Path(__file__).resolve().parents[1]
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=repo).strip(),'committed source required'
    plan=read(prospective)
    assert plan['roundoff_gate']==1e-12 and plan['normalized_accuracy_gate']==1e-10
    assert not any(plan[k] for k in ['actual_PD_qualified','full_restart_qualified','S3_complete','calibration_allowed'])
    budget=reserves(root.parent)
    available=int(next(line.split()[1] for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith('MemAvailable:')))*1024
    assert available>=2**30;budget['available_memory_bytes']=available
    root.mkdir();(root/'frozen').mkdir();artifacts=[]
    files=[
        'src/shared/peridynamics/pd_coulomb_contact.h','src/shared/peridynamics/pd_coulomb_disk.h',
        'src/shared/peridynamics/pd_many_contact.h','src/shared/peridynamics/pd_exact_contact_balance.h',
        'tests/3d_examples/extra_3d_PD_framework_tests/coulomb_contact_probe.cpp',
        'tests/3d_examples/extra_3d_PD_framework_tests/coulomb_contact_probe_io.h',
        'tests/3d_examples/extra_3d_PD_framework_tests/bounded_gel_probe_io.h',
        'tests/3d_examples/extra_3d_PD_framework_tests/coulomb_contact_tests.cpp',
        'tests/3d_examples/extra_3d_PD_framework_tests/CMakeLists.txt',
        'tools/peri_coulomb_global.py','tools/peri_coulomb_global_controls.py','tools/peri_coulomb_global_audit.py',
        'tools/peri_coulomb_contact_audit.py','tools/test_peri_coulomb_contact.py',
        'tools/peri_pd_reference_audit.py','tools/peri_normal_projection_audit.py','tools/peri_dense_contact_audit.py',
        'tools/peri_many_contact_audit.py','tools/peri_elastic_normal.py',
        'docs/peri/COUPLED_GLOBAL_OBSERVATIONS_20261105.md']
    for file in files:
        source=repo/file;target=root/'frozen'/source.name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
    binary=root/'frozen'/executable.name;shutil.copy2(executable,binary)
    artifacts.extend([record(executable),record(binary),record(prospective)])
    build=repo/'build/peri-latest-20261103'
    obj='extra_3d_PD_framework_tests/CMakeFiles/extra_3d_PD_coulomb_contact_probe.dir/coulomb_contact_probe.cpp.o'
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
    print(json.dumps(dict(status='prepared_global_algebra_only_not_passed',cases=len(plan['controls']),source=plan['source_commit'],compiler_dependencies=len(paths))))


def limits():
    resource.setrlimit(resource.RLIMIT_CPU,(90,90));resource.setrlimit(resource.RLIMIT_AS,(2**30,2**30))
    resource.setrlimit(resource.RLIMIT_FSIZE,(64*2**20,64*2**20));resource.setrlimit(resource.RLIMIT_CORE,(0,0))


def records(path):
    with Path(path).open() as stream:return [json.loads(line) for line in stream if line.strip()]


def live_identity(pid):
    try:
        base=Path('/proc')/str(pid);fields=(base/'stat').read_text().rsplit(')',1)[1].split()
        return dict(pid=pid,ppid=int(fields[1]),start_ticks=fields[19],boot_id=identity()['boot_id'],
                    argv=(base/'cmdline').read_bytes().decode().strip('\0').split('\0'))
    except FileNotFoundError:return None


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
                observed=live_identity(child.pid)
                try:code=child.wait(timeout=90)
                except subprocess.TimeoutExpired:child.kill();code=child.wait();timeout=True
            attempt=dict(name=case['name'],command=command,return_code=code,timed_out=timeout,wall_seconds=time.monotonic()-start,
                stdout=record(out),stderr=record(err),live_observation=observed,
                max_child_RSS_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,last_accepted_time=0,accepted_physical_steps=0)
            attempts.append(attempt)
            assert not timeout and attempt['max_child_RSS_KiB']<=1024*1024
            raw=records(out);meta=raw[0]
            assert meta['argv']==command and meta['pid']==child.pid and meta['ppid']==identity()['pid']
            assert meta['boot_id']==identity()['boot_id'] and meta['start_ticks'].isdigit()
            if observed and observed['argv']==command:
                assert {k:meta[k] for k in observed}==observed
            assert err.stat().st_size==0
            attempt['runtime_gate']=inspect(case,raw,code)
        verify(prep['artifacts']+prep['inputs'])
    except Exception as error:failure=repr(error);failure_traceback=traceback.format_exc()
    save(root/'terminal.json',dict(status='runtime_completed_pending_independent_review' if failure is None else 'stopped_for_review',
        finished_at=stamp(),attempts=attempts,failure=failure,failure_traceback=failure_traceback))
    print(json.dumps(dict(status=read(root/'terminal.json')['status'],attempts=len(attempts),failure=failure)))
    return int(failure is not None)


def comparisons(raw):
    from decimal import localcontext
    with localcontext() as ctx:
        ctx.prec=95
        a=raw['physical_slide_two_components'][2]['solve'];b=raw['redundant_rows'][2]['solve']
        close_tree(a['velocity'],b['velocity']);close_tree(a['impulses'],[b['impulses'][k]+b['impulses'][3+k] for k in range(3)])
        for label in ['normalization_preserves_physical_mass','tangent_basis_covariance']:
            b=raw[label][2]['solve'];close_tree(a['velocity'],b['velocity'])
        b=raw['tangent_basis_covariance'][2]['solve'];close_tree(b['impulses'],[a['impulses'][0],a['impulses'][2],-a['impulses'][1]])
        b=raw['coordinate_frame_covariance'][2]['solve'];close_tree(b['velocity'],[a['velocity'][1],a['velocity'][2],a['velocity'][0]])
        close_tree(b['impulses'],a['impulses'])
    return True


def review(root):
    assert not (root/'independent_review.json').exists(),'immutable review already exists'
    resource.setrlimit(resource.RLIMIT_AS,(2*2**30,2*2**30));resource.setrlimit(resource.RLIMIT_CPU,(600,600))
    prep=read(root/'preparation.json');plan=read(root/'plan.json');terminal=read(root/'terminal.json')
    inputs=[record(root/name) for name in ['state.json','terminal.json','preparation.json']]
    result=dict(status='coupled_global_review_failed',recorded_at=stamp(),identity=identity());assessed=[]
    try:
        verify(inputs);verify(prep['artifacts']+prep['inputs'])
        assert terminal['status']=='runtime_completed_pending_independent_review' and len(terminal['attempts'])==len(plan['controls'])
        raw={}
        for case,attempt in zip(plan['controls'],terminal['attempts']):
            assert case['name']==attempt['name']
            verify([attempt['stdout'],attempt['stderr']]);data=records(attempt['stdout']['path']);raw[case['name']]=data
            assessed.append(inspect(case,data,attempt['return_code']))
        comparisons(raw);negative=[]
        controls_=[
            ('wrong_signed_transpose','full_joint_original_Jacobian'),
            ('symmetrized_original_tangent','full_joint_original_Jacobian'),
            ('dropped_normal_radius_derivative','full_joint_original_Jacobian'),
            ('square_box_instead_of_disk','full_joint_original_Jacobian'),
            ('friction_loss_wrong_sign','physical_slide_two_components'),
            ('premature_float_gather','exact_original_gather_cancellation'),
            ('forged_precision_acceptance','precision_rejected_before_iteration')]
        for label,name in controls_:
            data=copy.deepcopy(raw[name]);case=next(c for c in plan['controls'] if c['name']==name)
            r=data[2].get('original')
            if r is None:
                solve=data[2]['solve']
                r=solve['trace'][-1]['original'] if solve['trace'] else solve['precision_rejection']['original']
            if label=='wrong_signed_transpose':r['jacobian'][0][3]=1.
            elif label=='symmetrized_original_tangent':r['jacobian'][0][1]=r['jacobian'][1][0]=.125
            elif label=='dropped_normal_radius_derivative':r['jacobian'][4][3]=0.
            elif label=='square_box_instead_of_disk':
                r['residual'][4:6]=[0.,0.];r['disk'][0]['residual']=[0.,0.]
            elif label=='friction_loss_wrong_sign':r['disk'][0]['friction_loss']=-r['disk'][0]['friction_loss']
            elif label=='premature_float_gather':r['normalized_relative_velocity'][0][0]=0.
            else:r['precision_budget_ok']=True
            path=root/(label+'.negative.jsonl')
            with path.open('x') as stream:
                for row in data:stream.write(json.dumps(row,allow_nan=False)+'\n')
            try:inspect(case,data,case['expected_code'])
            except AssertionError:negative.append(dict(name=label,rejected=True,artifact=record(path)))
            else:raise AssertionError(('negative accepted',label))
        verify(inputs);verify(prep['artifacts']+prep['inputs'])
        result.update(status='original_simultaneous_global_normal_disk_coupling_independently_passed',cases=assessed,
            maximum_normalized_error=max(c.get('maximum_normalized_error',0) for c in assessed),negative_controls=negative,
            covariance_and_redundant_scatter_passed=True,accepted_physical_steps=0,
            scope='Bounded manufactured simultaneous global equations and analytic one-particle mass impacts; no actual PD trajectory or complete restart qualification.',
            actual_PD_friction_qualified=False,full_restart_qualified=False,S3_complete=False,calibration_allowed=False)
    except Exception as error:result.update(failure=repr(error),failure_traceback=traceback.format_exc(),completed_assessments=assessed)
    result.update(state=record(root/'state.json'),terminal=record(root/'terminal.json'),preparation=record(root/'preparation.json'))
    save(root/'independent_review.json',result)
    print(json.dumps(dict(status=result['status'],cases=len(assessed),failure=result.get('failure'))))
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
