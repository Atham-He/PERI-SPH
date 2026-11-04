"""Freeze/run/review full actual-PD normal observations; no old-state rerun."""
import argparse
import copy
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time
import traceback
from peri_elastic_normal import identity,read,record,reserves,save,stamp,verify
from peri_bounded_gel_audit import inspect,read_jsonl

def controls():
    def base(name,mode=1):
        return dict(name=name,mode=mode,steps=16,dt=.01,stroke=.01,rotation=0,azimuths=16,rho0=1,parameters=[5,1,.5,.3],
            eta=.1,support=1.25,cap_width=1/64,relaxed_fraction=0,perturb_node=16,invalid_injection=0,
            initial_affine=[1,0,0,0,1,0,0,0,1],velocity_rate=[0]*9,local_perturbation=[0]*3,local_velocity_perturbation=[0]*3,
            translation=[0]*3,expected_code=0)
    operator=base('mixed_nonaffine_operator',0)
    operator.update(steps=1,dt=.1,initial_affine=[1.08,.2,-.03,.01,.93,.04,0,-.02,1.02],
        velocity_rate=[.1,1.,-.3,.2,-.2,.4,-.1,.1,.1],relaxed_fraction=.25,
        local_perturbation=[.01,-.008,.006],local_velocity_perturbation=[.2,-.1,.15])
    free=base('actual_sphere_free_translation',2);free.update(steps=1,dt=.02,translation=[.1,-.2,.3])
    normal=base('normal_compression');rotated=base('normal_compression_rotated');rotated['rotation']=1
    cases=[operator,free,normal,rotated]
    for label,flag in [('invalid_mass',1),('invalid_material_cursor',2),('invalid_reference_cache',3),('invalid_signed_surface_cache',4)]:
        case=base(label);case.update(invalid_injection=flag,expected_code=1);cases.append(case)
    return cases

def encode(c):
    values=[c['mode'],c['steps'],c['dt'],c['stroke'],c['rotation'],c['azimuths'],c['rho0'],*c['parameters'],c['eta'],c['support'],
        c['cap_width'],c['relaxed_fraction'],c['perturb_node'],c['invalid_injection']]
    for field in ['initial_affine','velocity_rate','local_perturbation','local_velocity_perturbation','translation']:values.extend(c[field])
    return ' '.join(map(str,values))+'\n'

def prepare(root,executable,prospective,reuse_root=None):
    assert not root.exists(),'new unused study directory required'
    repo=Path(__file__).resolve().parents[1]
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=repo).strip(),'source must be committed'
    budget=reserves(root.parent);plan=read(prospective)
    memory={line.split(':')[0]:int(line.split()[1])*1024 for line in Path('/proc/meminfo').read_text().splitlines() if line.startswith(('MemTotal:','MemAvailable:'))}
    assert memory['MemAvailable']>=2*2**30,('memory_reserve',memory)
    budget.update(memory,affinity_cpus=sorted(os.sched_getaffinity(0)),native_threads=1)
    assert plan['normalized_physical_gate']==1e-10 and plan['dense_penetration_stroke_gate']==.01
    assert plan['steps']==16 and plan['duration']==.16 and plan['dt']==.01
    root.mkdir();(root/'frozen').mkdir();artifacts=[]
    files=['src/shared/peridynamics/pd_bounded_normal_scene.h','src/shared/peridynamics/pd_bounded_gel.h','src/shared/peridynamics/pd_fixed_reference.h',
        'src/shared/peridynamics/pd_surface_quadrature.h','src/shared/peridynamics/pd_gel_affine_path.h','src/shared/peridynamics/pd_gel_material.h',
        'src/shared/peridynamics/pd_gel_material.cpp','src/shared/peridynamics/pd_many_contact.h','src/shared/peridynamics/pd_exact_contact_balance.h',
        'src/shared/peridynamics/pd_nonsymmetric_constraints.h','src/shared/peridynamics/pd_dense_metric_contact.h',
        'tests/3d_examples/extra_3d_PD_framework_tests/bounded_gel_probe.cpp','tests/3d_examples/extra_3d_PD_framework_tests/bounded_gel_probe_io.h',
        'tests/3d_examples/extra_3d_PD_framework_tests/fixed_pd_reference_tests.cpp','tests/3d_examples/extra_3d_PD_framework_tests/bounded_normal_scene_tests.cpp',
        'tests/3d_examples/extra_3d_PD_framework_tests/CMakeLists.txt',
        'tools/peri_bounded_gel.py','tools/peri_bounded_gel_audit.py','tools/peri_pd_reference_audit.py','tools/peri_pd_material_reference.py',
        'tools/test_peri_bounded_gel.py','tools/peri_gel_affine_audit.py','tools/peri_normal_projection_audit.py',
        'tools/peri_dense_contact_audit.py','tools/peri_many_contact_audit.py','tools/peri_elastic_normal.py','docs/peri/BOUNDED_NORMAL_QUALIFICATION_20261104.md']
    for file in files:
        source=repo/file;target=root/'frozen'/source.name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
    frozen=root/'frozen'/executable.name;shutil.copy2(executable,frozen);artifacts.extend([record(executable),record(frozen),record(prospective)])
    build=repo/'build/peri-latest-20261103';obj='extra_3d_PD_framework_tests/CMakeFiles/extra_3d_PD_bounded_gel_probe.dir/bounded_gel_probe.cpp.o'
    commands=subprocess.check_output(['ninja','-C',str(build),'-t','commands',obj],text=True)
    assert not any(x in commands for x in ['-ffast-math','-Ofast','-ffinite-math-only'])
    (root/'frozen/compiler_command.txt').write_text(commands);artifacts.append(record(root/'frozen/compiler_command.txt'))
    dependencies=subprocess.check_output(['ninja','-C',str(build),'-t','deps',obj],text=True);paths=[]
    for line in dependencies.splitlines():
        if line.startswith('    '):
            p=Path(line.strip());p=p if p.is_absolute() else build/p
            if p.is_file():paths.append(p.resolve())
    paths=sorted(set(paths));assert paths and sum(p.stat().st_size for p in paths)<512*2**20
    manifest=[]
    for index,source in enumerate(paths):
        target=root/'frozen/compiler_headers'/str(index)/source.name;target.parent.mkdir(parents=True);shutil.copy2(source,target)
        entry=dict(original=record(source),frozen=record(target));manifest.append(entry);artifacts.extend(entry.values())
    save(root/'frozen/compiler_dependencies.json',manifest);artifacts.append(record(root/'frozen/compiler_dependencies.json'))
    import _decimal,math
    artifacts.append(record(Path(sys.executable).resolve()))
    modules={}
    for module in [_decimal,math]:
        origin=module.__spec__.origin;modules[module.__name__]=origin
        if getattr(module,'__file__',None):artifacts.append(record(Path(module.__file__)))
        else:assert origin in ['built-in','frozen'],('unidentified interpreter module',module.__name__,origin)
    save(root/'frozen/interpreter_modules.json',modules);artifacts.append(record(root/'frozen/interpreter_modules.json'))
    for line in subprocess.check_output(['ldd',str(executable)],text=True).splitlines():
        for token in line.split():
            if token.startswith('/') and Path(token).is_file():artifacts.append(record(token))
    if plan.get('interval_preflight'):artifacts.append(record(plan['interval_preflight']))
    plan.update(controls=controls(),prepared_at=stamp(),identity=identity(),resources=budget,python_version=sys.version,
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip())
    if reuse_root is not None:
        # New assessor qualification may use existing complete raw. Never
        # execute an already attempted physical input or alter its old state.
        oldprep=read(reuse_root/'preparation.json');oldplan=read(reuse_root/'plan.json')
        oldterminal=read(reuse_root/'terminal.json');oldstate=read(reuse_root/'state.json')
        assert oldterminal['status']=='stopped_for_review' and oldterminal['failure']=="AssertionError('invalid gel material history')"
        assert record(frozen)['sha256']==record(oldprep['executable'])['sha256']
        historical=[entry for entry in oldprep['artifacts'] if Path(entry['path']).is_relative_to(reuse_root)]
        historical.extend(oldprep['inputs']);historical.extend(record(reuse_root/name) for name in ['preparation.json','state.json','plan.json','terminal.json'])
        reuse={}
        for oldattempt in oldterminal['attempts']:
            case=next(c for c in plan['controls'] if c['name']==oldattempt['name'])
            original=next(c for c in oldplan['controls'] if c['name']==oldattempt['name'])
            assert case==original and (reuse_root/(case['name']+'.input')).read_text()==encode(case)
            assert not oldattempt['timed_out'] and oldattempt['return_code']==case['expected_code']
            historical.extend([oldattempt['stdout'],oldattempt['stderr']]);meta=read_jsonl(oldattempt['stdout']['path'])[0]
            assert meta['argv']==oldattempt['command'] and meta['ppid']==oldstate['identity']['pid']
            assert meta['boot_id']==oldstate['identity']['boot_id'] and meta['start_ticks'].isdigit()
            reuse[case['name']]=dict(oldattempt,reused_from=str(reuse_root),historical_source_commit=oldplan['source_commit'],reused_native_dispatch=False)
        assert len(reuse)==6 and len(plan['controls'])-len(reuse)==2
        verify(historical);artifacts.extend(historical);plan['reused_controls']=reuse
    save(root/'plan.json',plan);artifacts.append(record(root/'plan.json'));inputs=[]
    for case in plan['controls']:
        path=root/(case['name']+'.input');path.write_text(encode(case));inputs.append(record(path))
    save(root/'preparation.json',dict(artifacts=artifacts,inputs=inputs,executable=str(frozen),compiler_headers=len(paths)))
    print(json.dumps(dict(status='prepared_not_passed',controls=len(plan['controls']),source_commit=plan['source_commit'],compiler_headers=len(paths))))

def limits():
    resource.setrlimit(resource.RLIMIT_CPU,(120,120));resource.setrlimit(resource.RLIMIT_AS,(2**30,2**30))
    resource.setrlimit(resource.RLIMIT_FSIZE,(256*2**20,256*2**20));resource.setrlimit(resource.RLIMIT_CORE,(0,0))

def run(root):
    assert not (root/'state.json').exists(),'old state cannot be redispatched'
    prep=read(root/'preparation.json');plan=read(root/'plan.json');attempts=[];failure=None;failure_traceback=None
    save(root/'state.json',dict(status='attempted',started_at=stamp(),identity=identity(),preparation=record(root/'preparation.json')))
    try:
        verify(prep['artifacts']+prep['inputs'])
        for case in plan['controls']:
            reserves(root);out=root/(case['name']+'.jsonl');err=root/(case['name']+'.stderr');start=time.monotonic();timed_out=False
            command=[prep['executable'],str(root/(case['name']+'.input'))]
            if case['name'] in plan.get('reused_controls',{}):
                attempt=copy.deepcopy(plan['reused_controls'][case['name']]);out=Path(attempt['stdout']['path']);code=attempt['return_code']
                verify([attempt['stdout'],attempt['stderr']])
            else:
                with out.open('xb') as stdout,err.open('xb') as stderr:
                    child=subprocess.Popen(command,stdout=stdout,stderr=stderr,preexec_fn=limits,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
                    try:code=child.wait(timeout=120)
                    except subprocess.TimeoutExpired:child.kill();code=child.wait();timed_out=True
                attempt=dict(name=case['name'],command=command,return_code=code,timed_out=timed_out,wall_seconds=time.monotonic()-start,
                    stdout=record(out),stderr=record(err),cumulative_child_maxRSS_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
            attempts.append(attempt);assert not timed_out and attempt['cumulative_child_maxRSS_KiB']<=1024*1024
            # Preserve last COMPLETE committed cursor even when a later record
            # or gate is invalid. Never count a candidate/rejected step.
            prefix=[]
            with out.open() as stream:
                for line in stream:
                    try:entry=json.loads(line)
                    except (ValueError,TypeError):break
                    prefix.append(entry)
            accepted=[r for r in prefix if r['type']=='step' and r['committed']]
            attempt['last_complete_accepted_time']=accepted[-1]['state']['time'] if accepted else 0
            attempt['last_complete_accepted_steps']=len(accepted);attempt['last_complete_record']=prefix[-1]['type'] if prefix else None
            records=read_jsonl(out);meta=records[0]
            assert meta['argv']==attempt['command'] and meta['start_ticks'].isdigit()
            if 'reused_from' not in attempt:
                assert meta['pid']==child.pid and meta['ppid']==identity()['pid'] and meta['boot_id']==identity()['boot_id']
            attempt['runtime_gate']=inspect(case,records,code,False)
        verify(prep['artifacts']+prep['inputs'])
    except Exception as e:failure=repr(e);failure_traceback=traceback.format_exc()
    save(root/'terminal.json',dict(status='runtime_completed_pending_independent_review' if failure is None else 'stopped_for_review',
        finished_at=stamp(),attempts=attempts,failure=failure,failure_traceback=failure_traceback))
    print(json.dumps(dict(status=read(root/'terminal.json')['status'],attempts=len(attempts),failure=failure)))
    return int(failure is not None)

def review(root):
    assert not (root/'independent_review.json').exists(),'immutable review exists'
    resource.setrlimit(resource.RLIMIT_AS,(8*2**30,8*2**30));resource.setrlimit(resource.RLIMIT_CPU,(3600,3600))
    prep=read(root/'preparation.json');plan=read(root/'plan.json');terminal=read(root/'terminal.json');assessed=[]
    terminal_inputs=[record(root/name) for name in ['state.json','terminal.json','preparation.json']]
    result=dict(status='bounded_PD_normal_review_failed',recorded_at=stamp(),identity=identity())
    try:
        verify([read(root/'state.json')['preparation']]);verify(prep['artifacts']+prep['inputs'])
        assert terminal['status']=='runtime_completed_pending_independent_review' and len(terminal['attempts'])==len(plan['controls'])
        for case,a in zip(plan['controls'],terminal['attempts']):
            verify([a['stdout'],a['stderr']]);assessed.append(inspect(case,read_jsonl(a['stdout']['path']),a['return_code'],True))
        raw_paths={a['name']:a['stdout']['path'] for a in terminal['attempts']}
        base=read_jsonl(raw_paths['normal_compression']);rotated=read_jsonl(raw_paths['normal_compression_rotated'])
        a=[r for r in base if r['type']=='step'];b=[r for r in rotated if r['type']=='step']
        assert len(a)==len(b)==16
        for x,y in zip(a,b):
            for field in ['position','velocity']:
                expected=[[-v[1],v[0],v[2]] for v in x['state'][field]]
                assert max(abs(p-q) for row,other in zip(expected,y['state'][field]) for p,q in zip(row,other))<=1e-10
            for field in ['particle_plate_impulse','plate_reaction_interval']:
                expected=[[-v[1],v[0],v[2]] for v in x[field]]
                assert max(abs(p-q) for row,other in zip(expected,y[field]) for p,q in zip(row,other))<=1e-10
            for field in ['kinetic_new','algorithmic_loss','normal_contact_loss','energy_residual']:
                assert abs(x[field]-y[field])<=1e-10
        negative=[]
        for label,name in [('collapsed_global_reference','mixed_nonaffine_operator'),('clipped_signed_MLS','mixed_nonaffine_operator'),
            ('symmetrized_true_Jacobian','mixed_nonaffine_operator'),('effective_J_as_kinetic_mass','mixed_nonaffine_operator'),
            ('swapped_material_algorithmic_loss','normal_compression')]:
            records=copy.deepcopy(read_jsonl(raw_paths[name]));case=next(c for c in plan['controls'] if c['name']==name)
            if label=='collapsed_global_reference':records[1]['reference']['shape_gradient']=[[0,0,0] for _ in range(33*33)]
            elif label=='clipped_signed_MLS':records[1]['surface']['weights']=[max(0,w) for w in records[1]['surface']['weights']]
            elif label=='symmetrized_true_Jacobian':
                j=records[2]['evaluation']['jacobian'];records[2]['evaluation']['jacobian']=[[(j[i][k]+j[k][i])/2 for k in range(99)] for i in range(99)]
            elif label=='effective_J_as_kinetic_mass':
                j=records[2]['evaluation']['jacobian'];v=records[2]['velocity'];records[2]['trial_kinetic_actual_mass']=sum(v[i]*j[i][k]*v[k]/2 for i in range(99) for k in range(99))
            else:
                row=next(r for r in records if r['type']=='step');row['algorithmic_loss'],row['evaluation']['physical_material_loss']=row['evaluation']['physical_material_loss'],row['algorithmic_loss']
            path=root/(label+'.negative.jsonl')
            with path.open('x') as stream:
                for row in records:stream.write(json.dumps(row,allow_nan=False)+'\n')
            try:inspect(case,records,0,False)
            except AssertionError:negative.append(dict(name=label,rejected=True,artifact=record(path)))
            else:raise AssertionError(('negative control accepted',label))
        verify(terminal_inputs);verify(prep['artifacts']+prep['inputs'])
        result.update(status='actual_bounded_PD_normal_reference_material_operator_trajectory_independently_passed',cases=assessed,negative_controls=negative,
            frame_trajectory_comparison=True,individual_redundant_dual_equality_required=False,normal_only=True,
            scope='one fixed33node/114sample228row normal sphere with actual all-point histories and true nonsymmetric operator; not time/space/horizon or full friction/restart qualification',
            S3_complete=False,full3D_friction_qualified=False,full_restart_qualified=False,calibration_allowed=False)
    except Exception as e:result.update(failure=repr(e),failure_traceback=traceback.format_exc(),completed_assessments=assessed)
    result.update(state=record(root/'state.json'),terminal=record(root/'terminal.json'),preparation=record(root/'preparation.json'))
    save(root/'independent_review.json',result)
    print(json.dumps(dict(status=result['status'],cases=len(assessed),failure=result.get('failure'))))
    return int(not result['status'].endswith('_passed'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run','review']);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--executable',type=Path);p.add_argument('--prospective',type=Path);p.add_argument('--reuse-root',type=Path);args=p.parse_args()
    if args.action=='prepare':
        fresh=not args.root.exists()
        try:prepare(args.root,args.executable,args.prospective,args.reuse_root)
        except Exception as error:
            if fresh and args.root.exists():save(args.root/'preparation_error.json',dict(status='preparation_failed_no_native_dispatch',
                identity=identity(),recorded_at=stamp(),failure=repr(error),failure_traceback=traceback.format_exc()))
            raise
    else:sys.exit(globals()[args.action](args.root))
