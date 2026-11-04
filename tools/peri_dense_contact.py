"""Bounded S3a full-matrix controls. New freeze/state; no old solver reruns."""
import argparse
import copy
import json
from pathlib import Path
import random
import resource
import shutil
import subprocess
import sys
import time
from peri_elastic_normal import identity, read, record, reserves, save, stamp, verify
from peri_dense_contact_audit import check, exact


def cases():
    result=[]
    def add(name,matrix,rows=(),mode=0,mass=None,position=None,velocity=None,dt=.1,expected_code=0,sweeps=10000,V=1):
        n=len(matrix)
        result.append(dict(name=name,matrix=matrix,rows=list(rows),mode=mode,mass=mass or [1]*n,
            position=position or [0]*n,velocity=velocity or [0]*n,dt=dt,expected_code=expected_code,
            max_sweeps=sweeps,velocity_scale=V,tolerance=1e-12))
    row=lambda C,new,old=0:dict(C=C,new=new,old=old)
    add('offdiagonal_response',[[2,1],[1,2]],[row([1,0],1)],dt=1)
    add('two_active_rows',[[2,1],[1,2]],[row([1,0],1),row([0,1],.5)],dt=1)
    add('signed_MLS',[[2,.5],[.5,1]],[row([1.5,-.5],.02)])
    add('redundant_rows',[[2,.5],[.5,1]],[row([1,0],.01),row([1,0],.01)])
    a=[[0.0]*6 for _ in range(6)]
    for i in range(6):
        a[i][i]=2
        if i: a[i][i-1]=a[i-1][i]=-.25
    add('flattened_two_node_3D',a,[row([0,.75,0,0,.25,0],.02),row([.24,0,.32,.36,0,.48],.01)])
    add('small_SI_metric',[[2e-9,1e-9],[1e-9,2e-9]],[row([1,0],1e-6)],dt=1e-3,V=.01)
    add('coupled_elastic',[[2,-1],[-1,2]],[row([1,0],.01)],mode=1,mass=[.01,.02])
    add('free_elastic',[[2,-1],[-1,2]],mode=1,mass=[.01,.02],position=[.1,-.05],velocity=[.2,-.1],dt=.05)
    add('rigid_mode',[[1,-1],[-1,1]],mode=1,mass=[.01,.02],velocity=[.1,.1],dt=.05)
    add('positive_old_gap',[[2,-1],[-1,2]],[row([1,0],.03)],mode=1,mass=[.01,.02],position=[.02,0])
    add('signed_elastic',[[2,-1],[-1,2]],[row([1.5,-.5],.01)],mode=1,mass=[.01,.02])
    generator=random.Random(20261104)
    for index in range(10):
        b=[[generator.randint(-2,2) for _ in range(4)] for _ in range(4)]
        k=[[sum(b[s][i]*b[s][j] for s in range(4))+(1 if i==j else 0) for j in range(4)] for i in range(4)]
        c=[[1 if i==j else generator.choice([-.25,0,.25]) for i in range(4)] for j in range(3)]
        add(f'deterministic_coupled_{index}',k,[row(coefficients,.001*(i+1)) for i,coefficients in enumerate(c)],
            mode=1,mass=[.01,.02,.03,.04],dt=.025)
    add('infeasible',[[1]],[row([1],.1),row([-1],.1)],expected_code=2,sweeps=64)
    add('nonSPD',[[1,2],[2,1]],[row([1,0],.01)],expected_code=1)
    add('singular',[[.1,.1],[.1,.1]],[row([1,0],.01)],expected_code=1)
    add('asymmetric',[[2,1],[0,2]],[row([1,0],.01)],expected_code=1)
    add('nonPSD_stiffness',[[1,2],[2,1]],mode=1,mass=[.01,.02],expected_code=1)
    add('zero_row',[[1]],[row([0],0)],expected_code=1)
    add('initial_penetration',[[1]],[row([1],.02,old=.01)],expected_code=1)
    return result


def prepare(root,executable,prospective):
    assert not root.exists(),'never reuse old study directory'
    repo=Path(__file__).resolve().parents[1]
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=repo).strip()
    resources=reserves(root.parent); plan=read(prospective); new_cases=cases()
    # Classify controls with exact arithmetic before dispatch. No native output
    # is supplied as a reference state or adjusted expected result.
    for case in new_cases:
        try: reference=exact(case); code=0 if reference is not None else 2
        except ValueError: code=1
        assert code==case['expected_code'],(case['name'],'prospective classification')
    root.mkdir(); (root/'frozen').mkdir(); artifacts=[]
    sources=['src/shared/peridynamics/pd_dense_metric_contact.h','src/shared/peridynamics/pd_normal_projection.h',
        'tests/3d_examples/extra_3d_PD_framework_tests/dense_metric_probe.cpp',
        'tests/3d_examples/extra_3d_PD_framework_tests/dense_metric_contact_tests.cpp',
        'tests/3d_examples/extra_3d_PD_framework_tests/CMakeLists.txt',
        'tools/peri_dense_contact.py','tools/peri_dense_contact_audit.py','tools/peri_normal_projection_audit.py',
        'tools/peri_elastic_normal.py','tools/test_peri_dense_contact.py','docs/peri/DENSE_METRIC_STAGE_S3A_20261104.md']
    for source in sources:
        original=repo/source; target=root/'frozen'/original.name
        shutil.copy2(original,target); artifacts.extend([record(original),record(target)])
    target=root/'frozen'/executable.name; shutil.copy2(executable,target)
    artifacts.extend([record(executable),record(target),record(prospective)])
    for line in subprocess.check_output(['ldd',str(executable)],text=True).splitlines():
        for token in line.split():
            if token.startswith('/') and Path(token).is_file(): artifacts.append(record(token))
    plan.update(cases=new_cases,prepared_at=stamp(),identity=identity(),resources=resources,
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip())
    save(root/'plan.json',plan); artifacts.append(record(root/'plan.json')); inputs=[]
    for case in new_cases:
        values=[case['mode'],len(case['mass']),len(case['rows']),case['dt'],case['velocity_scale'],case['tolerance'],case['max_sweeps']]
        for row in case['matrix']: values.extend(row)
        for name in ['mass','position','velocity']: values.extend(case[name])
        for row in case['rows']: values.extend([row['old'],row['new'],*row['C']])
        path=root/(case['name']+'.input'); path.write_text(' '.join(map(str,values))+'\n'); inputs.append(record(path))
    save(root/'preparation.json',dict(artifacts=artifacts,inputs=inputs,executable=str(target)))
    print(json.dumps(dict(status='prepared_not_passed',cases=len(new_cases),source_commit=plan['source_commit'])))


def limits():
    resource.setrlimit(resource.RLIMIT_CPU,(5,5)); resource.setrlimit(resource.RLIMIT_AS,(256*2**20,256*2**20))
    resource.setrlimit(resource.RLIMIT_CORE,(0,0)); resource.setrlimit(resource.RLIMIT_FSIZE,(32*2**20,32*2**20))


def run(root):
    assert not (root/'state.json').exists(),'attempted state cannot be redispatched'
    prep=read(root/'preparation.json'); plan=read(root/'plan.json'); attempts=[]; failure=None
    save(root/'state.json',dict(started_at=stamp(),status='attempted',identity=identity(),preparation=record(root/'preparation.json')))
    try:
        verify(prep['artifacts']+prep['inputs'])
        for case in plan['cases']:
            reserves(root); out=root/(case['name']+'.json'); err=root/(case['name']+'.stderr')
            command=[prep['executable'],str(root/(case['name']+'.input'))]; started=time.monotonic(); timed_out=False
            with out.open('xb') as stdout,err.open('xb') as stderr:
                child=subprocess.Popen(command,stdout=stdout,stderr=stderr,preexec_fn=limits)
                try: code=child.wait(timeout=5)
                except subprocess.TimeoutExpired: child.kill(); code=child.wait(); timed_out=True
            attempt=dict(name=case['name'],command=command,return_code=code,stdout=record(out),stderr=record(err),
                wall_seconds=time.monotonic()-started,timed_out=timed_out,
                cumulative_child_maxRSS_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
            attempts.append(attempt)
            assert not timed_out and attempt['cumulative_child_maxRSS_KiB']<=256*1024
            native=read(out)
            assert native['identity']['pid']==child.pid and native['identity']['ppid']==identity()['pid']
            assert native['identity']['argv']==command and native['identity']['boot_id']==identity()['boot_id']
            attempt['runtime_gate']=check(case,native,code)
        verify(prep['artifacts']+prep['inputs'])
    except Exception as error: failure=repr(error)
    save(root/'terminal.json',dict(status='runtime_completed_pending_independent_review' if failure is None else 'stopped_for_review',
        finished_at=stamp(),attempts=attempts,failure=failure))
    print(json.dumps(dict(status=read(root/'terminal.json')['status'],attempts=len(attempts),failure=failure)))
    return int(failure is not None)


def review(root):
    assert not (root/'independent_review.json').exists(),'receipt already exists'
    prep=read(root/'preparation.json'); plan=read(root/'plan.json'); terminal=read(root/'terminal.json')
    result=dict(status='S3a_independent_review_failed',recorded_at=stamp(),identity=identity()); assessed=[]
    try:
        verify([read(root/'state.json')['preparation']]); verify(prep['artifacts']+prep['inputs'])
        assert terminal['status']=='runtime_completed_pending_independent_review'
        assert len(terminal['attempts'])==len(plan['cases'])
        for case,attempt in zip(plan['cases'],terminal['attempts']):
            verify([attempt['stdout'],attempt['stderr']]); assert case['name']==attempt['name']
            assessed.append(check(case,read(attempt['stdout']['path']),attempt['return_code']))
        negative=[]
        case=plan['cases'][0]; bad=copy.deepcopy(read(root/(case['name']+'.json')))
        diagonal=copy.deepcopy(case); diagonal['matrix']=[[a if i==j else 0 for j,a in enumerate(row)] for i,row in enumerate(case['matrix'])]
        wrong=exact(diagonal); bad['velocity']=list(map(float,wrong['velocity'])); bad['normal_impulse']=list(map(float,wrong['normal_impulse']))
        for name,example,raw in [('diagonal_operator_substitution',case,bad)]:
            save(root/(name+'.negative.json'),raw)
            try: check(example,raw,0)
            except AssertionError: negative.append(dict(name=name,rejected=True,artifact=record(root/(name+'.negative.json'))))
            else: raise AssertionError('diagonal negative control accepted')
        case=next(c for c in plan['cases'] if c['name']=='coupled_elastic'); bad=copy.deepcopy(read(root/'coupled_elastic.json'))
        metric=exact(case)['metric']; velocity=bad['velocity']
        bad['kinetic_new']=sum(float(metric[i][j])*velocity[i]*velocity[j]/2 for i in range(len(velocity)) for j in range(len(velocity)))
        save(root/'effective_mass.negative.json',bad)
        try: check(case,bad,0)
        except AssertionError: negative.append(dict(name='effective_metric_as_physical_mass',rejected=True,artifact=record(root/'effective_mass.negative.json')))
        else: raise AssertionError('physical-mass negative control accepted')
        result.update(status='S3a_full_metric_algebra_independently_passed',cases=assessed,negative_controls=negative,
            scope='bounded dense coupled-node normal-only quadratic research controls',
            S3_complete=False,nonlinear_gel_integration_qualified=False,production_restart_qualified=False,calibration_allowed=False)
    except Exception as error: result.update(failure=repr(error),completed_assessments=assessed)
    result['terminal']=record(root/'terminal.json'); result['state']=record(root/'state.json'); result['preparation']=record(root/'preparation.json')
    save(root/'independent_review.json',result)
    print(json.dumps(dict(status=result['status'],assessed_cases=len(assessed),failure=result.get('failure'))))
    return int(not result['status'].endswith('_passed'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('action',choices=['prepare','run','review']); parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--executable',type=Path); parser.add_argument('--prospective',type=Path); args=parser.parse_args()
    if args.action=='prepare': prepare(args.root,args.executable,args.prospective)
    else: sys.exit(globals()[args.action](args.root))
