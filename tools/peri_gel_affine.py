"""Frozen material-point path/tangent controls; S3 and global solver remain gated."""
import argparse
import copy
import json
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time
from decimal import Decimal as D
from peri_elastic_normal import identity, read, record, reserves, save, stamp, verify
from peri_gel_affine_audit import check,green,multiply,sub,scale


def cases():
    identity_f=[1,0,0,0,1,0,0,0,1]
    shear=[1,.2,0,0,1,0,0,0,1]
    mixed=[1.08,.2,-.03,.01,.93,.04,0,-.02,1.02]
    rotation=[0,-1,0,1,0,0,0,0,1]
    prehistory=[.01,.002,0,.002,-.01,0,0,0,0]
    result=[]
    def add(name,final=None,dt=.15,initial=None,viscous=None,branch=500,old_loss=0,steps=None,expected=0):
        f0=initial if initial is not None else identity_f
        result.append(dict(name=name,parameters=[5000,1000,branch,.3],initial=f0,
            green=list(map(float,green(list(map(lambda v:D(str(v)),f0))))),
            viscous=viscous if viscous is not None else [0]*9,old_dissipation=old_loss,
            steps=steps if steps is not None else [dict(dt=dt,F=final if final is not None else f0)],expected_code=expected))
    add('identity_hold'); add('isochoric_shear',shear); add('mixed',mixed)
    add('small_relaxation_increment',mixed,dt=1e-10); add('large_relaxation_increment',mixed,dt=.6)
    add('volumetric_compression',[.9,0,0,0,.9,0,0,0,.9])
    add('volumetric_extension',[1.1,0,0,0,1.1,0,0,0,1.1])
    add('relaxation',mixed,initial=mixed)
    add('prehistory',mixed,viscous=prehistory,old_loss=3.25)
    half=[1,.1,0,0,1,0,0,0,1]
    add('ramp_hold_release',steps=[dict(dt=dt,F=f) for dt,f in [(.1,half),(.1,shear),(.2,shear),(.1,half),(.1,identity_f),(.2,identity_f)]])
    add('zero_branch',mixed,branch=0)
    rotated=list(map(float,multiply(list(map(D,rotation)),list(map(lambda v:D(str(v)),mixed)))))
    add('rotated_mixed',rotated,initial=rotation)
    add('rigid_hold',initial=rotation)
    add('no_advance',dt=0,viscous=prehistory,old_loss=3.25)
    add('interior_determinant_failure',[-1,0,0,0,-1,0,0,0,1],expected=1)
    add('history_mismatch',expected=1); result[-1]['green']=[.01,0,0,0,-.01,0,0,0,0]
    add('negative_dt',dt=-.1,expected=1)
    add('invalid_tau',expected=1); result[-1]['parameters'][3]=0
    return result


def prepare(root,executable,prospective):
    assert not root.exists(),'no study directory reuse'
    repo=Path(__file__).resolve().parents[1]
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=repo).strip(),'freeze committed source'
    resources=reserves(root.parent); plan=read(prospective); controls=cases()
    root.mkdir(); (root/'frozen').mkdir(); artifacts=[]
    for source in ['src/shared/peridynamics/pd_gel_affine_path.h','src/shared/peridynamics/pd_gel_material.h',
        'src/shared/peridynamics/pd_gel_material.cpp','tests/3d_examples/extra_3d_PD_framework_tests/gel_affine_path_probe.cpp',
        'tests/3d_examples/extra_3d_PD_framework_tests/gel_affine_path_tests.cpp',
        'tests/3d_examples/extra_3d_PD_framework_tests/CMakeLists.txt','tools/peri_gel_affine.py',
        'tools/peri_gel_affine_audit.py','tools/peri_elastic_normal.py','tools/test_peri_gel_affine.py',
        'docs/peri/GEL_AFFINE_PATH_STAGE_S3B_20261104.md']:
        original=repo/source; frozen=root/'frozen'/original.name
        shutil.copy2(original,frozen); artifacts.extend([record(original),record(frozen)])
    frozen_exe=root/'frozen'/executable.name; shutil.copy2(executable,frozen_exe)
    artifacts.extend([record(executable),record(frozen_exe),record(prospective)])
    for line in subprocess.check_output(['ldd',str(executable)],text=True).splitlines():
        for token in line.split():
            if token.startswith('/') and Path(token).is_file(): artifacts.append(record(token))
    plan.update(cases=controls,prepared_at=stamp(),identity=identity(),resources=resources,
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip())
    save(root/'plan.json',plan); artifacts.append(record(root/'plan.json')); inputs=[]
    for case in controls:
        values=[*case['parameters'],*case['initial'],*case['green'],*case['viscous'],case['old_dissipation'],len(case['steps'])]
        for step in case['steps']: values.extend([step['dt'],*step['F']])
        path=root/(case['name']+'.input'); path.write_text(' '.join(map(str,values))+'\n'); inputs.append(record(path))
    save(root/'preparation.json',dict(artifacts=artifacts,inputs=inputs,executable=str(frozen_exe)))
    print(json.dumps(dict(status='prepared_not_passed',cases=len(controls),source_commit=plan['source_commit'])))


def limits():
    resource.setrlimit(resource.RLIMIT_CPU,(5,5)); resource.setrlimit(resource.RLIMIT_AS,(256*2**20,256*2**20))
    resource.setrlimit(resource.RLIMIT_CORE,(0,0)); resource.setrlimit(resource.RLIMIT_FSIZE,(32*2**20,32*2**20))


def parse(path): return [json.loads(line) for line in Path(path).read_text().splitlines()]


def run(root):
    assert not (root/'state.json').exists(),'attempted state cannot be redispatched'
    prep=read(root/'preparation.json'); plan=read(root/'plan.json'); attempts=[]; failure=None
    save(root/'state.json',dict(status='attempted',started_at=stamp(),identity=identity(),preparation=record(root/'preparation.json')))
    try:
        verify(prep['artifacts']+prep['inputs'])
        for case in plan['cases']:
            reserves(root); out=root/(case['name']+'.jsonl'); err=root/(case['name']+'.stderr')
            command=[prep['executable'],str(root/(case['name']+'.input'))]; started=time.monotonic(); timeout=False
            with out.open('xb') as stdout,err.open('xb') as stderr:
                child=subprocess.Popen(command,stdout=stdout,stderr=stderr,preexec_fn=limits)
                try: code=child.wait(timeout=5)
                except subprocess.TimeoutExpired: child.kill(); code=child.wait(); timeout=True
            attempt=dict(name=case['name'],command=command,return_code=code,timed_out=timeout,wall_seconds=time.monotonic()-started,
                stdout=record(out),stderr=record(err),cumulative_child_maxRSS_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
            attempts.append(attempt)
            assert not timeout and attempt['cumulative_child_maxRSS_KiB']<=256*1024
            raw=parse(out); meta=raw[0]
            assert meta['type']=='identity' and meta['pid']==child.pid and meta['ppid']==identity()['pid']
            assert meta['argv']==command and meta['boot_id']==identity()['boot_id']
            attempt['runtime_gate']=check(case,raw,code)
        verify(prep['artifacts']+prep['inputs'])
    except Exception as error: failure=repr(error)
    save(root/'terminal.json',dict(status='runtime_completed_pending_independent_review' if failure is None else 'stopped_for_review',
        finished_at=stamp(),attempts=attempts,failure=failure))
    print(json.dumps(dict(status=read(root/'terminal.json')['status'],attempts=len(attempts),failure=failure)))
    return int(failure is not None)


def review(root):
    assert not (root/'independent_review.json').exists(),'do not repeat receipt'
    prep=read(root/'preparation.json'); plan=read(root/'plan.json'); terminal=read(root/'terminal.json')
    result=dict(status='S3b_material_path_independent_review_failed',recorded_at=stamp(),identity=identity()); assessed=[]
    try:
        verify([read(root/'state.json')['preparation']]); verify(prep['artifacts']+prep['inputs'])
        assert terminal['status']=='runtime_completed_pending_independent_review' and len(terminal['attempts'])==len(plan['cases'])
        for case,attempt in zip(plan['cases'],terminal['attempts']):
            verify([attempt['stdout'],attempt['stderr']]); assessed.append(check(case,parse(attempt['stdout']['path']),attempt['return_code']))
        by_name={c['name']:c for c in assessed}
        base=next(r for r in parse(root/'mixed.jsonl') if r['type']=='step')
        rotated=next(r for r in parse(root/'rotated_mixed.jsonl') if r['type']=='step')
        R=[0,-1,0,1,0,0,0,0,1]
        wanted=[float(v) for v in multiply(list(map(D,R)),list(map(lambda x:D(str(x)),base['average_pk1'])))]
        assert max(abs(a-b) for a,b in zip(wanted,rotated['average_pk1']))/5000<1e-10
        for field in ['green','viscous']:
            assert max(abs(a-b) for a,b in zip(base[field],rotated[field]))<1e-10
        for field in ['work','physical_dissipation','stored_new']:
            assert abs(base[field]-rotated[field])/5000<1e-10
        assert by_name['zero_branch']['maximum_consistent_tangent_asymmetry']<1e-10
        # Nonsymmetry is a measured admissibility result, not a failed material
        # accuracy gate or permission to symmetrize the true global Jacobian.
        asymmetry=by_name['mixed']['maximum_consistent_tangent_asymmetry']
        assert asymmetry>1e-8
        controls=[]; case=next(c for c in plan['cases'] if c['name']=='isochoric_shear')
        bad=copy.deepcopy(parse(root/'isochoric_shear.jsonl')); row=bad[1]
        # Wrong linear-Green update at the same F endpoints must be rejected.
        import math
        h=case['steps'][0]['dt']/case['parameters'][3]; alpha=-math.expm1(-h)/h
        endpoint=list(map(D,map(str,case['steps'][0]['F'])))
        wrong=list(map(float,scale(green(endpoint),D(str(1-alpha)))))
        row['viscous']=wrong
        save(root/'linear_Green_path.negative.json',bad)
        try: check(case,bad,0)
        except AssertionError: controls.append(dict(name='linear_Green_substitution',rejected=True,artifact=record(root/'linear_Green_path.negative.json')))
        else: raise AssertionError('wrong kinematic path accepted')
        case=next(c for c in plan['cases'] if c['name']=='mixed'); bad=copy.deepcopy(parse(root/'mixed.jsonl'))
        tangent=bad[1]['tangent']
        bad[1]['tangent']=[[(tangent[j][i]+tangent[i][j])/2 for i in range(9)] for j in range(9)]
        save(root/'symmetrized_tangent.negative.json',bad)
        try: check(case,bad,0)
        except AssertionError: controls.append(dict(name='consistent_tangent_symmetrization',rejected=True,artifact=record(root/'symmetrized_tangent.negative.json')))
        else: raise AssertionError('wrong symmetric tangent accepted')
        result.update(status='S3b_material_path_and_local_tangent_independently_passed',cases=assessed,negative_controls=controls,
            measured_mixed_tangent_asymmetry=asymmetry,SPD_global_solver_admissible_for_this_exact_tangent=False,
            next_requirement='nonlinear constrained residual/global work-conjugate PD operator preserving the nonsymmetric consistent tangent',
            S3b_global_operator_qualified=False,S3_complete=False,production_integrated=False,calibration_allowed=False)
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
