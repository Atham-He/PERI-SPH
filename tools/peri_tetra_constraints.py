"""Freeze/run/single independent review of bounded nonsymmetric global controls."""
import argparse
import copy
import json
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time
import traceback
from decimal import Decimal as D
from peri_elastic_normal import identity,read,record,reserves,save,stamp,verify
from peri_gel_affine_audit import green
from peri_tetra_constraints_audit import inspect,gradient,mv


def cases():
    reference=[0,0,0,1,0,0,0,1,0,0,0,1]
    mixed=[0,0,0,1.08,.01,0,.2,.93,-.02,-.03,.04,1.02]
    lower=[0]*12;lower[1]=1
    upper=[0]*12;upper[7]=-1
    def walls(height,stroke): return [dict(C=lower,old=0,new=0),dict(C=upper,old=-height,new=-height+stroke)]
    result=[]
    def add(name,position=None,velocity=None,rows=(),mode=0,branch=.5,dt=.025,expected=0,status=None,newton=32,invalid=None,sample=None):
        q=list(reference if position is None else position)
        f0=mv(gradient(),list(map(lambda x:D(str(x)),q)))
        result.append(dict(name=name,mode=mode,dt=dt,volume=1/6,parameters=[5,1,branch,.3],velocity_scale=1,tolerance=1e-12,
            max_newton=newton,max_line_search=24,node_mass=[.01]*4,position=q,velocity=list(velocity or [0]*12),
            green=list(map(float,green(f0))),viscous=[0]*9,old_dissipation=0,rows=copy.deepcopy(list(rows)),
            expected_code=expected,expected_status=status,invalid_reason=invalid,sample_velocity=sample))
    add('free_identity')
    add('two_wall_compression',rows=walls(1,.01))
    add('mixed_prestress',position=mixed,rows=walls(.93,.005))
    shear_v=[0]*12;shear_v[6]=.2
    add('free_shear_velocity',velocity=shear_v,dt=.05)
    add('rigid_translation',velocity=[.1,-.2,.3]*4)
    rotated=[]
    for node in range(4): x,y,z=mixed[3*node:3*node+3];rotated.extend([-y,x,z])
    rotrows=walls(.93,.005)
    for row in rotrows:
        c=row['C'];new=[]
        for node in range(4): x,y,z=c[3*node:3*node+3];new.extend([-y,x,z])
        row['C']=new
    add('rotated_mixed',position=rotated,rows=rotrows)
    add('zero_branch',position=mixed,rows=walls(.93,.005),branch=0)
    sample=[(a-b)/.15 for a,b in zip(mixed,reference)]
    add('operator_sample',mode=2,dt=.15,sample=sample)
    result.append(dict(name='quadratic_limit',mode=1,dt=.1,velocity_scale=1,tolerance=1e-12,max_newton=32,max_line_search=24,
        mass=[.01,.02],position=[0,0],velocity=[0,0],stiffness=[[2,-1],[-1,2]],rows=[dict(C=[1,0],old=0,new=.01)],expected_code=0))
    opposite=[-v for v in lower]
    add('infeasible',rows=[dict(C=lower,old=0,new=.01),dict(C=opposite,old=0,new=.01)],expected=2,status='no_KKT_feasible_active_set')
    add('Newton_limit',position=mixed,rows=walls(.93,.005),expected=2,status='Newton_iteration_limit',newton=1)
    add('invalid_mass',expected=1,invalid='nonpositive_mass');result[-1]['node_mass'][0]=0
    add('zero_row',rows=[dict(C=[0]*12,old=0,new=0)],expected=1,invalid='zero_constraint_row')
    singular=[0,0,0,-1,0,0,0,-1,0,0,0,1]
    add('singular_trial_path',mode=2,dt=.1,sample=[(a-b)/.1 for a,b in zip(singular,reference)],expected=1,invalid='invalid_trial_path')
    return result


def encode(case):
    if case['mode']==1:
        values=[1,len(case['mass']),case['dt'],case['velocity_scale'],case['tolerance'],case['max_newton'],case['max_line_search']]
        for key in ['mass','position','velocity']:values.extend(case[key])
        for row in case['stiffness']:values.extend(row)
    else:
        values=[case['mode'],case['dt'],case['volume'],*case['parameters'],case['velocity_scale'],case['tolerance'],case['max_newton'],case['max_line_search']]
        for key in ['node_mass','position','velocity','green','viscous']:values.extend(case[key])
        values.append(case['old_dissipation'])
    values.append(len(case['rows']))
    for row in case['rows']:values.extend([row['old'],row['new'],*row['C']])
    if case['mode']==2:values.extend(case['sample_velocity'])
    return ' '.join(map(str,values))+'\n'


def prepare(root,executable,prospective):
    assert not root.exists(),'new study directory required'
    repo=Path(__file__).resolve().parents[1]
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=repo).strip()
    plan=read(prospective); resources=reserves(root.parent); controls=cases()
    root.mkdir();(root/'frozen').mkdir();artifacts=[]
    files=['src/shared/peridynamics/pd_nonsymmetric_constraints.h','src/shared/peridynamics/pd_tetra_gel_constraints.h',
        'src/shared/peridynamics/pd_gel_affine_path.h','src/shared/peridynamics/pd_gel_material.h','src/shared/peridynamics/pd_gel_material.cpp',
        'src/shared/peridynamics/pd_dense_metric_contact.h','src/shared/peridynamics/pd_normal_projection.h',
        'tests/3d_examples/extra_3d_PD_framework_tests/tetra_constraints_probe.cpp',
        'tests/3d_examples/extra_3d_PD_framework_tests/tetra_constraints_tests.cpp',
        'tests/3d_examples/extra_3d_PD_framework_tests/CMakeLists.txt',
        'tools/peri_tetra_constraints.py','tools/peri_tetra_constraints_audit.py','tools/test_peri_tetra_constraints.py',
        'tools/peri_gel_affine_audit.py','tools/peri_dense_contact_audit.py','tools/peri_normal_projection_audit.py',
        'tools/peri_elastic_normal.py','docs/peri/TETRA_GLOBAL_CONSTRAINTS_20261104.md']
    for source in files:
        original=repo/source;frozen=root/'frozen'/original.name;shutil.copy2(original,frozen);artifacts.extend([record(original),record(frozen)])
    frozen_exe=root/'frozen'/executable.name;shutil.copy2(executable,frozen_exe)
    artifacts.extend([record(executable),record(frozen_exe),record(prospective),record(plan['quadratic_reference_raw']),record(plan['quadratic_reference_review'])])
    for line in subprocess.check_output(['ldd',str(executable)],text=True).splitlines():
        for token in line.split():
            if token.startswith('/') and Path(token).is_file():artifacts.append(record(token))
    plan.update(cases=controls,resources=resources,prepared_at=stamp(),identity=identity(),
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip())
    save(root/'plan.json',plan);artifacts.append(record(root/'plan.json'));inputs=[]
    for case in controls:
        path=root/(case['name']+'.input');path.write_text(encode(case));inputs.append(record(path))
    save(root/'preparation.json',dict(artifacts=artifacts,inputs=inputs,executable=str(frozen_exe)))
    print(json.dumps(dict(status='prepared_not_passed',cases=len(controls),source_commit=plan['source_commit'])))


def limits():
    resource.setrlimit(resource.RLIMIT_CPU,(5,5));resource.setrlimit(resource.RLIMIT_AS,(256*2**20,256*2**20))
    resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_FSIZE,(32*2**20,32*2**20))


def run(root):
    assert not (root/'state.json').exists(),'attempted state cannot be redispatched'
    prep=read(root/'preparation.json');plan=read(root/'plan.json');attempts=[];failure=None;failure_traceback=None
    save(root/'state.json',dict(status='attempted',started_at=stamp(),identity=identity(),preparation=record(root/'preparation.json')))
    try:
        verify(prep['artifacts']+prep['inputs'])
        for case in plan['cases']:
            reserves(root);out=root/(case['name']+'.json');err=root/(case['name']+'.stderr')
            command=[prep['executable'],str(root/(case['name']+'.input'))];started=time.monotonic();timeout=False
            with out.open('xb') as stdout,err.open('xb') as stderr:
                child=subprocess.Popen(command,stdout=stdout,stderr=stderr,preexec_fn=limits)
                try:code=child.wait(timeout=5)
                except subprocess.TimeoutExpired:child.kill();code=child.wait();timeout=True
            attempt=dict(name=case['name'],command=command,return_code=code,timed_out=timeout,wall_seconds=time.monotonic()-started,
                stdout=record(out),stderr=record(err),cumulative_child_maxRSS_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
            attempts.append(attempt)
            assert not timeout and attempt['cumulative_child_maxRSS_KiB']<=256*1024
            raw=read(out);meta=raw['identity']
            assert meta['pid']==child.pid and meta['ppid']==identity()['pid'] and meta['argv']==command and meta['boot_id']==identity()['boot_id']
            attempt['runtime_gate']=inspect(case,raw,code)
        verify(prep['artifacts']+prep['inputs'])
    except Exception as error:failure=repr(error);failure_traceback=traceback.format_exc()
    save(root/'terminal.json',dict(status='runtime_completed_pending_independent_review' if failure is None else 'stopped_for_review',
        finished_at=stamp(),attempts=attempts,failure=failure,failure_traceback=failure_traceback))
    print(json.dumps(dict(status=read(root/'terminal.json')['status'],attempts=len(attempts),failure=failure)))
    return int(failure is not None)


def review(root):
    assert not (root/'independent_review.json').exists(),'do not repeat review receipt'
    prep=read(root/'preparation.json');plan=read(root/'plan.json');terminal=read(root/'terminal.json')
    result=dict(status='global_tetra_independent_review_failed',recorded_at=stamp(),identity=identity());assessed=[]
    try:
        verify([read(root/'state.json')['preparation']]);verify(prep['artifacts']+prep['inputs'])
        assert terminal['status']=='runtime_completed_pending_independent_review' and len(terminal['attempts'])==len(plan['cases'])
        for case,attempt in zip(plan['cases'],terminal['attempts']):
            verify([attempt['stdout'],attempt['stderr']]);assessed.append(inspect(case,read(attempt['stdout']['path']),attempt['return_code']))
        base=read(root/'mixed_prestress.json');rotated=read(root/'rotated_mixed.json');R=[[0,-1,0],[1,0,0],[0,0,1]]
        for field in ['position','velocity','internal_force','contact_impulse']:
            for a in range(4):
                expected=[sum(R[i][j]*base[field][3*a+j] for j in range(3)) for i in range(3)]
                assert max(abs(expected[i]-rotated[field][3*a+i]) for i in range(3))<=1e-10
        for field in ['wall_work','kinetic_new','algorithmic_loss','material_loss']:
            assert abs(base[field]-rotated[field])<=1e-10
        for i in range(12):
            for j in range(12):
                value=sum(R[i%3][r]*base['jacobian'][3*(i//3)+r][3*(j//3)+c]*R[j%3][c] for r in range(3) for c in range(3))
                assert abs(value-rotated['jacobian'][i][j])<=1e-10
        old=read(plan['quadratic_reference_raw']);new=read(root/'quadratic_limit.json')
        for field in ['velocity','normal_impulse']:
            assert max(abs(a-b) for a,b in zip(new[field],old[field]))<=1e-10
        assert read(root/'zero_branch.json')['material_loss']==0
        translation=read(root/'rigid_translation.json')
        assert max(map(abs,translation['internal_force']))<=1e-10
        operator=read(root/'operator_sample.json');assert operator['jacobian_asymmetry']>1e-8
        controls=[]
        for label in ['wrong_transpose_force','symmetrized_global_jacobian','effective_jacobian_as_kinetic_mass']:
            if label=='effective_jacobian_as_kinetic_mass':
                case=next(c for c in plan['cases'] if c['name']=='two_wall_compression');bad=copy.deepcopy(read(root/'two_wall_compression.json'))
                v=bad['velocity'];j=bad['jacobian'];bad['kinetic_new']=sum(v[a]*j[a][b]*v[b]/2 for a in range(12) for b in range(12))
            else:
                case=next(c for c in plan['cases'] if c['name']=='operator_sample');bad=copy.deepcopy(operator)
                if label=='wrong_transpose_force':bad['internal_force'][0]*=1.01
                else:
                    j=bad['jacobian'];bad['jacobian']=[[(j[a][b]+j[b][a])/2 for b in range(12)] for a in range(12)]
            path=root/(label+'.negative.json');save(path,bad)
            try:inspect(case,bad,0)
            except AssertionError:controls.append(dict(name=label,rejected=True,artifact=record(path)))
            else:raise AssertionError(('negative control accepted',label))
        result.update(status='bounded_nonsymmetric_global_tetra_constraints_independently_passed',cases=assessed,negative_controls=controls,
            global_frame_rotation_checked=True,quadratic_limit_compared_to_existing_raw_no_old_solver_rerun=True,
            scope='one-step4node3Dvariational/control algebra; not a PD sphere or trajectory convergence qualification',
            angular_conservation_qualified=False,trajectory_accuracy_qualified=False,production_integrated=False,
            production_restart_qualified=False,S3_complete=False,calibration_allowed=False)
    except Exception as error:result.update(failure=repr(error),completed_assessments=assessed)
    result['terminal']=record(root/'terminal.json');result['state']=record(root/'state.json');result['preparation']=record(root/'preparation.json')
    save(root/'independent_review.json',result)
    print(json.dumps(dict(status=result['status'],assessed_cases=len(assessed),failure=result.get('failure'))))
    return int(not result['status'].endswith('_passed'))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run','review']);parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--executable',type=Path);parser.add_argument('--prospective',type=Path);args=parser.parse_args()
    if args.action=='prepare':prepare(args.root,args.executable,args.prospective)
    else:sys.exit(globals()[args.action](args.root))
