"""Freeze/run/independent-review manufactured many-contact algebra controls."""
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
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from peri_elastic_normal import identity, read, record, reserves, save, stamp, verify
from peri_many_contact_audit import inspect

def exact_decimal(value):
    value=Q(value)
    with localcontext() as ctx:
        ctx.prec=95
        x=Decimal(value.numerator)/Decimal(value.denominator)
        assert Q(x)==value,'control inputs must have exact finite decimals'
        return str(x)

def manufactured(name,n=3,rows=(),beta=Q(1,4),root=None):
    a=[[Q(0)]*n for _ in range(n)]
    for i in range(n):
        a[i][i]=Q(2)
        if i+1<n:a[i][i+1]=Q(3,16);a[i+1][i]=Q(-1,16)
    x=root or [Q(1 if i%2==0 else -1,2**(2+i%3)) for i in range(n)]
    c=[row[0] for row in rows];gaps=[row[1] for row in rows];lam=[row[2] for row in rows]
    load=[sum((a[i][j]*x[j] for j in range(n)),Q(0))+beta*x[i]**3-sum((c[s][i]*lam[s] for s in range(len(rows))),Q(0)) for i in range(n)]
    case=dict(name=name,mode=0,A=[[exact_decimal(v) for v in row] for row in a],beta=[exact_decimal(beta)]*n,
        load=list(map(exact_decimal,load)),position=['0']*n,dt='1',mass_scale='1',speed='1',tolerance='1e-12',
        max_newton=48,max_line_search=32,domain_radius='10',root_reference=list(map(exact_decimal,x)),
        dual_witness=list(map(exact_decimal,lam)),gap_reference=list(map(exact_decimal,gaps)),expected_code=0,
        rows=[dict(C=list(map(exact_decimal,row)),old='0',new=exact_decimal(sum((row[i]*x[i] for i in range(n)),Q(0))-g)) for row,g,_ in rows])
    return case

def controls():
    axis=lambda n,i:[Q(int(j==i)) for j in range(n)]
    rows=[(axis(3,0),Q(0),Q(1,8)) for _ in range(8)]+[(axis(3,1),Q(1,4),Q(0)),(axis(3,2),Q(0),Q(0))]
    mixed=manufactured('cubic_mixed',rows=rows)
    cases=[manufactured('free_linear',beta=Q(0)),mixed]
    perm=copy.deepcopy(mixed);perm['name']='rows_permuted'
    for key in ['rows','dual_witness','gap_reference']:perm[key].reverse()
    cases.append(perm)
    rot=copy.deepcopy(mixed);rot['name']='frame_axis_rotation'
    order=[1,0,2];sign=[-1,1,1]
    def vector(v):return [exact_decimal(sign[i]*Q(v[order[i]])) for i in range(3)]
    rot['A']=[[exact_decimal(sign[i]*sign[j]*Q(mixed['A'][order[i]][order[j]])) for j in range(3)] for i in range(3)]
    for key in ['position','load','root_reference']:rot[key]=vector(mixed[key])
    for row in rot['rows']:row['C']=vector(row['C'])
    cases.append(rot)
    scaled=copy.deepcopy(mixed);scaled['name']='positive_row_rescaling'
    for i,row in enumerate(scaled['rows']):
        factor=Q(2) if i%2 else Q(1,2)
        for key in ['old','new']:row[key]=exact_decimal(Q(row[key])*factor)
        row['C']=[exact_decimal(Q(x)*factor) for x in row['C']]
        scaled['dual_witness'][i]=exact_decimal(Q(scaled['dual_witness'][i])/factor)
        scaled['gap_reference'][i]=exact_decimal(Q(scaled['gap_reference'][i])*factor)
    cases.append(scaled)
    normalizers=copy.deepcopy(mixed);normalizers['name']='same_physical_equation_other_normalizers';normalizers['mass_scale']='2';normalizers['speed']='0.5'
    cases.append(normalizers)
    many=[]
    for i in range(256):many.append((axis(3,i%3),Q(0) if i<128 or i>=192 else Q(1,4),Q(1,128) if i<128 else Q(0)))
    cases.append(manufactured('max_rows_256',rows=many))
    cases.append(manufactured('max_dofs_128',n=128))
    cases.append(manufactured('biactive_origin',rows=[(axis(3,i),Q(0),Q(0)) for i in range(3)],root=[Q(0)]*3))
    cases.append(dict(name='rank_deficient_consistent_linear',mode=1,A=[['1','1'],['2','2']],rhs=['3','6'],minimum_norm_reference=['1.5','1.5'],expected_code=0))
    cases.append(dict(name='inconsistent_original_linear',mode=1,A=[['1','1'],['2','2']],rhs=['3','7'],expected_code=1))
    limit=copy.deepcopy(mixed);limit.update(name='Newton_limit',max_newton=1,expected_code=2,rejection='Newton_limit');cases.append(limit)
    domain=manufactured('domain_blocked',n=1,beta=Q(0),root=[Q(1)])
    domain.update(domain_radius='1e-20',expected_code=2,rejection='domain_blocked');cases.append(domain)
    infeasible=manufactured('infeasible_opposite_rows',n=1,rows=[([Q(1)],Q(0),Q(0)),([Q(-1)],Q(0),Q(0))])
    for row in infeasible['rows']:row['new']='0.125'
    infeasible.update(expected_code=2,rejection='infeasible');cases.append(infeasible)
    for name,label in [('zero_row','zero_row'),('initial_penetration','initial_penetration'),('negative_cubic','negative_cubic'),('zero_mass_scale','zero_mass_scale')]:
        case=copy.deepcopy(mixed);case.update(name=name,expected_code=1,invalid=label)
        if label=='zero_row':case['rows'][0]['C']=['0']*3
        elif label=='initial_penetration':case['rows'][0]['old']='0.125'
        elif label=='negative_cubic':case['beta'][0]='-0.25'
        else:case['mass_scale']='0'
        cases.append(case)
    return cases

def encode(case):
    values=[case['mode'],len(case['A'])]
    if case['mode']==1:
        for row in case['A']:values.extend(row)
        values.extend(case['rhs'])
    else:
        values.extend([case['dt'],case['mass_scale'],case['speed'],case['tolerance'],case['max_newton'],case['max_line_search'],case['domain_radius']])
        for key in ['position','beta','load']:values.extend(case[key])
        for row in case['A']:values.extend(row)
        values.append(len(case['rows']))
        for row in case['rows']:values.extend([row['old'],row['new'],*row['C']])
    return ' '.join(map(str,values))+'\n'

def prepare(root,executable,prospective):
    assert not root.exists(),'new unused state directory required'
    repo=Path(__file__).resolve().parents[1]
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=repo).strip(),'commit first'
    plan=read(prospective);budget=reserves(root.parent)
    assert plan['nonlinear_normalized_root_gate']==1e-12 and plan['independent_residual_Jac_KKT_gate']==1e-10
    root.mkdir();(root/'frozen').mkdir();artifacts=[]
    files=['src/shared/peridynamics/pd_many_contact.h','src/shared/peridynamics/pd_nonsymmetric_constraints.h',
        'src/shared/peridynamics/pd_dense_metric_contact.h','src/shared/peridynamics/pd_normal_projection.h',
        'tests/3d_examples/extra_3d_PD_framework_tests/many_contact_probe.cpp',
        'tests/3d_examples/extra_3d_PD_framework_tests/many_contact_tests.cpp',
        'tests/3d_examples/extra_3d_PD_framework_tests/CMakeLists.txt','tools/peri_many_contact.py',
        'tools/peri_many_contact_audit.py','tools/test_peri_many_contact.py','tools/peri_elastic_normal.py',
        'docs/peri/MANY_CONTACT_CONTROL_20261104.md']
    for file in files:
        source=repo/file;target=root/'frozen'/source.name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
    frozen=root/'frozen'/executable.name;shutil.copy2(executable,frozen);artifacts.extend([record(executable),record(frozen),record(prospective)])
    for filename in ['/usr/include/eigen3/Eigen/src/SVD/JacobiSVD.h','/usr/include/eigen3/Eigen/src/SVD/SVDBase.h']:
        source=Path(filename);target=root/'frozen'/source.name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
    linked=subprocess.check_output(['ldd',str(executable)],text=True);(root/'frozen/ldd.txt').write_text(linked);artifacts.append(record(root/'frozen/ldd.txt'))
    for line in linked.splitlines():
        for token in line.split():
            if token.startswith('/') and Path(token).is_file():artifacts.append(record(token))
    plan.update(controls=controls(),resources=budget,prepared_at=stamp(),identity=identity(),
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip())
    save(root/'plan.json',plan);artifacts.append(record(root/'plan.json'));inputs=[]
    for case in plan['controls']:
        path=root/(case['name']+'.input');path.write_text(encode(case));inputs.append(record(path))
    save(root/'preparation.json',dict(artifacts=artifacts,inputs=inputs,executable=str(frozen)))
    print(json.dumps(dict(status='prepared_not_passed',cases=len(plan['controls']),source_commit=plan['source_commit'])))

def limits():
    resource.setrlimit(resource.RLIMIT_CPU,(5,5));resource.setrlimit(resource.RLIMIT_AS,(256*2**20,256*2**20))
    resource.setrlimit(resource.RLIMIT_FSIZE,(32*2**20,32*2**20));resource.setrlimit(resource.RLIMIT_CORE,(0,0))

def run(root):
    assert not (root/'state.json').exists(),'attempted state never redispatched'
    prep=read(root/'preparation.json');plan=read(root/'plan.json');attempts=[];failure=None;failure_traceback=None
    save(root/'state.json',dict(status='attempted',started_at=stamp(),identity=identity(),preparation=record(root/'preparation.json')))
    try:
        verify(prep['artifacts']+prep['inputs'])
        for case in plan['controls']:
            reserves(root);out=root/(case['name']+'.json');err=root/(case['name']+'.stderr');started=time.monotonic();timeout=False
            command=[prep['executable'],str(root/(case['name']+'.input'))]
            env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
            with out.open('xb') as stdout,err.open('xb') as stderr:
                child=subprocess.Popen(command,stdout=stdout,stderr=stderr,env=env,preexec_fn=limits)
                try:code=child.wait(timeout=5)
                except subprocess.TimeoutExpired:child.kill();code=child.wait();timeout=True
            attempt=dict(name=case['name'],command=command,return_code=code,timed_out=timeout,wall_seconds=time.monotonic()-started,
                stdout=record(out),stderr=record(err),cumulative_child_maxRSS_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
            attempts.append(attempt)
            assert not timeout and attempt['cumulative_child_maxRSS_KiB']<=256*1024
            raw=read(out);meta=raw['identity']
            assert meta['argv']==command and meta['pid']==child.pid and meta['ppid']==identity()['pid']
            assert meta['boot_id']==identity()['boot_id'] and meta['start_ticks'].isdigit()
            attempt['runtime_gate']=inspect(case,raw,code)
        verify(prep['artifacts']+prep['inputs'])
    except Exception as e:failure=repr(e);failure_traceback=traceback.format_exc()
    save(root/'terminal.json',dict(status='runtime_completed_pending_independent_review' if failure is None else 'stopped_for_review',
        finished_at=stamp(),attempts=attempts,failure=failure,failure_traceback=failure_traceback))
    print(json.dumps(dict(status=read(root/'terminal.json')['status'],attempts=len(attempts),failure=failure)))
    return int(failure is not None)

def review(root):
    assert not (root/'independent_review.json').exists(),'immutable review already exists'
    prep=read(root/'preparation.json');plan=read(root/'plan.json');terminal=read(root/'terminal.json');assessed=[]
    result=dict(status='many_contact_independent_review_failed',recorded_at=stamp(),identity=identity())
    try:
        verify([read(root/'state.json')['preparation']]);verify(prep['artifacts']+prep['inputs'])
        assert terminal['status']=='runtime_completed_pending_independent_review' and len(terminal['attempts'])==len(plan['controls'])
        for case,attempt in zip(plan['controls'],terminal['attempts']):
            verify([attempt['stdout'],attempt['stderr']]);assessed.append(inspect(case,read(attempt['stdout']['path']),attempt['return_code']))
        base=next(x for x in assessed if x['name']=='cubic_mixed')
        for name in ['rows_permuted','positive_row_rescaling','same_physical_equation_other_normalizers']:
            item=next(x for x in assessed if x['name']==name)
            for field in ['velocity','scatter']:assert max(abs(x-y) for x,y in zip(base[field],item[field]))<=1e-10
        rotated=next(x for x in assessed if x['name']=='frame_axis_rotation')
        for field in ['velocity','scatter']:
            expected=[-base[field][1],base[field][0],base[field][2]]
            assert max(abs(x-y) for x,y in zip(expected,rotated[field]))<=1e-10
        negative=[];case=next(x for x in plan['controls'] if x['name']=='cubic_mixed')
        for label in ['symmetrized_operator_Jacobian','corrupted_nonzero_transpose_impulse','corrupted_generalized_Newton_matrix']:
            bad=copy.deepcopy(read(root/'cubic_mixed.json'))
            if label=='symmetrized_operator_Jacobian':
                j=bad['jacobian'];bad['jacobian']=[[(j[i][k]+j[k][i])/2 for k in range(len(j))] for i in range(len(j))]
            elif label=='corrupted_nonzero_transpose_impulse':bad['normal_impulse'][0]*=1.01
            else:bad['trace'][0]['generalized_jacobian'][0][len(case['position'])]*=1.01
            path=root/(label+'.negative.json');save(path,bad)
            try:inspect(case,bad,0)
            except AssertionError:negative.append(dict(name=label,rejected=True,artifact=record(path)))
            else:raise AssertionError(('negative control accepted',label))
        result.update(status='manufactured_many_row_nonsymmetric_contact_independently_passed',cases=assessed,negative_controls=negative,
            exact_monotonicity_certificates=True,unique_primal_and_scatter_checked=True,redundant_dual_uniqueness_claimed=False,
            tested_max_DOF=128,tested_max_rows=256,bounds_tested_together=False,
            scope='manufactured algebra, not actual PD/material trajectory, sphere, angular or energy qualification',
            production_integrated=False,production_restart_qualified=False,S3_complete=False,calibration_allowed=False)
    except Exception as e:result.update(failure=repr(e),failure_traceback=traceback.format_exc(),completed_assessments=assessed)
    result.update(terminal=record(root/'terminal.json'),state=record(root/'state.json'),preparation=record(root/'preparation.json'))
    save(root/'independent_review.json',result)
    print(json.dumps(dict(status=result['status'],cases=len(assessed),failure=result.get('failure'))))
    return int(not result['status'].endswith('_passed'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run','review']);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--executable',type=Path);p.add_argument('--prospective',type=Path);args=p.parse_args()
    if args.action=='prepare':prepare(args.root,args.executable,args.prospective)
    else:sys.exit(globals()[args.action](args.root))
