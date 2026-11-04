"""Versioned exact contact balance/new guarded solver qualification controls."""
import argparse
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import random
import resource
import shutil
import subprocess
import sys
import time
import traceback
from fractions import Fraction as Q
from peri_elastic_normal import identity, read, record, reserves, save, stamp, verify
from peri_many_contact import manufactured, controls as old_controls, encode as encode_many, exact_decimal
from peri_guarded_contact_audit import bits, inspect, oracle

def controls():
    result=[]
    def arithmetic(name,r,c,lam,scale=1.,rounding=0,invalid=None):
        if invalid:expected=1
        else:expected=0 if oracle(r,c,lam,scale)[2]<=Q(float(1e-12)) else 2
        result.append(dict(name=name,mode=2,momentum_bits=bits(r),coefficient_bits=list(map(bits,c)),impulse_bits=list(map(bits,lam)),
            scale_bits=bits(scale),rounding=rounding,expected_code=expected,invalid=invalid))
    arithmetic('small_momentum_between_huge_impulses',.0625,[1,1,1],[2.**70,.125,-2.**70])
    arithmetic('product_low_bits',1.,[1+2.**-27],[1-2.**-27])
    arithmetic('finite_result_despite_individual_product_overflow',.25,[2.**800,-2.**800],[2.**800,2.**800])
    arithmetic('subnormal_product',0.,[2.**-1074],[2.**-1074])
    arithmetic('normal_boundary_minus_subnormal_product',2.**-1022,[2.**-1074],[2.**-1074])
    arithmetic('tie_even_positive',1.,[-2.**-53],[1.])
    arithmetic('tie_odd_positive',1+2.**-52,[-2.**-53],[1.])
    arithmetic('tie_even_negative',-1.,[2.**-53],[1.])
    arithmetic('rounded_normalized_third',1.,[],[],3.)
    arithmetic('same_third_upward_environment',1.,[],[],3.,1)
    arithmetic('absolute_rounding_budget_rejection',1e9,[],[],3.)
    randomizer=random.Random(20261104);pairs=[]
    for i in range(128):
        coefficient=math.ldexp(1.+randomizer.randrange(16)/16,randomizer.randrange(-900,901))
        impulse=math.ldexp(1.+randomizer.randrange(16)/16,randomizer.randrange(-100,101))
        pairs.extend([(coefficient,impulse),(coefficient,-impulse)])
    randomizer.shuffle(pairs)
    arithmetic('all256_permuted_extreme_product_pairs',.1875,[x for x,y in pairs],[y for x,y in pairs],3.)
    arithmetic('invalid_scale',1.,[],[],0.,invalid='invalid_scale')
    arithmetic('nonfinite_supplied_momentum',math.inf,[],[],1.,invalid='nonfinite')
    arithmetic('nonfinite_contact',1.,[math.nan],[1.],1.,invalid='nonfinite')
    arithmetic('normalized_result_overflow',float.fromhex('0x1.fffffffffffffp1023'),[],[],2.**-1074,invalid='overflow')
    axis=lambda n,i:[Q(int(i==j)) for j in range(n)]
    rows=[(axis(6,i%6),Q(0),Q(1,16)) for i in range(24)]
    rows += [(axis(6,i%6),Q(1,8),Q(0)) for i in range(6)]
    rows += [(axis(6,i%6),Q(0),Q(0)) for i in range(6)]
    coupled=manufactured('new_six_coordinate_many_row_root',n=6,rows=rows,beta=Q(1,8),root=[Q(3 if i%2==0 else -3,2**(3+i%3)) for i in range(6)])
    result.append(coupled)
    permutation=copy.deepcopy(coupled);permutation['name']='new_root_reversed_rows'
    for key in ['rows','dual_witness','gap_reference']:permutation[key].reverse()
    result.append(permutation)
    free=manufactured('new_twelve_coordinate_cubic_free',n=12,beta=Q(3,8))
    result.append(free)
    infeasible=manufactured('new_two_coordinate_infeasible_geometry',n=2,rows=[([Q(1),Q(0)],Q(0),Q(0)),([Q(-1),Q(0)],Q(0),Q(0))],beta=Q(1,8))
    for row in infeasible['rows']:row['new']='0.0625'
    infeasible.update(expected_code=2,rejection='infeasible');result.append(infeasible)
    initial=manufactured('initial_precision_rejection',n=1,beta=Q(0))
    initial.update(A=[['1']],load=['-1000000000'],mass_scale='3',expected_code=2,rejection='precision_initial');result.append(initial)
    trial=manufactured('trial_precision_rejection',n=1,beta=Q(1),root=[Q(-1,2)])
    trial.update(A=[[exact_decimal(Q(2)**-40)]],load=[exact_decimal(-Q(1,8)-Q(2)**-41)],mass_scale='3',domain_radius='1e15',expected_code=2,rejection='precision_trial');result.append(trial)
    # The original state stopped before these four inputs were attempted.
    for case in old_controls()[14:]:
        case=copy.deepcopy(case);case['name']='previously_unattempted_'+case['name'];result.append(case)
    return result

def encode(case):
    if case['mode']!=2:return encode_many(case)
    values=[2,len(case['coefficient_bits']),case['rounding'],case['momentum_bits'],case['scale_bits'],*case['coefficient_bits'],*case['impulse_bits']]
    return ' '.join(map(str,values))+'\n'

def prepare(root,executable,prospective):
    assert not root.exists(),'new study required'
    repo=Path(__file__).resolve().parents[1];assert not subprocess.check_output(['git','status','--porcelain'],cwd=repo).strip()
    plan=read(prospective);budget=reserves(root.parent)
    assert plan['native_seconds_max']==5 and plan['accumulation_roundoff_budget']==1e-12
    root.mkdir();(root/'frozen').mkdir();artifacts=[]
    files=['src/shared/peridynamics/pd_exact_contact_balance.h','src/shared/peridynamics/pd_many_contact.h',
        'src/shared/peridynamics/pd_nonsymmetric_constraints.h','src/shared/peridynamics/pd_dense_metric_contact.h','src/shared/peridynamics/pd_normal_projection.h',
        'tests/3d_examples/extra_3d_PD_framework_tests/many_contact_probe.cpp','tests/3d_examples/extra_3d_PD_framework_tests/many_contact_tests.cpp',
        'tests/3d_examples/extra_3d_PD_framework_tests/exact_contact_balance_tests.cpp','tests/3d_examples/extra_3d_PD_framework_tests/CMakeLists.txt',
        'tools/peri_guarded_contact.py','tools/peri_guarded_contact_audit.py','tools/test_peri_guarded_contact.py',
        'tools/peri_many_contact.py','tools/peri_many_contact_audit.py','tools/peri_elastic_normal.py','docs/peri/EXACT_CONTACT_BALANCE_20261104.md']
    for file in files:
        source=repo/file;target=root/'frozen'/source.name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
    frozen=root/'frozen'/executable.name;shutil.copy2(executable,frozen);artifacts.extend([record(executable),record(frozen),record(prospective)])
    # Freeze every compiler-reported header, including the actual Boost/Eigen
    # implementation, instead of assuming an installed version is immutable.
    build=repo/'build/peri-latest-20261103';obj='extra_3d_PD_framework_tests/CMakeFiles/extra_3d_PD_guarded_contact_probe.dir/many_contact_probe.cpp.o'
    command_text=subprocess.check_output(['ninja','-C',str(build),'-t','commands',obj],text=True)
    assert not any(flag in command_text for flag in ['-ffast-math','-Ofast','-ffinite-math-only'])
    (root/'frozen/compiler_command.txt').write_text(command_text);artifacts.append(record(root/'frozen/compiler_command.txt'))
    dependency_text=subprocess.check_output(['ninja','-C',str(build),'-t','deps',obj],text=True)
    dependencies=[]
    for line in dependency_text.splitlines():
        if line.startswith('    '):
            source=Path(line.strip());source=source if source.is_absolute() else build/source
            if source.is_file():dependencies.append(source.resolve())
    dependencies=sorted(set(dependencies));assert dependencies
    assert sum(p.stat().st_size for p in dependencies)<512*2**20,'preparation source-copy budget'
    manifest=[]
    for index,source in enumerate(dependencies):
        target=root/'frozen/compiler_headers'/str(index)/source.name;target.parent.mkdir(parents=True);shutil.copy2(source,target)
        entry=dict(original=record(source),frozen=record(target));manifest.append(entry);artifacts.extend(entry.values())
    save(root/'frozen/compiler_dependencies.json',manifest);artifacts.append(record(root/'frozen/compiler_dependencies.json'))
    for line in subprocess.check_output(['ldd',str(executable)],text=True).splitlines():
        for token in line.split():
            if token.startswith('/') and Path(token).is_file():artifacts.append(record(token))
    proof=read(plan['preserved_v12_adjudication'])
    assert proof['old_study_marked_passed'] is False and proof['native_solver_called'] is False
    for key in ['preserved_state','preserved_terminal','preserved_preparation','preserved_native_ELF','frozen_assessor']:
        verify([proof[key]]);artifacts.append(proof[key])
    artifacts.append(record(plan['preserved_v12_adjudication']))
    old=read(plan['preserved_v12_terminal']);assert old['status']=='stopped_for_review'
    for attempt in old['attempts']:
        verify([attempt['stdout'],attempt['stderr']]);artifacts.extend([attempt['stdout'],attempt['stderr']])
    attempted_hashes={record(a['command'][1])['sha256'] for a in old['attempts']}
    plan.update(controls=controls(),resources=budget,prepared_at=stamp(),identity=identity(),
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip())
    save(root/'plan.json',plan);artifacts.extend([record(root/'plan.json'),record(plan['preserved_v12_terminal'])]);inputs=[]
    for case in plan['controls']:
        path=root/(case['name']+'.input');path.write_text(encode(case));item=record(path)
        assert item['sha256'] not in attempted_hashes,'old attempted solver input must not be rerun'
        inputs.append(item)
    save(root/'preparation.json',dict(artifacts=artifacts,inputs=inputs,executable=str(frozen),compiler_header_count=len(dependencies)))
    print(json.dumps(dict(status='prepared_not_passed',controls=len(plan['controls']),compiler_headers=len(dependencies),source_commit=plan['source_commit'])))

def limits():
    resource.setrlimit(resource.RLIMIT_CPU,(5,5));resource.setrlimit(resource.RLIMIT_AS,(256*2**20,256*2**20))
    resource.setrlimit(resource.RLIMIT_FSIZE,(32*2**20,32*2**20));resource.setrlimit(resource.RLIMIT_CORE,(0,0))

def run(root):
    assert not (root/'state.json').exists(),'old state cannot be redispatched'
    prep=read(root/'preparation.json');plan=read(root/'plan.json');attempts=[];failure=None;failure_traceback=None
    save(root/'state.json',dict(status='attempted',started_at=stamp(),identity=identity(),preparation=record(root/'preparation.json')))
    try:
        verify(prep['artifacts']+prep['inputs'])
        for case in plan['controls']:
            reserves(root);out=root/(case['name']+'.json');err=root/(case['name']+'.stderr');start=time.monotonic();timed_out=False
            command=[prep['executable'],str(root/(case['name']+'.input'))]
            with out.open('xb') as stdout,err.open('xb') as stderr:
                child=subprocess.Popen(command,stdout=stdout,stderr=stderr,preexec_fn=limits,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
                try:code=child.wait(timeout=5)
                except subprocess.TimeoutExpired:child.kill();code=child.wait();timed_out=True
            a=dict(name=case['name'],command=command,return_code=code,timed_out=timed_out,wall_seconds=time.monotonic()-start,
                stdout=record(out),stderr=record(err),cumulative_child_maxRSS_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
            attempts.append(a);assert not timed_out and a['cumulative_child_maxRSS_KiB']<=256*1024
            raw=read(out);meta=raw['identity'];assert meta['argv']==command and meta['pid']==child.pid and meta['ppid']==identity()['pid']
            assert meta['boot_id']==identity()['boot_id'] and meta['start_ticks'].isdigit()
            a['runtime_gate']=inspect(case,raw,code)
        verify(prep['artifacts']+prep['inputs'])
    except Exception as e:failure=repr(e);failure_traceback=traceback.format_exc()
    save(root/'terminal.json',dict(status='runtime_completed_pending_independent_review' if failure is None else 'stopped_for_review',
        finished_at=stamp(),attempts=attempts,failure=failure,failure_traceback=failure_traceback))
    print(json.dumps(dict(status=read(root/'terminal.json')['status'],attempts=len(attempts),failure=failure)))
    return int(failure is not None)

def review(root):
    assert not (root/'independent_review.json').exists(),'immutable review exists'
    prep=read(root/'preparation.json');plan=read(root/'plan.json');terminal=read(root/'terminal.json');assessed=[]
    result=dict(status='guarded_contact_review_failed',recorded_at=stamp(),identity=identity())
    try:
        verify([read(root/'state.json')['preparation']]);verify(prep['artifacts']+prep['inputs'])
        assert terminal['status']=='runtime_completed_pending_independent_review' and len(terminal['attempts'])==len(plan['controls'])
        for case,a in zip(plan['controls'],terminal['attempts']):
            verify([a['stdout'],a['stderr']]);assessed.append(inspect(case,read(a['stdout']['path']),a['return_code']))
        base=next(x for x in assessed if x['name']=='new_six_coordinate_many_row_root')
        perm=next(x for x in assessed if x['name']=='new_root_reversed_rows')
        for field in ['velocity','scatter']:assert max(abs(x-y) for x,y in zip(base[field],perm[field]))<=1e-10
        negatives=[]
        for label,name in [('hidden_roundoff_enclosure','rounded_normalized_third'),('false_precision_acceptance','absolute_rounding_budget_rejection'),('lost_small_momentum','small_momentum_between_huge_impulses'),('symmetrized_true_operator','new_six_coordinate_many_row_root')]:
            case=next(x for x in plan['controls'] if x['name']==name);bad=copy.deepcopy(read(root/(name+'.json')))
            if label=='hidden_roundoff_enclosure':bad.update(roundoff_bound=0,bound_bits=bits(0))
            elif label=='false_precision_acceptance':bad['within_budget']=True
            elif label=='lost_small_momentum':bad.update(normalized_value=0,value_bits=bits(0))
            else:
                j=bad['jacobian'];bad['jacobian']=[[(j[i][k]+j[k][i])/2 for k in range(len(j))] for i in range(len(j))]
            path=root/(label+'.negative.json');save(path,bad)
            try:inspect(case,bad,case['expected_code'])
            except AssertionError:negatives.append(dict(name=label,rejected=True,artifact=record(path)))
            else:raise AssertionError(('negative control accepted',label))
        result.update(status='exact_contact_balance_and_new_guarded_many_contact_independently_passed',cases=assessed,negative_controls=negatives,
            no_old_attempted_solver_input_rerun=True,preserved_v12_stopped=True,accumulation_scope='exact supplied floating products/sum and correctly rounded normalization only',
            production_integrated=False,S3_complete=False,production_restart_qualified=False,calibration_allowed=False)
    except Exception as e:result.update(failure=repr(e),failure_traceback=traceback.format_exc(),completed_assessments=assessed)
    result.update(state=record(root/'state.json'),terminal=record(root/'terminal.json'),preparation=record(root/'preparation.json'))
    save(root/'independent_review.json',result)
    print(json.dumps(dict(status=result['status'],cases=len(assessed),failure=result.get('failure'))))
    return int(not result['status'].endswith('_passed'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run','review']);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--executable',type=Path);p.add_argument('--prospective',type=Path);args=p.parse_args()
    if args.action=='prepare':prepare(args.root,args.executable,args.prospective)
    else:sys.exit(globals()[args.action](args.root))
