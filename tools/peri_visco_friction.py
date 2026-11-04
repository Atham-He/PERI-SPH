"""Freeze, run and separately review bounded S2 controls; never redispatch state."""
import argparse
import copy
import json
import math
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time

from peri_elastic_normal import identity, read, record, reserves, save, stamp, verify
from peri_visco_friction_audit import assess, check_checkpoint, rows


def prepare(root,executable,prospective):
    assert not root.exists(),'never reuse a study directory'
    repo=Path(__file__).resolve().parents[1]
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=repo).strip(),'commit source before freeze'
    plan=read(prospective); resources=reserves(root.parent)
    assert plan['stage']=='S2' and plan['status']=='prospective_frozen_before_dispatch'
    for case in plan['cases']:
        dt,total=map(float,case['config'][:2])
        assert dt>0 and total/dt<=10000
    root.mkdir(); (root/'frozen').mkdir(); artifacts=[]
    paths=['src/shared/peridynamics/pd_normal_projection.h','src/shared/peridynamics/pd_visco_friction_step.h',
        'tests/3d_examples/extra_3d_PD_framework_tests/visco_friction_probe.cpp',
        'tests/3d_examples/extra_3d_PD_framework_tests/visco_friction_tests.cpp',
        'tests/3d_examples/extra_3d_PD_framework_tests/CMakeLists.txt',
        'tools/peri_visco_friction.py','tools/peri_visco_friction_audit.py','tools/peri_elastic_normal.py',
        'tools/test_peri_visco_friction.py','docs/peri/VISCO_FRICTION_STAGE_S2_20261104.md']
    for relative in paths:
        original=repo/relative; frozen=root/'frozen'/original.name
        shutil.copy2(original,frozen); artifacts.extend([record(original),record(frozen)])
    frozen_exe=root/'frozen'/executable.name
    shutil.copy2(executable,frozen_exe); artifacts.extend([record(executable),record(frozen_exe),record(prospective)])
    for line in subprocess.check_output(['ldd',str(executable)],text=True).splitlines():
        for token in line.split():
            if token.startswith('/') and Path(token).is_file(): artifacts.append(record(token))
    for name in ['s1_raw_reference','s1_review_reference']:
        artifacts.append(record(plan[name]))
    plan.update(source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
                resources=resources,prepared_at=stamp(),preparer_identity=identity())
    save(root/'plan.json',plan); artifacts.append(record(root/'plan.json'))
    inputs=[]
    for case in plan['cases']:
        path=root/(case['name']+'.input')
        path.write_text(' '.join([str(case['kind'])]+list(map(str,case['config'])))+'\n')
        inputs.append(record(path))
    save(root/'preparation.json',dict(artifacts=artifacts,inputs=inputs,executable=str(frozen_exe)))
    print(json.dumps(dict(status='prepared_not_passed',cases=len(inputs),source_commit=plan['source_commit'])))


def limits():
    resource.setrlimit(resource.RLIMIT_CPU,(5,5))
    resource.setrlimit(resource.RLIMIT_AS,(256*2**20,256*2**20))
    resource.setrlimit(resource.RLIMIT_FSIZE,(32*2**20,32*2**20))
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))


def gate(case,assessment):
    if case.get('finest'):
        for group in ['force_relative_L2','velocity_relative_L2']:
            for value in assessment[group]:
                if value is not None: assert value<=.02,(case['name'],group,value)
        if assessment['algorithmic_fraction'] is not None:
            assert assessment['algorithmic_fraction']<=.01,(case['name'],'algorithmic_fraction',assessment['algorithmic_fraction'])


def run(root):
    assert not (root/'state.json').exists(),'attempted state cannot be redispatched'
    prep=read(root/'preparation.json'); plan=read(root/'plan.json')
    save(root/'state.json',dict(status='attempted',started_at=stamp(),identity=identity(),preparation=record(root/'preparation.json')))
    attempts=[]; checkpoints=[]; assessed=[]; failure=None
    def native(name,case,extra=(),trace_source=None,tail=False,expected_return=0):
        reserves(root)
        directory=root/name; directory.mkdir()
        trace=directory/'trace.csv'; cp=directory/'history.cp'; out=directory/'native.json'; err=directory/'native.stderr'
        event=None
        if trace_source:
            shutil.copy2(trace_source,trace)
            if tail:
                with trace.open('ab') as stream: stream.write(b'UNCOMMITTED_FAULT_TAIL,not_an_accepted_step\n')
            event=dict(source=record(trace_source),before_restore=record(trace),injected_tail=tail)
        command=[prep['executable'],str(root/(case['name']+'.input')),str(trace),str(cp),*map(str,extra)]
        started=time.monotonic()
        with out.open('xb') as stdout,err.open('xb') as stderr:
            child=subprocess.Popen(command,stdout=stdout,stderr=stderr,preexec_fn=limits)
            timed_out=False
            try: code=child.wait(timeout=5)
            except subprocess.TimeoutExpired:
                child.kill(); code=child.wait(); timed_out=True
        entry=dict(name=name,command=command,return_code=code,wall_seconds=time.monotonic()-started,
            stdout=record(out),stderr=record(err),trace=record(trace) if trace.exists() else None,
            checkpoint=record(cp) if cp.exists() else None,restore_event=event,
            child_max_RSS_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,timed_out=timed_out)
        attempts.append(entry)
        assert not timed_out and code==expected_return,entry
        assert entry['child_max_RSS_KiB']<=256*1024
        metadata=read(out)
        assert metadata['identity']['pid']==child.pid and metadata['identity']['ppid']==identity()['pid']
        assert metadata['identity']['argv']==command and metadata['identity']['boot_id']==identity()['boot_id']
        if code==0:
            assert metadata['completed'] and metadata['committed_bytes']==trace.stat().st_size
        return trace,cp
    try:
        verify(prep['artifacts']+prep['inputs'])
        for case in plan['cases']:
            trace,cp=native(case['name'],case)
            assessment=assess(case,trace); gate(case,assessment); assessed.append(assessment)
            checkpoints.append(dict(name=case['name'],**check_checkpoint(cp,case,trace)))
            assert sum(path.stat().st_size for path in root.rglob('trace.csv'))<=32*2**20
        # Full restart is a distinct control using new traces copied from one
        # prefix. No completed trajectory or old checkpoint is modified.
        case=next(c for c in plan['cases'] if c['name']==plan['restart_case'])
        count=round(float(case['config'][1])/float(case['config'][0])); half=count//2
        prefix,half_cp=native('restart_half',case,['--stop',half])
        for name,stop,tail in [('restart_no_advance',half,False),('restart_split',count,False),('restart_tail_rollback',count,True)]:
            trace,cp=native(name,case,['--restore',half_cp,'--stop',stop],prefix,tail)
            target_prefix=prefix if stop==half else root/case['name']/'trace.csv'
            target_cp=half_cp if stop==half else root/case['name']/'history.cp'
            assert trace.read_bytes()==target_prefix.read_bytes(),(name,'trace mismatch')
            assert cp.read_bytes()==target_cp.read_bytes(),(name,'full checkpoint mismatch')
        verify(prep['artifacts']+prep['inputs'])
    except Exception as error:
        failure=repr(error)
    save(root/'terminal.json',dict(status='runtime_completed_pending_independent_review' if failure is None else 'stopped_for_review',
        finished_at=stamp(),attempts=attempts,runtime_assessments=assessed,checkpoint_checks=checkpoints,failure=failure))
    print(json.dumps(dict(status=read(root/'terminal.json')['status'],attempts=len(attempts),failure=failure)))
    return int(failure is not None)


def review(root):
    assert not (root/'independent_review.json').exists(),'do not repeat receipt'
    prep=read(root/'preparation.json'); plan=read(root/'plan.json'); terminal=read(root/'terminal.json')
    result=dict(status='S2_independent_review_failed',recorded_at=stamp(),identity=identity())
    assessments=[]
    try:
        verify([read(root/'state.json')['preparation']]); verify(prep['artifacts']+prep['inputs'])
        assert terminal['status']=='runtime_completed_pending_independent_review'
        for attempt in terminal['attempts']:
            verify([attempt[name] for name in ['stdout','stderr','trace','checkpoint'] if attempt[name] is not None])
        for case in plan['cases']:
            trace=root/case['name']/'trace.csv'; cp=root/case['name']/'history.cp'
            assessment=assess(case,trace); gate(case,assessment)
            assessment['checkpoint']=check_checkpoint(cp,case,trace); assessments.append(assessment)
        by_name={a['name']:a for a in assessments}
        for group in ['normal','combined','friction']:
            sequence=[a for a in assessments if a['name'].startswith(group+'_dt')]
            assert len(sequence)>=3
            # First-order normal endpoint velocities must improve. Ideal
            # constant-load friction has an exact endpoint/interval solution.
            if group!='friction':
                assert all(sequence[i+1]['velocity_relative_L2'][0]<sequence[i]['velocity_relative_L2'][0] for i in range(len(sequence)-1))
        case=next(c for c in plan['cases'] if c['name']==plan['restart_case'])
        for name,target in [('restart_no_advance','restart_half'),('restart_split',case['name']),('restart_tail_rollback',case['name'])]:
            for file in ['trace.csv','history.cp']:
                assert (root/name/file).read_bytes()==(root/target/file).read_bytes()
        base=rows(root/plan['restart_case']/'trace.csv'); boost=rows(root/'boost_control'/'trace.csv')
        for a,b in zip(base,boost):
            for field in ['q1','vn1','x1','jn','jt','Dmaterial','Dfriction','Dnormal','Dalgorithm']:
                assert abs(a[field]-b[field])<=1e-11,(field,'Galilean control')
            assert abs((b['vt1']-a['vt1'])-3.25)<=1e-11
            assert abs((b['qt1']-a['qt1'])-3.25*a['t1'])<=1e-11
        zero=rows(root/'zero_branch_control'/'trace.csv'); old=read(plan['s1_raw_reference'])['steps']
        assert len(zero)==len(old)
        for a,b in zip(zero,old):
            for field,old_field in [('q1','q_new'),('vn1','v_new'),('jn','impulse')]:
                assert abs(a[field]-b[old_field])<=1e-12
        normal=rows(root/'normal_dt4'/'trace.csv'); zero_mu=rows(root/'zero_friction_control'/'trace.csv')
        for a,b in zip(normal,zero_mu):
            for field in ['q1','vn1','x1','jn','Dmaterial','Dalgorithm']:
                assert a[field]==b[field]
            assert b['jt']==0
        separation=rows(root/'separation_control'/'trace.csv')
        assert all(row['jn']>0 for row in separation[:8])
        assert all(row['jn']==0 and row['jt']==0 and row['mode']==0 for row in separation[8:])
        assert {1,2,3}.issubset(set(by_name['friction_dt3']['modes']))
        # Negative controls retain total dissipation but exchange its physical
        # and algorithmic labels, and separately corrupt nonzero friction force.
        controls=[]; target=root/'negative_control_data'
        target.mkdir()
        for field in ['dissipation_substitution','friction_force']:
            trace=target/(field+'.csv')
            lines=(root/plan['restart_case']/'trace.csv').read_text().splitlines()
            names=lines[0].split(','); cells=lines[80].split(',')
            if field=='dissipation_substitution':
                physical=names.index('Dmaterial'); algorithm=names.index('Dalgorithm')
                cells[physical]=str(float(cells[physical])+float(cells[algorithm])); cells[algorithm]='0'
            else:
                index=next(i for i,line in enumerate(lines[1:],1) if float(line.split(',')[names.index('ft')])!=0)
                cells=lines[index].split(','); cells[names.index('ft')]=str(float(cells[names.index('ft')])*1.01)
            if field=='dissipation_substitution': index=80
            lines[index]=','.join(cells); trace.write_text('\n'.join(lines)+'\n')
            try: assess(case,trace)
            except AssertionError: controls.append(dict(name=field,rejected=True,artifact=record(trace)))
            else: raise AssertionError(('negative control accepted',field))
        result.update(status='S2_reduced_viscoelastic_friction_independently_passed',cases=assessments,
            negative_controls=controls,restart='full_reduced_history_no_advance_split_and_tail_rollback_byte_exact',
            scope='two-coordinate SLS normal constraint and scalar Coulomb prototype only',
            production_PD_integration_qualified=False,calibration_allowed=False)
    except Exception as error:
        result.update(failure=repr(error),completed_assessments=assessments)
    result['terminal']=record(root/'terminal.json'); result['preparation']=record(root/'preparation.json')
    result['state']=record(root/'state.json')
    save(root/'independent_review.json',result)
    print(json.dumps(dict(status=result['status'],failure=result.get('failure'),assessed_cases=len(assessments))))
    return int(not result['status'].endswith('_passed'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('action',choices=['prepare','run','review'])
    parser.add_argument('--root',type=Path,required=True); parser.add_argument('--executable',type=Path)
    parser.add_argument('--prospective',type=Path); args=parser.parse_args()
    if args.action=='prepare': prepare(args.root,args.executable,args.prospective)
    else: sys.exit(globals()[args.action](args.root))
