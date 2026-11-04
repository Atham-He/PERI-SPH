"""Versioned, bounded S1 runner; review runs in a separate process.

Prepare freezes sources, executable, linked libraries, plan, inputs and assessor.
Neither run nor review may redispatch an existing state or rewrite a receipt.
"""
import argparse
import copy
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time


def stamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def record(path):
    path=Path(path)
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):
            digest.update(block)
    return dict(path=str(path),size=path.stat().st_size,sha256=digest.hexdigest())


def save(path,value):
    path=Path(path)
    with path.open('x',encoding='utf8') as stream:
        json.dump(value,stream,ensure_ascii=False,indent=2,allow_nan=False)


def read(path):
    return json.loads(Path(path).read_text())


def verify(entries):
    for expected in entries:
        assert record(expected['path'])==expected, ('frozen_integrity',expected['path'])


def reserves(root):
    # New studies use F. No relaxation of the original sealed E-volume guard.
    host=shutil.disk_usage('/mnt/f').free/2**30
    linux=shutil.disk_usage(root).free/2**30
    assert host>=120 and linux>=150, ('disk_reserve',host,linux)
    return dict(host_F_free_GiB=host,output_free_GiB=linux)


def identity():
    fields=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()
    return dict(pid=os.getpid(),ppid=os.getppid(),start_ticks=fields[19],
                boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),
                argv=Path('/proc/self/cmdline').read_bytes().decode().strip('\0').split('\0'))


def prepare(root,executable,prospective):
    assert not root.exists(), 'existing study is immutable'
    repo=Path(__file__).resolve().parents[1]
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=repo).strip(), 'source must be committed'
    source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()
    parent_resources=reserves(root.parent)
    plan=read(prospective)
    assert plan['time_steps']==[.05,.025,.0125]
    model=plan['physical_model']
    physical=[str(model[name]) for name in ['mass','spring_stiffness','initial_position','initial_velocity',
                                           'plate_total_displacement','duration']]
    root.mkdir()
    (root/'frozen').mkdir()
    source_paths=['src/shared/peridynamics/pd_normal_projection.h',
                  'src/shared/peridynamics/pd_elastic_normal_step.h',
                  'tests/3d_examples/extra_3d_PD_framework_tests/elastic_normal_probe.cpp',
                  'tests/3d_examples/extra_3d_PD_framework_tests/elastic_normal_tests.cpp',
                  'tests/3d_examples/extra_3d_PD_framework_tests/CMakeLists.txt',
                  'tools/peri_elastic_normal.py','tools/peri_elastic_normal_audit.py',
                  'tools/test_peri_elastic_normal.py',
                  'docs/peri/ELASTIC_CONSTRAINT_STAGE_S1_20261104.md']
    artifacts=[]
    for index,relative in enumerate(source_paths):
        original=repo/relative
        target=root/'frozen'/original.name
        shutil.copy2(original,target)
        artifacts.extend([record(original),record(target)])
    frozen_exe=root/'frozen'/executable.name
    shutil.copy2(executable,frozen_exe)
    artifacts.extend([record(executable),record(frozen_exe),record(prospective)])
    linked=subprocess.check_output(['ldd',str(executable)],text=True)
    (root/'frozen'/'ldd.txt').write_text(linked)
    for line in linked.splitlines():
        for token in line.split():
            if token.startswith('/') and Path(token).is_file():
                artifacts.append(record(token))
    artifacts.append(record(root/'frozen'/'ldd.txt'))
    save(root/'plan.json',dict(prospective=plan,physical=physical,source_commit=source_commit,
        controller_identity=identity(),resources=parent_resources,
        velocity_cursor='accepted_endpoint',force_cursor='interval_average_compared_at_midpoint',
        limits=dict(native_wall_seconds=5,native_RSS_MiB=256,output_MiB=32,steps=10000),
        metric_is_physical_mass=False,production_integration=False))
    artifacts.append(record(root/'plan.json'))
    inputs=[]
    for i,dt in enumerate(plan['time_steps']):
        path=root/f'input_{i}.txt'
        path.write_text(' '.join(physical+[str(dt),'0.1'])+'\n')
        inputs.append(record(path))
    save(root/'preparation.json',dict(recorded_at=stamp(),artifacts=artifacts,inputs=inputs,
                                    executable=str(frozen_exe)))
    print(json.dumps(dict(status='prepared_not_passed',root=str(root),source_commit=source_commit)))


def native_limits():
    resource.setrlimit(resource.RLIMIT_CPU,(5,5))
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    resource.setrlimit(resource.RLIMIT_FSIZE,(32*2**20,32*2**20))


def run(root):
    assert not (root/'state.json').exists(), 'never redispatch an attempted state'
    prep=read(root/'preparation.json')
    save(root/'state.json',dict(status='attempted',started_at=stamp(),controller_identity=identity(),
                                preparation=record(root/'preparation.json')))
    attempts=[]; failure=None
    try:
        verify(prep['artifacts']+prep['inputs'])
        for i,item in enumerate(prep['inputs']):
            reserves(root)
            out=root/f'raw_{i}.json'; err=root/f'raw_{i}.stderr'
            started=time.monotonic()
            resource_failure=None
            with out.open('xb') as stdout,err.open('xb') as stderr:
                child=subprocess.Popen([prep['executable'],item['path']],stdout=stdout,stderr=stderr,
                                       preexec_fn=native_limits)
                peak=0
                while child.poll() is None:
                    try:
                        status=Path(f'/proc/{child.pid}/status').read_text()
                        rss=next(int(x.split()[1])*1024 for x in status.splitlines() if x.startswith('VmRSS:'))
                        peak=max(peak,rss)
                    except (FileNotFoundError,StopIteration):
                        pass
                    if peak>256*2**20 or time.monotonic()-started>5:
                        child.kill(); child.wait()
                        resource_failure='native memory/wall-time resource guard'
                        break
                    time.sleep(.005)
                return_code=child.wait()
            attempt=dict(index=i,return_code=return_code,elapsed_seconds=time.monotonic()-started,
                         peak_observed_RSS_bytes=peak,stdout=record(out),stderr=record(err))
            attempts.append(attempt)
            assert resource_failure is None,resource_failure
            assert return_code==0,attempt
            assert out.stat().st_size<=32*2**20
            raw=read(out)
            assert raw['completed'] and raw['identity']['argv']==[prep['executable'],item['path']]
            assert raw['identity']['boot_id']==identity()['boot_id']
            assert raw['identity']['pid']==child.pid and raw['identity']['ppid']==os.getpid()
            assert raw['identity']['start_ticks'].isdigit()
            # Runtime arithmetic gates stop subsequent native cases. Full exact
            # recurrence/continuum assessment remains the separate review's job.
            for row in raw['steps']:
                assert all(math.isfinite(v) for v in row.values()),'nonfinite accepted field'
                assert row['maximum_momentum_residual']<=1e-12
                assert abs(row['energy_residual'])<=1e-12
                assert row['wall_new']-row['q_new']<=2e-11
        verify(prep['artifacts']+prep['inputs'])
    except Exception as error:
        failure=repr(error)
    save(root/'terminal.json',dict(status='runtime_completed_pending_independent_review' if failure is None else 'stopped_for_review',
                                 finished_at=stamp(),attempts=attempts,failure=failure))
    print(json.dumps(read(root/'terminal.json')))
    return 0 if failure is None else 1


def review(root):
    target=root/'independent_review.json'
    assert not target.exists(), 'receipt exists; do not repeat review'
    from peri_elastic_normal_audit import assess
    prep=read(root/'preparation.json'); terminal=read(root/'terminal.json'); plan=read(root/'plan.json')
    receipts=[]; result={'status':'independent_review_failed','recorded_at':stamp(),'identity':identity()}
    try:
        assert terminal['status']=='runtime_completed_pending_independent_review'
        verify([read(root/'state.json')['preparation']])
        verify(prep['artifacts']+prep['inputs'])
        for entry in terminal['attempts']:
            verify([entry['stdout'],entry['stderr']])
        assert len(terminal['attempts'])==3
        for i,dt in enumerate(plan['prospective']['time_steps']):
            raw=read(root/f'raw_{i}.json')
            assessed=assess(raw,plan['physical'],dt)
            receipts.append(assessed)
        finest=receipts[-1]
        assert finest['force_midpoint_relative_L2']<=.02
        assert finest['velocity_endpoint_relative_L2']<=.02
        assert finest['loss_fraction']<=.01
        assert all(receipts[i+1]['force_midpoint_relative_L2']<receipts[i]['force_midpoint_relative_L2'] and
                   receipts[i+1]['velocity_endpoint_relative_L2']<receipts[i]['velocity_endpoint_relative_L2']
                   for i in range(2))
        # Negative controls: replacing physical mass by A in kinetic energy and
        # perturbing a nonzero force must be rejected by the frozen assessor.
        original=read(root/'raw_2.json')
        controls=[]
        for field in ['kinetic_new','reaction_interval']:
            bad=copy.deepcopy(original)
            if field=='kinetic_new':
                a=float(plan['physical'][0])+.5*.0125**2*float(plan['physical'][1])
                bad['steps'][60][field]=.5*a*bad['steps'][60]['v_new']**2
            else:
                bad['steps'][60][field]*=1.01
            try:
                assess(bad,plan['physical'],.0125)
            except AssertionError:
                controls.append(dict(field=field,rejected=True))
            else:
                raise AssertionError(('negative control accepted',field))
        result.update(status='S1_independent_elastic_constraint_review_passed',cases=receipts,
                      negative_controls=controls,scientific_scope='scalar diagonal elastic normal-only research control',
                      production_PD_integration_qualified=False,material_calibration_allowed=False)
    except Exception as error:
        result.update(failure=repr(error),completed_assessments=receipts)
    result['sealed_terminal']=record(root/'terminal.json')
    result['sealed_preparation']=record(root/'preparation.json')
    result['sealed_state']=record(root/'state.json')
    save(target,result)
    print(json.dumps(result))
    return 0 if result['status'].endswith('_passed') else 1


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('action',choices=['prepare','run','review'])
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--executable',type=Path)
    parser.add_argument('--prospective',type=Path)
    args=parser.parse_args()
    if args.action=='prepare':
        prepare(args.root,args.executable,args.prospective)
    else:
        sys.exit(globals()[args.action](args.root))
