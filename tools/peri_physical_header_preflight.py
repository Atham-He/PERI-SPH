"""New producer identity for actual initial-header path; ZERO PD advances."""
import argparse,hashlib,json,os,shutil,subprocess,sys,sysconfig,traceback
from pathlib import Path
from peri_archive_component_build import prepare,compile_and_freeze
from peri_physical_session_preflight import assess
from peri_actual_checkpoint_campaign import Execution,verified_execution_component
from peri_elastic_normal import read,save,verify,record,reserves,stamp,identity

def main(root):
    repo=Path(__file__).resolve().parents[1];base=root.parent;prior=base/'engineering-s4-physical-session-v2-preflight-v1'
    assert base==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists(),'unique new header state only'
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    reserves(base);component=verified_execution_component(base)
    old_review=read(prior/'independent_preflight.json');verify([read(prior/'terminal.json')['receipt']])
    assert old_review['zero_step_owned_install_resave_bytes_identical'] and old_review['new_PD_steps']==0
    assert len(old_review['clock_semantic_rejections'])==6
    root.mkdir(mode=0o700)
    try:
        stage='tests/3d_examples/extra_3d_PD_framework_tests/'
        names=[stage+x for x in ['bounded_gel_physical_checkpoint_probe.cpp','bounded_gel_physical_time_observer.h','coulomb_contact_probe_io.h','bounded_gel_probe_io.h']]
        names += ['tools/'+x for x in ['peri_physical_header_preflight.py','peri_physical_header_audit.py','peri_physical_checkpoint_support.py',
            'peri_physical_session_preflight.py','peri_archive_component_build.py','peri_run_body_archive_engineering.py']]
        prepared=prepare(root,repo,names,stage+'bounded_gel_physical_checkpoint_probe.cpp','physical_header')
        reference=root/'frozen/reference';reference.mkdir()
        software=[]
        for source in sorted((repo/'tools').glob('*.py')):
            target=reference/source.name;shutil.copy2(source,target);software += [record(source),record(target)]
        for name in ['peri_physical_header_audit.py','peri_fixed_python_raw_audit.py']:
            source=repo/'tools'/name;target=root/'frozen'/name;shutil.copy2(source,target);software += [record(source),record(target)]
        # Backend selection happens in a FRESH child before Decimal/Fraction or
        # reference imports. The actual _pydecimal entry comes from frozen/.
        old_backend=base/'cross-runtime-v19-raw/frozen/_pydecimal.py'
        target=root/'frozen/_pydecimal.py';shutil.copy2(old_backend,target);software += [record(old_backend),record(target)]
        stdlib=Path(sysconfig.get_path('stdlib'))
        for source in sorted(stdlib.rglob('*')):
            if source.is_file() and source.suffix in ['.py','.so'] and not {'__pycache__','test','dist-packages'}.intersection(source.relative_to(stdlib).parts):
                target=root/'python_support'/source.relative_to(stdlib);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
                software += [record(source),record(target)]
        software += [record(sys.executable)]
        for line in Path('/proc/self/maps').read_text().splitlines():
            p=Path(line.split()[-1])
            if str(p).startswith('/') and p.is_file() and (p.suffix=='.so' or '.so.' in p.name):software.append(record(p))
        production=repo/'build/peri-latest-20261103/extra_3d_PD_gel_compression/bin/extra_3d_PD_gel_compression'
        assert record(production)['sha256']=='bf398b1415cbded594c4164760989dc91c99f7935ec3c4f047eb7c8bbee88b69'
        case=dict(name='v2_actual_header_ZERO_steps',mode=1,steps=80,protocol_kind=2,dt=.005,stroke=.01,rotation=0,azimuths=8,rho0=1.,
                  parameters=[5.,1.,.5,.3],eta=.1,support=1.25,cap_width=1/64,tangent_stroke=[.0006,-.0004],mu=[.22,.22],
                  checkpoint_steps=[16,40,64,80],physical_time=dict(denominator=1000000,nominal_ticks=5000,
                      normal_event_ticks=[160000,240000,400000],tangent_event_ticks=[20000,60000,80000]))
        anchors=[record(prior/'independent_preflight.json'),record(prior/'terminal.json'),record(production)]
        plan=dict(prepared_at=stamp(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
                  source_snapshots=prepared['source_snapshots'],compile_command=prepared['compile_command'],link_command=prepared['link_command'],
                  software=software,anchors=anchors,execution_component=component,case=case,
                  native_modes=['--initialize(required new producer identity)','--describe-initial(NEW actual header path)'],
                  reference_backend='fixed pure Python95 selected before Fraction/reference imports',
                  compile_budget=dict(CPU_seconds=180,wall_seconds=240,AS_GiB=2),native_budget=dict(CPU_seconds=30,wall_seconds=45,AS_GiB=1,threads=1),
                  initial_header_reference_budget=dict(CPU_seconds=180,wall_seconds=240,AS_GiB=2,threads=1),
                  normalized_gate=1e-10,zero_PD_steps=True,all_trial_material_T_qualification=False,
                  declared_header_negative_controls=['old_format','wrong_denominator','wrong_nominal_ticks','changed_physical_end','wrong_total_steps','wrong_checkpoint_selection'],
                  failure_policy='Preserve exact failure/raw/returncode and stop dependent modes;no redispatch/fallback/reset/gate changes.',
                  total_bound_GiB=1,host_F_reserve_GiB=120,output_reserve_GiB=150,S4_complete=False,calibration_allowed=False)
        verify(software+anchors);save(root/'plan.json',plan);save(root/'state.json',dict(status='new_actual_header_preflight_attempted',identity=identity(),started_at=stamp()))
        build=compile_and_freeze(prepared);shutil.copy2(prior/'original.input',root/'original.input')
        frozen=dict(plan=record(root/'plan.json'),build=build,software=software,inputs=[record(root/'original.input')],anchors=anchors)
        save(root/'frozen_execution.json',frozen);executor=Execution(root,plan);phases=[]
        def check():verify(software+anchors+build['support']+build['source_snapshots']+[build['ELF'],build['source_manifest'],build['build_runtime_manifest']])
        check();r=executor.run('initialize',[build['ELF']['path'],'--initialize',str(root)],wall=45,cpu=30,gib=1)
        phases.append(record(root/'initialize.attempt.json'));rows=[json.loads(x) for x in (root/'initialize.stdout').read_text().splitlines()]
        assert rows[-1]['type']=='initialization';initial=assess(root,base/'engineering-actual-checkpoint-driver-preflight-v5',build,rows[-1]);trusted=initial['contract_sha256']
        check();r=executor.run('describe_initial',[build['ELF']['path'],'--describe-initial',str(root),trusted],wall=45,cpu=30,gib=1)
        phases.append(record(root/'describe_initial.attempt.json'));assert r['stderr']['size']==0
        inputs=[record(root/name) for name in ['original.input','initial_candidate.bin','initial.phy','contract.bin']]
        admission=dict(case=case,raw=r['stdout'],artifacts=software+anchors+build['support']+build['source_snapshots']+inputs+
                       [record(root/'describe_initial.attempt.json'),build['ELF'],build['source_manifest'],build['build_runtime_manifest']],
                       new_PD_steps=0,CPU_seconds=180,wall_seconds=240,AS_GiB=2)
        save(root/'header_reference_admission.json',admission)
        check();reference_attempt=executor.run('header_reference',[sys.executable,str(root/'frozen/peri_physical_header_audit.py'),'--root',str(root)],wall=240,cpu=180,gib=2)
        phases.append(record(root/'header_reference.attempt.json'));q=read(root/'header_independent_review.json')
        assert 'failure' not in q and q['new_PD_steps']==q['independent_trial_material_points']==0
        assert q['maximum_normalized_error']<=1e-10 and len(q['clock_header_negative_controls'])==6
        verify(admission['artifacts']);check()
        receipt=dict(status='actual_v2_header_initial_geometry_and_complete_producer_identity_independently_qualified_ZERO_steps',
                     finished_at=stamp(),source_commit=plan['source_commit'],plan=record(root/'plan.json'),frozen_execution=record(root/'frozen_execution.json'),
                     actual_header_reference=record(root/'header_independent_review.json'),initial_complete_bytes=initial,phases=phases,
                     new_PD_steps=0,full_nonzero_restart_qualified=False,time_convergence_qualified=False,S4_complete=False,calibration_allowed=False)
        save(root/'independent_preflight.json',receipt);save(root/'terminal.json',dict(status=receipt['status'],return_code=0,finished_at=stamp(),receipt=record(root/'independent_preflight.json'),new_PD_steps=0))
        print(json.dumps(dict(status=receipt['status'],header_maximum_error=q['maximum_normalized_error'],new_PD_steps=0)))
    except BaseException as e:
        save(root/'failure.json',dict(recorded_at=stamp(),failure=repr(e),traceback=traceback.format_exc(),new_PD_steps=0));raise

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);a=p.parse_args();main(a.root)
