"""One specimen: automatic dynamics, complete-state restart, Gc work/contact.
No parameter search, physical state reset or notification/automation service.
The prescribed-opening phase is a spatial virtual-work check, not extra time.
"""
import argparse,hashlib,json,os,shutil,subprocess,traceback
from pathlib import Path
from datetime import datetime,timezone
import xml.etree.ElementTree as ET

def utc():return datetime.now(timezone.utc).isoformat()
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def put(path,value):
    path=Path(path);temp=path.with_name(path.name+'.tmp');temp.write_text(json.dumps(value,indent=2)+chr(10));temp.replace(path)
def identity(pid):
    p=Path('/proc')/str(pid);stat=(p/'stat').read_text().rsplit(')',1)[1].split()
    return {'pid':pid,'start_ticks':stat[19],'boot_id':Path('/proc/sys/kernel/random/boot_id').read_text().strip()}

def run(a):
    root=Path(a.output_root).resolve();root.mkdir(parents=True,exist_ok=True);source=Path(a.source).resolve();build=Path(a.build).resolve();binary=build/'pd_physical_crack_benchmark'
    if subprocess.check_output(['git','status','--porcelain'],cwd=source,text=True).strip():raise RuntimeError('Commit the canonical source before formal verification')
    dep=(build/'CMakeFiles/pd_physical_crack_benchmark.dir/pd_physical_crack_benchmark.cpp.o.d').read_text().replace(chr(92)+chr(10),' ')
    files={Path(p) for p in dep.split() if p.startswith(str(source)+'/') and Path(p).is_file()}
    hashes={str(p.relative_to(source)):sha(p) for p in sorted(files)}
    physics={p:d for p,d in hashes.items() if not p.endswith('/pd_automatic_growth_benchmark.h')}
    fingerprint=hashlib.sha256(json.dumps(physics,sort_keys=True).encode()).hexdigest();put(root/'benchmark_physics_identity.json',{'fingerprint':fingerprint,'files':physics})
    manifest={'created_UTC':utc(),'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=source,text=True).strip(),'binary_sha256':sha(binary),'diagnostic_replay':False,'verified_source_files_sha256':hashes,'physics_fingerprint':fingerprint,'protocol':'same irreversible specimen: automatic dynamics, complete-state restart, prescribed spatial opening work, closure and sliding','capsule_running':False,'stages':[]}
    env=os.environ.copy()
    for k in list(env):
        if k.startswith('PERI_GROWTH_'):del env[k]
    env['PERI_GROWTH_PHYSICS_FINGERPRINT']=fingerprint
    def native(name,output,changes,expected):
        if output.exists() and any(output.iterdir()):raise RuntimeError('Verification output is not fresh: '+str(output))
        if sha(binary)!=manifest['binary_sha256'] or any(sha(source/p)!=d for p,d in hashes.items()):raise RuntimeError('Pinned verification source or executable changed')
        runtime=root/('runtime_'+name);runtime.mkdir(exist_ok=True);log_path=root/(name+'.log');settings=env.copy();settings.update(changes);settings['PERI_COHESIVE_AUDIT_FILE']=str(root/(name+'_cohesive_failure.txt'));settings['PERI_TRANSFER_AUDIT_FILE']=str(root/(name+'_transfer_failure.txt'))
        with log_path.open('w') as log:
            process=subprocess.Popen([str(binary),str(output)],cwd=runtime,env=settings,stdout=log,stderr=subprocess.STDOUT)
            stage={'name':name,'output':str(output),'log':str(log_path),'started_UTC':utc(),**identity(process.pid)};manifest['stages'].append(stage);manifest.update(status='RUNNING_'+name.upper(),pid=stage['pid'],start_ticks=stage['start_ticks'],started_UTC=stage['started_UTC']);put(root/'benchmark_process.json',manifest);print('STARTED',name,process.pid,flush=True)
            try:code=process.wait()
            except BaseException:
                if process.poll() is None:
                    process.terminate()
                    try:process.wait(timeout=10)
                    except subprocess.TimeoutExpired:process.kill();process.wait()
                raise
        stage.update(exit_code=code,ended_UTC=utc());put(output/'process_identity.json',dict(manifest,exit_code=code));put(root/'benchmark_process.json',manifest);print('EXIT',name,code,flush=True)
        if code!=expected:raise RuntimeError(name+' native exit '+str(code)+'; see '+str(log_path))
        return stage
    if a.completed_state:
        completed=Path(a.completed_state).resolve();origin=completed.parent;old=json.loads((origin/'process_identity.json').read_text())
        if old.get('exit_code')!=86:raise RuntimeError('Completed dynamic origin did not finish its dynamic stage')
        if old.get('physics_fingerprint')!=fingerprint:
            previous_receipt=(a.recheck_cohesive_history or a.recheck_cohesive_kinematics or
                              a.recheck_platen_step_bound or a.recheck_layered_gradient_transfer)
            if not a.recheck_platen_transfer and not previous_receipt:
                raise RuntimeError('Completed dynamic origin is not compatible and accepted')
            # For a scoped reconstruction change, first verify the exact
            # saved forces, energy and cohesive state, then run work/contact.
            previous_path=Path(previous_receipt).resolve() if previous_receipt else origin.parent/'benchmark_process.json'
            previous=json.loads(previous_path.read_text())
            unchanged_except={'tests/3d_examples/extra_3d_PD_solid_core_capsule/pd_physical_crack_benchmark.cpp'}
            if a.recheck_cohesive_history:
                unchanged_except.update({'src/shared/peridynamics/pd_cohesive_interface_imex.h',
                                         'src/shared/peridynamics/pd_cohesive_surface_evolution.h'})
            elif a.recheck_layered_gradient_transfer:
                unchanged_except={'src/shared/peridynamics/pd_layered_fracture_reference.h'}
            elif a.recheck_platen_step_bound:
                unchanged_except={'src/shared/peridynamics/pd_extrinsic_crack_growth.h',
                                  'src/shared/peridynamics/pd_solid_reference_coupling.h'}
            elif a.recheck_cohesive_kinematics:
                unchanged_except={'src/shared/peridynamics/pd_cohesive_surface_evolution.h',
                                  'src/shared/peridynamics/pd_extrinsic_crack_growth.h',
                                  'tests/3d_examples/extra_3d_PD_solid_core_capsule/recovery_checkpoint.h'}
            else:
                unchanged_except.update({'src/shared/peridynamics/pd_consistent_reconstruction.h',
                                         'src/shared/peridynamics/pd_fracture_visible_surface.h'})
            old_hashes=previous['verified_source_files_sha256']
            changed={p for p in old_hashes.keys()|hashes.keys() if old_hashes.get(p)!=hashes.get(p)}
            if (previous.get('status')!='PASSED' or previous.get('exit_code')!=0 or
                previous.get('loaded_dynamic_fingerprint',previous.get('physics_fingerprint'))!=old.get('physics_fingerprint') or
                Path(previous['completed_dynamic_checkpoint']).resolve()!=completed or
                previous.get('completed_dynamic_checkpoint_sha256')!=sha(completed) or
                not changed<=unchanged_except):
                raise RuntimeError('Prior published dynamics are outside the requested recheck scope')
            env['PERI_GROWTH_PHYSICS_FINGERPRINT']=old['physics_fingerprint']
            manifest.update(loaded_dynamic_fingerprint=old['physics_fingerprint'],
                            prior_dynamic_receipt_sha256=sha(previous_path),
                            rechecked_platen_transfer=a.recheck_platen_transfer,
                            rechecked_cohesive_history=bool(a.recheck_cohesive_history),
                            rechecked_cohesive_kinematics=bool(a.recheck_cohesive_kinematics),
                            rechecked_platen_step_bound=bool(a.recheck_platen_step_bound),
                            rechecked_layered_gradient_transfer=bool(a.recheck_layered_gradient_transfer))
        manifest['resumed_dynamic_origin']=str(origin)
    else:
        origin=root/'benchmark_dynamic';stage=native('benchmark_dynamic',origin,{'PERI_GROWTH_STOP_AFTER_COMPLETED_DYNAMIC':'1'},86)
        if json.loads((origin/'precrack_stage_acceptance.json').read_text())['status']!='PASSED':raise RuntimeError('Precrack gate failed')
        completed=origin/'automatic_completed_state.bin'
    if not completed.is_file():raise RuntimeError('No complete dynamic checkpoint')
    manifest.update(completed_dynamic_checkpoint=str(completed),completed_dynamic_checkpoint_sha256=sha(completed))
    output=root/'benchmark';native('benchmark_work_contact',output,{'PERI_GROWTH_RESUME_COMPLETED_DYNAMIC':str(completed)},0)
    acceptance=json.loads((output/'acceptance.json').read_text());growth=json.loads((output/'automatic_growth_acceptance.json').read_text());protocol=json.loads((output/'opening_protocol.json').read_text())
    ownership=json.loads((output/'surface_ownership_acceptance.json').read_text())
    platen=json.loads((output/'platen_transfer_acceptance.json').read_text())
    history=json.loads((output/'cohesive_history_roundoff_acceptance.json').read_text())
    if acceptance.get('overall_status')!='PASSED' or any(x.get('status')!='PASSED' for x in [growth,ownership,platen,history]) or not protocol.get('completed_dynamic_state_restart_passed'):raise RuntimeError('Formal fracture, history and transfer acceptance is not PASSED')
    # Preserve both native segments; expose one full dynamic curve and a PVD
    # referring to its original frames. Opening has a separate load-parameter PVD.
    for name in ['automatic_growth.csv','automatic_growth_events.csv']:
        (output/name).rename(output/(Path(name).stem+'_resumed_no_time_advance.csv'));shutil.copyfile(origin/name,output/name)
    for p in origin.glob('preinsert_*.csv'):shutil.copyfile(p,output/p.name)
    (output/'automatic_growth.pvd').rename(output/'automatic_growth_resumed_end.pvd')
    tree=ET.parse(origin/'automatic_growth.pvd')
    for dataset in tree.findall('.//DataSet'):dataset.set('file',os.path.relpath((origin/dataset.get('file')).resolve(),output).replace(os.sep,'/'))
    tree.write(output/'automatic_growth.pvd',encoding='utf-8',xml_declaration=True)
    put(output/'dynamic_curve_origin.json',{'directory':str(origin),'checkpoint_sha256':sha(completed),'same_complete_state':True,'no_time_or_material_history_reset':True,'prescribed_opening_is_separate_spatial_work_test':True})
    manifest.update(status='PASSED',exit_code=0,ended_UTC=utc(),opening_protocol=protocol,surface_ownership=ownership,platen_transfer=platen,cohesive_history_roundoff=history);put(root/'benchmark_process.json',manifest);print('FORMAL_FRACTURE_ACCEPTANCE_PASSED',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--build',required=True);p.add_argument('--output-root',required=True);p.add_argument('--completed-state');scope=p.add_mutually_exclusive_group();scope.add_argument('--recheck-platen-transfer',action='store_true');scope.add_argument('--recheck-cohesive-history',metavar='PRIOR_ACCEPTANCE_MANIFEST');scope.add_argument('--recheck-cohesive-kinematics',metavar='PRIOR_ACCEPTANCE_MANIFEST');scope.add_argument('--recheck-platen-step-bound',metavar='PRIOR_ACCEPTANCE_MANIFEST')
    scope.add_argument('--recheck-layered-gradient-transfer',metavar='PRIOR_ACCEPTANCE_MANIFEST')
    a=p.parse_args()
    try:run(a)
    except BaseException as error:
        root=Path(a.output_root);root.mkdir(parents=True,exist_ok=True);path=root/'benchmark_process.json';state=json.loads(path.read_text()) if path.exists() else {};state.update(status='FAILED',exit_code=1,error=str(error),traceback=traceback.format_exc(),recorded_UTC=utc());put(path,state);raise
