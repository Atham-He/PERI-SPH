"""One physical capsule history, with a verified full-checkpoint handoff.
No parameter search, retry, memory cap, notification service or scheduler.
The first accepted load point closes segment 0; segment 1 continues its EXACT
full state. A failure is saved/reported and never converted into a new case.
"""
import argparse, csv, hashlib, json, math, os, struct, subprocess, time, traceback
from pathlib import Path
from datetime import datetime, timezone
import xml.etree.ElementTree as ET

def utc():return datetime.now(timezone.utc).isoformat()
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for data in iter(lambda:f.read(1024*1024),b''):h.update(data)
    return h.hexdigest()
def put(path,value):
    path=Path(path);tmp=path.with_name(path.name+'.tmp');tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)
def process_identity(pid):
    p=Path('/proc')/str(pid);stat=(p/'stat').read_text().rsplit(')',1)[1].split()
    init=Path('/proc/1/stat').read_text().rsplit(')',1)[1].split()
    return {'pid':pid,'start_ticks':stat[19],'state':stat[0],'boot_id':Path('/proc/sys/kernel/random/boot_id').read_text().strip(),'distro_init_start_ticks':init[19],'argv':[s.decode(errors='replace') for s in (p/'cmdline').read_bytes().split(b'\0') if s]}
def row(path):
    path=Path(path)
    if not path.exists():return None
    with path.open('rb') as f:
        header=f.readline().decode().strip().split(',');f.seek(0,2);size=f.tell();f.seek(max(0,size-32768));lines=f.read().decode(errors='replace').splitlines()
    for line in reversed(lines):
        values=line.strip().split(',')
        if len(values)==len(header) and values and values[0].isdigit():return dict(zip(header,values))
    return None

def index(root,segments):
    frames={}
    for entry in segments:
        output=Path(entry['output'])
        p=output/'capsule.pvd'
        if not p.exists():continue
        try:
            for a in ET.parse(p).findall('.//DataSet'):
                f=(p.parent/a.attrib['file']).resolve()
                t=float(a.attrib['timestep'])
                if f.is_file() and t<=entry.get('trajectory_end_s',float('inf')):frames[t]=f
        except (ET.ParseError,OSError):continue
    xml=ET.Element('VTKFile',type='Collection',version='0.1');collection=ET.SubElement(xml,'Collection')
    for t,f in sorted(frames.items()):ET.SubElement(collection,'DataSet',timestep=format(t,'.17g'),file=os.path.relpath(f,root).replace(os.sep,'/'))
    tmp=root/'capsule.pvd.tmp'
    ET.ElementTree(xml).write(tmp,encoding='utf-8',xml_declaration=True)
    try:
        tmp.replace(root/'capsule.pvd')
    except PermissionError:
        # Windows viewers can hold the old index open without delete sharing.
        # Leave its readable frame list in place and retry on the next refresh.
        print('Trajectory index is locked; the capsule continues running',flush=True)

def accepted_checkpoint(path, fingerprint, step):
    """Identify a retained slot at a completed stage without altering its bytes.

    The native reader still validates the full checksum and restored forces.
    """
    with path.open('rb') as f:
        def read(fmt):
            return struct.unpack(fmt,f.read(struct.calcsize(fmt)))[0]
        def string():
            length=read('<Q')
            if length>1048576:raise RuntimeError('Invalid checkpoint metadata string')
            return f.read(length).decode()
        magic,version,identity,folder=string(),read('<Q'),string(),Path(string())
        if (magic,version,identity)!=('PERI_SNSNI_RESTART',4,fingerprint) or folder.resolve()!=path.parent:
            raise RuntimeError('Accepted checkpoint does not match its original history')
        count=read('<Q')
        if count!=3:raise RuntimeError('Unexpected accepted checkpoint CSV layout')
        names=[string() for _ in range(count)]
        if names!=['progress.csv','equilibria.csv','fracture_events.csv'] or read('<Q')!=count:
            raise RuntimeError('Unexpected accepted checkpoint CSV layout')
        sizes=[read('<Q') for _ in range(count)]
        stage,stage_time,elapsed,checks=read('<Q'),read('<d'),read('<d'),read('<Q')
    with (folder/'equilibria.csv').open('rb') as f:
        rows=list(csv.DictReader(f.read(sizes[1]).decode().splitlines()))
    if not rows or int(rows[-1]['step'])!=step or int(rows[-1]['stage'])+1!=stage or stage_time!=0 or checks!=0 or float(rows[-1]['elapsed_s'])!=elapsed:
        raise RuntimeError('Retained slot does not end at the requested accepted load point')
    return {'version':version,'fingerprint':identity,'step':step,'elapsed_s':elapsed}

def run(a):
    root=Path(a.output_root).resolve();root.mkdir(parents=True,exist_ok=True)
    source=Path(a.source).resolve();binary=Path(a.binary).resolve();cfg=Path(a.config).resolve();cache=Path(a.cells).resolve();bank=Path(a.bank).resolve()
    if bool(a.resume_checkpoint)!=bool(a.resume_manifest):
        raise RuntimeError('Continuation requires both checkpoint and its original production manifest')
    origin=Path(a.resume_manifest).resolve() if a.resume_manifest else None
    previous=json.loads(origin.read_text()) if origin else None
    previous_execution=previous['identity'].get('execution','cpu') if previous else 'cpu'
    execution=a.execution or previous_execution
    old_identity=previous['identity'] if previous else {}
    spectral_mode='mass_normalized_gradient_gram'
    tighter_shell_step=a.tighter_shell_step_bound or old_identity.get('shell_step_bound')==spectral_mode
    reaction_window=a.reaction_settling_window_s
    if reaction_window is None:
        reaction_window=old_identity.get('reaction_settling_window_s')
    if reaction_window is not None and (not math.isfinite(reaction_window) or reaction_window<=0):
        raise ValueError('Reaction settling window must be finite and positive')

    acceptance=json.loads((Path(a.benchmark_root)/'acceptance.json').read_text())
    preflight=json.loads(Path(a.preflight).read_text())
    benchmark_manifest=json.loads(Path(a.benchmark_manifest).read_text());preflight_manifest=json.loads(Path(a.preflight_manifest).read_text())
    if acceptance.get('overall_status')!='PASSED' or preflight.get('status')!='PASSED':
        raise RuntimeError('Physical fracture or coupled-driver acceptance failed')
    if a.fracture_contact_only:
        if (acceptance.get('surface_law')!='contact_only' or benchmark_manifest.get('fracture_contact_only') is not True
                or preflight.get('fracture_contact_only') is not True
                or '--fracture-contact-only' not in preflight_manifest.get('argv',[])):
            raise RuntimeError('Contact-only fracture requires a matched benchmark and preflight')
    else:
        for name in ('automatic_growth','surface_ownership','platen_transfer','cohesive_history_roundoff'):
            gate=json.loads((Path(a.benchmark_root)/(name+'_acceptance.json')).read_text())
            if gate.get('status')!='PASSED':raise RuntimeError('Physical fracture gate failed: '+name)
        if a.column_tied_interface:
            column=json.loads((Path(a.benchmark_root)/'column_interface_acceptance.json').read_text())
            if column.get('status')!='PASSED':raise RuntimeError('Fixed-column fracture benchmark failed')
    # Older CPU receipts predate the execution option. CUDA always needs
    # acceptance from the same backend for both fracture and coupled mechanics.
    def recorded_execution(receipt):
        return next((arg.partition('=')[2] for arg in reversed(receipt.get('argv',[]))
                     if arg.startswith('--execution=')), 'cpu')

    if recorded_execution(benchmark_manifest)!=execution or recorded_execution(preflight_manifest)!=execution:
        raise RuntimeError('Fracture benchmark and coupled preflight must use the requested execution backend')
    if execution=='gpu':
        coupled_parameters=json.loads(Path(a.preflight).with_name('parameters.json').read_text())
        if coupled_parameters.get('execution_backend')!='cuda':
            raise RuntimeError('Coupled preflight did not run with CUDA')
        if a.fracture_contact_only and acceptance.get('execution_backend')!='cuda':
            raise RuntimeError('Fracture benchmark did not run with CUDA')
    if benchmark_manifest.get('exit_code')!=0 or benchmark_manifest.get('diagnostic_replay') is not False:
        raise RuntimeError('Formal benchmark did not exit successfully or was a diagnostic replay')
    if a.column_tied_interface and '--column-tied-interface' not in preflight_manifest.get('argv',[]):
        raise RuntimeError('Fixed-column interface requires a matched coupled preflight')
    if subprocess.check_output(['git','status','--porcelain'],cwd=source,text=True).strip():raise RuntimeError('Pin a clean canonical source before the production history')
    values=dict(line.split('=',1) for line in cfg.read_text().splitlines() if line.strip() and not line.lstrip().startswith('#'))
    material=Path(values['dpc-table'].strip());identity={'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=source,text=True).strip(),'binary_sha256':sha(binary),'config_sha256':sha(cfg),'DPC_table_sha256':sha(material),'core_cache_sha256':sha(cache),'crack_bank_sha256':sha(bank),'checkpoint_schema':(5 if a.column_tied_interface else 4)+(2 if a.fracture_contact_only else 0)}
    identity['execution']=execution
    if a.fracture_contact_only:identity['fracture_surface_law']='contact_only'
    if a.affine_crack_contact:
        if (not a.fracture_contact_only or preflight.get('affine_crack_contact') is not True
                or '--affine-crack-contact' not in preflight_manifest.get('argv',[])):
            raise RuntimeError('Affine crack contact requires a matched contact-only preflight')
        identity['crack_contact_interpolation']='affine'
        identity['checkpoint_schema']+=2
    if a.affine_inner_interface:
        if (not a.affine_crack_contact or a.column_tied_interface or preflight.get('affine_inner_interface') is not True
                or '--affine-inner-interface' not in preflight_manifest.get('argv',[])):
            raise RuntimeError('Affine inner interface requires a matched affine contact preflight')
        identity['inner_shell_interpolation']='affine'
        identity['checkpoint_schema']=10
    legacy_checkpoint_schema = identity['checkpoint_schema']
    if a.fracture_contact_only:
        identity['checkpoint_schema'] = 11
        identity['legacy_contact_checkpoint_schema'] = legacy_checkpoint_schema
    settling_option=None
    preflight_args=preflight_manifest.get('argv',[])
    if tighter_shell_step:
        if '--tighter-shell-step-bound' not in preflight_args:
            raise RuntimeError('The shell time-step bound requires a matching coupled preflight')
        identity['shell_step_bound']=spectral_mode
    if reaction_window is not None:
        settling_option=f'--reaction-settling-window-s={reaction_window:.17g}'
        if settling_option not in preflight_args:
            raise RuntimeError('The plate-reaction criterion requires a matching coupled preflight')
        identity['reaction_settling_window_s']=reaction_window
    if previous and execution!=previous_execution:
        for option in ('shell_step_bound','reaction_settling_window_s'):
            if identity.get(option)!=old_identity.get(option):
                raise RuntimeError('A backend handoff must keep the time-step and reaction acceptance settings')
    # A PASSED file from an earlier executable/input cannot authorize a new run.
    for name,path in [('binary_sha256',binary),('config_sha256',cfg),('DPC_table_sha256',material),('core_cache_sha256',cache),('crack_bank_sha256',bank)]:
        if preflight_manifest.get(name)!=sha(path):raise RuntimeError('Coupled preflight identity mismatch: '+name)
    if preflight_manifest.get('exit_code')!=0:raise RuntimeError('Coupled preflight did not exit successfully')
    if a.local_platen_step_bound and '--local-platen-step-bound' not in preflight_manifest.get('argv',[]):
        raise RuntimeError('Local platen step restriction requires a matched coupled preflight')
    for receipt in [benchmark_manifest,preflight_manifest]:
        hashes=receipt.get('verified_source_files_sha256',{})
        if not hashes:raise RuntimeError('Verification receipt has no source identity')
        for relative,digest in hashes.items():
            if sha(source/relative)!=digest:raise RuntimeError('Verified source changed after acceptance: '+relative)
    benchmark_binary=Path(a.benchmark_binary)
    if benchmark_manifest.get('binary_sha256')!=sha(benchmark_binary):raise RuntimeError('Accepted benchmark executable has changed')
    execution_fingerprint=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest();fingerprint=execution_fingerprint;owner=process_identity(os.getpid())
    inherited=[];checkpoint=None
    if a.resume_accepted_step and (not a.column_tied_interface or not a.resume_checkpoint):
        raise RuntimeError('An explicit retained equilibrium is supported for a fixed-column handoff')
    if a.resume_checkpoint:
        checkpoint=Path(a.resume_checkpoint).resolve()
        if a.resume_accepted_step:
            latest=accepted_checkpoint(checkpoint,previous['fingerprint'],a.resume_accepted_step)
        else:
            latest=json.loads((checkpoint.parent/'checkpoint.latest.json').read_text())
            if latest['fingerprint']!=previous['fingerprint'] or Path(latest['file']).resolve()!=checkpoint:raise RuntimeError('Continuation checkpoint differs from its original manifest')
        allowed_versions={identity['checkpoint_schema']}
        if a.resume_accepted_step and not a.fracture_contact_only:allowed_versions.add(4)
        if a.fracture_contact_only:
            allowed_versions.update({legacy_checkpoint_schema, legacy_checkpoint_schema - 2})
        if latest['version']!=previous['identity']['checkpoint_schema'] or latest['version'] not in allowed_versions:
            raise RuntimeError('Continuation checkpoint schema differs from its permitted source')
        for name in ['config_sha256','DPC_table_sha256','core_cache_sha256','crack_bank_sha256']:
            if identity[name]!=previous['identity'][name]:raise RuntimeError('Continuation would change a frozen physical input: '+name)
        fingerprint=previous['fingerprint']
        inherited=previous['segments']
        if not inherited or root==origin.parent:raise RuntimeError('Continuation needs an original history and a fresh output root')
        if a.resume_accepted_step:
            entries=[entry for entry in inherited if Path(entry['output']).resolve()==checkpoint.parent]
            if len(entries)!=1 or entries[0] is not inherited[-1]:raise RuntimeError('Accepted checkpoint must belong to the last original segment')
            entries[0]['trajectory_end_s']=latest['elapsed_s']
        elif execution!=previous_execution:
            # A live backend handoff may leave later samples in the origin.
            # Join the trajectory at the common committed checkpoint.
            with (checkpoint.parent/'progress.csv').open(newline='') as stream:
                saved=next((sample for sample in csv.DictReader(stream)
                            if int(sample['step'])==latest['step']),None)
            if saved:
                inherited[-1]['trajectory_end_s']=float(saved['elapsed_s'])
        receipt={'recorded_UTC':utc(),'status':'INPUTS_VERIFIED_NATIVE_STATE_CHECK_PENDING','origin_manifest':str(origin),'origin_manifest_sha256':sha(origin),'origin_source_commit':previous['identity']['source_commit'],'source_commit':identity['source_commit'],'checkpoint':str(checkpoint),'checkpoint_sha256':sha(checkpoint),'checkpoint_step':latest['step'],'history_fingerprint':fingerprint,'execution_fingerprint':execution_fingerprint,'frozen_inputs_unchanged':True,'checkpoint_bytes_modified':False,'loading_history_reset':False,'native_force_and_checksum_verification_required':True}
        receipt.update(origin_execution=previous_execution,execution=execution)
        put(root/'resume_compatibility.json',receipt)
    manifest={'created_UTC':utc(),'source':str(source),'binary':str(binary),'input_config':str(cfg),'identity':identity,'fingerprint':fingerprint,'execution_fingerprint':execution_fingerprint,'owner':owner,'same_history_checkpoint_handoff':True,'material_parameters_changed':False,'capsule_kinematics':'NOSB_PD_shell_plus_SNSNI_solid_SPH_core','benchmark_root':str(Path(a.benchmark_root).resolve()),'preflight':str(Path(a.preflight).resolve()),'segments':list(inherited)}
    if (root/'production_manifest.json').exists():raise RuntimeError('Output root already has a production manifest')
    put(root/'production_manifest.json',manifest)
    sequence=[len(inherited)] if checkpoint else range(2)
    for segment in sequence:
        if sha(binary)!=identity['binary_sha256'] or sha(cfg)!=identity['config_sha256'] or sha(material)!=identity['DPC_table_sha256'] or sha(bank)!=identity['crack_bank_sha256'] or sha(cache)!=identity['core_cache_sha256']:raise RuntimeError('Pinned source/input/binary changed before continuation')
        output=root/('mechanics_segment_'+str(segment));runtime=root/('runtime_segment_'+str(segment));runtime.mkdir(exist_ok=True)
        if output.exists() and any(output.iterdir()):raise RuntimeError('Production output segment already contains results')
        cmd=[str(binary),'--config='+str(cfg),'--snsni-cells='+str(cache),'--fracture-facets='+str(bank),'--output='+str(output),'--checkpoint-fingerprint='+fingerprint,'--execution='+execution]
        if a.column_tied_interface:cmd.append('--column-tied-interface')
        if a.fracture_contact_only:cmd.append('--fracture-contact-only')
        if a.affine_crack_contact:cmd.append('--affine-crack-contact')
        if a.affine_inner_interface:cmd.append('--affine-inner-interface')
        if a.local_platen_step_bound:cmd.append('--local-platen-step-bound')
        if tighter_shell_step:cmd.append('--tighter-shell-step-bound')
        if settling_option:cmd.append(settling_option)
        if a.checkpoint_wall_seconds is not None:cmd.append('--checkpoint-wall-seconds='+str(a.checkpoint_wall_seconds))
        if checkpoint is None:cmd.append('--checkpoint-stop-after-equilibrium')
        else:cmd.append('--resume-state='+str(checkpoint))
        log_path=root/('native_segment_'+str(segment)+'.log')
        with log_path.open('w') as log:
            process=subprocess.Popen(cmd,cwd=runtime,stdout=log,stderr=subprocess.STDOUT)
            try:native=process_identity(process.pid)
            except OSError:native={'pid':process.pid,'argv':cmd}
            entry={'segment':segment,'output':str(output),'log':str(log_path),'native':native,'started_UTC':utc(),'argv':cmd};manifest['segments'].append(entry);put(root/'production_manifest.json',manifest)
            print('STARTED',segment,process.pid,flush=True)
            try:
                while process.poll() is None:
                    state={'recorded_UTC':utc(),'status':'RUNNING','owner':owner,'native':native,'segment':segment,'output':str(output),'latest_progress':row(output/'progress.csv'),'latest_equilibrium':row(output/'equilibria.csv'),'source_commit':identity['source_commit'],'fingerprint':fingerprint,'capsule_running':True}
                    put(root/'owner_progress.json',state)
                    validation=output/'resume_validation.json'
                    handoff=output/('affine_inner_interface_handoff.json' if a.affine_inner_interface else 'affine_crack_contact_handoff.json')
                    if a.affine_crack_contact and inherited and validation.exists() and handoff.exists():
                        change=json.loads(handoff.read_text())
                        if change['status']=='PASSED' and 'trajectory_end_s' not in inherited[-1]:
                            inherited[-1]['trajectory_end_s']=change['elapsed_s']
                            put(root/'production_manifest.json',manifest)
                    index(root,manifest['segments'])
                    if a.resume_checkpoint and validation.exists():
                        native_check=json.loads(validation.read_text());receipt.update(status='PASSED' if native_check.get('status')=='PASSED' else 'FAILED',native_state_validation=native_check);put(root/'resume_compatibility.json',receipt)
                    if a.status_file:
                        target=Path(a.status_file);review=json.loads(target.read_text()) if target.exists() else {};review.update(status='PHYSICAL_CAPSULE_RUNNING',physical_capsule_running=True,capsule_restart_allowed=True,source_commit=identity['source_commit'],source_worktree_clean=True,owner_progress=str(root/'owner_progress.json'),production_manifest=str(root/'production_manifest.json'));put(target,review)
                    time.sleep(30)
            except BaseException:
                if process.poll() is None:
                    process.terminate()
                    try:process.wait(timeout=10)
                    except subprocess.TimeoutExpired:process.kill();process.wait()
                raise
            code=process.wait();entry.update(exit_code=code,ended_UTC=utc());put(root/'production_manifest.json',manifest);index(root,manifest['segments'])
        if segment==0 and code==86:
            latest=json.loads((output/'checkpoint.latest.json').read_text());checkpoint=Path(latest['file']);accepted=row(output/'equilibria.csv')
            if latest['version']!=identity['checkpoint_schema'] or latest['fingerprint']!=fingerprint or not accepted:raise RuntimeError('First equilibrium did not publish its compatible complete physical checkpoint')
            receipt={'recorded_UTC':utc(),'checkpoint':str(checkpoint),'checkpoint_sha256':sha(checkpoint),'accepted_equilibrium':accepted,'fingerprint':fingerprint,'native_FNV_verification_required_on_load':True,'no_loading_or_material_history_reset':True}
            put(root/'checkpoint_handoff.json',receipt);print('CONTINUING_COMPLETE_STATE',checkpoint,flush=True);continue
        if code!=0:raise RuntimeError('Physical capsule native exit '+str(code)+'; see '+str(log_path))
        if not (output/'summary.json').exists():raise RuntimeError('Native exit zero without completed-curve summary')
        state={'recorded_UTC':utc(),'status':'COMPLETE','capsule_running':False,'latest_equilibrium':row(output/'equilibria.csv'),'summary':json.loads((output/'summary.json').read_text()),'source_commit':identity['source_commit']};put(root/'owner_progress.json',state)
        if a.status_file:
            target=Path(a.status_file);review=json.loads(target.read_text());review.update(status='PHYSICAL_CAPSULE_COMPLETE',physical_capsule_running=False,owner_progress=str(root/'owner_progress.json'));put(target,review)
        return
    raise RuntimeError('Continuation protocol ended without completed history')

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for option in ['source','binary','config','cells','bank','output-root','benchmark-root','preflight','benchmark-manifest','benchmark-binary','preflight-manifest']:p.add_argument('--'+option,required=True)
    p.add_argument('--status-file');p.add_argument('--resume-checkpoint');p.add_argument('--resume-manifest');p.add_argument('--resume-accepted-step',type=int);p.add_argument('--column-tied-interface',action='store_true')
    p.add_argument('--checkpoint-wall-seconds',type=int,help='minimum wall seconds between periodic checkpoints (default: 600)')
    p.add_argument('--fracture-contact-only',action='store_true')
    p.add_argument('--affine-crack-contact',action='store_true')
    p.add_argument('--affine-inner-interface',action='store_true')
    p.add_argument('--local-platen-step-bound',action='store_true',help='check the platen step limit using the potentially contacting traces')
    p.add_argument('--tighter-shell-step-bound',action='store_true',
                   help='use the mass-normalized bound on the active shell gradient')
    p.add_argument('--reaction-settling-window-s',type=float,
                   help='accept a held load using the platen reaction stability window')
    p.add_argument('--execution',choices=('cpu','gpu'),
                   help='solver backend (default: original backend on continuation, CPU for a new run)')
    a=p.parse_args()
    if a.checkpoint_wall_seconds is not None and a.checkpoint_wall_seconds<=0:
        p.error('--checkpoint-wall-seconds requires a positive integer')
    try:run(a)
    except BaseException as error:
        root=Path(a.output_root);root.mkdir(parents=True,exist_ok=True);put(root/'owner_progress.json',{'recorded_UTC':utc(),'status':'FAILED','capsule_running':False,'error':str(error),'traceback':traceback.format_exc()})
        if a.status_file:
            target=Path(a.status_file);review=json.loads(target.read_text()) if target.exists() else {};review.update(status='PHYSICAL_CAPSULE_STOPPED',physical_capsule_running=False,capsule_restart_allowed=False,owner_progress=str(root/'owner_progress.json'));put(target,review)
        raise
