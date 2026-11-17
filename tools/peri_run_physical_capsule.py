"""One physical capsule history, with a verified full-checkpoint handoff.
No parameter search, retry, memory cap, notification service or scheduler.
The first accepted load point closes segment 0; segment 1 continues its EXACT
full state. A failure is saved/reported and never converted into a new case.
"""
import argparse, csv, hashlib, json, os, subprocess, time, traceback
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
    for output in segments:
        p=output/'capsule.pvd'
        if not p.exists():continue
        try:
            for a in ET.parse(p).findall('.//DataSet'):
                f=(p.parent/a.attrib['file']).resolve()
                if f.is_file():frames[float(a.attrib['timestep'])]=f
        except (ET.ParseError,OSError):continue
    xml=ET.Element('VTKFile',type='Collection',version='0.1');collection=ET.SubElement(xml,'Collection')
    for t,f in sorted(frames.items()):ET.SubElement(collection,'DataSet',timestep=format(t,'.17g'),file=os.path.relpath(f,root).replace(os.sep,'/'))
    tmp=root/'capsule.pvd.tmp';ET.ElementTree(xml).write(tmp,encoding='utf-8',xml_declaration=True);tmp.replace(root/'capsule.pvd')

def run(a):
    root=Path(a.output_root).resolve();root.mkdir(parents=True,exist_ok=True)
    source=Path(a.source).resolve();binary=Path(a.binary).resolve();cfg=Path(a.config).resolve();cache=Path(a.cells).resolve();bank=Path(a.bank).resolve()
    acceptance=json.loads((Path(a.benchmark_root)/'acceptance.json').read_text());growth=json.loads((Path(a.benchmark_root)/'automatic_growth_acceptance.json').read_text());preflight=json.loads(Path(a.preflight).read_text())
    benchmark_manifest=json.loads(Path(a.benchmark_manifest).read_text());preflight_manifest=json.loads(Path(a.preflight_manifest).read_text())
    if acceptance.get('overall_status')!='PASSED' or growth.get('status')!='PASSED' or preflight.get('status')!='PASSED':raise RuntimeError('Physical growth/energy or coupled-driver gate is not PASSED')
    if benchmark_manifest.get('exit_code')!=0 or benchmark_manifest.get('diagnostic_replay') is not False:raise RuntimeError('Formal benchmark did not exit successfully or was a diagnostic replay')
    if subprocess.check_output(['git','status','--porcelain'],cwd=source,text=True).strip():raise RuntimeError('Pin a clean canonical source before the production history')
    values=dict(line.split('=',1) for line in cfg.read_text().splitlines() if line.strip() and not line.lstrip().startswith('#'))
    material=Path(values['dpc-table'].strip());identity={'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=source,text=True).strip(),'binary_sha256':sha(binary),'config_sha256':sha(cfg),'DPC_table_sha256':sha(material),'core_cache_sha256':sha(cache),'crack_bank_sha256':sha(bank),'checkpoint_schema':4}
    # A PASSED file from an earlier executable/input cannot authorize a new run.
    for name,path in [('binary_sha256',binary),('config_sha256',cfg),('DPC_table_sha256',material),('core_cache_sha256',cache),('crack_bank_sha256',bank)]:
        if preflight_manifest.get(name)!=sha(path):raise RuntimeError('Coupled preflight identity mismatch: '+name)
    if preflight_manifest.get('exit_code')!=0:raise RuntimeError('Coupled preflight did not exit successfully')
    for receipt in [benchmark_manifest,preflight_manifest]:
        hashes=receipt.get('verified_source_files_sha256',{})
        if not hashes:raise RuntimeError('Verification receipt has no source identity')
        for relative,digest in hashes.items():
            if sha(source/relative)!=digest:raise RuntimeError('Verified source changed after acceptance: '+relative)
    benchmark_binary=Path(a.benchmark_binary)
    if benchmark_manifest.get('binary_sha256')!=sha(benchmark_binary):raise RuntimeError('Accepted benchmark executable has changed')
    fingerprint=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest();owner=process_identity(os.getpid())
    manifest={'created_UTC':utc(),'source':str(source),'binary':str(binary),'input_config':str(cfg),'identity':identity,'fingerprint':fingerprint,'owner':owner,'same_history_checkpoint_handoff':True,'material_parameters_changed':False,'capsule_kinematics':'NOSB_PD_shell_plus_SNSNI_solid_SPH_core','benchmark_root':str(Path(a.benchmark_root).resolve()),'preflight':str(Path(a.preflight).resolve()),'segments':[]}
    put(root/'production_manifest.json',manifest);segments=[];checkpoint=None
    for segment in range(2):
        if sha(binary)!=identity['binary_sha256'] or sha(cfg)!=identity['config_sha256'] or sha(material)!=identity['DPC_table_sha256'] or sha(bank)!=identity['crack_bank_sha256'] or sha(cache)!=identity['core_cache_sha256']:raise RuntimeError('Pinned source/input/binary changed before continuation')
        output=root/('mechanics_segment_'+str(segment));runtime=root/('runtime_segment_'+str(segment));runtime.mkdir(exist_ok=True)
        if output.exists() and any(output.iterdir()):raise RuntimeError('Production output segment already contains results')
        cmd=[str(binary),'--config='+str(cfg),'--snsni-cells='+str(cache),'--fracture-facets='+str(bank),'--output='+str(output),'--checkpoint-fingerprint='+fingerprint]
        if checkpoint is None:cmd.append('--checkpoint-stop-after-equilibrium')
        else:cmd.append('--resume-state='+str(checkpoint))
        log_path=root/('native_segment_'+str(segment)+'.log')
        with log_path.open('w') as log:
            process=subprocess.Popen(cmd,cwd=runtime,stdout=log,stderr=subprocess.STDOUT)
            try:native=process_identity(process.pid)
            except OSError:native={'pid':process.pid,'argv':cmd}
            entry={'segment':segment,'output':str(output),'log':str(log_path),'native':native,'started_UTC':utc(),'argv':cmd};manifest['segments'].append(entry);segments.append(output);put(root/'production_manifest.json',manifest)
            print('STARTED',segment,process.pid,flush=True)
            try:
                while process.poll() is None:
                    state={'recorded_UTC':utc(),'status':'RUNNING','owner':owner,'native':native,'segment':segment,'output':str(output),'latest_progress':row(output/'progress.csv'),'latest_equilibrium':row(output/'equilibria.csv'),'source_commit':identity['source_commit'],'fingerprint':fingerprint,'capsule_running':True}
                    put(root/'owner_progress.json',state);index(root,segments)
                    if a.status_file:
                        target=Path(a.status_file);review=json.loads(target.read_text()) if target.exists() else {};review.update(status='PHYSICAL_CAPSULE_RUNNING',physical_capsule_running=True,capsule_restart_allowed=True,source_commit=identity['source_commit'],source_worktree_clean=True,owner_progress=str(root/'owner_progress.json'),production_manifest=str(root/'production_manifest.json'));put(target,review)
                    time.sleep(30)
            except BaseException:
                if process.poll() is None:
                    process.terminate()
                    try:process.wait(timeout=10)
                    except subprocess.TimeoutExpired:process.kill();process.wait()
                raise
            code=process.wait();entry.update(exit_code=code,ended_UTC=utc());put(root/'production_manifest.json',manifest);index(root,segments)
        if segment==0 and code==86:
            latest=json.loads((output/'checkpoint.latest.json').read_text());checkpoint=Path(latest['file']);accepted=row(output/'equilibria.csv')
            if latest['version']!=4 or latest['fingerprint']!=fingerprint or not accepted:raise RuntimeError('First equilibrium did not publish its compatible complete physical checkpoint')
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
    p.add_argument('--status-file');a=p.parse_args()
    try:run(a)
    except BaseException as error:
        root=Path(a.output_root);root.mkdir(parents=True,exist_ok=True);put(root/'owner_progress.json',{'recorded_UTC':utc(),'status':'FAILED','capsule_running':False,'error':str(error),'traceback':traceback.format_exc()})
        if a.status_file:
            target=Path(a.status_file);review=json.loads(target.read_text()) if target.exists() else {};review.update(status='PHYSICAL_CAPSULE_STOPPED',physical_capsule_running=False,capsule_restart_allowed=False,owner_progress=str(root/'owner_progress.json'));put(target,review)
        raise
