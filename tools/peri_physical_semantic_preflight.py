"""One-shot v2 full semantic guard preflight, no physical solver entry."""
import argparse,json,shutil,subprocess,sys,traceback
from pathlib import Path
from peri_actual_checkpoint_campaign import Execution,verified_execution_component
from peri_physical_checkpoint_negative import specification,run_controls
from peri_physical_checkpoint_support import initial_contract
from peri_elastic_normal import read,record,verify,save,reserves,identity,stamp

def main(root):
    repo=Path(__file__).resolve().parents[1];base=root.parent;native=base/'engineering-s4-physical-header-v1'
    header=base/'analysis-s4-header-v3-admission-path';binding=base/'engineering-s4-v2-reference-binding-v1'
    assert base==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists(),'unique new semantic state required'
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    q=read(header/'independent_review.json');assert q['status']=='existing_actual_v2_header_independently_qualified_by_NEW_path_admitted_reader_ZERO_steps'
    assert read(header/'terminal.json')['return_code']==0 and read(binding/'terminal.json')['return_code']==0
    build=read(native/'frozen_execution.json')['build'];spec=specification();assert spec['schema']==2 and len(spec['controls'])==22
    component=verified_execution_component(base);reserves(base);root.mkdir(mode=0o700)
    try:
        reference=root/'frozen/reference';reference.mkdir(parents=True);artifacts=[]
        for source in sorted((repo/'tools').glob('*.py')):
            target=reference/source.name;shutil.copy2(source,target);artifacts += [record(source),record(target)]
        artifacts += read(native/'plan.json')['software']+build['support']+build['source_snapshots']+[build['ELF']]
        for anchor in [header,binding]:
            t=read(anchor/'terminal.json');verify([t['receipt']]);artifacts += [record(anchor/'terminal.json'),t['receipt']]
        for name in ['source_manifest.json','build_runtime_manifest.json','original.input','initial_candidate.bin','initial.phy','contract.bin']:
            source=native/name;target=root/name;shutil.copy2(source,target);artifacts += [record(source),record(target)]
        exe=root/'frozen/physical_PD';shutil.copy2(build['ELF']['path'],exe);artifacts.append(record(exe))
        case=read(native/'plan.json')['case'];trust=q['initial_complete_bytes']['contract_sha256']
        initial=initial_contract(root,case,bytes.fromhex(trust));fields=initial['fields'];archive=(root/'initial_candidate.bin').read_bytes()
        derived=dict(checkpoints={0:dict(archive=archive,physical=initial['physical'])})
        plan=dict(prepared_at=stamp(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
                  artifacts=artifacts,executable=str(exe),trusted_contract=trust,case=case,negative_specification=spec,
                  execution_component=component,CPU_seconds=30,wall_seconds=45,AS_GiB=1,threads=1,
                  physical_solver_calls=0,physical_dispatch_allowed=False,metadata_QA_methods=4,metadata_QA_passed=True,
                  full_restart_qualified=False,S3_prior_qualification_preserved=True,S4_complete=False,calibration_allowed=False,
                  failure_policy='Stop dependent controls on any failure,preserve exact causes/raw/returncodes/no old state redispatch.')
        plan['critical_runtime']=[record(exe)]+[record(root/n) for n in ['source_manifest.json','build_runtime_manifest.json','original.input']]
        plan['critical_runtime'] += [x for x in build['support'] if '.so' in Path(x['path']).name]
        verify(artifacts);save(root/'plan.json',plan);save(root/'state.json',dict(status='new_v2_semantic_preflight_attempted',identity=identity(),started_at=stamp()))
        result=run_controls(root,plan,Execution(root,plan),derived,step=0,generation=1)
        q=read(result['path']);assert q['semantic_rejections']==22 and q['positive_binding_controls']==1 and not q['qualified_physical_nonzero_seed']
        verify(artifacts);save(root/'independent_review.json',dict(status='new_v2_complete_semantic_admission_22reject1valid_ZERO_steps_qualified',
            finished_at=stamp(),review=result,plan=record(root/'plan.json'),physical_solver_calls=0,full_nonzero_restart_qualified=False,S4_complete=False,calibration_allowed=False))
        save(root/'terminal.json',dict(status='new_v2_complete_semantic_admission_22reject1valid_ZERO_steps_qualified',return_code=0,
            finished_at=stamp(),receipt=record(root/'independent_review.json')))
        print(json.dumps(dict(status='new_v2_complete_semantic_admission_22reject1valid_ZERO_steps_qualified',physical_solver_calls=0)))
    except BaseException as error:
        save(root/'failure.json',dict(recorded_at=stamp(),failure=repr(error),traceback=traceback.format_exc(),physical_solver_calls=0));raise

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path);a=p.parse_args();main(a.root)
