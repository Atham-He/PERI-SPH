"""One prospectively frozen synthetic-only crash-localization attempt."""
from pathlib import Path
import json,shutil,subprocess,sys,traceback
from peri_elastic_normal import read,record,verify,save,reserves,stamp,identity
from peri_actual_checkpoint_campaign import Execution,verified_execution_component


def main(root):
    repo=Path(__file__).resolve().parents[1];base=root.parent
    component=base/'engineering-s4-durable-coordinator-v1'
    assert base==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists()
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    reserves(base);execution=verified_execution_component(base);root.mkdir(mode=0o700)
    try:
        frozen=root/'frozen';shutil.copytree(component/'frozen',frozen)
        artifacts=[record(p) for p in sorted(frozen.rglob('*')) if p.is_file()]
        for name in ['peri_synthetic_json_heap_diagnostic.py','peri_synthetic_json_heap_preflight.py']:
            source=repo/'tools'/name;target=frozen/name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
        artifacts+=read(component/'plan.json')['artifacts']
        I=Path('/mnt/e/PERI-SPH/PER_sph_casule/framework_work/2026-11-05_S4_convergence_stage/incidents/20261105T113813Z-actual-index-sigsegv')
        artifacts += [record(I/'incident_receipt.json'),record(I/'core_static_inventory.json')]
        plan=dict(prepared_at=stamp(),source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
            artifacts=artifacts,execution_component=execution,synthetic_input=dict(unknowns=447,nodes=33,trace_records=4,repetitions=8),
            CPU_seconds=180,wall_seconds=240,AS_GiB=2,threads=1,
            scope='Synthetic originalpack/canonical large nested JSON temporary allocation and release only,all phases fsynced/faulthandler.Not a replay/minimization of the old physical raw.',
            physical_solver_calls=0,material_calculations=0,old_raw_reads=0,physical_or_actual_index_dispatch_admitted=False,
            failure_policy='Retain exact attempt/stage/stack/returncode;stop dependencies.No same-root retry or old state/raw rerun.')
        verify(artifacts);save(root/'plan.json',plan);save(root/'state.json',dict(status='new_synthetic_heap_diagnostic_attempted',identity=identity(),started_at=stamp()))
        Execution(root,plan).run('diagnostic',[sys.executable,'-B',str(frozen/'peri_synthetic_json_heap_diagnostic.py'),str(root)],wall=240,cpu=180,gib=2)
        q=read(root/'diagnostic_result.json');assert q['repetitions']==8 and q['precision']==95 and not q['C_decimal_loaded']
        assert q['physical_solver_calls']==q['material_calculations']==q['old_raw_reads']==0
        verify(artifacts);save(root/'independent_diagnostic.json',dict(status=q['status'],plan=record(root/'plan.json'),
            attempt=record(root/'diagnostic.attempt.json'),result=record(root/'diagnostic_result.json'),
            stage_witness=record(root/'stage_witness.jsonl'),faulthandler=record(root/'faulthandler.txt'),root_cause_resolved=False,
            originalfailedindex_retried=False,physical_solver_calls=0,material_calculations=0,old_raw_reads=0))
        save(root/'terminal.json',dict(return_code=0,finished_at=stamp(),receipt=record(root/'independent_diagnostic.json')))
        print(json.dumps(q))
    except BaseException as error:
        save(root/'failure.json',dict(recorded_at=stamp(),failure=repr(error),traceback=traceback.format_exc()));raise


if __name__=='__main__':main(Path(sys.argv[1]))
