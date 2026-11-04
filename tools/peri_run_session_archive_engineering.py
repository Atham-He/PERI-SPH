"""Freeze one new detached-contract component test; never call a PD solver."""
import argparse
import json
import os
import shlex
import shutil
import subprocess
import traceback
from pathlib import Path
from peri_run_body_archive_engineering import run, reserves, save, identity, stamp
from peri_bounded_body_archive_audit import record
from peri_bounded_session_archive_audit import make_fields, encode, assess


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',required=True);args=parser.parse_args()
    root=Path(args.root).resolve()
    assert root.parent==Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material') and not root.exists(), 'no redispatch'
    reserves(root);repo=Path(__file__).resolve().parents[1]
    assert subprocess.check_output(['git','status','--porcelain'],cwd=repo,text=True)==''
    root.mkdir(mode=0o700)
    try:
        source=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()
        body_root=root.parent/'engineering-nonzero-body-archive'
        body_receipt=json.loads((body_root/'engineering_receipt.json').read_text())
        body_fixture=body_root/'nonzero_body.bin'
        assert record(body_fixture)==body_receipt['fixture']
        original_input=root.parent/'study-v19-early-friction-cycle/actual_PD_early_biaxial_slip_cycle.input'
        anchor=[record(body_fixture),record(body_root/'engineering_receipt.json'),record(original_input)]
        files=['tools/bounded_session_archive_checks.cpp','src/shared/peridynamics/pd_bounded_session_archive.h',
            'tools/peri_bounded_session_archive_audit.py','tools/peri_run_session_archive_engineering.py',
            'tools/peri_run_body_archive_engineering.py','tools/peri_bounded_body_archive_audit.py',
            'tools/peri_canonical_archive_audit.py','docs/peri/BOUNDED_SESSION_CONTRACT_20261105.md']
        frozen=root/'frozen';snapshot=[]
        for name in files:
            target=frozen/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(repo/name,target)
            snapshot += [record(repo/name),record(target)]
        before={str(p.resolve()):record(p) for p in (repo/'src').rglob('*') if p.suffix in ['.h','.hpp','.cpp']}
        builddir=repo/'build/peri-latest-20261103'
        old_object='extra_3d_PD_framework_tests/CMakeFiles/extra_3d_PD_bounded_gel_friction_probe.dir/bounded_gel_friction_probe.cpp.o'
        old_binary='extra_3d_PD_framework_tests/bin/extra_3d_PD_bounded_gel_friction_probe'
        compiler=shlex.split(subprocess.check_output(['ninja','-t','commands','-s',old_object],cwd=builddir,text=True))
        for flag,value in [('-o',root/'session.o'),('-MF',root/'session.d'),('-MT',root/'session.o'),
                           ('-c',frozen/'tools/bounded_session_archive_checks.cpp')]:
            compiler[compiler.index(flag)+1]=str(value)
        compiler += ['-Wall','-Wextra','-Wno-misleading-indentation','-Wno-unused-parameter','-Werror']
        text=subprocess.check_output(['ninja','-t','commands','-s',old_binary],cwd=builddir,text=True).strip()
        assert text.startswith(': && ') and text.endswith(' && :')
        linker=shlex.split(text[5:-5]);assert linker.count(old_object)==1
        linker[linker.index(old_object)]=str(root/'session.o');linker[linker.index('-o')+1]=str(root/'session_checks')
        plan=dict(prepared_at=stamp(),source_commit=source,physical_advances=0,source_snapshots=snapshot,anchors=anchor,
            compile_command=compiler,link_command=linker,build_cwd=str(builddir),
            compile_budget=dict(CPU_seconds=180,wall_seconds=240,AS_GiB=2),
            native_budget=dict(CPU_seconds=30,wall_seconds=90,AS_GiB=1,threads=1),
            file_bound_MiB=256,total_bound_GiB=1,host_reserve_GiB=120,output_reserve_GiB=150,
            fixture_recipe='Preserved nonzero body+original17input; actual frozen source/ELF/build-runtime manifests; external fingerprint generated independently with struct/hashlib. Cursor fixture hashes are actual NON-STREAM files, never a physical output proof.',
            gates=['140 fields, entire canonical file bitwise equal','20 valid-hash specific semantic refusals',
                   'external trust/generation/current-rounding/FTZ four explicit refusals','baseline unchanged',
                   '96 protocol accepted time/actual-dt bit pairs equal independent Python','zero physical advance/live commit'],
            failure_policy='Preserve exact exception/returncode/artifacts and stop this state, no solver call/retry/fallback/reset/gate relaxation.',
            output_prefix_verified=False,full_restart_qualified=False,S3_complete=False,S4_complete=False,calibration_allowed=False)
        save(root/'plan.json',plan);save(root/'state.json',dict(status='single_contract_component_attempted',identity=identity(os.getpid()),started_at=stamp()))
        compilation=run(root,'compile',compiler,builddir,180,240,2)
        link=run(root,'link',linker,builddir,30,60,2)
        deptext=(root/'session.d').read_text().replace('\\\n',' ')
        dependencies=[Path(p).resolve() for p in shlex.split(deptext.split(':',1)[1])]
        libraries=[builddir/p for p in linker if p.endswith(('.a','.so')) or '.so.' in p]
        dynamic=subprocess.check_output(['ldd',str(root/'session_checks')],text=True)
        (root/'dynamic_libraries.txt').write_text(dynamic)
        for line in dynamic.splitlines():
            assert 'not found' not in line
            libraries += [Path(x) for x in line.strip().split() if x.startswith('/')]
        inputs=sorted(set(dependencies+libraries+[Path(compiler[0]),Path(linker[0])]),key=str)
        support=[]
        for i,p in enumerate(inputs):
            original=record(p)
            if str(p) in before:assert original==before[str(p)], ('compile input changed',p)
            target=root/'support'/str(i)/p.name;target.parent.mkdir(parents=True,exist_ok=False);shutil.copy2(p,target)
            support += [original,record(target)]
        save(root/'source_manifest.json',dict(source_commit=source,source_files=snapshot,
            actual_compile_sources=[record(p) for p in dependencies]))
        save(root/'build_runtime_manifest.json',dict(compiler_command=compiler,linker_command=linker,support=support,ELF=record(root/'session_checks')))
        identities=dict(SourceManifest=record(root/'source_manifest.json')['sha256'],
            BuildRuntimeManifest=record(root/'build_runtime_manifest.json')['sha256'],ELF=record(root/'session_checks')['sha256'])
        fields,trusted=make_fields(body_fixture.read_bytes(),original_input,identities)
        fixture=root/'session_candidate.bin'
        with fixture.open('xb') as stream:stream.write(encode(fields,37));stream.flush();os.fsync(stream.fileno())
        save(root/'external_trust.json',dict(contract_SHA256=trusted,fixture=record(fixture),identities=identities,
            body_fixture=record(body_fixture),original_input=record(original_input),physical_output_verified=False))
        execution=dict(recorded_at=stamp(),plan=record(root/'plan.json'),ELF=record(root/'session_checks'),
            fixture=record(fixture),external_trust=record(root/'external_trust.json'),support=support,
            source_snapshots=snapshot,anchors=anchor,source_manifest=record(root/'source_manifest.json'),
            build_runtime_manifest=record(root/'build_runtime_manifest.json'),compile=compilation,link=link)
        for item in support+snapshot+anchor:assert record(item['path'])==item
        save(root/'frozen_execution.json',execution)
        evidence=root/'evidence';evidence.mkdir()
        command=[str(root/'session_checks'),str(fixture),str(evidence),trusted]
        native=run(root,'native',command,root,30,90,1)
        for item in support+snapshot+anchor+[execution[k] for k in ['plan','ELF','fixture','external_trust','source_manifest','build_runtime_manifest']]:
            assert record(item['path'])==item
        report=json.loads((root/'native.stdout').read_text())
        if native['identity']:assert report['native_identity']==native['identity']
        else:
            assert report['native_identity']['argv']==command and report['native_identity']['ppid']==os.getpid()
            assert report['native_identity']['boot_id']==identity(os.getpid())['boot_id']
        result=assess(fields,fixture,evidence,report)
        result.update(recorded_at=stamp(),source_commit=source,plan=record(root/'plan.json'),
            frozen_execution=record(root/'frozen_execution.json'),native=native,external_trust=record(root/'external_trust.json'),
            actual_compiler_dependency_count=len(dependencies),support_original_and_copy_records=len(support),
            evidence=[record(p) for p in sorted(evidence.iterdir())])
        save(root/'engineering_receipt.json',result)
        save(root/'terminal.json',dict(status=result['status'],finished_at=stamp(),return_code=0,
            receipt=record(root/'engineering_receipt.json'),physical_advances=0,output_prefix_verified=False,full_restart_qualified=False))
        print(json.dumps(dict(status=result['status'],fields=140,refusals=24)))
    except BaseException as error:
        save(root/'failure.json',dict(recorded_at=stamp(),error=repr(error),traceback=traceback.format_exc(),
            physical_advances=0,output_prefix_verified=False,full_restart_qualified=False,S3_complete=False))
        raise


if __name__=='__main__':main()
