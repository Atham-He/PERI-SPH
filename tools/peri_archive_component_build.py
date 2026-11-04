"""Isolated component build; freeze actual includes/libs, never rebuild production."""
import shlex
import shutil
import subprocess
from pathlib import Path
from peri_run_body_archive_engineering import run,save
from peri_bounded_body_archive_audit import record


def prepare(root,repo,files,cpp,stem):
    root=Path(root);repo=Path(repo);frozen=root/'frozen';snapshots=[]
    for name in files:
        target=frozen/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(repo/name,target)
        snapshots += [record(repo/name),record(target)]
    before={str(p.resolve()):record(p) for p in (repo/'src').rglob('*') if p.suffix in ['.h','.hpp','.cpp']}
    builddir=repo/'build/peri-latest-20261103'
    old_object='extra_3d_PD_framework_tests/CMakeFiles/extra_3d_PD_bounded_gel_friction_probe.dir/bounded_gel_friction_probe.cpp.o'
    old_binary='extra_3d_PD_framework_tests/bin/extra_3d_PD_bounded_gel_friction_probe'
    compiler=shlex.split(subprocess.check_output(['ninja','-t','commands','-s',old_object],cwd=builddir,text=True))
    for flag,value in [('-o',root/(stem+'.o')),('-MF',root/(stem+'.d')),('-MT',root/(stem+'.o')),('-c',frozen/cpp)]:
        compiler[compiler.index(flag)+1]=str(value)
    compiler += ['-Wall','-Wextra','-Wno-misleading-indentation','-Wno-unused-parameter','-Werror']
    text=subprocess.check_output(['ninja','-t','commands','-s',old_binary],cwd=builddir,text=True).strip()
    assert text.startswith(': && ') and text.endswith(' && :')
    linker=shlex.split(text[5:-5]);assert linker.count(old_object)==1
    linker[linker.index(old_object)]=str(root/(stem+'.o'));linker[linker.index('-o')+1]=str(root/(stem+'_checks'))
    return dict(root=root,repo=repo,before=before,source_snapshots=snapshots,stem=stem,builddir=builddir,
                compile_command=compiler,link_command=linker)


def compile_and_freeze(prepared):
    p=prepared;root=p['root'];builddir=p['builddir'];compiler=p['compile_command'];linker=p['link_command']
    compilation=run(root,'compile',compiler,builddir,180,240,2)
    link=run(root,'link',linker,builddir,30,60,2)
    text=(root/(p['stem']+'.d')).read_text().replace('\\\n',' ')
    dependencies=[Path(x).resolve() for x in shlex.split(text.split(':',1)[1])]
    libraries=[builddir/x for x in linker if x.endswith(('.a','.so')) or '.so.' in x]
    elf=root/(p['stem']+'_checks');dynamic=subprocess.check_output(['ldd',str(elf)],text=True)
    with (root/'dynamic_libraries.txt').open('x') as stream:stream.write(dynamic)
    for line in dynamic.splitlines():
        assert 'not found' not in line
        libraries += [Path(x) for x in line.strip().split() if x.startswith('/')]
    inputs=sorted(set(dependencies+libraries+[Path(compiler[0]),Path(linker[0])]),key=str);support=[]
    for i,path in enumerate(inputs):
        original=record(path)
        if str(path) in p['before']:assert original==p['before'][str(path)], ('compile input changed',path)
        target=root/'support'/str(i)/path.name;target.parent.mkdir(parents=True,exist_ok=False);shutil.copy2(path,target)
        support += [original,record(target)]
    save(root/'source_manifest.json',dict(source_files=p['source_snapshots'],actual_compile_sources=[record(x) for x in dependencies]))
    save(root/'build_runtime_manifest.json',dict(compiler_command=compiler,linker_command=linker,support=support,ELF=record(elf)))
    for item in support+p['source_snapshots']:assert record(item['path'])==item
    return dict(compile=compilation,link=link,ELF=record(elf),support=support,source_snapshots=p['source_snapshots'],
        source_manifest=record(root/'source_manifest.json'),build_runtime_manifest=record(root/'build_runtime_manifest.json'),
        actual_compiler_dependency_count=len(dependencies),support_original_and_copy_records=len(support))
