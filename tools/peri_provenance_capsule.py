"""Prospective historical-source capsules; actual runtime paths remain strict.

No equation execution, implicit source fallback, checkout mutation, or old-plan
rewrite. A capsule qualifies provenance packaging only, never an old failed run.
"""
from pathlib import Path, PurePosixPath
import hashlib
import json
import os
import re
import stat
import subprocess

POLICY = 'exact_historical_git_source_and_strict_actual_artifacts_v1'


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True).encode()


def fail(message):
    raise ValueError(message)


def record(path):
    path = Path(path)
    before = path.stat()
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    after = path.stat()
    if (before.st_ino, before.st_size, before.st_mtime_ns) != (after.st_ino, after.st_size, after.st_mtime_ns):
        fail('file changed during metadata read')
    return dict(path=str(path), size=after.st_size, sha256=digest)


def bound(rec):
    if set(rec) != {'path', 'size', 'sha256'}:
        fail('exact record schema')
    path = Path(rec['path'])
    regular_unlinked(path)
    with path.open('rb') as stream:
        before = os.fstat(stream.fileno())
        data = stream.read()
        after = os.fstat(stream.fileno())
    current = path.stat()
    if (before.st_ino, before.st_size, before.st_mtime_ns) != (after.st_ino, after.st_size, after.st_mtime_ns):
        fail('bound file changed during read')
    if (current.st_ino, current.st_size, current.st_mtime_ns) != (after.st_ino, after.st_size, after.st_mtime_ns):
        fail('bound file replaced during read')
    if len(data) != rec['size'] or hashlib.sha256(data).hexdigest() != rec['sha256']:
        fail('exact record binding')
    return data


def absolute(path):
    p = PurePosixPath(path)
    if not p.is_absolute() or str(p) != path or '..' in p.parts:
        fail('canonical absolute origin required')
    return path


def regular_unlinked(path):
    path = Path(path)
    for parent in [path, *path.parents]:
        if parent.is_symlink():
            fail('capsule symlink forbidden')
    if not stat.S_ISREG(path.lstat().st_mode) or path.lstat().st_nlink != 1:
        fail('capsule regular file required')


def topology(path):
    path = Path(path)
    return dict(link_text=os.readlink(path), resolved_target=str(path.resolve(strict=True))) if path.is_symlink() else None


def inputs(contract_rec):
    contract = json.loads(bound(contract_rec))
    keys = {'schema', 'policy', 'source_commit', 'repository', 'catalog',
            'historical_source_origins', 'required_actual_runtime_records', 'runtime_symlinks'}
    if set(contract) != keys or type(contract['schema']) is not int or contract['schema'] != 1 or contract['policy'] != POLICY:
        fail('contract schema or policy')
    repo = absolute(contract['repository'])
    if not re.fullmatch('[0-9a-f]{40}', contract['source_commit']):
        fail('exact source commit required')
    catalog = json.loads(bound(contract['catalog']))
    if catalog['source_commit'] != contract['source_commit']:
        fail('catalog source commit binding')
    records = {}
    for rec in catalog['artifacts']:
        if set(rec) != {'path', 'size', 'sha256'}:
            fail('catalog record schema')
        path = absolute(rec['path'])
        if type(rec['size']) is not int or rec['size'] < 0 or not re.fullmatch('[0-9a-f]{64}', rec['sha256']):
            fail('catalog size or SHA')
        if path in records and records[path] != rec:
            fail('conflicting duplicate origin')
        records[path] = rec
    source = contract['historical_source_origins']
    if source != sorted(set(source)) or not source:
        fail('unique sorted source origin list required')
    for path in source:
        absolute(path)
        if path not in records or not path.startswith(repo + '/tools/') or not path.endswith('.py'):
            fail('historical source policy boundary')
    source_set = set(source)
    runtime = {k: v for k, v in records.items() if k not in source_set}
    required = contract['required_actual_runtime_records']
    if not required:
        fail('explicit actual runtime binding required')
    for rec in required:
        if rec['path'] not in runtime or runtime[rec['path']] != rec:
            fail('actual runtime cannot be archived as source')
    if len({r['path'] for r in required}) != len(required):
        fail('duplicate actual runtime binding')
    links = contract['runtime_symlinks']
    if not set(links).issubset(runtime):
        fail('runtime topology contains non-runtime origin')
    for path, rec in runtime.items():
        if topology(path) != links.get(path):
            fail('actual runtime symlink topology changed')
        if record(path) != rec:
            fail('actual runtime or artifact changed')
    return contract, records, source, runtime


def object_bytes(root, digest):
    if not re.fullmatch('[0-9a-f]{64}', digest):
        fail('object digest')
    p = root / 'objects' / (digest + '.blob')
    regular_unlinked(p)
    return p.read_bytes()


def prepare(root, contract_rec):
    root = Path(absolute(str(root)))
    if any(p.is_symlink() for p in root.parents):
        fail('capsule parent symlink forbidden')
    if root.exists():
        fail('new capsule directory required; no retry')
    contract, records, source, runtime = inputs(contract_rec)
    exported = []
    # Export only the exact declared Git version. Never read an altered checkout
    # as a substitute, search for a convenient version, or modify HEAD.
    for origin in source:
        relative = str(Path(origin).relative_to(contract['repository']))
        spec = contract['source_commit'] + ':' + relative
        data = subprocess.check_output(['git', '-C', contract['repository'], 'show', spec])
        rec = records[origin]
        if len(data) != rec['size'] or hashlib.sha256(data).hexdigest() != rec['sha256']:
            fail('declared historical Git bytes mismatch')
        git_oid = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
        exported.append((rec, data, git_oid))
    root.mkdir(parents=False)
    (root / 'objects').mkdir()
    entries = []
    for rec, data, git_oid in exported:
        dest = root / 'objects' / (rec['sha256'] + '.blob')
        if not dest.exists():
            with dest.open('xb') as stream:
                stream.write(data)
        entries.append(dict(role='historical_source', origin=rec, git_blob_oid=git_oid,
                            object='objects/' + rec['sha256'] + '.blob'))
    manifest = dict(schema=1, policy=POLICY, contract=contract_rec,
                    source_entries=entries, actual_records=[runtime[k] for k in sorted(runtime)],
                    unique_original_record_count=len(records), numerical_or_restart_qualification=False)
    path = root / 'manifest.json'
    with path.open('xb') as stream:
        stream.write(canonical(manifest))
    verify(root, contract_rec, record(path))
    return record(path)


def verify(root, contract_rec, manifest_rec):
    root = Path(root)
    regular_unlinked(root / 'manifest.json')
    if {p.name for p in root.iterdir()} != {'objects', 'manifest.json'}:
        fail('unbound capsule root entry')
    if manifest_rec['path'] != str(root / 'manifest.json'):
        fail('manifest path binding')
    manifest = json.loads(bound(manifest_rec))
    contract, records, source, runtime = inputs(contract_rec)
    expected_keys = {'schema', 'policy', 'contract', 'source_entries', 'actual_records',
                     'unique_original_record_count', 'numerical_or_restart_qualification'}
    if set(manifest) != expected_keys or type(manifest['schema']) is not int or manifest['schema'] != 1 or manifest['policy'] != POLICY:
        fail('manifest schema or policy')
    if manifest['contract'] != contract_rec or manifest['numerical_or_restart_qualification'] is not False:
        fail('trusted contract or qualification scope')
    if manifest['unique_original_record_count'] != len(records):
        fail('complete original catalog required')
    if manifest['actual_records'] != [runtime[k] for k in sorted(runtime)]:
        fail('actual artifact roles or completeness changed')
    entries = manifest['source_entries']
    if [e.get('origin', {}).get('path') for e in entries] != source:
        fail('historical source completeness or ordering')
    for entry in entries:
        if set(entry) != {'role', 'origin', 'git_blob_oid', 'object'}:
            fail('source entry schema')
        expected = records[entry['origin']['path']]
        if entry['role'] != 'historical_source' or entry['origin'] != expected:
            fail('exact historical role or origin binding')
        digest = expected['sha256']
        if entry['object'] != 'objects/' + digest + '.blob':
            fail('object path binding')
        data = object_bytes(root, digest)
        if len(data) != expected['size'] or hashlib.sha256(data).hexdigest() != digest:
            fail('historical capsule bytes changed')
        oid = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
        if entry['git_blob_oid'] != oid:
            fail('Git blob identity')
    allowed = {e['object'] for e in entries}
    actual = {str(p.relative_to(root)) for p in (root / 'objects').iterdir()}
    if actual != allowed:
        fail('unbound capsule object')
    return dict(status='exact_provenance_capsule_only_qualified',
                historical_source_count=len(source), actual_artifact_count=len(runtime),
                unique_original_record_count=len(records), source_commit=contract['source_commit'],
                actual_runtime_paths_preserved=True, numerical_or_restart_qualification=False)
