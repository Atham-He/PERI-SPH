"""Independent byte/role reader for prospective source provenance capsules.

Does not import the producer, any PD/gel arithmetic, or historical reviewers.
"""
from pathlib import Path, PurePosixPath
import hashlib
import json
import os
import stat


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest_file(path):
    path = Path(path)
    with path.open('rb') as f:
        return dict(path=str(path), size=os.fstat(f.fileno()).st_size,
                    sha256=hashlib.file_digest(f, 'sha256').hexdigest())


def checked(rec, regular=False):
    p = Path(rec['path'])
    if regular:
        require(not any(q.is_symlink() for q in [p, *p.parents]), 'independent symlink rejection')
        info = p.lstat()
        require(stat.S_ISREG(info.st_mode) and info.st_nlink == 1, 'independent regular private file')
    before = p.stat()
    data = p.read_bytes()
    after = p.stat()
    require((before.st_ino, before.st_size, before.st_mtime_ns) ==
            (after.st_ino, after.st_size, after.st_mtime_ns), 'independent stable read')
    require(set(rec) == {'path', 'size', 'sha256'}, 'independent record schema')
    require(len(data) == rec['size'] and hashlib.sha256(data).hexdigest() == rec['sha256'], 'independent byte binding')
    return data


def audit(root, contract_record, manifest_record):
    root = Path(root)
    require(manifest_record['path'] == str(root / 'manifest.json'), 'independent manifest path')
    c = json.loads(checked(contract_record, True))
    m = json.loads(checked(manifest_record, True))
    require(set(c) == {'schema', 'policy', 'source_commit', 'repository', 'catalog',
                     'historical_source_origins', 'required_actual_runtime_records', 'runtime_symlinks'}, 'independent contract schema')
    require(type(c['schema']) is int and c['schema'] == 1, 'independent contract version')
    require(c['policy'] == 'exact_historical_git_source_and_strict_actual_artifacts_v1', 'independent policy')
    catalog = json.loads(checked(c['catalog'], True))
    require(c['source_commit'] == catalog['source_commit'] and len(c['source_commit']) == 40 and
            set(c['source_commit']) <= set('0123456789abcdef'), 'independent source commit')
    pins = {}
    for rec in catalog['artifacts']:
        path = rec['path']
        require(str(PurePosixPath(path)) == path and path.startswith('/') and '..' not in PurePosixPath(path).parts, 'independent absolute origin')
        require(set(rec) == {'path', 'size', 'sha256'} and type(rec['size']) is int and rec['size'] >= 0 and
                len(rec['sha256']) == 64 and set(rec['sha256']) <= set('0123456789abcdef'), 'independent catalog pin')
        require(path not in pins or pins[path] == rec, 'independent conflicting duplicate')
        pins[path] = rec
    sources = c['historical_source_origins']
    require(sources == sorted(set(sources)) and bool(sources), 'independent source partition')
    require(str(PurePosixPath(c['repository'])) == c['repository'] and c['repository'].startswith('/') and
            '..' not in PurePosixPath(c['repository']).parts, 'independent repository origin')
    for path in sources:
        require(path in pins and path.startswith(c['repository'] + '/tools/') and path.endswith('.py'), 'independent source boundary')
    runtime = {p: r for p, r in pins.items() if p not in set(sources)}
    required = c['required_actual_runtime_records']
    require(bool(required) and len({r['path'] for r in required}) == len(required), 'independent required runtime set')
    for rec in required:
        require(runtime.get(rec['path']) == rec, 'independent runtime cannot be historical')
    require(set(c['runtime_symlinks']).issubset(runtime), 'independent runtime topology set')
    for path, rec in runtime.items():
        p = Path(path)
        expected_topology = c['runtime_symlinks'].get(path)
        actual_topology = dict(link_text=os.readlink(p), resolved_target=str(p.resolve(strict=True))) if p.is_symlink() else None
        require(actual_topology == expected_topology, 'independent runtime topology binding')
        before = p.stat()
        require(digest_file(p) == rec, 'independent actual runtime byte binding')
        after = p.stat()
        require((before.st_ino, before.st_size, before.st_mtime_ns) ==
                (after.st_ino, after.st_size, after.st_mtime_ns), 'independent stable runtime hash')
    require(set(m) == {'schema', 'policy', 'contract', 'source_entries', 'actual_records',
                      'unique_original_record_count', 'numerical_or_restart_qualification'}, 'independent manifest schema')
    require(type(m['schema']) is int and m['schema'] == 1 and m['policy'] == c['policy'], 'independent manifest identity')
    require(m['contract'] == contract_record and m['numerical_or_restart_qualification'] is False, 'independent scope and contract')
    require(type(m['unique_original_record_count']) is int and m['unique_original_record_count'] == len(pins), 'independent catalog completeness')
    require(m['actual_records'] == [runtime[k] for k in sorted(runtime)], 'independent actual artifact completeness')
    require([e.get('origin', {}).get('path') for e in m['source_entries']] == sources, 'independent source completeness')
    require({p.name for p in root.iterdir()} == {'manifest.json', 'objects'}, 'independent root inventory')
    expected_objects = set()
    for entry in m['source_entries']:
        require(set(entry) == {'role', 'origin', 'git_blob_oid', 'object'}, 'independent source entry schema')
        require(entry['role'] == 'historical_source' and entry['origin'] == pins[entry['origin']['path']], 'independent origin role')
        rec = entry['origin']
        relative = 'objects/' + rec['sha256'] + '.blob'
        require(entry['object'] == relative, 'independent object name')
        p = root / relative
        data = checked(dict(path=str(p), size=rec['size'], sha256=rec['sha256']), True)
        oid = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
        require(entry['git_blob_oid'] == oid, 'independent Git blob')
        expected_objects.add(relative)
    require({str(p.relative_to(root)) for p in (root / 'objects').iterdir()} == expected_objects, 'independent object inventory')
    return dict(status='independently_qualified_exact_provenance_capsule_only',
                source_count=len(sources), actual_count=len(runtime), unique_count=len(pins),
                numerical_or_restart_qualification=False)
