"""New finite metadata controls, including correct-digest semantic mutations.

Fixtures contain a tiny local Git repository, not an upstream clone or a PD run.
"""
from pathlib import Path
import copy
import hashlib
import json
import os
import subprocess
import sys
import time

import peri_provenance_capsule as producer
import peri_provenance_capsule_audit as independent


def dump(path, value):
    Path(path).write_text(json.dumps(value, sort_keys=True, separators=(',', ':')))


def fixture(work, label, link_runtime=False, checkout_runtime=False):
    case = work / label
    case.mkdir()
    repo = case / 'repo'
    (repo / 'tools').mkdir(parents=True)
    src = repo / 'tools/source.py'
    src.write_bytes(b'value = 17\n')
    run = (repo / 'tools/runtime.py') if checkout_runtime else (case / 'runtime.bin')
    run.write_bytes(b'ACTUAL-RUNTIME')
    env = dict(os.environ, GIT_AUTHOR_DATE='2026-11-05T00:00:00Z',
               GIT_COMMITTER_DATE='2026-11-05T00:00:00Z')
    for args in [['init', '--quiet'], ['add', '.'],
                 ['-c', 'user.name=CapsuleFixture', '-c', 'user.email=fixture@invalid', 'commit', '--quiet', '-m', 'exact fixture']]:
        subprocess.run(['git', '-C', str(repo), *args], env=env, check=True, capture_output=True)
    commit = subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD']).decode().strip()
    if link_runtime:
        target = case / 'runtime_alias'
        target.symlink_to(run.name)
        run = target
    catalog = case / 'catalog.json'
    rows = [producer.record(src), producer.record(run)]
    dump(catalog, dict(source_commit=commit, artifacts=[*rows, rows[0]]))
    contract = dict(schema=1, policy=producer.POLICY, repository=str(repo),
                    source_commit=commit, catalog=producer.record(catalog),
                    historical_source_origins=[str(src)],
                    required_actual_runtime_records=[producer.record(run)],
                    runtime_symlinks={str(run): producer.topology(run)} if run.is_symlink() else {})
    path = case / 'contract.json'
    dump(path, contract)
    return dict(case=case, repo=repo, src=src, runtime=run, catalog=catalog,
                contract_path=path, contract=contract, contract_record=producer.record(path),
                capsule=case / 'capsule')


def refresh_contract(f):
    dump(f['contract_path'], f['contract'])
    f['contract_record'] = producer.record(f['contract_path'])


def package(f):
    f['manifest_record'] = producer.prepare(f['capsule'], f['contract_record'])


def mutate_manifest(f, transform):
    p = f['capsule'] / 'manifest.json'
    m = json.loads(p.read_text())
    transform(m)
    dump(p, m)
    f['manifest_record'] = producer.record(p)


def accepted(f):
    a = producer.verify(f['capsule'], f['contract_record'], f['manifest_record'])
    b = independent.audit(f['capsule'], f['contract_record'], f['manifest_record'])
    assert a['unique_original_record_count'] == b['unique_count'] == 2
    assert a['historical_source_count'] == b['source_count'] == 1
    assert a['actual_artifact_count'] == b['actual_count'] == 1
    assert not a['numerical_or_restart_qualification'] and not b['numerical_or_restart_qualification']


def rejected(f):
    failures = []
    for fn in [producer.verify, independent.audit]:
        try:
            fn(f['capsule'], f['contract_record'], f['manifest_record'])
        except (ValueError, OSError, KeyError, TypeError) as error:
            failures.append(type(error).__name__ + ':' + str(error))
        else:
            raise AssertionError('semantic mutation accepted by ' + fn.__module__)
    return failures


def prepare_rejected(f):
    try:
        package(f)
    except (ValueError, OSError, KeyError, TypeError, subprocess.CalledProcessError) as error:
        assert not f['capsule'].exists(), 'admission rejection must precede output'
        return [type(error).__name__ + ':' + str(error)]
    raise AssertionError('invalid prospective contract accepted')


def run(work):
    work = Path(work)
    work.mkdir()
    results = []

    def valid(label, before=None, after=None, **kwargs):
        f = fixture(work, label, **kwargs)
        if before: before(f)
        package(f)
        if after: after(f)
        accepted(f)
        results.append(dict(label=label, expected='both_independent_readers_accept', passed=True))

    def negative(label, mutation, **kwargs):
        f = fixture(work, label, **kwargs)
        package(f)
        mutation(f)
        results.append(dict(label=label, expected='both_readers_reject', passed=True, reasons=rejected(f)))

    def admission(label, mutation, **kwargs):
        f = fixture(work, label, **kwargs)
        mutation(f)
        results.append(dict(label=label, expected='reject_before_capsule_creation', passed=True, reasons=prepare_rejected(f)))

    valid('exact_complete_duplicate_deduplicated')
    valid('checkout_changed_after_seal', after=lambda f: f['src'].write_bytes(b'value = 99\n'))
    valid('exact_declared_git_bytes_despite_checkout_change', before=lambda f: f['src'].write_bytes(b'newer uncommitted source\n'))

    def corrupt_object(f):
        p = next((f['capsule'] / 'objects').iterdir())
        p.write_bytes(p.read_bytes() + b'!')
    negative('historical_object_corruption', corrupt_object)
    negative('correct_digest_missing_source', lambda f: mutate_manifest(f, lambda m: m.update(source_entries=[])))
    negative('correct_digest_extra_source', lambda f: mutate_manifest(f, lambda m: m['source_entries'].append(copy.deepcopy(m['source_entries'][0]))))
    negative('correct_digest_origin_changed', lambda f: mutate_manifest(f, lambda m: m['source_entries'][0]['origin'].update(path=str(f['repo'] / 'tools/other.py'))))
    negative('correct_digest_role_changed', lambda f: mutate_manifest(f, lambda m: m['source_entries'][0].update(role='actual_runtime')))
    negative('correct_digest_runtime_omitted', lambda f: mutate_manifest(f, lambda m: m.update(actual_records=[])))
    negative('correct_digest_object_escape', lambda f: mutate_manifest(f, lambda m: m['source_entries'][0].update(object='../outside.bin')))
    negative('correct_digest_wrong_git_oid', lambda f: mutate_manifest(f, lambda m: m['source_entries'][0].update(git_blob_oid='0' * 40)))
    negative('correct_digest_false_physics_qualification', lambda f: mutate_manifest(f, lambda m: m.update(numerical_or_restart_qualification=True)))

    def changed_runtime(f):
        f['runtime'].write_bytes(b'CHANGED-RUNTIME')
        mutate_manifest(f, lambda m: m.update(actual_records=[producer.record(f['runtime'])]))
    negative('runtime_changed_even_with_updated_manifest_digest', changed_runtime)

    def runtime_alias(f):
        alternate = f['case'] / 'same_runtime_bytes'
        alternate.write_bytes(f['runtime'].read_bytes())
        f['runtime'].unlink()
        f['runtime'].symlink_to(alternate.name)
    negative('same_bytes_actual_runtime_symlink_retarget', runtime_alias, link_runtime=True)

    def object_alias(f):
        p = next((f['capsule'] / 'objects').iterdir())
        external = f['case'] / 'same_object.bin'
        external.write_bytes(p.read_bytes())
        p.unlink()
        p.symlink_to(external)
    negative('same_bytes_capsule_symlink', object_alias)

    def object_hardlink(f):
        p = next((f['capsule'] / 'objects').iterdir())
        os.link(p, f['case'] / 'external_hardlink')
    negative('same_bytes_capsule_hardlink', object_hardlink)
    negative('extra_unbound_object', lambda f: (f['capsule'] / 'objects/unbound.blob').write_bytes(b'foreign'))
    negative('extra_unbound_root_file', lambda f: (f['capsule'] / 'unbound.bin').write_bytes(b'foreign'))

    def runtime_relabel(f):
        f['contract']['historical_source_origins'] = sorted([str(f['src']), str(f['runtime'])])
        refresh_contract(f)
    admission('required_actual_runtime_cannot_be_source', runtime_relabel, checkout_runtime=True)

    def conflicting(f):
        cat = json.loads(f['catalog'].read_text())
        altered = copy.deepcopy(cat['artifacts'][0])
        altered['sha256'] = '0' * 64
        cat['artifacts'].append(altered)
        dump(f['catalog'], cat)
        f['contract']['catalog'] = producer.record(f['catalog'])
        refresh_contract(f)
    admission('correct_digest_conflicting_duplicate_origin', conflicting)

    def wrong_commit(f):
        f['contract']['source_commit'] = '0' * 40
        refresh_contract(f)
    admission('correct_digest_catalog_commit_mismatch', wrong_commit)

    def wrong_git_content(f):
        cat = json.loads(f['catalog'].read_text())
        for rec in cat['artifacts']:
            if rec['path'] == str(f['src']): rec['sha256'] = '0' * 64
        dump(f['catalog'], cat)
        f['contract']['catalog'] = producer.record(f['catalog'])
        refresh_contract(f)
    admission('expected_content_not_in_declared_git_version', wrong_git_content)

    def parent_link(f):
        target = f['case'] / 'outside'
        target.mkdir()
        link = f['case'] / 'linked'
        link.symlink_to(target, target_is_directory=True)
        f['capsule'] = link / 'capsule'
    admission('capsule_parent_symlink_before_writes', parent_link)

    def changed_contract(f):
        f['contract']['historical_source_origins'] = []
        # Deliberately retain the externally trusted original contract record.
        dump(f['contract_path'], f['contract'])
    negative('externally_pinned_contract_change', changed_contract)
    assert len(results) == 24
    return dict(status='new_provenance_capsule_24_metadata_controls_passed', controls=results,
                producer_and_independent_readers=True, solver_calls=0, material_calculations=0,
                old_equation_reviews_repeated=False, actual_old_state_resumed=False,
                actual_full80_qualified=False, physical_restart_qualified=False, S4_complete=False)


if __name__ == '__main__':
    started = time.process_time()
    result = run(Path(sys.argv[1]))
    result['CPU_seconds'] = time.process_time() - started
    result['C_decimal_loaded'] = '_decimal' in sys.modules
    print(json.dumps(result, indent=2))
