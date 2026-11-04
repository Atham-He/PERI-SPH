"""Independent immutable-selection and filesystem-byte oracle, no PD solve."""
import hashlib,re,struct
from pathlib import Path
from peri_canonical_archive_audit import parse
from peri_bounded_session_archive_audit import encode

TAIL=b'synthetic_uncommitted_later_tail\n'
REFUSALS={'duplicate_owner_lock','unsafe_branch','stale_expected','wrong_explicit_generation',
    'wrong_selected_archive_digest','wrong_contract_selection','symlink_branch','existing_directory',
    'semantic_bad_candidate','changed_after_proof','injected_before_stage','injected_after_files',
    'injected_before_rename','existing_pending_pointer','symlink_selected_archive','symlink_HEAD',
    'selected_prefix_identity','branch_path_shape','publication_uncertain_block'}


def selection(data,generation):
    actual,fields=parse(data);assert actual==generation
    assert set(fields)=={'SelectionSchema','Branch','ArchiveSHA256','PhysicalPrefixSHA256','ContractSHA256'}
    assert fields['SelectionSchema']==dict(tag=2,shape=[],raw=struct.pack('<Q',1))
    branch=fields['Branch'];assert branch['tag']==1 and len(branch['shape'])==1
    assert 3<=branch['shape'][0]<=96 and len(branch['raw'])==branch['shape'][0]
    name=branch['raw'].decode('ascii');assert re.fullmatch(r'g_[A-Za-z0-9_-]+',name)
    for key in ['ArchiveSHA256','PhysicalPrefixSHA256','ContractSHA256']:
        assert fields[key]['tag']==1 and fields[key]['shape']==[32] and len(fields[key]['raw'])==32
    return name,fields


def verify_branch(root,name,generation,source,physical,pointer=None):
    root=Path(root);branch=root/'store'/name;archive=(branch/'state.bin').read_bytes()
    g,fields=parse(archive);source_g,expected=parse(source)
    assert source_g==37 and g==generation and fields==expected and archive==encode(expected,generation)
    assert (branch/'prefix.phy').read_bytes()==physical
    assert (branch/'active.phy').read_bytes()==physical+(TAIL if name=='g_37_initial' else b'')
    if pointer is not None:
        selected,packet=selection(pointer,generation);assert selected==name
        assert packet['ArchiveSHA256']['raw']==hashlib.sha256(archive).digest()
        assert packet['PhysicalPrefixSHA256']['raw']==hashlib.sha256(physical).digest()
        assert packet['ContractSHA256']['raw']==physical[16:48]
    return dict(branch=name,generation=generation,archive_bytes=len(archive),prefix_bytes=len(physical),archive_SHA256=hashlib.sha256(archive).hexdigest())


def assess_publish(root,native):
    root=Path(root);source=(root/'session_candidate.bin').read_bytes();physical=(root/'physical.phy').read_bytes()
    assert native['status']=='generation_publish_fault_controls_passed_NOT_physical_restart'
    assert native['durable_publishes']==2 and native['injected_post_publish_uncertainty'] is True
    assert native['uncertain_owner_known'] and native['further_publish_refused']
    assert native['source_generation']==37 and native['terminal_selected_generation']==38
    labels=[label for label,reason in native['refusals']]
    assert len(labels)==len(set(labels))==19 and set(labels)==REFUSALS
    assert native['nonlinear_solver_calls']==native['physical_advances']==0 and native['full_restart_qualified'] is False
    initial=(root/'selection_initial.bin').read_bytes();following=(root/'selection_next.bin').read_bytes();uncertain=(root/'selection_uncertain.bin').read_bytes()
    branches=[verify_branch(root,'g_37_initial',37,source,physical,initial),
              verify_branch(root,'g_37_next',37,source,physical,following),
              verify_branch(root,'g_38_uncertain',38,source,physical,uncertain)]
    assert (root/'store/HEAD').read_bytes()==uncertain
    for name in ['g_fault_files','g_fault_rename','g_pending']:
        assert (root/'store'/name/'state.bin').read_bytes()==source
        assert (root/'store'/name/'prefix.phy').read_bytes()==physical
        assert (root/'store'/name/'active.phy').read_bytes()==physical
    assert (root/'store/previous_g_37_next.bin').read_bytes()==initial
    assert (root/'store/previous_g_38_uncertain.bin').read_bytes()==following
    assert (root/'store/previous_g_fault_rename.bin').read_bytes()==initial
    assert selection((root/'store/pending_g_fault_rename.bin').read_bytes(),37)[0]=='g_fault_rename'
    assert (root/'store/pending_g_pending.bin').read_bytes()==b'keep_protected'
    assert (root/'protected.bin').read_bytes()==b'keep_protected'
    assert (root/'store/g_symlink').is_symlink()
    assert (root/'store/g_badarchive/state.bin').is_symlink()
    assert (root/'symlink_head_store/HEAD').is_symlink()
    assert not (root/'symlink_head_store/g_no_follow').exists()
    _,shape=parse((root/'selection_badshape.bin').read_bytes())
    assert shape['Branch']['tag']==1 and shape['Branch']['shape']==[2,3]  # Valid SHA, semantic shape invalid.
    try:selection((root/'selection_badshape.bin').read_bytes(),37)
    except AssertionError:pass
    else:raise AssertionError('independent branch-shape rejection failed')
    for name in ['g_stale','g_bad_semantic','g_bad_binding','g_fault_stage','g_forbidden_after_uncertainty']:
        assert not (root/'store'/name).exists()
    assert not (root/'escape').exists()
    return dict(status='generation_publication_and_expected_fault_evidence_independently_matched',branches=branches,
        actual_durable_publishes=2,controlled_uncertainty_publishes=1,expected_refusals=19,
        synthetic_tail_bytes=len(TAIL),physical_advances=0,full_restart_qualified=False)


def assess_restore(root,publish,native):
    root=Path(root);source=(root/'session_candidate.bin').read_bytes();physical=(root/'physical.phy').read_bytes()
    assert native['status']=='explicit_older_generation_separate_process_restore_passed_NOT_physical_restart'
    assert native['explicit_generation']==37 and native['ignored_current_generation']==38 and native['durable_publishes']==1
    assert native['archive_bytes']==len(source) and native['physical_prefix_bytes']==len(physical)
    assert native['nonlinear_solver_calls']==native['physical_advances']==0 and not native['full_restart_qualified']
    restored=(root/'selection_restored.bin').read_bytes();uncertain=(root/'selection_uncertain.bin').read_bytes()
    branch=verify_branch(root,'g_37_restored',37,source,physical,restored)
    assert (root/'store/HEAD').read_bytes()==restored
    assert (root/'store/previous_g_37_restored.bin').read_bytes()==uncertain
    for info in publish['branches']:
        pointer=(root/({'g_37_initial':'selection_initial.bin','g_37_next':'selection_next.bin','g_38_uncertain':'selection_uncertain.bin'}[info['branch']])).read_bytes()
        verify_branch(root,info['branch'],info['generation'],source,physical,pointer)
    return dict(status='F_generation_transactions_owned_swaps_and_explicit_separate_process_restore_independently_matched_NOT_physical_restart',
        durable_publishes=3,controlled_uncertain_publish=1,expected_refusals=19,explicit_restore_generation=37,
        retained_newer_generation=38,restored_branch=branch,all_prior_branches_and_synthetic_tail_preserved=True,
        actual_file_and_directory_fsync_and_rename_success=True,power_loss_experiment_performed=False,
        new_physical_solver_calls=0,physical_advances=0,accepted_step_writer_qualified=False,
        running_solver_ELF_identity_bound=False,full_restart_qualified=False,S3_complete=False,S4_complete=False,
        limits='Converted existing payload and synthetic uncommitted tail; filesystem/owned checkpoint component, not a new accepted PD trajectory or complete physical restart qualification.')
