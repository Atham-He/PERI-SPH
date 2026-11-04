"""Independent byte oracle for saved-record durable writer replay (not Newton)."""
import hashlib,json,struct
from pathlib import Path
from peri_canonical_archive_audit import parse
from peri_bounded_session_archive_audit import encode,make_fields
from peri_bounded_mechanical_stream_audit import parse_stream
from peri_bounded_generation_store_audit import selection

LABELS={'wrong_ELF','wrong_source','wrong_build','wrong_input','wrong_bound_contract','wrong_FP',
        'invalid_reaction','wrong_old_hash','extra_old_tail','before_append','partial_write','full_write_no_sync','sync_error'}


def create(root,old,original_input,identities):
    root=Path(root);old=Path(old)
    _,old_candidate=parse((old/'session_candidate.bin').read_bytes())
    body={k.removeprefix('body/'):v for k,v in old_candidate.items() if k.startswith('body/')}
    fields,trusted=make_fields(encode(body,37),original_input,identities)
    records=parse_stream((old/'physical.phy').read_bytes());assert len(records)==41
    header=b'PERIPHY1'+struct.pack('<Q',1)+bytes.fromhex(trusted)
    pieces=[header];prefixes=[]
    for i,record in enumerate(records):
        data=encode(record,i+1);pieces += [struct.pack('<Q',len(data)),data]
        prefixes.append(b''.join(pieces))
    candidates=[]
    def field(tag,shape,values):
        if not isinstance(values,list):values=[values]
        return dict(tag=tag,shape=shape,raw=struct.pack('<'+{1:'B',2:'Q',3:'d'}[tag]*len(values),*values))
    for i,record in enumerate(records):
        candidate=dict(fields)
        for key,value in record.items():
            if key.startswith('body/'):candidate[key]=value
        candidate['output/Bytes']=field(2,[],len(prefixes[i]));candidate['output/Records']=field(2,[],i+1)
        candidate['output/LastStep']=field(2,[],i)
        candidate['output/LastTime']=record['body/accepted/AcceptedTime']
        candidate['output/HeaderSHA256']=field(1,[32],list(hashlib.sha256(header).digest()))
        candidate['output/PrefixSHA256']=field(1,[32],list(hashlib.sha256(prefixes[i]).digest()))
        candidate['accepted_algorithm/NextPhase']=field(2,[],3 if i==40 else (0 if i<16 else (1 if i<24 else 2)))
        candidates.append(candidate)
    for name,data in [('reference_stream.phy',prefixes[40]),('initial.phy',prefixes[0]),('initial_candidate.bin',encode(candidates[0],1))]:
        with (root/name).open('xb') as stream:stream.write(data)
    return dict(trusted=trusted,records=records,prefixes=prefixes,candidates=candidates,
                inputs=[root/'reference_stream.phy',root/'initial.phy',root/'initial_candidate.bin'])


def branch(root,sub,name,step,generation,prepared):
    directory=Path(root)/sub/name
    assert (directory/'state.bin').read_bytes()==encode(prepared['candidates'][step],generation)
    assert (directory/'prefix.phy').read_bytes()==prepared['prefixes'][step]
    return directory


def assess_phase(root,prepared,mode,native):
    root=Path(root)
    assert native['nonlinear_solver_calls']==native['new_physical_steps']==0 and not native['full_restart_qualified']
    if mode=='--faults':
        assert native['status']=='accepted_writer_identity_prepare_and_partial_tail_controls_passed_NOT_physical_restart'
        assert native['primary_controls']==13 and native['blocked_continuation_controls']==18
        assert {x['label'] for x in native['controls']}==LABELS and len(native['controls'])==13
        initial=prepared['prefixes'][0];frame=struct.pack('<Q',len(encode(prepared['records'][1],2)))+encode(prepared['records'][1],2)
        for control in native['controls']:
            label=control['label'];directory=branch(root,'faults/'+label,'g_1_initial',0,1,prepared)
            active=(directory/'active.phy').read_bytes();w=control['witness']
            assert not w['committed'] and control['writer_stopped']==(label not in ['wrong_ELF','wrong_source','wrong_build','wrong_input'])
            assert selection((root/'faults'/label/'HEAD').read_bytes(),1)[0]=='g_1_initial'
            assert not (root/'faults'/label/'g_forbidden').exists()
            if label in ['partial_write','full_write_no_sync','sync_error']:
                written=len(frame)//2 if label=='partial_write' else len(frame)
                assert active==initial+frame[:written] and w['written_bytes']==written and w['frame_bytes']==len(frame)
                assert w['old_step']==0 and w['proposed_step']==1 and w['old_committed_bytes']==len(initial)
                assert w['old_time']==0 and w['proposed_time']==.01 and w['error_number']==5 and w['injected']
                assert w['failure_stage']==(4 if label=='sync_error' else 2)
                assert control['status']==(3 if label=='sync_error' else 1)
            elif label=='wrong_old_hash':assert active==b'Q'+initial[1:] and w['written_bytes']==0
            elif label=='extra_old_tail':assert active==initial+b'UNSELECTED\n' and w['written_bytes']==0
            else:assert active==initial and w['written_bytes']==0
        return dict(status='writer_faults_independently_matched',primary_controls=13,blocked_controls=18)
    assert native['status']=='saved_record_durable_replay_phase_passed_NOT_new_physical_trajectory'
    assert native['mode']==mode and native['running_replay_ELF_identity_bound'] is True
    expected_count={'--continuous':40,'--split-first':16,'--resume':24,'--rollback':24}[mode]
    assert native['accepted_record_replays']==expected_count
    if mode=='--split-first':
        assert (root/'split_first_final.bin').read_bytes()==encode(prepared['candidates'][16],16)
        directory=branch(root,'split','g_16_mid',16,16,prepared)
        assert (directory/'active.phy').read_bytes()==prepared['prefixes'][16]
        assert (root/'split/g_1_initial/active.phy').read_bytes()==prepared['prefixes'][16]
        assert selection((root/'split_selection_16.bin').read_bytes(),16)[0]=='g_16_mid'
    else:
        name={'--continuous':'continuous','--resume':'split_resumed','--rollback':'rollback'}[mode]
        assert (root/(name+'_final.bin')).read_bytes()==encode(prepared['candidates'][40],100)
        sub='continuous' if mode=='--continuous' else 'split'
        directory=branch(root,sub,'g_100_rollback_final' if mode=='--rollback' else 'g_100_complete',40,100,prepared)
        assert (directory/'active.phy').read_bytes()==prepared['prefixes'][40]
        if mode in ['--resume','--rollback']:
            assert (root/('rollback_noadvance.bin' if mode=='--rollback' else 'split_noadvance.bin')).read_bytes()==encode(prepared['candidates'][16],16)
            restored=branch(root,'split','g_16_rollback' if mode=='--rollback' else 'g_16_resumed',16,16,prepared)
            assert (restored/'active.phy').read_bytes()==prepared['prefixes'][40]
            branch(root,'split','g_16_mid',16,16,prepared)
        if mode=='--rollback':
            branch(root,'split','g_100_complete',40,100,prepared)
            assert (root/'continuous_final.bin').read_bytes()==(root/'split_resumed_final.bin').read_bytes()==(root/'rollback_final.bin').read_bytes()
    return dict(status='saved_record_phase_all_archive_prefix_bytes_independently_matched',mode=mode,replayed_record_commits=expected_count)


def final_assessment(root,prepared,phases):
    assert len(phases)==5
    assert sum(x.get('replayed_record_commits',0) for x in phases)==104
    root=Path(root)
    assert (root/'continuous_final.bin').read_bytes()==(root/'split_resumed_final.bin').read_bytes()==(root/'rollback_final.bin').read_bytes()
    return dict(status='actual_running_replay_identity_and_durable_record_writer_full_bytes_and_failure_tails_independently_matched_NOT_new_physical_restart',
        saved_record_commits=104,new_physical_solver_calls=0,new_physical_steps=0,
        running_replay_ELF_identity_bound=True,primary_fault_controls=13,blocked_continuation_controls=18,
        continuous_split_replay_and_explicit_rollback_archive_bytes_identical=True,
        all_replay_physical_prefix_bytes_identical=True,saved_step16_separate_process_noadvance_bytes_identical=True,
        stopped_partial_full_uncommitted_tails_preserved=True,frames=41,final_saved_step=40,
        final_archive_bytes=len((root/'continuous_final.bin').read_bytes()),final_stream_bytes=len(prepared['prefixes'][40]),
        new_physical_native_restart_qualified=False,full_restart_qualified=False,S3_complete=False,S4_complete=False,
        limits='Real append/fsync/owned commits using existing saved-record replays; not a newly solved PD trajectory,independent physical qualification,or whole S3 completion.')
