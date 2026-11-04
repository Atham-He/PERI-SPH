"""Independent physical-byte and complete Newton-packet restart comparisons.

No native launch and no claim of material/trajectory convergence. A successful
whole reference receipt is a separate mandatory prerequisite for importing
mathematical evidence by exact full-step byte identity.
"""
import hashlib,json,re,struct
from pathlib import Path
from peri_canonical_archive_audit import parse
from peri_bounded_session_archive_audit import encode
from peri_bounded_mechanical_stream_audit import body_fields,mechanical_fields,loads,packed


def record(path):
    path=Path(path)
    with path.open('rb') as stream:digest=hashlib.file_digest(stream,'sha256').hexdigest()
    return dict(path=str(path),size=path.stat().st_size,sha256=digest)


def trace_packets(path,first,last,*,resume_generation=None):
    """Require complete original rows, exact clocks and no hidden extra rows."""
    assert 1<=first<=last<=40
    with Path(path).open('rb') as stream:
        def line():
            raw=next(stream);assert raw.endswith(b'\n'),'partial Newton record'
            return loads(raw),raw
        identity,raw=line();assert identity['type']=='identity';yield identity,raw
        header,raw=line()
        if first==1:
            assert resume_generation is None and header['type']=='header'
            old_time=header['initial_state']['time'];assert old_time==0
        else:
            assert resume_generation is not None and header['type']=='resume_header'
            assert header['selected_generation']==resume_generation and header['accepted_step']==first-1
            assert header['new_initial_state'] is False;old_time=header['accepted_time']
        yield header,raw
        for step in range(first,last+1):
            row,raw=line()
            assert row['type']=='step' and row['committed'] is True and row['local_step_gates_passed'] is True and row['candidate_ready'] is True
            assert row['reason']=='all_applicable_step_gates_passed'
            state=row['state'];interval=row['interval']
            assert state['step']==step and interval['step']==step-1 and interval['old_time']==old_time
            assert state['time']==interval['new_time']>old_time and interval['actual_dt']>0
            old_time=state['time'];yield row,raw
        end,raw=line();assert end['type']=='end' and end['completed'] is True
        assert end['accepted_steps']==last and end['last_accepted_time']==old_time and end['actual_solver_calls']==last-first+1
        yield end,raw
        assert stream.read()==b'','extra records after native end'


def trace_index(path,first,last,*,resume_generation=None):
    rows=[];header=None;identity=None;header_sha256=None;header_bytes=0
    for row,raw in trace_packets(path,first,last,resume_generation=resume_generation):
        if row['type']=='identity':identity=row
        elif row['type'] in ['header','resume_header']:
            header=row;header_sha256=hashlib.sha256(raw).hexdigest();header_bytes=len(raw)
        elif row['type']=='step':
            rows.append(dict(step=row['state']['step'],time=row['state']['time'],bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
    assert len(rows)==last-first+1
    return dict(raw=record(path),identity=identity,header=header,header_sha256=header_sha256,header_bytes=header_bytes,steps=rows,first_step=first,last_step=last,
        scope='SHA of each ENTIRE native step line,including every original Newton/trial/point tangent packet;not endpoint-only comparison')


def compare_new_branch(qualified_continuous,new_branch):
    """Identity is evidence import only; caller authenticates the full review.

If a new trial line differs, return its exact step list for a prospectively
frozen independent new-context assessment. Never downgrade full-state/output
byte equality to a tolerance merely because the final force looks close.
"""
    assert qualified_continuous['first_step']==1 and qualified_continuous['last_step']==40
    first,last=new_branch['first_step'],new_branch['last_step'];assert first<=last<=40
    if first==1:
        assert (new_branch['header_sha256'],new_branch['header_bytes'])==(qualified_continuous['header_sha256'],qualified_continuous['header_bytes'])
    if first>1:
        previous=qualified_continuous['steps'][first-2]
        assert new_branch['header']['accepted_step']==previous['step']
        assert struct.pack('<d',new_branch['header']['accepted_time'])==struct.pack('<d',previous['time'])
    old={x['step']:x for x in qualified_continuous['steps']};different=[]
    for row in new_branch['steps']:
        previous=old[row['step']]
        assert struct.pack('<d',row['time'])==struct.pack('<d',previous['time'])
        if (row['bytes'],row['sha256'])!=(previous['bytes'],previous['sha256']):different.append(row['step'])
    return dict(status='all_new_native_step_packets_byte_identical_to_qualified_continuous' if not different else 'new_trial_packets_require_independent_assessment',
        first_step=first,last_step=last,compared_complete_steps=last-first+1,different_steps=different,
        original_full_reference_not_repeated=not different,full_restart_qualified=False)


def initial_contract(root):
    root=Path(root);data=(root/'initial_candidate.bin').read_bytes();generation,fields=parse(data)
    assert generation==1 and len(fields)==140
    contract={k.removeprefix('contract/'):v for k,v in fields.items() if k.startswith('contract/')}
    trusted=hashlib.sha256(encode(contract,1)).digest()
    body={k.removeprefix('body/'):v for k,v in fields.items() if k.startswith('body/')};assert len(body)==98
    expected=encode(mechanical_fields(body),1);header=b'PERIPHY1'+struct.pack('<Q',1)+trusted
    physical=header+struct.pack('<Q',len(expected))+expected
    assert (root/'initial.phy').read_bytes()==physical
    assert fields['output/PrefixSHA256']['raw']==hashlib.sha256(physical).digest()
    assert fields['output/HeaderSHA256']['raw']==hashlib.sha256(header).digest()
    return dict(fields=fields,body=body,trusted=trusted,header=header,physical=physical)


def candidate_fields(initial,body,prefix,step,time):
    fields=dict(initial['fields']);fields.update({'body/'+name:value for name,value in body.items()})
    def set_(name,value,tag=2,shape=()):
        expected=packed(tag,shape,value);assert fields[name]['tag']==tag and fields[name]['shape']==list(shape);fields[name]=expected
    set_('accepted_algorithm/NextPhase',3 if step==40 else (0 if step<16 else (1 if step<24 else 2)))
    for name,value in [('Records',step+1),('Bytes',len(prefix)),('LastStep',step)]:set_('output/'+name,value)
    set_('output/LastTime',time,3);set_('output/PrefixSHA256',list(hashlib.sha256(prefix).digest()),1,(32,))
    return fields


def derive_continuous(root,path):
    """Rebuild ALL accepted physical fields from actual new native JSON bits.

The initial packet is externally anchored to the independent zero-advance
preflight. This function does not substitute the output for a whole independent
material/reference/geometry review; callers require BOTH proofs.
"""
    initial=initial_contract(root);prefix=initial['physical'];checkpoints={};index=[];header=None;header_sha256=None;header_bytes=0
    for row,raw in trace_packets(path,1,40):
        if row['type']=='header':
            header=row;header_sha256=hashlib.sha256(raw).hexdigest();header_bytes=len(raw)
        elif row['type']=='step':
            body=body_fields(initial['body'],row['state'],row['friction_history'])
            step=row['state']['step'];physical=encode(mechanical_fields(body,row),step+1)
            prefix+=struct.pack('<Q',len(physical))+physical
            index.append(dict(step=step,time=row['state']['time'],bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
            if step in [8,20,32,40]:
                fields=candidate_fields(initial,body,prefix,step,row['state']['time'])
                checkpoints[step]=dict(archive=encode(fields,step),physical=prefix,body=body,
                    full_native_state=row['state'],full_native_friction_history=row['friction_history'])
    assert len(index)==40 and set(checkpoints)=={8,20,32,40}
    return dict(initial=initial,checkpoints=checkpoints,index=dict(raw=record(path),header=header,header_sha256=header_sha256,header_bytes=header_bytes,steps=index,first_step=1,last_step=40))


def selected(store,selection_leaf,generation,expected,trusted):
    store=Path(store);selection=(store/selection_leaf).read_bytes();actual,fields=parse(selection)
    assert actual==generation and set(fields)=={'SelectionSchema','Branch','ArchiveSHA256','PhysicalPrefixSHA256','ContractSHA256'}
    assert fields['SelectionSchema']==packed(2,(),1)
    branch=fields['Branch']['raw'].decode('ascii');assert re.fullmatch(r'g_[A-Za-z0-9_-]{1,94}',branch)
    assert fields['ContractSHA256']==dict(tag=1,shape=[32],raw=trusted)
    archive=store/branch/'state.bin';physical=store/branch/'prefix.phy'
    assert not archive.is_symlink() and not physical.is_symlink()
    archive_bytes=archive.read_bytes();prefix_bytes=physical.read_bytes()
    archived_generation,candidate=parse(archive_bytes);assert archived_generation==generation and len(candidate)==140
    assert prefix_bytes[:48]==b'PERIPHY1'+struct.pack('<Q',1)+trusted
    assert fields['ArchiveSHA256']==dict(tag=1,shape=[32],raw=hashlib.sha256(archive_bytes).digest())
    assert fields['PhysicalPrefixSHA256']==dict(tag=1,shape=[32],raw=hashlib.sha256(prefix_bytes).digest())
    assert archive_bytes==expected['archive'],'complete candidate byte mismatch'
    assert prefix_bytes==expected['physical'],'complete committed physical prefix byte mismatch'
    return dict(generation=generation,branch=branch,selection=record(store/selection_leaf),archive=record(archive),physical=record(physical),
        full_candidate_byte_equal=True,full_physical_prefix_byte_equal=True)


def file_tree(root):
    root=Path(root);result={}
    for path in sorted(root.rglob('*')):
        assert not path.is_symlink(),'unexpected symlink in checkpoint input/output tree'
        if path.is_file():result[str(path.relative_to(root))]=record(path)
    return result


def assert_preserved(before,after,*,allowed_changes=frozenset()):
    """Older branch/prefix/tails survive; only explicitly declared pointers move."""
    assert set(before)<=set(after),'prior checkpoint input/output file disappeared'
    for name,value in before.items():
        if name not in allowed_changes:assert after[name]==value,('prior data changed',name)


def restored_noadvance_metadata(path,*,generation,step,time,archive_sha256):
    """Metadata is only ONE part of proof; selected() must check full bytes.

Read-only --validate is deliberately not accepted as positive owned-state
noadvance evidence. Logical generation and physical step are independent.
"""
    raw=Path(path).read_bytes();lines=raw.splitlines(keepends=True)
    assert len(lines)==3 and all(line.endswith(b'\n') for line in lines)
    identity,header,end=[loads(line) for line in lines]
    assert identity['type']=='identity' and header['type']=='resume_header' and end['type']=='restored_noadvance'
    for value in [header,end]:
        assert type(value['selected_generation']) is int and value['selected_generation']==generation
        assert type(value['accepted_step']) is int and value['accepted_step']==step
        assert struct.pack('<d',value['accepted_time'])==struct.pack('<d',time)
    assert header['new_initial_state'] is False
    assert type(end['actual_solver_calls']) is int and end['actual_solver_calls']==0
    assert type(end['new_physical_steps']) is int and end['new_physical_steps']==0
    assert end['owned_state_installed_and_resaved'] is True and end['full_restart_qualified'] is False
    assert end['archive_sha256']==archive_sha256
    return dict(status='native_reports_owned_install_resave_without_advance_full_bytes_still_required',
        generation=generation,accepted_step=step,actual_solver_calls=0,identity=identity,raw=record(path),full_restart_qualified=False)
