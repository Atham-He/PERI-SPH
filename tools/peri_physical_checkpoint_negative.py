"""Correct-SHA semantic inputs for the ACTUAL readonly admission boundary.

All mutations and expected cause classes are declared before dispatch. A valid
binding control precedes corrupt inputs. Positive owned noadvance is separate.
"""
import hashlib,json,struct
from pathlib import Path
from peri_canonical_archive_audit import parse
from peri_bounded_session_archive_audit import encode
from peri_bounded_mechanical_stream_audit import packed
from peri_physical_checkpoint_support import file_tree,record


def specification():
    controls=[
        ('density','bounded archive body: accepted material cursor'),
        ('material_viscous','invalid gel material history'),
        ('neighbor','fixed PD reference: invalid directed packet'),
        ('active_mask','fixed PD reference: invalid directed packet'),
        ('nonaffine_cache','fixed PD reference: asymmetric quadratic cache'),
        ('contact_order','bounded archive body: contact order'),
        ('normal_dual','bounded archive body: stored dual cone'),
        ('tangent_dual','bounded archive body: stored dual cone'),
        ('ledger_correction','mechanical stream: final archive/physical output body mismatch'),
        ('travel_correction','mechanical stream: final archive/physical output body mismatch'),
        ('clock','bounded session: exact accepted protocol clock'),
        ('explicit_generation','canonical archive: schema or explicit generation mismatch'),
        ('external_trust','generation store: external contract selection mismatch'),
        ('expected_HEAD','checkpoint driver: expected HEAD digest mismatch'),
        ('interior_physical_record','physical observation byte mismatch'),
        ('actual_ELF','running contract: actual executable digest mismatch')]
    controls += [('accepted_tick','bounded session: accepted rational tick/history mismatch'),
        ('next_phase','bounded session: exact accepted protocol clock'),
        ('output_time','bounded session: structural output cursor'),
        ('time_denominator','external trusted contract mismatch'),
        ('normal_amplitude','external trusted contract mismatch'),
        ('tangent_amplitude','external trusted contract mismatch')]
    return dict(schema=2,positive_binding_control=True,controls=[dict(name=name,expected_nonzero_cause=cause) for name,cause in controls],
        zero_step_expected_causes=dict(travel_correction='bounded archive body: initial contact history',
            clock='bounded archive body: basic accepted clock',interior_physical_record='initial physical observation'),
        physical_solver_calls=0,live_or_output_mutation_allowed=False)


def changed(fields,name,index,transform):
    result=dict(fields);field=fields[name];code={1:'B',2:'Q',3:'d'}[field['tag']]
    values=list(struct.unpack('<'+code*(len(field['raw'])//struct.calcsize(code)),field['raw']));values[index]=transform(values[index])
    result[name]=dict(field,raw=struct.pack('<'+code*len(values),*values));return result


def corrupt_record(prefix,step):
    header=prefix[:48];offset=48;chunks=[];selected_step=1 if step>0 else 0
    for index in range(step+1):
        length=struct.unpack('<Q',prefix[offset:offset+8])[0];offset+=8;chunk=prefix[offset:offset+length];offset+=length
        generation,fields=parse(chunk);assert generation==index+1
        if index==selected_step:
            fields=changed(fields,'observation/lower/NormalWork',0,lambda value:value+.01);chunk=encode(fields,generation)
        chunks.append(struct.pack('<Q',len(chunk))+chunk)
    assert offset==len(prefix);return header+b''.join(chunks)


def mutation(archive,prefix,generation,name):
    actual,fields=parse(archive);assert actual==generation
    step=struct.unpack('<Q',fields['body/accepted/AcceptedStep']['raw'])[0]
    contract={k.removeprefix('contract/'):v for k,v in fields.items() if k.startswith('contract/')}
    trusted=hashlib.sha256(encode(contract,1)).digest();request_generation=generation;head_override=None
    if name=='density':fields=changed(fields,'body/accepted/GelAcceptedDensity',0,lambda x:x*1.01)
    elif name=='material_viscous':fields=changed(fields,'body/accepted/GelViscous',1,lambda x:x+.01)
    elif name=='neighbor':
        n=struct.unpack('<Q',fields['body/model/PDReferencePosition/count']['raw'])[0]
        fields=changed(fields,'body/model/PDReferenceTarget',0,lambda _:n)
    elif name=='active_mask':fields=changed(fields,'body/model/PDReferenceActive',0,lambda _:2)
    elif name=='nonaffine_cache':fields=changed(fields,'body/model/PDNonaffineUnitCache',1,lambda x:x+.01)
    elif name=='contact_order':fields=changed(fields,'body/contact/PlateIDs',0,lambda _:1)
    elif name=='normal_dual':fields=changed(fields,'body/accepted/LastNormalImpulse',0,lambda x:-2*(1+abs(x)))
    elif name=='tangent_dual':
        normal=struct.unpack('<d',fields['body/accepted/LastNormalImpulse']['raw'][:8])[0]
        fields=changed(fields,'body/accepted/LastTwoTangentImpulse',0,lambda x:10*(1+abs(x)+abs(normal)))
    elif name=='ledger_correction':fields=changed(fields,'body/accepted/ContactCumulativeLedgerCorrection',11,lambda x:x+.01)
    elif name=='travel_correction':fields=changed(fields,'body/accepted/TangentSlipCorrection',0,lambda x:x+.01)
    elif name=='clock':fields=changed(fields,'body/accepted/AcceptedTime',0,lambda x:x+.0001)
    elif name=='explicit_generation':request_generation+=1
    elif name=='external_trust':trusted=bytes([trusted[0]^1])+trusted[1:]
    elif name=='expected_HEAD':head_override='0'*64
    elif name=='interior_physical_record':prefix=corrupt_record(prefix,step)
    elif name=='accepted_tick':fields=changed(fields,'accepted_algorithm/AcceptedTick',0,lambda x:x+1)
    elif name=='next_phase':fields=changed(fields,'accepted_algorithm/NextPhase',0,lambda x:(x+1)%4)
    elif name=='output_time':fields=changed(fields,'output/LastTime',0,lambda x:x+.0001)
    elif name=='time_denominator':fields=changed(fields,'contract/TimeDenominator',0,lambda x:x*2)
    elif name=='normal_amplitude':fields=changed(fields,'contract/PhysicalNormalStroke',0,lambda x:x+.0001)
    elif name=='tangent_amplitude':fields=changed(fields,'contract/PhysicalTangentStroke0',0,lambda x:x+.0001)
    elif name=='actual_ELF':
        fields=changed(fields,'contract/ELFSHA256',0,lambda x:x^1)
        contract={k.removeprefix('contract/'):v for k,v in fields.items() if k.startswith('contract/')}
        trusted=hashlib.sha256(encode(contract,1)).digest();prefix=b'PERIPHY2'+struct.pack('<Q',2)+trusted+prefix[48:]
        fields['output/HeaderSHA256']=packed(1,(32,),list(hashlib.sha256(prefix[:48]).digest()))
    else:raise AssertionError(('undeclared mutation',name))
    if name in ['interior_physical_record','actual_ELF']:
        fields['output/PrefixSHA256']=packed(1,(32,),list(hashlib.sha256(prefix).digest()))
    result=encode(fields,generation);parse(result) # Correct envelope SHA and finite typed fields.
    return result,prefix,trusted,request_generation,head_override


def fixture(path,good_archive,good_prefix,archive,prefix,generation,trusted):
    path.mkdir();good=path/'g_good';bad=path/'g_bad';good.mkdir();bad.mkdir()
    for branch,a,p in [(good,good_archive,good_prefix),(bad,archive,prefix)]:
        with (branch/'state.bin').open('xb') as stream:stream.write(a)
        with (branch/'prefix.phy').open('xb') as stream:stream.write(p)
    def selection(branch,a,p,trust):
        return encode(dict(SelectionSchema=packed(2,(),1),Branch=packed(1,(len(branch),),list(branch.encode())),
            ArchiveSHA256=packed(1,(32,),list(hashlib.sha256(a).digest())),PhysicalPrefixSHA256=packed(1,(32,),list(hashlib.sha256(p).digest())),
            ContractSHA256=packed(1,(32,),list(trust))),generation)
    _,good_fields=parse(good_archive)
    good_contract={k.removeprefix('contract/'):v for k,v in good_fields.items() if k.startswith('contract/')}
    good_trust=hashlib.sha256(encode(good_contract,1)).digest()
    (path/'HEAD').write_bytes(selection('g_good',good_archive,good_prefix,good_trust))
    (path/'selected_bad.bin').write_bytes(selection('g_bad',archive,prefix,trusted))
    (path/'SESSION.lock').write_bytes(b'')
    return hashlib.sha256((path/'HEAD').read_bytes()).hexdigest()


def run_controls(root,plan,execution,derived,*,step=16,generation=None):
    from peri_elastic_normal import save,read,stamp,verify
    generation=generation if generation is not None else step
    baseline=derived['checkpoints'][step];archive=baseline['archive'];prefix=baseline['physical'];trusted=bytes.fromhex(plan['trusted_contract'])
    spec=specification();assert spec==plan['negative_specification']
    frozen_plan=record(root/'plan.json');verify([frozen_plan]+plan['artifacts'])
    if step>0:
        qualified=read(root/'continuous.independent_review.json')
        from peri_physical_checkpoint_reference import require_full_reference
        require_full_reference(qualified,plan['case'])
        assert qualified['status']=='new_physical_time_CP_continuous_full_reference_passed_NOT_full_restart' and qualified['raw']==derived['index']['raw']
        _,fields=parse(archive);scale=derived['index']['header']['normalization_impulse_scale']
        def doubles(name):return struct.unpack('<'+'d'*(len(fields[name]['raw'])//8),fields[name]['raw'])
        assert max(map(abs,doubles('body/accepted/GelViscous')))>1e-10,'nonzero material-history coverage missing'
        assert max(doubles('body/accepted/LastNormalImpulse'))>1e-8*scale,'meaningful normal dual missing'
        assert max(map(abs,doubles('body/accepted/LastTwoTangentImpulse')))>1e-10*scale,'meaningful tangent dual missing'
    cases=root/'semantic_inputs';cases.mkdir();results=[]
    # Bind the actual producer with a completely valid explicit input FIRST.
    valid=cases/'valid_binding';head=fixture(valid,archive,prefix,archive,prefix,generation,trusted)
    before=file_tree(valid)
    result=execution.run('semantic_valid_binding',[plan['executable'],'--validate',str(root),str(valid),trusted.hex(),'selected_bad.bin',str(generation),head],wall=45,cpu=30,gib=1)
    rows=[json.loads(x) for x in Path(result['stdout']['path']).read_text().splitlines()]
    assert len(rows)==2 and rows[1]['type']=='validated_noadvance' and rows[1]['new_physical_steps']==0 and rows[1]['archive_sha256']==hashlib.sha256(archive).hexdigest()
    assert all(rows[0][k]==v for k,v in result['birth_identity'].items()) and rows[0]['argv']==result['command']
    assert result['stderr']['size']==0 and file_tree(valid)==before
    results.append(record(root/'semantic_valid_binding.attempt.json'))
    for control in spec['controls']:
        verify([frozen_plan]+plan.get('critical_runtime',plan['artifacts']))
        name=control['name'];a,p,trust,requested,head_override=mutation(archive,prefix,generation,name)
        path=cases/name;expected_head=fixture(path,archive,prefix,a,p,generation,trusted if name=='external_trust' else trust)
        before=file_tree(path);save(root/('semantic_'+name+'.input_receipt.json'),dict(recorded_at=stamp(),inputs=before,mutation=name,
            expected_cause=spec['zero_step_expected_causes'].get(name,control['expected_nonzero_cause']) if step==0 else control['expected_nonzero_cause']))
        result=execution.run('semantic_'+name,[plan['executable'],'--validate',str(root),str(path),trust.hex(),'selected_bad.bin',str(requested),head_override or expected_head],
            wall=45,cpu=30,gib=1,expected_code=1)
        with Path(result['stdout']['path']).open() as stream:identity=json.loads(next(stream))
        assert all(identity[k]==v for k,v in result['birth_identity'].items()) and identity['argv']==result['command']
        witness=result['native_failure_witnesses'];assert len(witness)==1
        expected=spec['zero_step_expected_causes'].get(name,control['expected_nonzero_cause']) if step==0 else control['expected_nonzero_cause']
        assert expected in witness[0]['message'],(name,'wrong refusal cause',witness[0]['message'],expected)
        assert witness[0]['actual_solver_calls']==0 and witness[0]['append_witness'] is None
        assert witness[0]['mechanical_durable_step']==0 and file_tree(path)==before,'invalid input changed owned/output state'
        results.append(record(root/('semantic_'+name+'.attempt.json')))
        verify([frozen_plan]+plan.get('critical_runtime',plan['artifacts']))
    verify([frozen_plan]+plan['artifacts'])
    result=dict(status='physical_time_v2_readonly_CP_correct_SHA_semantic_refusals_and_valid_binding_passed',recorded_at=stamp(),
        accepted_step=step,logical_generation=generation,qualified_physical_nonzero_seed=step>0,
        positive_binding_controls=1,semantic_rejections=len(spec['controls']),results=results,specification=spec,
        no_private_or_live_or_output_changes=True,native_solver_calls=0,full_restart_qualified=False,S3_prior_qualification_preserved=True,S4_complete=False,calibration_allowed=False)
    save(root/'semantic_rejection_review.json',result);return record(root/'semantic_rejection_review.json')
