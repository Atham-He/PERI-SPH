"""Independent diagnosis of a stopped campaign; never calls a native solver."""
from decimal import Decimal as D,localcontext
from fractions import Fraction as F
from pathlib import Path
import json,math,sys,traceback
from peri_elastic_normal import read,record,save,stamp,identity,verify
from peri_coulomb_global_audit import inspect,check_header
from peri_coulomb_contact_audit import (
    frac,decimal,fdot,operator_reference,original_reference,close_tree,dec)


def rejected_gather(case,records,code):
    """Assess the preserved uncommitted sample, never label it a valid root."""
    with localcontext() as ctx:
        ctx.prec=95
        assert case['name']=='exact_original_gather_cancellation' and code==2
        assert [r['type'] for r in records]==['identity','header','global_sample','end']
        p=case['input'];h,r,end=records[1:];native=r['original'];n=p['n'];c=check_header(p,h)
        assert not r['accepted_physical_step'] and not end['completed']
        assert end['accepted_physical_steps']==end['last_accepted_time']==0
        assert not native['precision_budget_ok'] and len(native['jacobian'])==n+3*p['m']
        v=p['trial_velocity'];imp=p['trial_impulses'];I=h['impulse_scale'];V=h['velocity_scale']
        R,J=operator_reference(h,v)
        maximum=close_tree(native['supplied_momentum'],R,dec(I))
        maximum=max(maximum,close_tree(native['supplied_jacobian'],J,dec(I)/dec(V)))
        reference=original_reference(h['old_position'],h['contacts'],h['dt'],I,V,v,imp,R,J,c)
        maximum=max(maximum,close_tree(native['residual'],reference['residual']),close_tree(native['jacobian'],reference['jacobian']))
        actual_errors=[];certificates=[]
        for k,momentum in enumerate(native['supplied_momentum']):
            exact=(frac(momentum)-fdot([row[k] for row in c],imp))/frac(I)
            assert native['residual'][k]==float(exact)
            error=abs(frac(native['residual'][k])-exact)
            absolute=abs(exact)
            if not absolute:expected=F(0)
            else:
                # Exact binary exponent of the ORIGINAL rational result.
                e=absolute.numerator.bit_length()-absolute.denominator.bit_length()
                if absolute<F(2)**e:e-=1
                quantum=F(2)**max(e-52,-1074)
                expected=F(0) if (absolute/quantum).denominator==1 else F(2)**max(e-53,-1074)
            assert frac(native['residual_roundoff_bound'][k])==expected
            assert error<=expected
            actual_errors.append(error);certificates.append(expected)
        assert max(certificates)>frac(1e-12) and max(certificates)==1
        # The exact error is small but the declared pre-existing half-unit
        # certificate is conservative. No guard is bypassed or relaxed.
        assert max(actual_errors)<frac(1e-12)
        assert reference['relative_velocity'][0][0]==1
        assert native['normalized_relative_velocity'][0][0]==1
        assert abs(frac(native['normalized_relative_velocity'][0][0])-reference['relative_velocity'][0][0])<=frac(native['relative_velocity_roundoff_bound'][0][0])
        for a in range(2):
            assert abs(frac(native['normalized_relative_velocity'][0][a])-reference['relative_velocity'][0][a])<=frac(native['relative_velocity_roundoff_bound'][0][a])
        disk=native['disk'][0];want=reference['disks'][0]
        assert disk['representable'] and disk['within_budget'] and disk['branch']==want['branch']
        for field in ['residual','normal_derivative','impulse_derivative','velocity_derivative','friction_loss']:
            maximum=max(maximum,close_tree(disk[field],want[field]))
        assert 0<=native['maximum_jacobian_roundoff_bound']<=1e-12
        for actual,row in zip(native['jacobian'],reference['jacobian']):
            for a,b in zip(actual,row):
                assert abs(dec(a)-b)<=dec(native['maximum_jacobian_roundoff_bound'])+D('1e-80')*max(D(1),abs(b))
        assert all(x<=1e-12 for x in native['residual_roundoff_bound'][n:])
        assert native['merit']==native['squared_norm']==0  # not computed as a valid iterate
        return dict(status='original_half_unit_precision_certificate_rejection_independently_verified',
            original_prospective_expected_code=case['expected_code'],actual_return_code=code,original_study_passed=False,
            root_accepted=False,accepted_physical_steps=0,
            maximum_normalized_equation_error=maximum,maximum_exact_conversion_error=float(max(actual_errors)),
            maximum_preserved_certificate=float(max(certificates)),unchanged_precision_gate=1e-12,
            exact_relative_gather=1.,premature_float_gather=0.,
            limitation='Extreme signed weights are unsupported by the existing conservative accumulation certificate. This is not a material/PD failure and is not a passed valid sample. No relaxation, solver rerun or new certificate was used.')


def assess(root):
    assert not (root/'independent_review.json').exists(),'review already exists'
    prep=read(root/'preparation.json');plan=read(root/'plan.json');terminal=read(root/'terminal.json')
    receipt=dict(recorded_at=stamp(),identity=identity(),status='stopped_global_review_failed',study_passed=False,
        S3_complete=False,actual_PD_friction_qualified=False,calibration_allowed=False)
    inputs=[record(root/n) for n in ['state.json','terminal.json','preparation.json','plan.json']]
    try:
        verify(inputs);verify(prep['artifacts']+prep['inputs'])
        assert terminal['status']=='stopped_for_review' and len(terminal['attempts'])==17
        qualified=[]
        for case,attempt in zip(plan['controls'][:16],terminal['attempts'][:16]):
            verify([attempt['stdout'],attempt['stderr']]);assert case['name']==attempt['name']
            raw=[json.loads(line) for line in Path(attempt['stdout']['path']).read_text().splitlines()]
            qualified.append(inspect(case,raw,attempt['return_code']))
        case=plan['controls'][16];attempt=terminal['attempts'][16]
        verify([attempt['stdout'],attempt['stderr']]);assert case['name']==attempt['name']
        raw=[json.loads(line) for line in Path(attempt['stdout']['path']).read_text().splitlines()]
        stop=rejected_gather(case,raw,attempt['return_code'])
        assert attempt.get('runtime_gate') is None
        receipt.update(status='partial_original_global_and_conservative_precision_stop_independently_assessed',
            independently_qualified_original_controls=qualified,qualified_control_count=16,stop=stop,
            unattempted_controls=[c['name'] for c in plan['controls'][17:]],
            complete_global_qualification=False,full_restart_qualified=False,
            original_failure_preserved=terminal['failure'],
            next_allowed_finite_action='A separately frozen follow-up may import these immutable16 original controls and the correctly diagnosed uncommitted precision refusal, execute only13 NEVER-ATTEMPTED original inputs, and audit new negative controls. It must keep this original study failed and document that the extreme-weight valid-sample expectation was false. No existing native input is redispatched.')
        verify(inputs);verify(prep['artifacts']+prep['inputs'])
    except Exception as error:receipt.update(failure=repr(error),failure_traceback=traceback.format_exc())
    receipt.update(inputs=inputs,support=record(__file__))
    save(root/'independent_review.json',receipt)
    print(json.dumps(dict(status=receipt['status'],qualified=receipt.get('qualified_control_count'),failure=receipt.get('failure'))))
    return int('failure' in receipt)


if __name__=='__main__':sys.exit(assess(Path(sys.argv[1])))
