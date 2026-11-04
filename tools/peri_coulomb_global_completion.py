"""Only never-attempted inputs, after an immutable independent stop diagnosis."""
import argparse,copy,json,resource,shutil,subprocess,sys,traceback
from pathlib import Path
import peri_coulomb_global as original
from peri_coulomb_global_audit import inspect
from peri_coulomb_global_stop_audit import rejected_gather
from peri_elastic_normal import identity,read,record,reserves,save,stamp,verify


def prepare(root,prior,prospective):
    assert not root.exists(),'new unused directory required'
    repo=Path(__file__).resolve().parents[1]
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=repo).strip()
    prep=read(prior/'preparation.json');oldplan=read(prior/'plan.json');terminal=read(prior/'terminal.json')
    diagnosed=read(prior/'independent_review.json');plan=read(prospective);resources=reserves(root.parent)
    assert diagnosed['status']=='partial_original_global_and_conservative_precision_stop_independently_assessed'
    assert not diagnosed['study_passed'] and not diagnosed['complete_global_qualification']
    assert len(terminal['attempts'])==17 and diagnosed['qualified_control_count']==16
    assert plan['normalized_accuracy_gate']==1e-10 and plan['roundoff_gate']==1e-12
    verify(prep['artifacts']+prep['inputs']);verify(diagnosed['inputs'])
    inherited=[]
    for attempt in terminal['attempts']:
        for field in ['stdout','stderr']:verify([attempt[field]]);inherited.append(attempt[field])
    inherited.extend([record(prior/n) for n in ['preparation.json','state.json','terminal.json','plan.json','independent_review.json']])
    inherited.extend([record(prior/'frozen_stop_audit'/n) for n in ['support.json','peri_coulomb_global_stop_audit.py']])
    root.mkdir();(root/'frozen').mkdir();artifacts=prep['artifacts']+prep['inputs']+inherited+[record(prospective)]
    for file in ['tools/peri_coulomb_global_completion.py','tools/peri_coulomb_global_stop_audit.py',
                 'docs/peri/COUPLED_GLOBAL_STOP_AND_COMPLETION_20261105.md']:
        source=repo/file;target=root/'frozen'/source.name;shutil.copy2(source,target);artifacts.extend([record(source),record(target)])
    attempted={a['name'] for a in terminal['attempts']};remaining=[c for c in oldplan['controls'] if c['name'] not in attempted]
    assert len(remaining)==13 and [c['name'] for c in remaining]==diagnosed['unattempted_controls']
    plan.update(controls=remaining,prepared_at=stamp(),identity=identity(),resources=resources,prior_root=str(prior),
        compiled_numerical_source_commit=oldplan['source_commit'],
        controller_source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
        source_and_ELF_unchanged_from_prior=True,original_failed_study_remains_failed=True,
        imported_qualified_controls=diagnosed['independently_qualified_original_controls'],
        imported_conservative_precision_stop=diagnosed['stop'])
    save(root/'plan.json',plan);artifacts.append(record(root/'plan.json'));inputs=[]
    for case in remaining:
        path=root/(case['name']+'.input');shutil.copy2(prior/path.name,path)
        assert path.read_text().split()==case['tokens'];inputs.append(record(path))
    save(root/'preparation.json',dict(artifacts=artifacts,inputs=inputs,executable=prep['executable'],
        compiler_dependencies=prep['compiler_dependencies'],reused_original_ELF=True,no_attempted_native_input_repeated=True))
    print(json.dumps(dict(status='prepared_only13_never_attempted_algebra_inputs',native_cases=len(remaining),
        imported_qualified_controls=16,original_failed_study_preserved=True)))


def review(root):
    assert not (root/'independent_review.json').exists(),'immutable review already exists'
    resource.setrlimit(resource.RLIMIT_AS,(2*2**30,2*2**30));resource.setrlimit(resource.RLIMIT_CPU,(600,600))
    prep=read(root/'preparation.json');plan=read(root/'plan.json');terminal=read(root/'terminal.json');prior=Path(plan['prior_root'])
    oldplan=read(prior/'plan.json');oldterminal=read(prior/'terminal.json');oldreview=read(prior/'independent_review.json')
    inputs=[record(root/n) for n in ['state.json','terminal.json','preparation.json','plan.json']]
    result=dict(status='global_completion_review_failed',recorded_at=stamp(),identity=identity(),original_failed_study_preserved=True)
    assessed=[]
    try:
        verify(inputs);verify(prep['artifacts']+prep['inputs'])
        assert terminal['status']=='runtime_completed_pending_independent_review' and len(terminal['attempts'])==13
        assert oldreview['qualified_control_count']==16 and not oldreview['study_passed']
        raw={};all_cases={}
        # The16 qualified original assessments are imported by immutable SHA,
        # not re-run or re-reviewed. Only new cross-case covariance/scatter and
        # corrupted negative artifacts read their original physical values.
        for case,attempt in zip(oldplan['controls'][:17],oldterminal['attempts']):
            assert case['name']==attempt['name'];all_cases[case['name']]=case
            raw[case['name']]=original.records(attempt['stdout']['path'])
        for case,attempt in zip(plan['controls'],terminal['attempts']):
            assert case['name']==attempt['name'];verify([attempt['stdout'],attempt['stderr']])
            data=original.records(attempt['stdout']['path']);all_cases[case['name']]=case;raw[case['name']]=data
            assessed.append(inspect(case,data,attempt['return_code']))
        assert set(all_cases)=={c['name'] for c in oldplan['controls']}
        original.comparisons(raw);negative=[]
        for label,name in [
            ('wrong_signed_transpose','full_joint_original_Jacobian'),
            ('symmetrized_original_tangent','full_joint_original_Jacobian'),
            ('dropped_normal_radius_derivative','full_joint_original_Jacobian'),
            ('square_box_instead_of_disk','full_joint_original_Jacobian'),
            ('friction_loss_wrong_sign','physical_slide_two_components'),
            ('premature_float_gather','exact_original_gather_cancellation'),
            ('forged_precision_acceptance','precision_rejected_before_iteration')]:
            data=copy.deepcopy(raw[name]);case=all_cases[name];r=data[2].get('original')
            if r is None:
                s=data[2]['solve'];r=s['trace'][-1]['original'] if s['trace'] else s['precision_rejection']['original']
            if label=='wrong_signed_transpose':r['jacobian'][0][3]=1.
            elif label=='symmetrized_original_tangent':r['jacobian'][0][1]=r['jacobian'][1][0]=.125
            elif label=='dropped_normal_radius_derivative':r['jacobian'][4][3]=0.
            elif label=='square_box_instead_of_disk':r['residual'][4:6]=[0.,0.];r['disk'][0]['residual']=[0.,0.]
            elif label=='friction_loss_wrong_sign':r['disk'][0]['friction_loss']=-r['disk'][0]['friction_loss']
            elif label=='premature_float_gather':r['normalized_relative_velocity'][0][0]=0.
            else:r['precision_budget_ok']=True
            path=root/(label+'.negative.jsonl')
            with path.open('x') as stream:
                for row in data:stream.write(json.dumps(row,allow_nan=False)+'\n')
            try:
                if label=='premature_float_gather':rejected_gather(case,data,2)
                else:inspect(case,data,case['expected_code'])
            except AssertionError:negative.append(dict(name=label,rejected=True,artifact=record(path)))
            else:raise AssertionError(('negative accepted',label))
        verify(inputs);verify(prep['artifacts']+prep['inputs'])
        qualified=oldreview['independently_qualified_original_controls']
        result.update(status='bounded_original_global_normal_two_tangent_coupling_independently_qualified_with_explicit_precision_exclusions',
            imported_qualified_controls=qualified,new_assessments=assessed,
            imported_stop_assessment=oldreview['stop'],prior_independent_review=record(prior/'independent_review.json'),
            original_input_count=30,new_native_inputs=13,total_unique_native_inputs=30,supported_original_control_count=18,
            maximum_normalized_error=max([c.get('maximum_normalized_error',0) for c in qualified+assessed]+
                [oldreview['stop']['maximum_normalized_equation_error']]),
            negative_controls=negative,covariance_and_redundant_scatter_passed=True,
            precision_limitation='The originally intended valid extreme1e16 signed-weight sample was refused by a conservative half-unit certificate1>1e-12; its true conversion error is small. The original study remains failed, this sample is excluded from supported inputs, and no threshold/algorithm/source/ELF was changed.',
            accepted_physical_steps=0,actual_PD_friction_qualified=False,full_restart_qualified=False,S3_complete=False,S4_complete=False,calibration_allowed=False)
    except Exception as error:result.update(failure=repr(error),failure_traceback=traceback.format_exc(),completed_assessments=assessed)
    result.update(inputs=inputs);save(root/'independent_review.json',result)
    print(json.dumps(dict(status=result['status'],new_assessments=len(assessed),failure=result.get('failure'))))
    return int('failure' in result)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run','review']);p.add_argument('--root',required=True,type=Path)
    p.add_argument('--prior',type=Path);p.add_argument('--prospective',type=Path);args=p.parse_args()
    if args.action=='prepare':
        fresh=not args.root.exists()
        try:prepare(args.root,args.prior,args.prospective)
        except Exception as error:
            if fresh and args.root.exists():save(args.root/'preparation_error.json',dict(recorded_at=stamp(),failure=repr(error),native_dispatched=False))
            raise
    elif args.action=='run':sys.exit(original.run(args.root))
    else:sys.exit(review(args.root))
