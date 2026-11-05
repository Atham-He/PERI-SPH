"""NEW prospective CP admission tests; no native or old physical review calls."""
import copy,json,unittest
from pathlib import Path
from peri_actual_checkpoint_campaign import require_tail_qualification


class ExplicitTailAdmission(unittest.TestCase):
    def zero(self):
        anchor=dict(review=dict(path='synthetic-qualified-normal',size=1,sha256='ab'*32))
        q=dict(name='same114sample_zero_mu',accepted_steps=16,all_trial_point_tangents_checked=True,
            maximum_normalized_error=1e-15,independent_trial_material_points=33,
            geometry=dict(nodes=33,physical_DOF=99,normal_rows=228,exact_rigidity_rank=93),independent_stability=dict(rigid_modes=6))
        review=dict(status='never_attempted_original_tail_independently_passed',backend_selected_before_reference_imports=True,
            physical_solver_calls=0,original_failed_study_unchanged=True,assessment=q,
            normal_limit_comparison=dict(status='new_zero_mu_same114sample_normal_limit_passed',steps=16,
                maximum_normalized_error=1e-15,old_solver_dispatched=False,old_review_repeated=False,old_qualified_review_imported=anchor['review']))
        return review,anchor
    def test_explicit_zero_normal_limit_accepts_only_declared_scope(self):
        r,a=self.zero();require_tail_qualification('same114sample_zero_mu',r,a)
        changes=[('missing_comparison',lambda x:x.pop('normal_limit_comparison')),
            ('missing_scope',lambda x:x['normal_limit_comparison'].update(status='passed')),
            ('too_large_error',lambda x:x['normal_limit_comparison'].update(maximum_normalized_error=2e-10)),
            ('negative_error',lambda x:x['normal_limit_comparison'].update(maximum_normalized_error=-1)),
            ('nan_error',lambda x:x['normal_limit_comparison'].update(maximum_normalized_error=float('nan'))),
            ('partial_anchor_steps',lambda x:x['normal_limit_comparison'].update(steps=15)),
            ('old_native_repeated',lambda x:x['normal_limit_comparison'].update(old_solver_dispatched=True)),
            ('old_review_repeated',lambda x:x['normal_limit_comparison'].update(old_review_repeated=True)),
            ('different_anchor',lambda x:x['normal_limit_comparison'].update(old_qualified_review_imported=dict(path='other'))),
            ('different_mapping',lambda x:x['assessment']['geometry'].update(normal_rows=116)),
            ('different_nodes',lambda x:x['assessment']['geometry'].update(nodes=34))]
        for name,change in changes:
            bad=copy.deepcopy(r);change(bad)
            with self.subTest(name=name),self.assertRaises((AssertionError,KeyError)):
                require_tail_qualification('same114sample_zero_mu',bad,a)
        with self.assertRaises(AssertionError):require_tail_qualification('same114sample_zero_mu',r,None)
    def test_success_label_does_not_replace_complete_independent_evidence(self):
        r,a=self.zero()
        for field,value in [('all_trial_point_tangents_checked',False),('accepted_steps',15),('maximum_normalized_error',2e-10),('independent_trial_material_points',0)]:
            bad=copy.deepcopy(r);bad['assessment'][field]=value
            with self.subTest(field=field),self.assertRaises(AssertionError):require_tail_qualification('same114sample_zero_mu',bad,a)
        for field,value in [('backend_selected_before_reference_imports',False),('physical_solver_calls',1),('original_failed_study_unchanged',False)]:
            bad=copy.deepcopy(r);bad[field]=value
            with self.subTest(field=field),self.assertRaises(AssertionError):require_tail_qualification('same114sample_zero_mu',bad,a)
    def test_actual_completed_strict_receipt_is_imported_not_recomputed(self):
        p=Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material/study-v20-unattempted-friction-tails/actual_PD_strict_stick_hold_release.independent_review.json')
        if not p.exists():self.skipTest('completed fixture not mounted')
        r=json.loads(p.read_text());require_tail_qualification('actual_PD_strict_stick_hold_release',r)
        # A coarse algorithmic-loss ratio over1% is not retroactively a S3
        # failure; its future convergence target stays separate and unchanged.
        self.assertGreater(r['assessment']['algorithmic_loss_positive_input'],.01)
        bad=copy.deepcopy(r);bad['assessment']['observations']['strict_stick_rows']=0
        with self.assertRaises(AssertionError):require_tail_qualification('actual_PD_strict_stick_hold_release',bad)


if __name__=='__main__':unittest.main()
