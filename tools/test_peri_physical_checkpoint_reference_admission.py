"""Metadata admission controls; synthetic receipts confer NO physical pass."""
import copy,unittest
from peri_physical_checkpoint_reference import require_full_reference

def receipt():
    names=['wrong_physical_mass','altered_tangent_dual_history','kinetic_using_wrong_mass','swapped_material_algorithmic_losses',
           'wrong_friction_loss','wrong_wall_velocity_cursor','corrupt_compensation_history','endpoint_substituted_for_whole_interval',
           'wrong_physical_tick_cursor','wrong_protocol_denominator']
    return dict(status='new_physical_time_CP_continuous_full_reference_passed_NOT_full_restart',
        backend_selected_before_reference_imports=True,exception_triggered_backend_switch=False,native_solver_calls=0,
        assessment=dict(name='metadata_only',accepted_steps=80,all_trial_point_tangents_checked=True,independent_trial_material_points=1,
                        maximum_normalized_error=1e-12,geometry=dict(nodes=33,physical_DOF=99,normal_rows=116)),
        negative_controls=[dict(name=n,rejected=True) for n in names])

class NewV2ReferenceAdmission(unittest.TestCase):
    c=dict(name='metadata_only',steps=80)
    def test_full_named_scope_required(self):require_full_reference(receipt(),self.c)
    def test_old40step_label_does_not_pass_new80(self):
        r=receipt();r['assessment']['accepted_steps']=40
        with self.assertRaises(AssertionError):require_full_reference(r,self.c)
    def test_missing_one_trial_T_is_not_full(self):
        r=receipt();r['assessment']['all_trial_point_tangents_checked']=False
        with self.assertRaises(AssertionError):require_full_reference(r,self.c)
    def test_missing_clock_negative_is_not_complete(self):
        r=receipt();r['negative_controls'].pop()
        with self.assertRaises(AssertionError):require_full_reference(r,self.c)
    def test_nan_or_duplicate_control_refused(self):
        for kind in ['nan','duplicate']:
            r=receipt()
            if kind=='nan':r['assessment']['maximum_normalized_error']=float('nan')
            else:r['negative_controls'][-1]=r['negative_controls'][0].copy()
            with self.assertRaises(AssertionError):require_full_reference(r,self.c)
    def test_wrong_geometry_scope_refused(self):
        r=receipt();r['assessment']['geometry']['nodes']=4
        with self.assertRaises(AssertionError):require_full_reference(r,self.c)

if __name__=='__main__':unittest.main()
