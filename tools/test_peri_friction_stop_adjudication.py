"""Coverage must distinguish physical row identity and active-time history."""
import unittest
from peri_friction_stop_adjudication import contact_events, summarize


class StopCoverageChecks(unittest.TestCase):
    def step(self, row, sign, active=True):
        return dict(step=1, phase=0, upper_wall_velocity=[0,0], contacts=[dict(
            row=row, active=active, sliding=active, signs=[sign,sign])])

    def test_opposite_plates_are_not_temporal_reversal(self):
        self.assertEqual(summarize([self.step(0,1),self.step(1,-1)])['same_contact_temporal_reversals'],[0,0])

    def test_detached_reverse_is_not_frictional_reverse(self):
        self.assertEqual(summarize([self.step(0,1),self.step(0,-1,False)])['same_contact_temporal_reversals'],[0,0])

    def test_same_active_row_can_reverse(self):
        self.assertEqual(summarize([self.step(3,1),self.step(3,-1)])['same_contact_temporal_reversals'],[1,1])

    def test_exact_original_gather_and_threshold(self):
        h=dict(C_full_three_components=[[0,1,0],[1e16,1,-1e16],[0,0,1]],
               normalization_impulse_scale=1, plate_ids=[1])
        row=dict(solve=dict(velocity=[1,1,1],impulses=[1,0,0]),
                 contacts=[dict(wall_tangent_velocity=[0,1])])
        event=contact_events(h,row)[0]
        self.assertEqual(event['relative_velocity'],[1,0])
        self.assertEqual(event['signs'],[1,0])
        row['solve']['impulses'][0]=1e-8
        self.assertFalse(contact_events(h,row)[0]['active'])


if __name__=='__main__':unittest.main()
