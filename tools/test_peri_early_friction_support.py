import unittest
from peri_early_friction_support import continuation_controls
from peri_bounded_friction_support import controls, encode


class NewCoverageAdmission(unittest.TestCase):
    def test_new_input_and_only_never_attempted_tail(self):
        old=controls();new=continuation_controls([c['name'] for c in old[:3]])
        self.assertEqual(new[1:],old[3:])
        self.assertTrue(new[0]['require_two_way_slip'])
        self.assertNotIn(encode(new[0]),[encode(c) for c in old])
        for field in ['dt','stroke','rho0','parameters','mu','eta','support','cap_width','steps']:
            self.assertEqual(new[0][field],old[2][field])
        self.assertEqual(new[0]['protocol_kind'],2)

    def test_unexpected_or_additional_attempt_rejects(self):
        names=[c['name'] for c in controls()]
        for attempted in [names[:2],names[:4],list(reversed(names[:3])),names[:3]+names[:1]]:
            with self.assertRaises(AssertionError):continuation_controls(attempted)


if __name__=='__main__':unittest.main()
