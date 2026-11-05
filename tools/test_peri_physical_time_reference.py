import copy
import json
import unittest
from fractions import Fraction as F
import peri_physical_time_reference as ref

class ClockReferenceControls(unittest.TestCase):
    def complete_grid(self):
        t=0;records=[]
        while t<400000:
            r=ref.expected_interval(t,10000)
            r={k:([float(x) for x in v] if isinstance(v,list) else float(v) if isinstance(v,F) else v) for k,v in r.items()}
            r['actual_dt']=r['new_time']-r['old_time']
            records.append(r);t=r['new_tick']
        return records

    def test_exact_reference_clock(self):
        q=ref.assess_grid('\n'.join(map(json.dumps,self.complete_grid())),10000)
        self.assertEqual(q['intervals'],40)
        self.assertLess(q['maximum_normalized_error'],1e-12)

    def test_wrong_physical_speed_rejected(self):
        records=self.complete_grid()
        records[2]['new_wall'][0]+=.001
        with self.assertRaises(AssertionError):ref.assess_grid('\n'.join(map(json.dumps,records)),10000)

    def test_same_endpoints_wrong_interior_cubic_rejected(self):
        records=self.complete_grid()
        # Add .01*s*(1-s): unchanged endpoint values, wrong whole interval.
        records[0]['normal_cubic'][1]+=.01
        records[0]['normal_cubic'][2]-=.01
        with self.assertRaises(AssertionError):ref.assess_grid('\n'.join(map(json.dumps,records)),10000)

    def test_loss_of_tail_rejected(self):
        with self.assertRaises(AssertionError):ref.assess_grid('\n'.join(map(json.dumps,self.complete_grid()[:-1])),10000)

    def test_packet_layout(self):
        self.assertEqual(len(ref.packet(10000)),104)
        self.assertEqual(len(ref.negative_packets()),16)

if __name__=='__main__':unittest.main()
