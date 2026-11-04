import copy
import math
import unittest
from fractions import Fraction as Q
from peri_guarded_contact import controls
from peri_guarded_contact_audit import bits,check_balance,inspect,oracle

class ExactContactOracleTests(unittest.TestCase):
    def test_reference_recovers_low_product_and_cancellation_bits(self):
        exact,value,bound=oracle(1.,[1+2.**-27],[1-2.**-27],1.)
        self.assertEqual(exact,Q(2)**-54);self.assertEqual(Q(value),exact);self.assertEqual(bound,0)
        exact,value,bound=oracle(.0625,[1,1,1],[2.**70,.125,-2.**70],1.)
        self.assertEqual(exact,Q(-1,16));self.assertEqual(value,-.0625);self.assertEqual(bound,0)

    def test_outward_subnormal_bound_cannot_be_hidden_as_zero(self):
        exact,value,bound=oracle(0.,[2.**-1074],[2.**-1074],1.)
        self.assertLess(exact,0);self.assertEqual(value,0);self.assertEqual(bound,Q(2)**-1074)
        with self.assertRaises(AssertionError):check_balance(0.,[2.**-1074],[2.**-1074],1.,value,0.)

    def test_guard_is_absolute_normalized_accuracy_not_just_finiteness(self):
        exact,value,bound=oracle(1e9,[],[],3.)
        self.assertTrue(math.isfinite(value));self.assertGreater(bound,Q(float(1e-12)))
        allowed,_,_=check_balance(1e9,[],[],3.,value,float(bound));self.assertFalse(allowed)

    def test_halfway_rounding_uses_even_binary_significand(self):
        self.assertEqual(oracle(1.,[-2.**-53],[1.],1.)[1],1.)
        self.assertEqual(oracle(1.+2.**-52,[-2.**-53],[1.],1.)[1],1.+2.**-51)

    def test_false_primitive_acceptance_is_rejected(self):
        case=next(x for x in controls() if x['name']=='absolute_rounding_budget_rejection')
        _,value,bound=oracle(1e9,[],[],3.)
        raw=dict(manufactured_operator_only=True,arithmetic_only=True,rounding_mode_requested=0,
            value_bits=bits(value),bound_bits=bits(float(bound)),normalized_value=value,roundoff_bound=float(bound),within_budget=False)
        inspect(case,raw,2);raw['within_budget']=True
        with self.assertRaises(AssertionError):inspect(case,raw,2)

    def test_invalid_label_requires_actual_invalid_scale(self):
        case=next(x for x in controls() if x['name']=='invalid_scale')
        raw=dict(manufactured_operator_only=True,error='scale')
        inspect(case,raw,1);case=copy.deepcopy(case);case['scale_bits']=bits(1.)
        with self.assertRaises(AssertionError):inspect(case,raw,1)

if __name__=='__main__':unittest.main()
