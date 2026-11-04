import copy
import unittest
from decimal import Decimal as D, localcontext
from peri_many_contact import controls
from peri_many_contact_audit import exact_rank, inspect, strict_monotonicity_certificate, system

class ManyContactIndependentTests(unittest.TestCase):
    def test_every_declared_positive_root_has_an_exact_KKT_witness(self):
        with localcontext() as ctx:
            ctx.prec=95
            for case in controls():
                if case['mode'] or case['expected_code']:continue
                self.assertGreater(float(D(strict_monotonicity_certificate(case).split('/')[0])),0)
                f,_,kkt,_=system(case,list(map(D,case['root_reference'])),list(map(D,case['dual_witness'])))
                self.assertLessEqual(max(map(abs,f)),D('1e-80'),case['name'])
                self.assertLessEqual(kkt,D('1e-80'),case['name'])

    def test_dual_redistribution_is_not_a_different_physical_scatter(self):
        with localcontext() as ctx:
            ctx.prec=95
            case=next(x for x in controls() if x['name']=='cubic_mixed')
            v=list(map(D,case['root_reference']));lam=list(map(D,case['dual_witness']))
            lam[0]+=D('.03125');lam[1]-=D('.03125')
            f,_,kkt,_=system(case,v,lam)
            self.assertEqual(max(map(abs,f)),0);self.assertEqual(kkt,0)

    def test_exact_rank_distinguishes_redundancy_from_inconsistency(self):
        self.assertEqual(exact_rank([[1,1],[2,2]]),1)
        self.assertEqual(exact_rank([[1,1,3],[2,2,6]]),1)
        self.assertEqual(exact_rank([[1,1,3],[2,2,7]]),2)

    def test_actual_symmetric_part_must_have_a_positivity_certificate(self):
        case=copy.deepcopy(controls()[0]);case['A'][0][0]='-2'
        with self.assertRaises(AssertionError):strict_monotonicity_certificate(case)

    def test_invalidity_name_is_not_evidence(self):
        case=next(x for x in controls() if x['name']=='zero_row')
        raw=dict(manufactured_operator_only=True,error='zero row')
        inspect(case,raw,1)
        case['rows'][0]['C']=['1','0','0']
        with self.assertRaises(AssertionError):inspect(case,raw,1)

    def test_wrong_minimum_norm_solution_cannot_pass_consistent_equations(self):
        case=next(x for x in controls() if x['name']=='rank_deficient_consistent_linear')
        raw=dict(manufactured_operator_only=True,linear_rank=1,linear_residual=0,value=[1.5,1.5])
        inspect(case,raw,0)
        raw['value']=[1,2]
        with self.assertRaises(AssertionError):inspect(case,raw,0)

if __name__=='__main__':unittest.main()
