import copy
import unittest
from peri_dense_contact_audit import exact,check


class DenseOracleTests(unittest.TestCase):
    def setUp(self):
        self.case=dict(name='known',mode=0,matrix=[[2,1],[1,2]],mass=[1,1],position=[0,0],velocity=[0,0],
                       dt=1,rows=[dict(C=[1,0],old=0,new=1)])

    def native(self):
        return dict(converged=True,velocity=[1,-.5],normal_impulse=[1.5],sweeps=1,linear_relative_residual=0)

    def test_offdiagonal_exact_solution(self):
        result=exact(self.case)
        self.assertEqual(list(map(float,result['velocity'])),[1,-.5])
        self.assertEqual(list(map(float,result['normal_impulse'])),[1.5])
        self.assertEqual(check(self.case,self.native(),0)['maximum_normalized_error'],0)

    def test_diagonal_substitution_rejected(self):
        bad=self.native(); bad.update(velocity=[1,0],normal_impulse=[2])
        with self.assertRaises(AssertionError): check(self.case,bad,0)

    def test_exact_negative_principal_minor_rejected(self):
        self.case['matrix']=[[1,2],[2,1]]
        with self.assertRaises(ValueError): exact(self.case)
        self.case['mode']=1
        with self.assertRaises(ValueError): exact(self.case)

    def test_infeasible_geometry_has_no_primal(self):
        self.case.update(matrix=[[1]],mass=[1],position=[0],velocity=[0],
            rows=[dict(C=[1],old=0,new=1),dict(C=[-1],old=0,new=1)])
        self.assertIsNone(exact(self.case))


if __name__=='__main__': unittest.main()
