import unittest
from decimal import Decimal as D
from peri_tetra_constraints_audit import gradient,mv,inspect


class GlobalConstraintAssessorTests(unittest.TestCase):
    def test_reference_partition_translation_and_affine_exactness(self):
        b=gradient(); q=list(map(D,[0,0,0,1,0,0,0,1,0,0,0,1]))
        self.assertEqual(mv(b,q),list(map(D,[1,0,0,0,1,0,0,0,1])))
        for k in range(3):
            rigid=[D(int(i%3==k)) for i in range(12)]
            self.assertEqual(mv(b,rigid),[D(0)]*9)

    def test_quadratic_full_matrix_limit_and_wrong_neighbor_response(self):
        case=dict(name='quad',mode=1,mass=[.01,.02],position=[0,0],velocity=[0,0],stiffness=[[2,-1],[-1,2]],
                  dt=.1,rows=[dict(C=[1,0],old=0,new=.01)])
        raw=dict(converged=True,velocity=[.1,1/60],normal_impulse=[.0019166666666666667])
        inspect(case,raw,0)
        raw['velocity'][1]=0
        with self.assertRaises(AssertionError): inspect(case,raw,0)

    def test_infeasibility_needs_a_real_geometry_certificate(self):
        case=dict(name='infeasible',mode=0,expected_code=2,expected_status='no_KKT_feasible_active_set',
                  rows=[dict(C=[1],new=.1),dict(C=[-1],new=.1)])
        raw=dict(converged=False,status='no_KKT_feasible_active_set',old_history_unchanged=True)
        inspect(case,raw,2)
        case['rows'][1]['new']=-.1
        with self.assertRaises(AssertionError): inspect(case,raw,2)

    def test_invalidity_label_alone_cannot_pass(self):
        case=dict(name='mass',mode=0,expected_code=1,invalid_reason='nonpositive_mass',node_mass=[0,.01,.01,.01])
        inspect(case,dict(error='invalid mass'),1)
        case['node_mass'][0]=.01
        with self.assertRaises(AssertionError): inspect(case,dict(error='invalid mass'),1)


if __name__=='__main__': unittest.main()
