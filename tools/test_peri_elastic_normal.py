"""Assessor sensitivity using a hand-computed two-step mass/spring history."""
import copy
import unittest
from peri_elastic_normal_audit import assess


class ElasticAssessorTests(unittest.TestCase):
    def setUp(self):
        self.model=['0.01','1','0','0','0.2','2']
        self.raw={'completed':True,'steps':[]}
        for i,impulse,ke0,ke1,el0,el1,work,loss in [
            (1,.051,0,.00005,0,.005,.0051,.00005),
            (2,.15,.00005,.00005,.005,.02,.015,0)]:
            self.raw['steps'].append(dict(step=i,time_old=i-1,time_new=i,dt=1,
                q_old=.1*(i-1),v_old=0 if i==1 else .1,wall_old=.1*(i-1),wall_new=.1*i,
                q_new=.1*i,v_new=.1,impulse=impulse,reaction_interval=impulse,
                elastic_force=-.05 if i==1 else -.15,kinetic_old=ke0,kinetic_new=ke1,
                elastic_old=el0,elastic_new=el1,wall_work=work,algorithmic_loss=loss,
                energy_residual=0,maximum_momentum_residual=0,primal_velocity_error=0,
                natural_velocity_error=0,stationarity_impulse_error=0,complementarity_energy_error=0,sweeps=1))

    def test_hand_calculated_physical_balance(self):
        result=assess(self.raw,self.model,1)
        self.assertLess(result['maximum_normalized_exact_error'],1e-14)

    def test_effective_metric_cannot_be_kinetic_mass(self):
        self.raw['steps'][0]['kinetic_new']=.5*.51*.1**2
        with self.assertRaises(AssertionError):
            assess(self.raw,self.model,1)

    def test_corrupted_force_cannot_hide_in_balanced_net_force(self):
        self.raw['steps'][1]['reaction_interval']*=1.01
        with self.assertRaises(AssertionError):
            assess(self.raw,self.model,1)

    def test_partial_trajectory_is_not_complete(self):
        self.raw['steps'].pop()
        with self.assertRaises(AssertionError):
            assess(self.raw,self.model,1)


if __name__=='__main__':
    unittest.main()
