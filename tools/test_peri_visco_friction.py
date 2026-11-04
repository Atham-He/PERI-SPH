"""Hand-calculated zero-branch control tests the independent S2 reader."""
import csv
import math
from pathlib import Path
import tempfile
import unittest
from peri_visco_friction_audit import assess, check_checkpoint


class ViscoFrictionReaderTests(unittest.TestCase):
    def setUp(self):
        self.case=dict(name='hand_calculated',kind=2,config=['1','2','.01','.01','1','0','.3','.2','0','0','0','0','0','0','0'])
        self.data=[]; x=0.0; decay=math.exp(-1/.3); alpha=(1-decay)*.3
        for i,jn,ke0,ke1,el0,el1,work,loss in [
            (1,.051,0,.00005,0,.005,.0051,.00005),
            (2,.15,.00005,.00005,.005,.02,.015,0)]:
            xn=decay*x+alpha*.1; ut=.05 if i==1 else .1
            self.data.append(dict(step=i,t0=i-1,t1=i,dt=1,q0=.1*(i-1),vn0=0 if i==1 else .1,
                x0=x,qt0=0,vt0=0,q1=.1*i,vn1=.1,x1=xn,qt1=0,vt1=0,wall0=.1*(i-1),wall1=.1*i,
                wall_ut=ut,jn=jn,jt=0,fn=jn,ft=0,material_force=.05 if i==1 else .15,
                U0=el0,U1=el1,K0=ke0,K1=ke1,Wwall=work,Wbody=0,Dmaterial=0,Dfriction=0,Dnormal=0,
                Dalgorithm=loss,Eres=0,Pnres=0,Ptres=0,primal=0,natural=0,stationarity=0,
                complementarity=0,cone=0,friction_natural=0,slip=-ut,mode=3))
            x=xn
        self.temporary=tempfile.TemporaryDirectory(); self.path=Path(self.temporary.name)/'trace.csv'

    def tearDown(self): self.temporary.cleanup()

    def write(self):
        with self.path.open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(self.data[0]))
            writer.writeheader(); writer.writerows(self.data)

    def test_independent_physical_reference(self):
        self.write(); result=assess(self.case,self.path)
        self.assertLess(result['maximum_normalized_error'],1e-13)

    def test_energy_sum_cannot_hide_wrong_dissipation_labels(self):
        self.data[0]['Dmaterial']=self.data[0]['Dalgorithm']; self.data[0]['Dalgorithm']=0
        self.write()
        with self.assertRaises(AssertionError): assess(self.case,self.path)

    def test_friction_force_corruption_is_detected(self):
        self.data[0]['ft']=.001; self.write()
        with self.assertRaises(AssertionError): assess(self.case,self.path)

    def test_truncated_checkpoint_is_rejected(self):
        path=Path(self.temporary.name)/'bad.cp'; path.write_bytes(b'SVF2')
        self.write()
        with self.assertRaises(AssertionError): check_checkpoint(path,self.case,self.path)


if __name__=='__main__': unittest.main()
