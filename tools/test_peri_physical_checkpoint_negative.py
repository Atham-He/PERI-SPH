"""Real candidate bytes/semantic mutations, zero native or material calls."""
import hashlib,struct,unittest
from pathlib import Path
from peri_physical_checkpoint_negative import specification,mutation
from peri_canonical_archive_audit import parse

ROOT=Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material/engineering-s4-physical-header-v1')
class SemanticMetadata(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.a=(ROOT/'initial_candidate.bin').read_bytes();cls.p=(ROOT/'initial.phy').read_bytes()
    def test_declared_unique_complete_controls(self):
        spec=specification();self.assertEqual(spec['schema'],2)
        self.assertEqual(len(spec['controls']),22);self.assertEqual(len({x['name'] for x in spec['controls']}),22)
    def test_all_mutations_remain_correct_SHA_and_original_bytes(self):
        before=(hashlib.sha256(self.a).hexdigest(),hashlib.sha256(self.p).hexdigest())
        for item in specification()['controls']:
            a,p,t,g,h=mutation(self.a,self.p,1,item['name']);parse(a)
            self.assertTrue((a,p,t,g,h)!=(self.a,self.p,self.p[16:48],1,None))
        self.assertEqual(before,(hashlib.sha256(self.a).hexdigest(),hashlib.sha256(self.p).hexdigest()))
    def test_actual_ELF_rebind_preserves_explicit_v2_stream(self):
        a,p,t,g,h=mutation(self.a,self.p,1,'actual_ELF')
        self.assertEqual(p[:16],b'PERIPHY2'+struct.pack('<Q',2));self.assertEqual(p[16:48],t)
        _,f=parse(a);self.assertEqual(f['output/HeaderSHA256']['raw'],hashlib.sha256(p[:48]).digest())
        self.assertEqual(f['output/PrefixSHA256']['raw'],hashlib.sha256(p).digest())
    def test_six_clock_mutations_change_single_declared_field(self):
        _,old=parse(self.a)
        for name in ['accepted_tick','next_phase','output_time','time_denominator','normal_amplitude','tangent_amplitude']:
            a,p,t,g,h=mutation(self.a,self.p,1,name);_,new=parse(a)
            self.assertEqual(sum(new[k]!=v for k,v in old.items()),1);self.assertEqual(p,self.p)

if __name__=='__main__':unittest.main()
