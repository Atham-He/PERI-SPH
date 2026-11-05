"""Admission boundary for the prospective NEW80-step physical experiment."""
import copy,unittest
from pathlib import Path
from peri_elastic_normal import read
from peri_physical_checkpoint_campaign import campaign_case,require_ready

BASE=Path('/mnt/f/PERI-SPH-verification/2026-11-04-gel-material')
class CampaignAdmission(unittest.TestCase):
    def setUp(self):
        self.h=read(BASE/'analysis-s4-header-v3-admission-path/independent_review.json')
        self.initial=read(BASE/'analysis-s4-header-v3-admission-path/header_independent_review.json')
        self.sem=read(BASE/'engineering-s4-v2-complete-semantic-v1/semantic_rejection_review.json')
        self.binding=read(BASE/'engineering-s4-v2-reference-binding-v1/independent_review.json')
    def test_named_actual_prerequisites_admitted(self):require_ready(self.h,self.initial,self.sem,self.binding)
    def test_clock_grid_fixed_end_and_checkpoints(self):
        source=read(BASE/'engineering-s4-physical-header-v1/plan.json')['case'];old=copy.deepcopy(source)
        case=campaign_case(source);self.assertEqual(source,old);self.assertEqual(case['steps'],80)
        for field,value in [('steps',40),('dt',.0025),('checkpoint_steps',[8,40])]:
            bad=copy.deepcopy(source);bad[field]=value
            with self.assertRaises(AssertionError):campaign_case(bad)
    def test_changed_speed_not_admitted(self):
        bad=read(BASE/'engineering-s4-physical-header-v1/plan.json')['case'];bad['physical_time']['normal_event_ticks']=[80000,120000,200000]
        bad['physical_time']['tangent_event_ticks']=[10000,30000,40000]
        with self.assertRaises(AssertionError):campaign_case(bad)
    def test_missing_true_complete_byte_proof_rejected(self):
        self.initial['accepted_initial_body_all_bytes_match']=False
        with self.assertRaises(AssertionError):require_ready(self.h,self.initial,self.sem,self.binding)
    def test_partial_semantic_scope_rejected(self):
        self.sem['semantic_rejections']=16
        with self.assertRaises(AssertionError):require_ready(self.h,self.initial,self.sem,self.binding)
    def test_C_backend_or_material_selftest_not_admitted(self):
        for key,value in [('C_decimal_loaded',True),('precision',28),('material_calculations',1)]:
            bad=copy.deepcopy(self.binding);bad['observation'][key]=value
            with self.assertRaises(AssertionError):require_ready(self.h,self.initial,self.sem,bad)

if __name__=='__main__':unittest.main()
