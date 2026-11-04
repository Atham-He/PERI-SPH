import json,tempfile,unittest
from pathlib import Path
from peri_bounded_friction_support import complete_cursor,controls,encode,rotate_tensor,rotate_vector


class BoundedFrictionControlTests(unittest.TestCase):
    def cursor(self,rows,error='',code=0):
        with tempfile.TemporaryDirectory() as folder:
            out=Path(folder)/'out';err=Path(folder)/'err'
            out.write_text('\n'.join(json.dumps(r) for r in rows)+'\n');err.write_text(error)
            return complete_cursor(out,err,code)[0]

    def test_native_error_cursor_overrides_complete_looking_record(self):
        rows=[dict(type='identity'),dict(type='step',committed=True,state=dict(step=1,time=.01)),
              dict(type='step',committed=True,state=dict(step=2,time=.02))]
        error=json.dumps(dict(type='native_error',last_durable_accepted_steps=1,last_durable_accepted_time=.01))
        c=self.cursor(rows,error,1)
        self.assertEqual(c['complete_looking_steps'],2)
        self.assertEqual(c['last_durable_accepted_steps'],1)
        self.assertTrue(c['native_error_cursor_used'] and c['durability_confirmed'])

    def test_hard_kill_without_witness_does_not_invent_durable_cursor(self):
        c=self.cursor([dict(type='identity'),dict(type='step',committed=True,state=dict(step=1,time=.01))],code=-9)
        self.assertEqual(c['complete_looking_steps'],1)
        self.assertIsNone(c['last_durable_accepted_steps'])
        self.assertIsNone(c['last_durable_accepted_time'])
        self.assertFalse(c['durability_confirmed'])

    def test_successful_normal_exit_confirms_accepted_records(self):
        c=self.cursor([dict(type='identity'),dict(type='step',committed=True,state=dict(step=1,time=.01)),
                       dict(type='end',completed=True,accepted_steps=1,last_accepted_time=.01)])
        self.assertEqual(c['last_durable_accepted_steps'],1)
        self.assertTrue(c['durability_confirmed'])

    def test_new_cases_are_independent_protocols_with_same_map_zero_limit(self):
        cases=controls()
        self.assertEqual([c['steps'] for c in cases],[16,16,40,40,16])
        self.assertEqual(len({c['name'] for c in cases}),5)
        for c in cases:
            self.assertEqual(len(encode(c).split()),17)
            self.assertLessEqual(c['native_wall_seconds'],3600)
            self.assertLessEqual(c['output_GiB'],4)
        self.assertEqual(cases[-1]['azimuths'],16)
        self.assertEqual(cases[-1]['mu'],[0,0])
        self.assertEqual(cases[-1]['tangent_stroke'],[0,0])
        self.assertTrue(cases[2]['require_two_way_slip'])

    def test_tensor_and_vector_rotation_share_true_transpose(self):
        v=[1,2,3];t=[1,2,3,4,5,6,7,8,9]
        rt=rotate_tensor(t);rv=rotate_vector(v)
        tv=[sum(t[3*i+j]*v[j] for j in range(3)) for i in range(3)]
        self.assertEqual([sum(rt[3*i+j]*rv[j] for j in range(3)) for i in range(3)],rotate_vector(tv))

    def test_opposed_plates_are_not_temporal_reverse_slip(self):
        from peri_bounded_friction_audit import reversed_contact_counts
        signs={(0,0):{1},(1,0):{-1},(0,1):{1},(1,1):{-1}}
        self.assertEqual(reversed_contact_counts(signs),[0,0])
        signs[0,0].add(-1);signs[1,1].add(1)
        self.assertEqual(reversed_contact_counts(signs),[1,1])


if __name__=='__main__':unittest.main()
