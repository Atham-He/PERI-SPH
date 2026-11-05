import tempfile,unittest
from pathlib import Path
from peri_elastic_normal import record
from peri_physical_header_readonly_v2 import verify_one

class NewHeaderRecordAdmission(unittest.TestCase):
    def test_single_record_admitted(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'raw';p.write_bytes(b'NEW header record')
            verify_one(record(p))
    def test_changed_same_size_bytes_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'raw';p.write_bytes(b'original');r=record(p);p.write_bytes(b'changed!')
            with self.assertRaises(AssertionError):verify_one(r)
    def test_wrong_SHA_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'raw';p.write_bytes(b'NEW header record');r=record(p);r['sha256']='0'*64
            with self.assertRaises(AssertionError):verify_one(r)

if __name__=='__main__':unittest.main()
