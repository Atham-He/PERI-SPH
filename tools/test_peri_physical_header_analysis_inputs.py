import json,tempfile,unittest
from pathlib import Path
from peri_elastic_normal import record
from peri_physical_header_analysis_inputs import copy_admission

class NewReadOnlyInputTransport(unittest.TestCase):
    def fixture(self,d):
        native=Path(d)/'native';audit=Path(d)/'audit';native.mkdir();audit.mkdir()
        raw=native/'raw';raw.write_bytes(b'NEW metadata path fixture, no physical content')
        for n in ['initial_candidate.bin','initial.phy','contract.bin','describe_initial.attempt.json']:(native/n).write_bytes(b'path fixture only')
        (native/'header_reference_admission.json').write_text(json.dumps(dict(raw=record(raw))))
        return native,audit,raw
    def test_native_and_audit_roots_are_explicit(self):
        with tempfile.TemporaryDirectory() as d:
            native,audit,raw=self.fixture(d);records=copy_admission(native,audit)
            self.assertEqual(len(records),2);self.assertEqual(records[0]['sha256'],records[1]['sha256'])
            self.assertEqual((native/'header_reference_admission.json').read_bytes(),(audit/'header_reference_admission.json').read_bytes())
    def test_missing_admission_refused_before_execution(self):
        with tempfile.TemporaryDirectory() as d:
            native,audit,raw=self.fixture(d);(native/'header_reference_admission.json').unlink()
            with self.assertRaises(FileNotFoundError):copy_admission(native,audit)
    def test_missing_candidate_refused_before_execution(self):
        with tempfile.TemporaryDirectory() as d:
            native,audit,raw=self.fixture(d);(native/'initial_candidate.bin').unlink()
            with self.assertRaises(AssertionError):copy_admission(native,audit)
            self.assertFalse((audit/'header_reference_admission.json').exists())
    def test_altered_raw_refused_before_execution(self):
        with tempfile.TemporaryDirectory() as d:
            native,audit,raw=self.fixture(d);raw.write_bytes(b'changed')
            with self.assertRaises(AssertionError):copy_admission(native,audit)
            self.assertFalse((audit/'header_reference_admission.json').exists())

if __name__=='__main__':unittest.main()
