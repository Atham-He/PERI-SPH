import hashlib,json,tempfile,unittest
from pathlib import Path
from peri_buffered_metadata import BinaryReader,buffered_record,unique_records,verify_unique


class BufferedMetadataControls(unittest.TestCase):
    def test_identical_pins_each_unique_checked_once_without_omission(self):
        with tempfile.TemporaryDirectory() as d:
            pins=[]
            for k in range(3):
                p=Path(d)/str(k);p.write_bytes(bytes([k])*1234);pins.append(buffered_record(str(p)))
            called=[]
            def inspect(path):called.append(path);return buffered_record(path)
            q=verify_unique([pins[0],pins[1],pins[0],pins[2],pins[1]],inspect)
            self.assertEqual(called,[p['path'] for p in pins]);self.assertEqual(q['declared_records'],5)
            self.assertEqual(q['checked_unique_paths'],3);self.assertEqual(q['identical_redundancies'],2)

    def test_conflicting_pin_rejects_before_any_IO(self):
        a=dict(path='/synthetic',size=3,sha256='a'*64);b=dict(a,sha256='b'*64);called=[]
        with self.assertRaisesRegex(ValueError,'conflicting'):
            verify_unique([a,b],lambda p:called.append(p))
        self.assertEqual(called,[])

    def test_wrong_size_SHA_and_changed_file_reject(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'payload';p.write_bytes(b'original');pin=buffered_record(str(p))
            for wrong in [dict(pin,size=99),dict(pin,sha256='f'*64)]:
                with self.assertRaisesRegex(ValueError,'mismatch'):verify_unique([wrong])
            p.write_bytes(b'changed!')
            with self.assertRaisesRegex(ValueError,'mismatch'):verify_unique([pin])

    def test_multimegabyte_UTF8_lines_all_offsets_and_prefixes_identical(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'raw';lines=[json.dumps(dict(k=k,text='\u00e9\U0001f6e0'*(500000 if k==1 else 3)),ensure_ascii=False).encode()+b'\n' for k in range(3)]
            p.write_bytes(b''.join(lines));r=BinaryReader(p);prefix=b''
            try:
                for line in lines:
                    row,c=r.row();prefix+=line
                    self.assertEqual(row,json.loads(line));self.assertEqual(c['offset'],len(prefix))
                    self.assertEqual(c['prefix_sha256'],hashlib.sha256(prefix).hexdigest())
                    self.assertEqual(c['last_line_sha256'],hashlib.sha256(line).hexdigest())
                self.assertEqual(r.finish(buffered_record(str(p))),hashlib.sha256(prefix).hexdigest())
            finally:r.close()

    def test_partial_record_hidden_tail_and_wholeSHA_reject(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'raw';p.write_bytes(b'{"k":1}')
            r=BinaryReader(p)
            try:
                with self.assertRaisesRegex(ValueError,'newline'):r.row()
            finally:r.close()
            p.write_bytes(b'{"k":1}\n{"hidden":true}\n');r=BinaryReader(p)
            try:
                r.row()
                with self.assertRaisesRegex(ValueError,'hidden'):r.finish(buffered_record(str(p)))
            finally:r.close()
            p.write_bytes(b'{"k":1}\n');r=BinaryReader(p)
            try:
                r.row()
                with self.assertRaisesRegex(ValueError,'whole raw'):r.finish(dict(path=str(p),size=p.stat().st_size,sha256='0'*64))
            finally:r.close()

    def test_invalid_record_domain_and_distinct_literal_paths(self):
        for item in [dict(path='relative',size=0,sha256='a'*64),dict(path='/a',size=True,sha256='a'*64),dict(path='/a',size=0,sha256='x')]:
            with self.assertRaises(ValueError):unique_records([item])
        a=dict(path='/a/../b',size=0,sha256='a'*64);b=dict(path='/b',size=0,sha256='a'*64)
        self.assertEqual(len(unique_records([a,b])),2)
