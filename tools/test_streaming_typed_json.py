"""New exact-byte equivalence controls,not old component test reruns."""
import decimal as dec
from fractions import Fraction
import hashlib,io,json,math,resource,time,unittest
from peri_reference_checkpoint import pack,canonical,context_descriptor
from peri_streaming_typed_json import stream_typed,streamed_typed_sha


class StreamingTypedControls(unittest.TestCase):
    def equivalent(self,value):
        old=canonical(pack(value));out=io.BytesIO();size=stream_typed(value,out.write,buffer_bytes=127)
        self.assertEqual(out.getvalue(),old);self.assertEqual(size,len(old))
        self.assertEqual(streamed_typed_sha(value),hashlib.sha256(old).hexdigest())

    def test_scalars_and_binary_signed_zero(self):
        self.equivalent([None,True,False,0,-1234567890123456789012345,0.,-0.,1/7,b'\x00\xff',"a\n\"\\\u00e9\U0001f6e0"])

    def test_exact_decimal_quantum_flags_and_special_finite_domain(self):
        values=[dec.Decimal('-0E-95'),dec.Decimal('123.45000'),dec.Decimal('1E-999'),dec.Decimal('Infinity'),dec.Decimal('-Infinity')]
        before=context_descriptor();self.equivalent(values);self.assertEqual(before,context_descriptor())

    def test_rational_tuple_keys_and_heterogeneous_sets(self):
        self.equivalent({(Fraction(2,7),-0.):{None,False,1,'1'},'cache':frozenset({(1,2),(3,4)}),Fraction(-3,5):[Fraction(1,10)]})

    def test_unicode_key_order_is_canonical_encoded_order(self):
        self.equivalent({'\U0001f6e0':1,'\u00e9':2,'\\u00e9':3,'a':4,(-1,):5,2:6})

    def test_context_and_selected_capture_nested_shapes(self):
        self.equivalent(dict(context=context_descriptor(),capture=dict(score=1e-9,selected=dict(old=[[-0.,1.,2.]],
            z=[[dec.Decimal(1)/dec.Decimal(7)]*9],loss=[dec.Decimal('1E-90')],
            material_cache={(3,(.01,-0.),'.005',True):dict(T=[[Fraction(1,3),dec.Decimal('-0E-95')]])}))))

    def test_large_nonsymmetric_global_matrix_byte_and_hash_equivalence(self):
        n=447
        value=dict(original_nonsymmetric_J=[[float((13*i+7*j)%17)/17 for j in range(n)] for i in range(n)],
            point_T=[[[float((i+j+k)%13)/13 for k in range(81)] for j in range(33)] for i in range(2)])
        self.equivalent(value)

    def test_hash_sink_does_not_retain_emitted_payload(self):
        n=257;value=[list(range(n)) for _ in range(n)];chunks=[]
        size=stream_typed(value,lambda b:chunks.append(len(b)),buffer_bytes=4096)
        self.assertGreater(size,500000);self.assertLessEqual(max(chunks),8192)

    def test_invalid_values_nesting_and_byte_budget_reject(self):
        for value in [float('inf'),float('nan'),dec.Decimal('NaN'),object()]:
            with self.assertRaises(ValueError):streamed_typed_sha(value)
        value=0
        for _ in range(130):value=[value]
        with self.assertRaisesRegex(ValueError,'nesting'):streamed_typed_sha(value)
        with self.assertRaisesRegex(ValueError,'byte budget'):streamed_typed_sha(list(range(100)),byte_limit=8)


def run():
    import io,sys
    log=io.StringIO();start=time.process_time();result=unittest.TextTestRunner(stream=log,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(StreamingTypedControls))
    return dict(status='new_incremental_typed_canonical_bytes_equivalent' if result.wasSuccessful() else 'new_streaming_typed_controls_failed',
        tests_run=result.testsRun,errors=len(result.errors),failures=len(result.failures),log=log.getvalue(),CPU_seconds=time.process_time()-start,
        peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,root_cause_resolved=False,physical_solver_calls=0,material_calculations=0,old_raw_reads=0)
