import copy,decimal,struct,unittest
from fractions import Fraction as F
from peri_decimal_diagnostic import context_snapshot
from peri_decimal_packet_audit import exact_decimal,round_nearest_even,context_words,assess_operation,arithmetic


class MinimalArithmeticWitnessTests(unittest.TestCase):
    def test_integer_rounding_half_even_and_scales(self):
        for value,want in [(F(5,2),F(2)),(F(7,2),F(4)),(F(-5,2),F(-2)),(F(95,10),F(10)),(F(1,10**100),F(1,10**100))]:
            self.assertEqual(round_nearest_even(value,1),want)
        self.assertEqual(exact_decimal('-1.2500e-3'),F(-1,800))
        self.assertEqual(exact_decimal('-0E-100'),0)

    def test_normal_two_backend_witness_uses_fraction_reference(self):
        with decimal.localcontext() as c:
            c.prec=95;context=context_snapshot()
            operands={'a':{'decimal':'1.23456789012345678901234567890123456789'},'e':{'decimal':'9.87654321098765432109876543210987654321'}}
            report=assess_operation('a*e',operands,context)
            self.assertTrue(all(x['matches_exact_nearest_even'] for x in report['backends'].values()))
            self.assertFalse(report['study_passed'])

    def test_public_context_guards_and_unknown_bits_without_live_write(self):
        with decimal.localcontext() as c:
            c.prec=95;context=context_snapshot()
            self.assertEqual(context_words(context)['status'],'guarded_context_words_observed')
            copy_=copy.deepcopy(context);raw=bytearray.fromhex(copy_['owned_context_object_hex'])
            traps=struct.unpack_from('<I',raw,40)[0];struct.pack_into('<I',raw,40,traps|0x80000000)
            copy_['owned_context_object_hex']=raw.hex()
            self.assertEqual(context_words(copy_)['unknown_signal_bits']['traps'],0x80000000)
            struct.pack_into('<q',raw,16,28);copy_['owned_context_object_hex']=raw.hex()
            self.assertEqual(context_words(copy_)['status'],'unsupported_context_layout_public_fields_disagree')
            self.assertEqual(c.prec,95)

    def test_operation_rejects_code_and_nonfinite_inputs(self):
        with self.assertRaises(AssertionError):arithmetic('__import__("os")',{},exact_decimal)
        for text in ['NaN','Infinity','1E1000000000']:
            with self.assertRaises(AssertionError):exact_decimal(text)


if __name__=='__main__':unittest.main()
