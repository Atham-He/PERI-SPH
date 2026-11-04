import decimal,unittest
from peri_decimal_diagnostic import context_snapshot,simple,exception_operands,failed_context


class ReadOnlyDecimalObservations(unittest.TestCase):
    def test_context_observation_preserves_precision_traps_and_sticky_flags(self):
        with decimal.localcontext() as ctx:
            ctx.prec=95
            decimal.Decimal(1)/decimal.Decimal(3)
            flags=ctx.flags.copy();traps=ctx.traps.copy();observed=context_snapshot()
            self.assertEqual(ctx.prec,95);self.assertEqual(ctx.flags,flags);self.assertEqual(ctx.traps,traps)
            self.assertEqual(observed['precision'],95)
            self.assertEqual(len(bytes.fromhex(observed['owned_context_object_hex'])),type(ctx).__basicsize__)

    def test_operand_snapshot_keeps_exact_decimal_text_and_signed_zero(self):
        self.assertEqual(simple(decimal.Decimal('-0')),{'decimal':'-0'})
        self.assertEqual(simple(decimal.Decimal('1.2345678901234567890123456789E-100')),{'decimal':'1.2345678901234567890123456789E-100'})

    def test_traceback_context_survives_localcontext_exit(self):
        try:
            with decimal.localcontext() as ctx:
                ctx.prec=95
                raise RuntimeError('diagnostic witness')
        except RuntimeError as error:
            self.assertEqual(failed_context(error)['precision'],95)
            self.assertEqual(exception_operands(error),[])


if __name__=='__main__':unittest.main()
