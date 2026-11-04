"""Minimal arithmetic/context observations, never a solver or study pass.

Exact integer/Fraction rounding does not use either decimal implementation.
Context byte interpretation is allowed only after its public fields agree
with the explicit 64-bit little-endian CPython3.12 layout profile. No live
context memory is written, and unknown bits are never silently discarded.
"""
import ast
from fractions import Fraction as F
import re
import struct

SIGNALS=dict(InvalidOperation=0x3ba,FloatOperation=0x400,DivisionByZero=4,
    Overflow=0x800,Underflow=0x4000,Subnormal=0x2000,Inexact=0x40,Rounded=0x1000,Clamped=1)
ROUNDINGS=['ROUND_UP','ROUND_DOWN','ROUND_CEILING','ROUND_FLOOR','ROUND_HALF_UP','ROUND_HALF_DOWN','ROUND_HALF_EVEN','ROUND_05UP','ROUND_TRUNC']


def context_words(snapshot):
    raw=bytes.fromhex(snapshot['owned_context_object_hex'])
    if snapshot['context_type_basicsize']!=96 or len(raw)!=96:
        return dict(status='unsupported_context_layout_do_not_guess_offsets')
    names=['precision','Emax','Emin','traps','status','newtrap','round','clamp','allcr']
    values=dict(zip(names,struct.unpack_from('<qqqIIIIII',raw,16)))
    if any(values[k]!=snapshot[k] for k in ['precision','Emax','Emin','clamp']) or not 0<=values['round']<len(ROUNDINGS) or ROUNDINGS[values['round']]!=snapshot['rounding']:
        return dict(status='unsupported_context_layout_public_fields_disagree')
    public_agreement=all(bool(values[word]&mask)==bool(snapshot[public][name])
        for word,public in [('traps','traps'),('status','flags')] for name,mask in SIGNALS.items())
    if not public_agreement:return dict(status='context_layout_or_public_signal_mapping_disagrees',words=values)
    unknown={name:values[name]&~0x7fff for name in ['traps','status','newtrap']}
    return dict(status='guarded_context_words_observed',words=values,
                unknown_signal_bits=unknown,public_signal_mapping_agrees=True,
                runtime_cause_proven=False)


def exact_decimal(text):
    assert isinstance(text,str) and len(text)<=4096
    match=re.fullmatch(r'([+-]?)([0-9]+)(?:\.([0-9]*))?(?:[eE]([+-]?[0-9]+))?',text)
    assert match,'finite explicit decimal operand required'
    sign,whole,part,exponent=match.groups();part=part or '';power=int(exponent or 0)-len(part)
    assert abs(power)<=1000 and len(whole)+len(part)<=1024,'bounded operand profile'
    numerator=int(whole+part)*(-1 if sign=='-' else 1)
    return F(numerator*10**power) if power>=0 else F(numerator,10**(-power))


def power10(exponent):return F(10**exponent) if exponent>=0 else F(1,10**(-exponent))


def round_nearest_even(value,precision):
    assert isinstance(value,F) and 1<=precision<=256
    assert max(value.numerator.bit_length(),value.denominator.bit_length())<=12000,'bounded exact rounding witness'
    if not value:return value
    absolute=abs(value)
    exponent=len(str(absolute.numerator))-len(str(absolute.denominator))
    if absolute<power10(exponent):exponent-=1
    assert power10(exponent)<=absolute<power10(exponent+1)
    quantum=power10(exponent-precision+1);scaled=absolute/quantum
    quotient,remainder=divmod(scaled.numerator,scaled.denominator)
    if 2*remainder>scaled.denominator or (2*remainder==scaled.denominator and quotient%2):quotient+=1
    return (-1 if value<0 else 1)*quotient*quantum


def arithmetic(expression,operands,convert,round_result=lambda x:x):
    """An arithmetic AST whitelist; never eval source or follow attributes."""
    node=ast.parse(expression,mode='eval').body
    def visit(item):
        if isinstance(item,ast.Name):
            value=operands[item.id]
            assert isinstance(value,dict) and set(value)=={'decimal'}
            return convert(value['decimal'])
        if isinstance(item,ast.UnaryOp) and isinstance(item.op,(ast.UAdd,ast.USub)):
            v=visit(item.operand);return round_result(v if isinstance(item.op,ast.UAdd) else -v)
        assert isinstance(item,ast.BinOp) and isinstance(item.op,(ast.Add,ast.Sub,ast.Mult,ast.Div)),'arithmetic-only operation profile'
        left,right=visit(item.left),visit(item.right)
        if isinstance(item.op,ast.Add):answer=left+right
        elif isinstance(item.op,ast.Sub):answer=left-right
        elif isinstance(item.op,ast.Mult):answer=left*right
        else:answer=left/right
        return round_result(answer)
    return visit(node)


def backend_operation(backend,expression,operands,snapshot):
    assert snapshot['rounding']=='ROUND_HALF_EVEN','this witness does not alter unsupported rounding'
    # This is a NEW declared clean arithmetic witness, not repair or mutation
    # of the failed reader's private Context. Its difference is explicit.
    with backend.localcontext() as context:
        context.prec=snapshot['precision'];context.Emin=snapshot['Emin'];context.Emax=snapshot['Emax']
        context.rounding=backend.ROUND_HALF_EVEN;context.clamp=snapshot['clamp']
        for name,value in snapshot['traps'].items():context.traps[getattr(backend,name)]=value
        for name,value in snapshot['flags'].items():context.flags[getattr(backend,name)]=value
        before=dict(precision=context.prec,flags={k.__name__:v for k,v in context.flags.items()},traps={k.__name__:v for k,v in context.traps.items()})
        try:
            result=arithmetic(expression,operands,backend.Decimal)
            return dict(status='operation_returned',value=str(result),before=before,
                        flags={k.__name__:v for k,v in context.flags.items()})
        except Exception as error:return dict(status='operation_raised',failure=repr(error),before=before)


def assess_operation(expression,operands,snapshot):
    import decimal,_pydecimal
    assert snapshot['precision']==95 and snapshot['rounding']=='ROUND_HALF_EVEN'
    rounded=arithmetic(expression,operands,exact_decimal,lambda x:round_nearest_even(x,95))
    # The minimal witness handles normal finite magnitudes only. It refuses
    # exponent-boundary arithmetic instead of hiding overflow or underflow.
    if rounded:
        magnitude=len(str(abs(rounded.numerator)))-len(str(rounded.denominator))
        assert snapshot['Emin']+100<magnitude<snapshot['Emax']-100,'outside bounded normal-exponent witness'
    results={name:backend_operation(backend,expression,operands,snapshot)
             for name,backend in [('C_decimal',decimal),('pure_Python_decimal',_pydecimal)]}
    for value in results.values():
        value['matches_exact_nearest_even']=value['status']=='operation_returned' and exact_decimal(value['value'])==rounded
    return dict(expression=expression,exact_rounded_numerator=str(rounded.numerator),
        exact_rounded_denominator=str(rounded.denominator),backends=results,
        context_observation=context_words(snapshot),
        clean_witness_context_scope='New contexts copy declared public fields only; private unknown bits are reported, not recreated or cleared in the original reader.',
        original_context_repaired=False,original_failure_explained=False,study_passed=False)
