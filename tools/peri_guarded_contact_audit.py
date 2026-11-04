"""Independent exact-Fraction binary64 accumulation oracle and guard audit."""
import math
import struct
from decimal import Decimal as D, localcontext
from fractions import Fraction as Q
from peri_many_contact_audit import inspect as inspect_many, equation, system, decimal, dot, magnitude

def bits(value):return str(struct.unpack('Q',struct.pack('d',float(value)))[0])
def from_bits(value):return struct.unpack('d',struct.pack('Q',int(value)))[0]
def power(exponent):return Q(2)**exponent

def oracle(momentum,coefficients,impulses,scale):
    """Fraction constructs exact binary inputs; CPython rounds the ratio."""
    assert math.isfinite(scale) and scale>0
    assert all(math.isfinite(x) for x in [momentum,*coefficients,*impulses])
    value=(Q(momentum)-sum((Q(c)*Q(l) for c,l in zip(coefficients,impulses)),Q(0)))/Q(scale)
    rounded=float(value)  # may raise OverflowError, deliberately independent
    if not math.isfinite(rounded):raise OverflowError('oracle overflow')
    error=abs(value-Q(rounded))
    if not error:return value,rounded,Q(0)
    absolute=abs(value);exponent=absolute.numerator.bit_length()-absolute.denominator.bit_length()
    if absolute<power(exponent):exponent-=1
    unit=max(exponent-52,-1074)
    return value,rounded,power(max(unit-1,-1074))

def check_balance(momentum,coefficients,impulses,scale,value,bound):
    exact,rounded,enclosure=oracle(momentum,coefficients,impulses,scale)
    assert bits(value)==bits(rounded) or value==rounded==0
    assert Q(bound)==enclosure and abs(exact-Q(value))<=Q(bound)
    return enclosure<=Q(float(1e-12)),exact,enclosure

def precision_rejection(case,raw,code):
    assert code==2 and not raw['converged'] and raw['velocity']==raw['normal_impulse']==[]
    assert raw['status']=='contact_accumulation_precision_budget'
    p=raw['precision_rejection'];assert p['present']
    n=len(case['position']);v=p['velocity'];lam=p['normal_impulse'];scale=float(case['mass_scale'])*float(case['speed'])
    assert len(v)==len(p['supplied_momentum'])==len(p['normalized_values'])==len(p['roundoff_bounds'])==n
    assert len(lam)==len(case['rows'])
    decisions=[]
    for i in range(n):
        coefficients=[float(row['C'][i]) for row in case['rows']]
        decisions.append(check_balance(p['supplied_momentum'][i],coefficients,lam,scale,p['normalized_values'][i],p['roundoff_bounds'][i])[0])
    assert not all(decisions),'guard requires an exact over-budget rounding enclosure'
    if case['rejection']=='precision_initial':
        assert not p['during_line_search'] and p['iteration']==0 and raw['trace']==[]
        assert all(x==0 for x in v+lam)
    else:
        assert case['rejection']=='precision_trial' and p['during_line_search'] and p['iteration']==0
        assert len(raw['trace'])==1 and p['alpha']==1
        t=raw['trace'][0];assert t['iteration']==0 and t['alpha']==0 and t['domain_rejections']==0
        old=list(map(decimal,t['velocity']));old_lam=list(map(decimal,t['normal_impulse']))
        with localcontext() as ctx:
            ctx.prec=95
            f,h,_,_=system(case,old,old_lam);native_f=list(map(decimal,t['normalized_residual']))
            native_h=[list(map(decimal,row)) for row in t['generalized_jacobian']];d=list(map(decimal,t['direction']))
            assert len(d)==n+len(lam)
            for x,y in zip(f,native_f):assert abs(x-y)<=D('1e-10')
            for row,expected in zip(native_h,h):
                for x,y in zip(row,expected):assert abs(x-y)<=D('1e-10')*max(D(1),abs(y))
            norm=max(sum(map(abs,row),D(0)) for row in native_h)
            error=magnitude([dot(row,d)+b for row,b in zip(native_h,native_f)])
            assert error/(norm*magnitude(d)+magnitude(native_f))<=D('1e-12')
            for i in range(n):assert abs(decimal(v[i])-(old[i]+decimal(case['speed'])*d[i]))<=D('1e-10')*max(D(1),abs(decimal(v[i])))
    return dict(name=case['name'],status='exact_enclosure_proves_precision_rejection_without_commit',during_line_search=p['during_line_search'])

def inspect(case,raw,code):
    assert code==case['expected_code'] and raw['manufactured_operator_only']
    if case['mode']==2:
        r=from_bits(case['momentum_bits']);scale=from_bits(case['scale_bits'])
        c=list(map(from_bits,case['coefficient_bits']));lam=list(map(from_bits,case['impulse_bits']))
        if code==1:
            label=case['invalid']
            if label=='invalid_scale':assert not math.isfinite(scale) or scale<=0
            elif label=='nonfinite':assert any(not math.isfinite(x) for x in [r,*c,*lam])
            elif label=='overflow':
                try:oracle(r,c,lam,scale)
                except OverflowError:pass
                else:raise AssertionError('overflow certificate absent')
            else:raise AssertionError(('unproved invalidity',label))
            assert raw['error'];return dict(name=case['name'],status='exact_arithmetic_invalidity_proved',reason=label)
        assert raw['arithmetic_only'] and raw['rounding_mode_requested']==case['rounding']
        value=from_bits(raw['value_bits']);bound=from_bits(raw['bound_bits'])
        assert value==raw['normalized_value'] and bound==raw['roundoff_bound']
        allowed,exact,enclosure=check_balance(r,c,lam,scale,value,bound)
        assert raw['within_budget']==allowed and code==(0 if allowed else 2)
        return dict(name=case['name'],status='independent_exact_binary_product_sum_rounding_enclosure_passed',
            within_budget=allowed,exact_value=str(exact),roundoff_enclosure=str(enclosure),correct_rounding_bit_pattern=raw['value_bits'])
    if case.get('rejection') in ['precision_initial','precision_trial']:return precision_rejection(case,raw,code)
    assessed=inspect_many(case,raw,code)
    for t in raw.get('trace',[]):
        bounds=t['accumulation_roundoff_bounds'];assert len(bounds)==len(case['position'])
        assert all(0<=x<=1e-12 for x in bounds)
    if raw.get('converged'):assert not raw['precision_rejection']['present']
    return assessed
