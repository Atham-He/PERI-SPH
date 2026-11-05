"""Incremental exact bytes of canonical(pack(value)),without the full tree.

Distinct implementation for new admissions. Does not replace old frozen pack
or alter numerical values. Dictionary keys/set entries are sorted by their
canonical typed bytes,exactly as the original format specifies.
"""
import hashlib
import json
import math
import struct
from fractions import Fraction
import decimal as dec


def token(value):
    return json.dumps(value,ensure_ascii=True,sort_keys=True,separators=(',',':'),allow_nan=False).encode()


def typed_chunks(value,depth=0):
    if depth>128:raise ValueError('typed state nesting')
    kind=type(value)
    if value is None:yield b'["none"]';return
    if kind is bool:yield b'["bool",'+token(value)+b']';return
    if kind is int:yield b'["int",'+token(str(value))+b']';return
    if kind is float:
        if not math.isfinite(value):raise ValueError('nonfinite binary scalar')
        yield b'["float64",'+token(struct.pack('<d',value).hex())+b']';return
    if kind is str:yield b'["str",'+token(value)+b']';return
    if kind is bytes:yield b'["bytes",'+token(value.hex())+b']';return
    if isinstance(value,dec.Decimal):
        q=value.as_tuple()
        if value.is_nan():raise ValueError('NaN reference value')
        yield b'["decimal",'+token(q.sign)+b','+token(''.join(map(str,q.digits)))+b','+token(q.exponent)+b']';return
    if kind is Fraction:
        yield b'["rational",'+token(str(value.numerator))+b','+token(str(value.denominator))+b']';return
    if kind in (list,tuple):
        yield b'["'+(b'list' if kind is list else b'tuple')+b'",['
        for i,item in enumerate(value):
            if i:yield b','
            yield from typed_chunks(item,depth+1)
        yield b']]';return
    if kind in (set,frozenset):
        entries=sorted(b''.join(typed_chunks(item,depth+1)) for item in value)
        yield b'["'+(b'set' if kind is set else b'frozenset')+b'",['
        for i,item in enumerate(entries):
            if i:yield b','
            yield item
        yield b']]';return
    if kind is dict:
        entries=sorted((b''.join(typed_chunks(k,depth+1)),k) for k in value)
        yield b'["dict",['
        for i,(encoded,key) in enumerate(entries):
            if i:yield b','
            yield b'[';yield encoded;yield b','
            yield from typed_chunks(value[key],depth+1)
            yield b']'
        yield b']]';return
    raise ValueError('unsupported reference state type:'+kind.__name__)


def stream_typed(value,sink,byte_limit=2**31,buffer_bytes=65536):
    if type(byte_limit) is not int or byte_limit<1 or type(buffer_bytes) is not int or buffer_bytes<1:
        raise ValueError('explicit streaming byte budget')
    pending=bytearray();total=0
    for chunk in typed_chunks(value):
        total+=len(chunk)
        if total>byte_limit:raise ValueError('typed stream byte budget')
        pending.extend(chunk)
        if len(pending)>=buffer_bytes:sink(bytes(pending));pending.clear()
    if pending:sink(bytes(pending))
    return total


def streamed_typed_sha(value,byte_limit=2**31):
    h=hashlib.sha256();stream_typed(value,h.update,byte_limit)
    return h.hexdigest()
