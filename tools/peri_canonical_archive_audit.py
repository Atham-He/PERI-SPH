"""Independent struct/hashlib byte-format oracle; no PD restart claim."""
import hashlib,math,re,struct


def parse(data):
    assert 72<=len(data)<=16*1024*1024
    assert data[:8]==b'PERIAR01' and data[-8:]==b'PERIEND1'
    assert hashlib.sha256(data[:-40]).digest()==data[-40:-8]
    offset=8;limit=len(data)-40
    def take(n):
        nonlocal offset
        assert 0<=n<=limit-offset
        value=data[offset:offset+n];offset+=n;return value
    def u64():return struct.unpack('<Q',take(8))[0]
    assert u64()==1;generation=u64();assert generation>0
    count=u64();assert count<=256;result={};previous=''
    for _ in range(count):
        n=u64();assert 0<n<=128
        name=take(n).decode('ascii');assert re.fullmatch(r'[A-Za-z0-9_./-]+',name)
        assert name>previous;previous=name
        tag,rank=take(2);assert tag in [1,2,3] and rank<=3
        shape=[u64() for _ in range(rank)]
        assert all(n<=1024*1024 for n in shape) and math.prod(shape)<=1024*1024
        size=u64();width=1 if tag==1 else 8
        assert size==math.prod(shape)*width
        raw=take(size)
        if tag==3:assert all(math.isfinite(x) for x in struct.unpack('<'+'d'*(size//8),raw))
        result[name]=dict(tag=tag,shape=shape,raw=raw)
    assert offset==limit
    return generation,result


def assess_fixture(data,native):
    generation,fields=parse(data)
    assert generation==7 and set(fields)=={'accepted/step','accepted/values','reference/mask'}
    assert fields['accepted/step']==dict(tag=2,shape=[],raw=struct.pack('<Q',19))
    assert fields['reference/mask']==dict(tag=1,shape=[4],raw=b'\1\0\1\1')
    values=[0.,-0.,.1,float.fromhex('0x0.0000000000001p-1022'),
            float.fromhex('0x1.fffffffffffffp+1023'),float.fromhex('0x1p-1022')]
    assert fields['accepted/values']==dict(tag=3,shape=[2,3],raw=struct.pack('<6d',*values))
    assert native['fixture_bytes']==len(data) and native['physical_steps']==0 and not native['full_restart_qualified']
    assert native['large_SHA256']==hashlib.sha256(bytes((17*i+3)&255 for i in range(1024*1024))).hexdigest()
    return dict(status='canonical_byte_fixture_independently_matched',scalar_bit_patterns=6,
        signed_zero_preserved=True,subnormal_and_extreme_finite_preserved=True,
        kernel_SHA256_matches_independent_hashlib=True,physical_steps=0,full_restart_qualified=False)
