"""Lossless fixed-pure95 independent-auditor state and explicit commit store.

Import ONLY after fixed backend selection. This module never evaluates a
material/solver, repairs history, selects a latest checkpoint, or grants a
physical qualification. A separately qualified reader owns integration.
"""
import fcntl,hashlib,json,math,os,re,struct,sys
from pathlib import Path
import decimal as dec
from fractions import Fraction

MAX_BYTES=128*2**20
AUDIT_MUTABLE={'maximum','trial_points','statistics','ever_active','slip_signs','row_slip_signs'}
AUDIT_FIELDS=AUDIT_MUTABLE|{'case','h','full','ref','geometry','mass','I','E','C','ids','cache','G','capture'}
STAT_COUNTS={'sliding_rows','strict_stick_rows','previously_active_separation_rows','nonzero_tangent0_rows','nonzero_tangent1_rows'}
STAT_DECIMALS={'Maxwell_loss','normal_loss','friction_loss','algorithmic_loss','positive_input_work','minimum_observed_gap'}
SIGNALS=('Clamped','InvalidOperation','DivisionByZero','Inexact','Rounded','Subnormal','Overflow','Underflow','FloatOperation')
ROUNDING=('ROUND_CEILING','ROUND_DOWN','ROUND_FLOOR','ROUND_HALF_DOWN','ROUND_HALF_EVEN','ROUND_HALF_UP','ROUND_UP','ROUND_05UP')

def backend():
    if dec is not sys.modules.get('_pydecimal') or sys.modules.get('decimal') is not dec or '_decimal' in sys.modules or dec.getcontext().prec!=95:
        raise ValueError('fresh fixed purePython95 backend required,never fallback')

def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=True,allow_nan=False).encode()
def digest(value):return hashlib.sha256(value).hexdigest()
def sha(value):
    if not isinstance(value,str) or not re.fullmatch('[0-9a-f]{64}',value):raise ValueError('SHA256 format')
    return value
def integer(value,minimum=0):
    if type(value) is not int or value<minimum:raise ValueError('integer domain')
    return value

def pack(value,depth=0):
    if depth>128:raise ValueError('typed state nesting')
    child=lambda x:pack(x,depth+1)
    if value is None:return ['none']
    if type(value) is bool:return ['bool',value]
    if type(value) is int:return ['int',str(value)]
    if type(value) is float:
        if not math.isfinite(value):raise ValueError('nonfinite binary scalar')
        return ['float64',struct.pack('<d',value).hex()]
    if type(value) is str:return ['str',value]
    if type(value) is bytes:return ['bytes',value.hex()]
    if isinstance(value,dec.Decimal):
        t=value.as_tuple()
        if value.is_nan():raise ValueError('NaN reference value')
        return ['decimal',t.sign,''.join(map(str,t.digits)),t.exponent]
    if type(value) is Fraction:return ['rational',str(value.numerator),str(value.denominator)]
    if type(value) in (list,tuple):return ['list' if type(value) is list else 'tuple',[child(x) for x in value]]
    if type(value) in (set,frozenset):return ['set' if type(value) is set else 'frozenset',sorted((child(x) for x in value),key=canonical)]
    if type(value) is dict:
        entries=sorted(([child(k),child(v)] for k,v in value.items()),key=lambda x:canonical(x[0]))
        return ['dict',entries]
    raise ValueError('unsupported reference state type:'+type(value).__name__)

def unpack(value,depth=0):
    if depth>128 or type(value) is not list or not value or type(value[0]) is not str:raise ValueError('typed state syntax')
    kind=value[0];child=lambda x:unpack(x,depth+1)
    if kind=='none' and len(value)==1:return None
    if len(value)<2:raise ValueError('typed scalar arity')
    if kind=='bool' and len(value)==2 and type(value[1]) is bool:return value[1]
    if kind=='int' and len(value)==2 and type(value[1]) is str:
        n=int(value[1])
        if str(n)!=value[1]:raise ValueError('noncanonical integer')
        return n
    if kind=='float64' and len(value)==2 and re.fullmatch('[0-9a-f]{16}',value[1]):
        x=struct.unpack('<d',bytes.fromhex(value[1]))[0]
        if not math.isfinite(x):raise ValueError('nonfinite binary scalar')
        return x
    if kind=='str' and len(value)==2 and type(value[1]) is str:return value[1]
    if kind=='bytes' and len(value)==2 and type(value[1]) is str:
        b=bytes.fromhex(value[1])
        if b.hex()!=value[1]:raise ValueError('noncanonical bytes')
        return b
    if kind=='decimal' and len(value)==4:
        _,sign,digits,exponent=value
        if type(sign) is not int or sign not in (0,1) or type(digits) is not str or not re.fullmatch('[0-9]+',digits):raise ValueError('Decimal tuple')
        if not (type(exponent) is int or exponent=='F'):raise ValueError('Decimal exponent')
        x=dec.Decimal((sign,tuple(map(int,digits)),exponent))
        if pack(x)!=value:raise ValueError('noncanonical Decimal tuple')
        return x
    if kind=='rational' and len(value)==3:
        x=Fraction(int(value[1]),int(value[2]))
        if pack(x)!=value:raise ValueError('noncanonical rational')
        return x
    if kind in ('list','tuple','set','frozenset') and len(value)==2 and type(value[1]) is list:
        values=list(map(child,value[1]));out={'list':list,'tuple':tuple,'set':set,'frozenset':frozenset}[kind](values)
        if pack(out)!=value:raise ValueError('noncanonical or duplicate container entries')
        return out
    if kind=='dict' and len(value)==2 and type(value[1]) is list:
        result={}
        for pair in value[1]:
            if type(pair) is not list or len(pair)!=2:raise ValueError('dict entry')
            k,v=map(child,pair)
            if k in result:raise ValueError('duplicate decoded key')
            result[k]=v
        if pack(result)!=value:raise ValueError('noncanonical dict ordering')
        return result
    raise ValueError('unsupported typed tag/arity')

def context_descriptor(context=None):
    backend();c=dec.getcontext() if context is None else context
    q={k:getattr(c,k) for k in ('prec','rounding','Emin','Emax','capitals','clamp')}
    for name in ('flags','traps'):q[name]={s:bool(getattr(c,name)[getattr(dec,s)]) for s in SIGNALS}
    return q

def stage_context(q,profile):
    backend()
    if set(q)!={'prec','rounding','Emin','Emax','capitals','clamp','flags','traps'}:raise ValueError('complete Decimal context required')
    if q['prec']!=95 or q['rounding'] not in ROUNDING:raise ValueError('reference precision/rounding')
    if {k:v for k,v in q.items() if k!='flags'}!={k:v for k,v in profile.items() if k!='flags'}:raise ValueError('frozen Decimal settings/traps changed')
    for name in ('flags','traps'):
        if set(q[name])!=set(SIGNALS) or any(type(v) is not bool for v in q[name].values()):raise ValueError('Decimal signal map')
    c=dec.Context(**{k:q[k] for k in ('prec','rounding','Emin','Emax','capitals','clamp')})
    for name in ('flags','traps'):
        for signal,value in q[name].items():getattr(c,name)[getattr(dec,signal)]=value
    return c

def walk(value,path=()):
    yield path,value
    if type(value) is dict:
        for k,v in value.items():
            yield from walk(k,path+('<key>',k))
            yield from walk(v,path+(k,))
    elif type(value) in (list,tuple,set,frozenset):
        for k,v in enumerate(value):yield from walk(v,path+(k,))

def validate(state,contract):
    backend()
    required={'schema','accepted_step','old','reference_viscous','reference_loss','friction_history','audit','capture','raw_cursor','decimal_context','cpu_completed'}
    if set(state)!=required or state['schema']!=1:raise ValueError('reference checkpoint schema/member inventory')
    step=integer(state['accepted_step']);total=integer(contract['total_steps'],1);n=integer(contract['nodes'],1);m=integer(contract['contacts'],1)
    if contract['schema']!=1 or n!=33 or m!=116 or total>98304:raise ValueError('unsupported bounded inspector scope')
    if step>total or len(contract['index'])!=total+1:raise ValueError('accepted audit cursor')
    entry=contract['index'][step]
    if state['raw_cursor']!=entry['raw_cursor'] or state['old']['step']!=step or float(state['old']['time']).hex()!=entry['time_hex']:raise ValueError('raw prefix/accepted history cursor mismatch')
    for name in ('position','velocity','green','viscous','dissipation','accepted_stress','density'):
        if len(state['old'][name])!=n:raise ValueError('native state shape')
    if len(state['old']['last_normal_impulse'])!=m:raise ValueError('normal dual shape')
    if digest(canonical(pack(state['old'])))!=entry['native_state_sha'] or digest(canonical(pack(state['friction_history'])))!=entry['friction_sha']:raise ValueError('accepted native observation differs from pinned prefix')
    z=state['reference_viscous'];loss=state['reference_loss']
    if len(z)!=n or len(loss)!=n or any(len(row)!=9 for row in z):raise ValueError('independent history shape')
    for row in z:
        if not all(isinstance(x,dec.Decimal) and x.is_finite() for x in row):raise ValueError('independent Decimal history required')
        scale=1+max(map(lambda x:abs(Fraction(x)),row));tol=Fraction(1,10**80)*scale
        if abs(sum((Fraction(row[k]) for k in (0,4,8)),Fraction(0)))>tol or any(abs(Fraction(row[3*a+b])-Fraction(row[3*b+a]))>tol for a in range(3) for b in range(3)):raise ValueError('invalid deviatoric/symmetric reference history')
    if any(not isinstance(x,dec.Decimal) or not x.is_finite() or Fraction(x)<-Fraction(1,10**80) for x in loss):raise ValueError('invalid independent cumulative loss')
    a=state['audit']
    if set(a)!=AUDIT_MUTABLE or a['trial_points']!=entry['trial_points']:raise ValueError('complete audit counters required')
    if type(a['maximum']) is not float or not math.isfinite(a['maximum']) or not 0<=a['maximum']<=1e-10:raise ValueError('saved accuracy gate')
    stats=a['statistics']
    if set(stats)!=STAT_COUNTS|STAT_DECIMALS or any(type(stats[k]) is not int or stats[k]<0 for k in STAT_COUNTS) or any(not isinstance(stats[k],dec.Decimal) for k in STAT_DECIMALS):raise ValueError('complete reference statistics inventory')
    if type(a['ever_active']) is not set or any(type(x) is not int or not 0<=x<m for x in a['ever_active']):raise ValueError('active contact history')
    if len(a['slip_signs'])!=2 or any(type(x) is not set or not x<={-1,1} for x in a['slip_signs']):raise ValueError('tangent sign history')
    for key,signs in a['row_slip_signs'].items():
        if type(key) is not tuple or len(key)!=2 or not 0<=key[0]<m or key[1] not in (0,1) or type(signs) is not set or not signs<={-1,1}:raise ValueError('per-contact tangent history')
    capture=state['capture']
    if set(capture)!={'score','selected'} or type(capture['score']) is not float or not math.isfinite(capture['score']):raise ValueError('negative-control sample inventory')
    if step==0:
        if capture['selected'] is not None or capture['score']!=-1.:raise ValueError('initial negative-control sample')
    else:
        selected=capture['selected']
        if type(selected) is not dict or set(selected)!={'case','header','old','z','loss','history','row','material_cache'} or capture['score']<0:raise ValueError('complete selected negative-control context required')
        selected_step=selected['row']['state']['step'];integer(selected_step,1)
        if selected_step>step or selected['old']['step']!=selected_step-1 or type(selected['material_cache']) is not dict:raise ValueError('selected negative sample clock/cache')
        if digest(canonical(pack(selected['case'])))!=contract['case_sha'] or digest(canonical(pack(selected['header'])))!=contract['header_sha']:raise ValueError('negative sample static binding')
        selected_index=contract['index'][selected_step]
        if digest(canonical(pack(selected['row'])))!=selected_index['native_row_sha'] or digest(canonical(pack(selected['old'])))!=contract['index'][selected_step-1]['native_state_sha'] or digest(canonical(pack(selected['history'])))!=contract['index'][selected_step-1]['friction_sha']:raise ValueError('negative sample raw context binding')
    if type(state['cpu_completed']) is not float or not math.isfinite(state['cpu_completed']) or not 0<=state['cpu_completed']<=contract['CPU_cap']:raise ValueError('saved CPU diagnostic')
    stage_context(state['decimal_context'],contract['decimal_profile'])
    for path,value in walk(state):
        if isinstance(value,dec.Decimal) and not value.is_finite():
            if not (step==0 and path==('audit','statistics','minimum_observed_gap') and value==dec.Decimal('Infinity')):raise ValueError('nonfinite reference state')
    return state

def encode(state,contract):
    validate(state,contract);body=pack(state)
    data=canonical(dict(schema=1,contract_sha256=digest(canonical(contract)),payload=body,payload_sha256=digest(canonical(body))))
    if len(data)>MAX_BYTES:raise ValueError('reference state byte budget')
    return data

def decode(data,contract):
    backend()
    if len(data)>MAX_BYTES:raise ValueError('reference state byte budget')
    envelope=json.loads(data)
    if canonical(envelope)!=data or set(envelope)!={'schema','contract_sha256','payload','payload_sha256'} or envelope['schema']!=1:raise ValueError('canonical reference envelope')
    if envelope['contract_sha256']!=digest(canonical(contract)) or envelope['payload_sha256']!=digest(canonical(envelope['payload'])):raise ValueError('reference contract/payload binding')
    state=unpack(envelope['payload']);validate(state,contract)
    if encode(state,contract)!=data:raise ValueError('lossless canonical roundtrip')
    return state

def snapshot(audit,old,z,loss,history,capture,cursor,contract,cpu):
    if set(audit.__dict__)!=AUDIT_FIELDS or audit.full is not True or audit.cache!={}:raise ValueError('unsupported audit inventory or non-boundary cache')
    if digest(canonical(pack(audit.case)))!=contract['case_sha'] or digest(canonical(pack(audit.h)))!=contract['header_sha']:raise ValueError('actual inspector static binding')
    state=dict(schema=1,accepted_step=old['step'],old=old,reference_viscous=z,reference_loss=loss,friction_history=history,
        audit={k:getattr(audit,k) for k in AUDIT_MUTABLE},capture=dict(score=capture.score,selected=capture.selected),raw_cursor=cursor,
        decimal_context=context_descriptor(),cpu_completed=cpu)
    # A private lossless image,not aliases to mutable live arrays.
    return decode(encode(state,contract),contract)

def install(state,audit,capture,contract):
    validate(state,contract)
    if set(audit.__dict__)!=AUDIT_FIELDS or audit.full is not True or audit.cache!={}:raise ValueError('restore destination inventory')
    if digest(canonical(pack(audit.case)))!=contract['case_sha'] or digest(canonical(pack(audit.h)))!=contract['header_sha']:raise ValueError('restore destination static binding')
    staged=decode(encode(state,contract),contract);context=stage_context(staged['decimal_context'],contract['decimal_profile'])
    audit.__dict__.update(staged['audit']);capture.__dict__.update(staged['capture']);audit.capture=capture
    dec.setcontext(context)
    return staged['old'],staged['reference_viscous'],staged['reference_loss'],staged['friction_history']

def file_bytes(path):
    p=Path(path)
    if p.is_symlink() or not p.is_file():raise ValueError('regular checkpoint file required')
    if p.stat().st_size>(MAX_BYTES if p.name=='state.json' else 16*2**10):raise ValueError('checkpoint file byte budget')
    return p.read_bytes()
def sync_directory(path):
    fd=os.open(path,os.O_RDONLY|os.O_DIRECTORY)
    try:os.fsync(fd)
    finally:os.close(fd)
def exclusive_file(path,data):
    with Path(path).open('xb') as f:f.write(data);f.flush();os.fsync(f.fileno())

class Store:
    """Expected-HEAD CAS and explicit committed ancestor;no latest fallback."""
    def __init__(self,path,contract):
        self.path=Path(path);self.contract=contract;self.binding=digest(canonical(contract))
        if self.path.is_symlink() or not self.path.is_dir():raise ValueError('existing regular store directory required')
    def head(self):
        if (self.path/'HEAD').is_symlink():raise ValueError('HEAD alias rejected')
        return digest(file_bytes(self.path/'HEAD')) if (self.path/'HEAD').exists() else None
    def _chain(self,expected):
        if self.head()!=expected:raise ValueError('expected HEAD mismatch')
        if expected is None:return []
        data=file_bytes(self.path/'HEAD');chain=[]
        while data:
            q=json.loads(data)
            if canonical(q)!=data or set(q)!={'schema','contract_sha256','sequence','leaf','snapshot_sha256','parent_leaf','parent_sha256'} or q['schema']!=1 or q['contract_sha256']!=self.binding:raise ValueError('commit metadata/schema')
            leaf=q['leaf']
            if not isinstance(leaf,str) or not re.fullmatch('g_[A-Za-z0-9_-]{1,80}',leaf):raise ValueError('explicit branch leaf')
            if (self.path/leaf).is_symlink() or not (self.path/leaf).is_dir():raise ValueError('regular checkpoint branch required')
            archive=file_bytes(self.path/leaf/'commit.json')
            if archive!=data:raise ValueError('commit chain bytes')
            integer(q['sequence'],1);sha(q['snapshot_sha256']);chain.append(q)
            if len(chain)>98304+128:raise ValueError('commit history budget')
            if q['parent_leaf'] is None:
                if q['sequence']!=1 or q['parent_sha256'] is not None:raise ValueError('initial commit link')
                break
            parent=q['parent_leaf']
            if not isinstance(parent,str) or not re.fullmatch('g_[A-Za-z0-9_-]{1,80}',parent):raise ValueError('parent branch leaf')
            if (self.path/parent).is_symlink():raise ValueError('parent branch alias')
            data=file_bytes(self.path/parent/'commit.json')
            if digest(data)!=sha(q['parent_sha256']) or json.loads(data)['sequence']!=q['sequence']-1:raise ValueError('commit ancestry')
        return chain
    def read(self,leaf,expected_head):
        chain=self._chain(expected_head);q=next((x for x in chain if x['leaf']==leaf),None)
        if q is None:raise ValueError('explicit selection is not a committed ancestor')
        data=file_bytes(self.path/leaf/'state.json')
        if digest(data)!=q['snapshot_sha256']:raise ValueError('selected snapshot SHA mismatch')
        return decode(data,self.contract),data
    def publish(self,state,label,expected_head,hook=lambda stage:None):
        backend()
        if not re.fullmatch('[A-Za-z0-9_-]{1,80}',label):raise ValueError('new branch label')
        data=encode(state,self.contract);leaf='g_'+label
        lock=self.path/'SESSION.lock'
        if lock.is_symlink() or (lock.exists() and not lock.is_file()):raise ValueError('regular session lock required')
        with lock.open('a+b') as owner:
            fcntl.flock(owner.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
            chain=self._chain(expected_head);parent=chain[0] if chain else None
            target=self.path/leaf;target.mkdir()
            q=dict(schema=1,contract_sha256=self.binding,sequence=1 if parent is None else parent['sequence']+1,leaf=leaf,
                snapshot_sha256=digest(data),parent_leaf=None if parent is None else parent['leaf'],parent_sha256=expected_head)
            commit=canonical(q);exclusive_file(target/'state.json',data);exclusive_file(target/'commit.json',commit);sync_directory(target)
            hook('before_HEAD')
            temp=self.path/('.HEAD_'+label+'.tmp');exclusive_file(temp,commit);os.replace(temp,self.path/'HEAD');sync_directory(self.path)
            hook('after_HEAD')
        return dict(leaf=leaf,HEAD_sha256=digest(commit),snapshot_sha256=digest(data),sequence=q['sequence'])
