"""Synthetic full-shaped auditor memory;no constitutive or solver calls."""
import copy
import decimal as dec
from fractions import Fraction
from types import SimpleNamespace
from peri_reference_checkpoint import AUDIT_FIELDS,STAT_COUNTS,STAT_DECIMALS,canonical,pack,digest,context_descriptor,stage_context

def native(step,n=33,m=116):
    return dict(time=step*.005,step=step,position=[[0.,0.,0.] for _ in range(n)],velocity=[[step*.001,0.,0.] for _ in range(n)],
        green=[[0.]*9 for _ in range(n)],viscous=[[0.]*9 for _ in range(n)],dissipation=[0.]*n,
        accepted_stress=[[0.]*9 for _ in range(n)],density=[1.]*n,last_normal_impulse=[0.]*m,last_dt=0. if step==0 else .005,last_KKT=0.,last_newton_iterations=0)
def friction(m=116):
    return dict(last_tangent_impulse=[[0.,0.] for _ in range(m)],relative_tangent_travel_sum=[[0.,0.] for _ in range(m)],
        relative_tangent_travel_correction=[[0.,0.] for _ in range(m)],ledger_sum=[0.]*12,ledger_correction=[0.]*12,
        ledger_order='lower_Ixyz,upper_Ixyz,lower_work,upper_work,Maxwell,normal,friction,algorithmic')
def row(step):return dict(type='synthetic_accepted_metadata_NOT_PHYSICAL',state=native(step),friction_history=friction())

def initial(total=8):
    c=dec.getcontext().copy();c.prec=95;c.clear_flags();profile=context_descriptor(c)
    case=dict(name='synthetic_memory_only',steps=total);header=dict(nodes=33,contacts=116,mode='synthetic_not_PD')
    index=[]
    for step in range(total+1):
        index.append(dict(time_hex=float(step*.005).hex(),trial_points=33*step,
            native_state_sha=digest(canonical(pack(native(step)))),friction_sha=digest(canonical(pack(friction()))),native_row_sha=digest(canonical(pack(row(step)))),
            raw_cursor=dict(offset=100+step*1000,prefix_sha256=digest(('synthetic_prefix_'+str(step)).encode()),last_line_sha256=digest(('synthetic_line_'+str(step)).encode()))))
    contract=dict(schema=1,nodes=33,contacts=116,total_steps=total,index=index,decimal_profile=profile,CPU_cap=300.,
        case_sha=digest(canonical(pack(case))),header_sha=digest(canonical(pack(header))),
        raw_sha256=digest(b'SYNTHETIC_NEVER_PHYSICAL'),software_sha256=digest(b'SYNTHETIC_EXPLICIT_SOFTWARE'))
    statistics={k:0 for k in STAT_COUNTS};statistics.update({k:dec.Decimal(0) for k in STAT_DECIMALS});statistics['minimum_observed_gap']=dec.Decimal('Infinity')
    state=dict(schema=1,accepted_step=0,old=native(0),reference_viscous=[[dec.Decimal(0)]*9 for _ in range(33)],reference_loss=[dec.Decimal(0)]*33,
        friction_history=friction(),audit=dict(maximum=0.,trial_points=0,statistics=statistics,ever_active=set(),slip_signs=[set(),set()],row_slip_signs={}),
        capture=dict(score=-1.,selected=None),raw_cursor=index[0]['raw_cursor'],decimal_context=profile,cpu_completed=0.)
    return contract,state,case,header

def advance(previous,contract,case,header):
    dec.setcontext(stage_context(previous['decimal_context'],contract['decimal_profile']));state=copy.deepcopy(previous);step=previous['accepted_step']+1
    increment=dec.Decimal(1)/dec.Decimal(7)
    for i,history in enumerate(state['reference_viscous']):
        a=history[0]+increment;history[:]=[a,dec.Decimal('-0E-95'),dec.Decimal(0),dec.Decimal(0),-a/2,dec.Decimal(0),dec.Decimal(0),dec.Decimal(0),-a/2]
        state['reference_loss'][i]+=dec.Decimal('1E-90')
    state.update(accepted_step=step,old=native(step),friction_history=friction(),raw_cursor=contract['index'][step]['raw_cursor'],cpu_completed=step*.01)
    audit=state['audit'];audit['maximum']=1e-15;audit['trial_points']=33*step
    audit['ever_active'].add(step%116);audit['slip_signs'][0].add(1 if step%2 else -1);audit['slip_signs'][1].add(-1)
    audit['row_slip_signs'].setdefault((step%116,step%2),set()).add(1 if step%2 else -1)
    for k in STAT_COUNTS:audit['statistics'][k]+=1
    for k in STAT_DECIMALS-{'minimum_observed_gap'}:audit['statistics'][k]+=increment
    audit['statistics']['minimum_observed_gap']=-dec.Decimal(step)/dec.Decimal(10**6)
    cache={(step-1,(0.,-0.),'.005',True):dict(point_T=[[increment,dec.Decimal('-0E-95')]],rational=Fraction(1,3),typed_bytes=b'\x00\xff')}
    state['capture']=dict(score=step*1e-9,selected=dict(case=case,header=header,old=previous['old'],z=previous['reference_viscous'],
        loss=previous['reference_loss'],history=previous['friction_history'],row=row(step),material_cache=cache))
    state['decimal_context']=context_descriptor();return state

def inspector(state,case,header):
    fields={k:None for k in AUDIT_FIELDS};fields.update(state['audit']);fields.update(case=case,h=header,full=True,cache={},capture=None)
    return SimpleNamespace(**fields),SimpleNamespace(**state['capture'])
