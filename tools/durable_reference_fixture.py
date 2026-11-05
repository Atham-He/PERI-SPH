"""Explicitly synthetic stream/owner for coordinator controls,zero PD math."""
import copy
import json
import decimal as dec
from reference_checkpoint_fixture import initial,advance,inspector,native,friction,row
from peri_reference_checkpoint import snapshot
from peri_durable_reference_stream import build_index,state_contract
from peri_physical_time_friction_audit import PhysicalTimeFrictionInspector


def fixture(root):
    contract,state,case,header=initial()
    clock=dict(denominator=1000000,nominal_ticks=5000,normal_event_ticks=[10000,20000,40000],tangent_event_ticks=[5000,10000,15000])
    case.update(physical_time=clock,protocol_kind=2,dt=.005)
    header.update(type='header',format=2,steps_requested=8,time_denominator=1000000,nominal_ticks=5000,
        normal_event_ticks=clock['normal_event_ticks'],tangent_event_ticks=clock['tangent_event_ticks'],
        initial_state=native(0),initial_friction_history=friction())
    meta=dict(type='identity',pid=123,ppid=12,start_ticks='456',boot_id='00000000-0000-0000-0000-000000000000',argv=['SYNTHETIC_NOT_NATIVE'])
    rows=[meta,header]
    for step in range(1,9):
        r=row(step);old=float((step-1)*5000/1000000);new=float(step*5000/1000000)
        r.update(type='step',committed=True,solve=dict(trace=[dict(iteration=0)]),
            interval=dict(old_tick=(step-1)*5000,new_tick=step*5000,time_denominator=1000000,old_time=old,new_time=new,actual_dt=new-old))
        r['state']['time']=new;rows.append(r)
    rows.append(dict(type='end',completed=True,accepted_steps=8,last_accepted_time=.04))
    raw=root/'synthetic.raw';raw.write_bytes(b''.join(json.dumps(r,allow_nan=False).encode()+b'\n' for r in rows))
    from peri_reference_checkpoint import digest
    record=dict(path=str(raw),size=raw.stat().st_size,sha256=digest(raw.read_bytes()))
    index=build_index(raw,case,record,{k:v for k,v in meta.items() if k!='type'})
    contract=state_contract(index,digest(b'SYNTHETIC_SOFTWARE'),10000.)
    return raw,index,contract,rows


class ToyInspector(PhysicalTimeFrictionInspector):
    """Actual wrapper type/subclass ABI,but overridden SYNTHETIC state steps."""
    contract=None
    calls=[]
    def __init__(self,case,header,full):
        _,state,_,_=initial();owner,_=inspector(state,case,header)
        self.__dict__.update(copy.deepcopy(owner.__dict__))

    def initial(self):
        type(self).calls.append('initial')
        _,state,_,_=initial()
        return state['old'],state['reference_viscous'],state['reference_loss'],state['friction_history']

    def step(self,r,old,z,loss,history):
        type(self).calls.append(old['step']+1)
        s=snapshot(self,old,z,loss,history,self.capture,self.contract['index'][old['step']]['raw_cursor'],self.contract,0.)
        new=advance(s,self.contract,self.case,self.h)
        new['capture']['selected']['row']=copy.deepcopy(r)
        new['old']=copy.deepcopy(r['state']);new['friction_history']=copy.deepcopy(r['friction_history'])
        self.__dict__.update(new['audit']);self.capture.__dict__.update(new['capture'])
        return new['old'],new['reference_viscous'],new['reference_loss'],new['friction_history']
