"""Binary accepted-frontier index and durable wrapper of unchanged PD audit.

Import only after fixed pure95 selection. Indexing verifies bytes and clocks,
NOT physical equations. Actual full auditing requires a separately frozen plan,
launch reservation and full original negative/final controls.
"""
import copy
import hashlib
import json
from pathlib import Path
import time
import decimal as dec
from peri_reference_checkpoint import (Store, backend, canonical, digest, pack,
    snapshot, install, context_descriptor, encode, validate)
from peri_physical_time_friction_audit import PhysicalTimeFrictionInspector, schedule_grid
from peri_physical_friction_negative import OneStepCapture


def typed_sha(value):
    return digest(canonical(pack(value)))


def entry(step, old, history, row, offset, prefix, line, points):
    return dict(time_hex=float(old['time']).hex(), trial_points=points,
        native_state_sha=typed_sha(old), friction_sha=typed_sha(history), native_row_sha=typed_sha(row),
        raw_cursor=dict(offset=offset, prefix_sha256=prefix, last_line_sha256=digest(line)))


def binary_row(stream, rolling):
    line=stream.readline()
    if not line or len(line)>128*2**20 or not line.endswith(b'\n'):
        raise ValueError('complete bounded newline record required')
    value=json.loads(line); rolling.update(line)
    if type(value) is not dict:raise ValueError('raw row inventory')
    return value,line


def build_index(path,case,expected_raw,expected_identity):
    """Return a byte/clock binding only;never call Inspector/material methods."""
    backend(); grid=schedule_grid(case); total=len(grid)-1; den=case['physical_time']['denominator']
    if total!=case['steps']:raise ValueError('requested/grid steps')
    path=Path(path)
    if path.is_symlink() or not path.is_file() or path.stat().st_size!=expected_raw['size']:
        raise ValueError('pinned regular raw file/size')
    rolling=hashlib.sha256(); entries=[]; points=0
    with path.open('rb') as stream:
        meta,_=binary_row(stream,rolling)
        if meta.get('type')!='identity' or any(meta.get(k)!=v for k,v in expected_identity.items()):
            raise ValueError('full native birth/argv binding')
        header,line=binary_row(stream,rolling)
        if header.get('type')!='header' or header.get('format')!=2 or header['steps_requested']!=total:
            raise ValueError('explicit physical-time header')
        for h,c in [('time_denominator','denominator'),('nominal_ticks','nominal_ticks'),
                    ('normal_event_ticks','normal_event_ticks'),('tangent_event_ticks','tangent_event_ticks')]:
            if header[h]!=case['physical_time'][c]:raise ValueError('header physical clock binding')
        initial=header['initial_state']
        if initial['step']!=0 or initial['time']!=0.:raise ValueError('initial raw clock')
        entries.append(entry(0,initial,header['initial_friction_history'],header,stream.tell(),rolling.hexdigest(),line,0))
        for step in range(1,total+1):
            row,line=binary_row(stream,rolling); state=row.get('state',{}); interval=row.get('interval',{})
            if row.get('type')!='step' or state.get('step')!=step or row.get('committed') is not True:
                raise ValueError('complete accepted sequence')
            old_time=float(grid[step-1]/den); new_time=float(grid[step]/den)
            # Division of small exact integer ticks gives correctly rounded
            # binary endpoints,the same endpoints as the rational protocol.
            if (interval.get('old_tick'),interval.get('new_tick'),interval.get('time_denominator'))!=(grid[step-1],grid[step],den):
                raise ValueError('accepted raw rational cursor')
            if (float(interval['old_time']).hex(),float(interval['new_time']).hex(),float(state['time']).hex())!=(old_time.hex(),new_time.hex(),new_time.hex()):
                raise ValueError('accepted raw binary time')
            if float(interval['actual_dt']).hex()!=(new_time-old_time).hex() or float(old_time+interval['actual_dt']).hex()!=new_time.hex():
                raise ValueError('accepted actual dt')
            trace=row['solve']['trace']
            if not trace or [t['iteration'] for t in trace]!=list(range(len(trace))):raise ValueError('complete native Newton trace inventory')
            points+=33*len(trace)
            entries.append(entry(step,state,row['friction_history'],row,stream.tell(),rolling.hexdigest(),line,points))
        end,line=binary_row(stream,rolling)
        if end.get('type')!='end' or end.get('completed') is not True or end.get('accepted_steps')!=total or float(end['last_accepted_time']).hex()!=entries[-1]['time_hex']:
            raise ValueError('complete accepted endpoint')
        tail=stream.read();rolling.update(tail)
        if tail.strip():raise ValueError('hidden raw tail')
    if rolling.hexdigest()!=expected_raw['sha256']:raise ValueError('whole raw SHA')
    return dict(schema=1,scope='binary_raw_sequence_clock_index_NOT_equation_review',raw=expected_raw,
        identity=meta,header=header,case=case,index=entries,end=end,trial_points=points,
        final_prefix_sha256=entries[-1]['raw_cursor']['prefix_sha256'],whole_raw_sha256=rolling.hexdigest(),
        source_case_sha256=typed_sha(case),header_sha256=typed_sha(header),equation_review_qualified=False)


def state_contract(index,software_sha,CPU_cap):
    backend()
    return dict(schema=1,nodes=33,contacts=116,total_steps=len(index['index'])-1,index=index['index'],
        decimal_profile=context_descriptor(),CPU_cap=CPU_cap,case_sha=index['source_case_sha256'],header_sha=index['header_sha256'],
        raw_sha256=index['raw']['sha256'],software_sha256=software_sha)


def verify_index_frontier(stream,index,step):
    cursor=index['index'][step]['raw_cursor'];left=cursor['offset'];rolling=hashlib.sha256()
    stream.seek(0)
    while left:
        block=stream.read(min(left,2**20))
        if not block:raise ValueError('raw shorter than committed frontier')
        rolling.update(block);left-=len(block)
    if rolling.hexdigest()!=cursor['prefix_sha256']:raise ValueError('selected raw prefix SHA')
    return rolling


def mathematical_image(state,contract):
    """Explicit comparison projection;original measured bytes stay preserved.

    cpu_completed is the ONLY excluded diagnostic. The original snapshot bytes
    and original CPU value are recorded separately;this is never described as
    identity of the original entire checkpoint files.
    """
    validate(state,contract)
    value={k:v for k,v in state.items() if k!='cpu_completed'}
    return dict(excluded_fields=['cpu_completed'],measured_CPU_diagnostic=state['cpu_completed'],
        original_state_sha256=digest(encode(state,contract)),mathematical_state_sha256=typed_sha(value))


class DurableSession:
    """Actual Inspector owner at fully verified accepted boundaries only.

    Stores independently computed z/loss/context/statistics/capture,never
    installs the native viscous/dissipation fields as the independent state.
    Factory injection is for explicitly synthetic structural controls only.
    """
    def __init__(self,path,index,contract,store,expected_head,selected_leaf=None,
                 inspector_type=PhysicalTimeFrictionInspector,capture_type=OneStepCapture):
        backend();self.path=Path(path);self.index=index;self.contract=contract;self.store=store
        self.expected_head=expected_head;self.case=index['case'];self.header=index['header']
        self.audit=inspector_type(self.case,self.header,True);self.capture=capture_type();self.audit.capture=self.capture
        self.restored_CPU_diagnostic=None
        if selected_leaf is None:
            if expected_head is not None:raise ValueError('explicit checkpoint selection required,never automatic latest')
            self.old,self.z,self.loss,self.history=self.audit.initial()
        else:
            state,_=store.read(selected_leaf,expected_head)
            self.old,self.z,self.loss,self.history=install(state,self.audit,self.capture,contract)
            self.restored_CPU_diagnostic=state['cpu_completed']
        self.step=self.old['step'];self.stream=self.path.open('rb')
        try:self.rolling=verify_index_frontier(self.stream,index,self.step)
        except BaseException:
            self.stream.close();raise

    def close(self):
        self.stream.close()

    def image(self,cpu=None):
        return snapshot(self.audit,self.old,self.z,self.loss,self.history,self.capture,
            self.index['index'][self.step]['raw_cursor'],self.contract,time.process_time() if cpu is None else cpu)

    def commit(self,label,cpu=None,hook=lambda point:None):
        state=self.image(cpu)
        result=self.store.publish(state,label,self.expected_head,hook);self.expected_head=result['HEAD_sha256']
        return dict(**result,comparison=mathematical_image(state,self.contract),accepted_step=self.step)

    def advance(self):
        if self.step>=self.contract['total_steps']:raise ValueError('cannot advance beyond frozen endpoint')
        row,line=binary_row(self.stream,self.rolling);step=self.step+1;wanted=self.index['index'][step]
        if (self.stream.tell()!=wanted['raw_cursor']['offset'] or self.rolling.hexdigest()!=wanted['raw_cursor']['prefix_sha256']
            or digest(line)!=wanted['raw_cursor']['last_line_sha256'] or typed_sha(row)!=wanted['native_row_sha']):
            raise ValueError('new raw line/prefix differs from externally frozen index')
        self.old,self.z,self.loss,self.history=self.audit.step(row,self.old,self.z,self.loss,self.history)
        self.step=step
        # Full snapshot validation also checks independent history and cumulative
        # original material-point count against the frozen metadata index.
        return self.image()

    def finish_raw(self):
        if self.step!=self.contract['total_steps']:raise ValueError('partial audit cannot finalize')
        row,_=binary_row(self.stream,self.rolling);tail=self.stream.read();self.rolling.update(tail)
        if row!=self.index['end'] or tail.strip() or self.rolling.hexdigest()!=self.index['whole_raw_sha256']:
            raise ValueError('complete raw endpoint/tail/SHA')
        return row
