"""Distinct byte-equivalent typed-hash index engine with durable stage hooks.

Does not adopt a partial old index or change numerical values/packet coverage.
Full original row typed bytes are hashed incrementally,not constructed as a
second complete temporary Python tree. Root cause of the old SIGSEGV unknown.
"""
import hashlib
from pathlib import Path
from peri_reference_checkpoint import backend,digest
from peri_durable_reference_stream import binary_row
from peri_physical_time_friction_audit import schedule_grid
from peri_streaming_typed_json import streamed_typed_sha


def indexed_entry(step,old,history,row,offset,prefix,line,points,observe):
    observe('before_state_hash',step)
    state=streamed_typed_sha(old);friction=streamed_typed_sha(history)
    observe('before_complete_row_typed_hash',step)
    complete=streamed_typed_sha(row)
    observe('after_complete_row_typed_hash',step)
    return dict(time_hex=float(old['time']).hex(),trial_points=points,native_state_sha=state,friction_sha=friction,
        native_row_sha=complete,raw_cursor=dict(offset=offset,prefix_sha256=prefix,last_line_sha256=digest(line)))


def build_streamed_index(path,case,expected_raw,expected_identity,observe=lambda stage,step:None,accept=lambda step,value:None):
    backend();grid=schedule_grid(case);total=len(grid)-1;den=case['physical_time']['denominator'];path=Path(path)
    if total!=case['steps']:raise ValueError('requested/grid steps')
    if path.is_symlink() or not path.is_file() or path.stat().st_size!=expected_raw['size']:raise ValueError('pinned regular raw file/size')
    rolling=hashlib.sha256();entries=[];points=0
    with path.open('rb') as stream:
        observe('before_identity_parse',0);meta,_=binary_row(stream,rolling)
        if meta.get('type')!='identity' or any(meta.get(k)!=v for k,v in expected_identity.items()):raise ValueError('full native birth/argv binding')
        observe('before_header_parse',0);header,line=binary_row(stream,rolling)
        if header.get('type')!='header' or header.get('format')!=2 or header['steps_requested']!=total:raise ValueError('explicit physical-time header')
        for h,c in [('time_denominator','denominator'),('nominal_ticks','nominal_ticks'),('normal_event_ticks','normal_event_ticks'),('tangent_event_ticks','tangent_event_ticks')]:
            if header[h]!=case['physical_time'][c]:raise ValueError('header physical clock binding')
        initial=header['initial_state']
        if initial['step']!=0 or initial['time']!=0.:raise ValueError('initial raw clock')
        entries.append(indexed_entry(0,initial,header['initial_friction_history'],header,stream.tell(),rolling.hexdigest(),line,0,observe));accept(0,entries[-1])
        for step in range(1,total+1):
            observe('before_step_parse',step);row,line=binary_row(stream,rolling);observe('after_step_parse',step)
            state=row.get('state',{});interval=row.get('interval',{})
            if row.get('type')!='step' or state.get('step')!=step or row.get('committed') is not True:raise ValueError('complete accepted sequence')
            a=float(grid[step-1]/den);b=float(grid[step]/den)
            if (interval.get('old_tick'),interval.get('new_tick'),interval.get('time_denominator'))!=(grid[step-1],grid[step],den):raise ValueError('accepted raw rational cursor')
            if (float(interval['old_time']).hex(),float(interval['new_time']).hex(),float(state['time']).hex())!=(a.hex(),b.hex(),b.hex()):raise ValueError('accepted raw binary time')
            if float(interval['actual_dt']).hex()!=(b-a).hex() or float(a+interval['actual_dt']).hex()!=b.hex():raise ValueError('accepted actual dt')
            trace=row['solve']['trace']
            if not trace or [t['iteration'] for t in trace]!=list(range(len(trace))):raise ValueError('complete native Newton trace inventory')
            points+=33*len(trace)
            entries.append(indexed_entry(step,state,row['friction_history'],row,stream.tell(),rolling.hexdigest(),line,points,observe));accept(step,entries[-1])
            observe('before_row_release',step);del row,state,interval,trace,line;observe('after_row_release',step)
        observe('before_end_parse',total);end,_=binary_row(stream,rolling)
        if end.get('type')!='end' or end.get('completed') is not True or end.get('accepted_steps')!=total or float(end['last_accepted_time']).hex()!=entries[-1]['time_hex']:raise ValueError('complete accepted endpoint')
        tail=stream.read();rolling.update(tail)
        if tail.strip():raise ValueError('hidden raw tail')
    if rolling.hexdigest()!=expected_raw['sha256']:raise ValueError('whole raw SHA')
    return dict(schema=1,scope='binary_raw_sequence_clock_index_NOT_equation_review',raw=expected_raw,identity=meta,header=header,case=case,
        index=entries,end=end,trial_points=points,final_prefix_sha256=entries[-1]['raw_cursor']['prefix_sha256'],whole_raw_sha256=rolling.hexdigest(),
        source_case_sha256=streamed_typed_sha(case),header_sha256=streamed_typed_sha(header),equation_review_qualified=False)
