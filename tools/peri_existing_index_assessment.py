"""Independent assessment of an EXISTING index artifact,never a builder.

Consumes old producer bytes as unqualified evidence. Does not call old/new
index builders,solve/material functions or import a failed execution as passed.
"""
from pathlib import Path
import json
from fractions import Fraction
from peri_buffered_metadata import BinaryReader
from peri_streaming_typed_json import streamed_typed_sha
from peri_physical_time_friction_audit import schedule_grid


def require(condition,reason):
    if not condition:raise ValueError(reason)


def assess(index,path,case,expected_identity,expected_raw,frontiers,observe=lambda stage,step:None):
    fields={'schema','scope','raw','identity','header','case','index','end','trial_points','final_prefix_sha256','whole_raw_sha256',
            'source_case_sha256','header_sha256','equation_review_qualified'}
    require(set(index)==fields and index['schema']==1 and index['scope']=='binary_raw_sequence_clock_index_NOT_equation_review','existing index schema')
    require(index['equation_review_qualified'] is False,'metadata cannot grant equation qualification')
    require(index['raw']==expected_raw and index['whole_raw_sha256']==expected_raw['sha256'],'original pinned raw identity')
    expected_case=streamed_typed_sha(case)
    require(index['source_case_sha256']==expected_case==streamed_typed_sha(index['case']),'full case typed bytes')
    grid=schedule_grid(case);total=len(grid)-1;den=case['physical_time']['denominator']
    require(total==case['steps'] and len(index['index'])==total+1 and len(frontiers)==total+1,'all accepted frontiers required')
    reader=BinaryReader(path);points=0;old_time=0.;observe('before_full_native_identity',0)
    try:
        meta,_=reader.row();require(meta.get('type')=='identity','native identity row')
        require(all(meta.get(k)==v for k,v in expected_identity.items()),'full native argv/birth identity')
        require(streamed_typed_sha(meta)==streamed_typed_sha(index['identity']),'complete native identity bytes')
        header,cursor=reader.row()
        require(header.get('type')=='header' and header['format']==2 and header['steps_requested']==total,'physical-time header')
        require(header['primal_dofs']==99 and header['contact_count']==116 and header['unknowns']==447 and len(header['reference']['position'])==33,'bounded 33/99/116/447 scope')
        require(index['header_sha256']==streamed_typed_sha(header)==streamed_typed_sha(index['header']),'full header typed bytes')
        for h,c in [('time_denominator','denominator'),('nominal_ticks','nominal_ticks'),('normal_event_ticks','normal_event_ticks'),('tangent_event_ticks','tangent_event_ticks')]:
            require(header[h]==case['physical_time'][c],'frozen physical-time header cursor')
        old=header['initial_state'];require(old['step']==0 and float(old['time']).hex()==float(0.).hex(),'initial accepted clock')
        def check_entry(step,state,history,row,cursor,points):
            want=index['index'][step]
            require(set(want)=={'time_hex','trial_points','native_state_sha','friction_sha','native_row_sha','raw_cursor'},'complete frontier inventory')
            require(want['time_hex']==float(state['time']).hex() and want['trial_points']==points,'accepted time/point counter')
            require(want['raw_cursor']==cursor,'all byte offset/prefix/lineSHA')
            require(want['native_state_sha']==streamed_typed_sha(state),'entire native accepted state')
            require(want['friction_sha']==streamed_typed_sha(history),'entire native friction history')
            require(want['native_row_sha']==streamed_typed_sha(row),'complete native row including all trial payload')
            require(streamed_typed_sha(frontiers[step])==streamed_typed_sha(want),'independent durable frontier bytes')
        check_entry(0,old,header['initial_friction_history'],header,cursor,0);observe('accepted_metadata_verified',0)
        for step in range(1,total+1):
            observe('before_accepted_metadata',step);row,cursor=reader.row()
            require(row.get('type')=='step' and row.get('committed') is True and row['local_step_gates_passed'] is True and row['candidate_ready'] is True,'accepted runtime record flags')
            state=row['state'];interval=row['interval'];a=float(Fraction(grid[step-1],den));b=float(Fraction(grid[step],den));dt=b-a
            require(type(state['step']) is int and state['step']==step,'ordered accepted sequence')
            require((interval['old_tick'],interval['new_tick'],interval['time_denominator'])==(grid[step-1],grid[step],den),'accepted exact rational clock')
            require((float(interval['old_time']).hex(),float(interval['new_time']).hex(),float(state['time']).hex(),float(interval['actual_dt']).hex(),float(state['last_dt']).hex())==(a.hex(),b.hex(),b.hex(),dt.hex(),dt.hex()),'accepted binary endpoints/actual dt')
            require(float(old_time).hex()==a.hex() and float(old_time+interval['actual_dt']).hex()==b.hex(),'whole accepted temporal chain')
            require(row['solve']['converged'] is True and row['solve']['trace'],'complete Newton trace')
            trace=row['solve']['trace'];require([t['iteration'] for t in trace]==list(range(len(trace))),'all ordered trace records')
            points+=33*len(trace)
            check_entry(step,state,row['friction_history'],row,cursor,points)
            old_time=b;observe('accepted_metadata_verified',step)
        end,_=reader.row()
        require(end.get('type')=='end' and end.get('completed') is True and end['accepted_steps']==total,'complete native end marker')
        require(float(end['last_accepted_time']).hex()==float(old_time).hex()==float(header['end_time']).hex(),'complete end time')
        require(streamed_typed_sha(end)==streamed_typed_sha(index['end']),'entire native end record')
        require(index['trial_points']==points and index['final_prefix_sha256']==index['index'][-1]['raw_cursor']['prefix_sha256'],'full native point count/final prefix')
        whole=reader.finish(expected_raw);observe('whole_raw_complete_verified',total)
        return dict(accepted_steps=total,frontiers=total+1,point_metadata=points,final_time=old_time,raw_sha256=whole,
            all_complete_row_payloads_checked=True,all81prefix_line_offset_state_friction_row_typedSHA_checked=True,
            material_equation_qualified=False)
    finally:reader.close()
