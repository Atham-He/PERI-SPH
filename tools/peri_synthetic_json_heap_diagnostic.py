"""Synthetic-only stress of the original typed JSON temporary/deallocation.

No original raw,material equations,solver or numerical thresholds are used.
This may localize a runtime crash;not reproducing it does not establish cause.
"""
from pathlib import Path
import faulthandler,gc,hashlib,json,os,resource,sys,time


def main(root):
    fault=(root/'faulthandler.txt').open('x');faulthandler.enable(file=fault,all_threads=True)
    from peri_fixed_python_raw_audit import select_backend
    backend=select_backend(Path(__file__).parent/'reference');backend.getcontext().prec=95
    from peri_reference_checkpoint import pack,canonical
    plan=json.loads((root/'plan.json').read_text());spec=plan['synthetic_input']
    def event(stage,iteration,**extra):
        q=dict(stage=stage,iteration=iteration,CPU_seconds=time.process_time(),RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,**extra)
        with (root/'stage_witness.jsonl').open('a') as f:f.write(json.dumps(q)+'\n');f.flush();os.fsync(f.fileno())
    n=spec['unknowns'];nodes=spec['nodes']
    event('before_synthetic_build',-1)
    matrix=[[0. if i!=j else float((i%17)+1)/17 for j in range(n)] for i in range(n)]
    tangent=[[float((i+3*j)%19)/19 for j in range(81)] for i in range(nodes)]
    trace=[dict(iteration=k,original=dict(jacobian=matrix,residual=[float(k)/7]*n),points=[dict(endpoint_tangent=t,old=[0.]*9,new=[float(k)/13]*9) for t in tangent]) for k in range(spec['trace_records'])]
    row=dict(type='SYNTHETIC_MEMORY_DIAGNOSTIC_NOT_PHYSICS',trace=trace,state=dict(step=1,position=[[0.,-0.,1.]]*nodes))
    raw=json.dumps(row,allow_nan=False,separators=(',',':')).encode();rawsha=hashlib.sha256(raw).hexdigest()
    del matrix,tangent,trace,row;gc.collect();event('synthetic_input_ready',-1,raw_bytes=len(raw),raw_sha256=rawsha)
    hashes=[];payload_bytes=[]
    for k in range(spec['repetitions']):
        event('before_json_load',k);value=json.loads(raw);event('after_json_load',k)
        event('before_original_pack',k);typed=pack(value);event('after_original_pack',k)
        event('before_original_canonical',k);payload=canonical(typed);event('after_original_canonical',k,payload_bytes=len(payload))
        hashes.append(hashlib.sha256(payload).hexdigest());payload_bytes.append(len(payload))
        event('before_temporary_release',k);del payload,typed,value;gc.collect();event('after_temporary_release',k)
    assert len(set(hashes))==len(set(payload_bytes))==1
    q=dict(status='new_synthetic_large_typed_JSON_diagnostic_completed_no_crash_reproduced',repetitions=len(hashes),
        synthetic_input_bytes=len(raw),typed_bytes=payload_bytes[0],typed_sha256=hashes[0],CPU_seconds=time.process_time(),
        peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,precision=backend.getcontext().prec,
        C_decimal_loaded='_decimal' in sys.modules,root_cause_resolved=False,physical_solver_calls=0,material_calculations=0,old_raw_reads=0)
    with (root/'diagnostic_result.json').open('x') as f:json.dump(q,f,indent=2)
    print(json.dumps(q));fault.close()


if __name__=='__main__':main(Path(sys.argv[1]))
