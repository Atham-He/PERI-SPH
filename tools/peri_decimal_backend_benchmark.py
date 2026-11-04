"""Declared scalar-material timing/cross-runtime control, zero solver calls.

The pure-Python backend is selected BEFORE reference imports, never in response
to an exception. This benchmark is not a whole-trajectory replacement review.
"""
import argparse,json,resource,sys,time
from pathlib import Path


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--backend',choices=['C','Python'],required=True)
    parser.add_argument('--reference',type=Path,required=True);parser.add_argument('--input',type=Path,required=True)
    args=parser.parse_args()
    if args.backend=='Python':
        import _pydecimal as arithmetic
        sys.modules['decimal']=arithmetic
    else:import decimal as arithmetic
    sys.path.insert(0,str(args.reference))
    from peri_pd_material_reference import material_step
    fixture=json.loads(args.input.read_text())
    def encode(value):
        if isinstance(value,arithmetic.Decimal):return str(value)
        if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
        if isinstance(value,(list,tuple)):return [encode(v) for v in value]
        return value
    with arithmetic.localcontext() as context:
        context.prec=95
        params=[arithmetic.Decimal(x) for x in fixture['parameters']]
        f0,f1,z=([arithmetic.Decimal(x) for x in fixture[k]] for k in ['F0','F1','viscous0'])
        dt=arithmetic.Decimal(fixture['dt']);loss=arithmetic.Decimal(fixture['loss0'])
        values=[];timings=[]
        for _ in range(fixture['repetitions']):
            start=time.process_time();values.append(encode(material_step(f0,f1,z,loss,dt,params,True)))
            timings.append(time.process_time()-start)
        assert values[0]==values[1]
    print(json.dumps(dict(backend=args.backend,arithmetic_module=arithmetic.__file__,
        full_point_tangent=True,independent_decimal_digits=95,CPU_seconds=timings,
        peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,result=values[0],
        physical_solver_calls=0,study_passed=False,full_restart_qualified=False),allow_nan=False))


if __name__=='__main__':main()
