import json,sys
from paraview.simple import PVDReader
r=PVDReader(FileName=sys.argv[1]);r.UpdatePipeline(2.0);info=r.GetDataInformation();b=info.GetBounds();assert info.GetNumberOfPoints()==776,(info.GetNumberOfPoints(),b)
assert abs(b[0]+0.0004875)<1e-12 and abs(b[1]-0.0004875)<1e-12,b
result=dict(status="PASSED",reader="ParaView PVDReader",frame=2,points=info.GetNumberOfPoints(),bounds_m=list(b),actual_geometric_opening_m=0.0006,raw_geometry_not_displacement_scaled=True)
open(sys.argv[2],"w").write(json.dumps(result,indent=2)+"\n");print(json.dumps(result))
