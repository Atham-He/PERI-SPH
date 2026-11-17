"""Geometric candidate bank for the existing antipodal-Fibonacci PD shell.
No material/fracture parameter is inferred or fitted. Actual annular-sector
areas integrate t*dell_mid once; barriers extend into the empty cavity solely
for nonlocal visibility and do not add fracture surface area.
"""
import argparse, json, math, hashlib
from pathlib import Path
import numpy as np
from scipy.spatial import SphericalVoronoi

def generate(inner, thickness, spacing, output):
    outer=inner+thickness
    half=max(32,math.ceil(2*math.pi*outer**2/spacing**2))
    directions=[]
    for i in range(half):
        z=1-i/half;phi=i*math.pi*(3-math.sqrt(5));r=math.sqrt(max(0,1-z*z))
        n=np.array([r*math.cos(phi),r*math.sin(phi),z])
        directions.extend([n,-n])
    directions=np.array(directions)
    voronoi=SphericalVoronoi(directions,radius=1,center=np.zeros(3))
    voronoi.sort_vertices_of_regions()
    edges={}
    for region,vertices in enumerate(voronoi.regions):
        for a,b in zip(vertices,vertices[1:]+vertices[:1]):
            edges.setdefault(tuple(sorted((a,b))),[]).append(region)
    xi=np.array([-math.sqrt(3/5),0,math.sqrt(3/5)])
    weight=np.array([5/9,8/9,5/9])
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    areas=[]
    with output.open('w') as f:
        f.write(f'PD_COHESIVE_BANK_V1 {len(edges)} {len(directions)}\n')
        def point(v): f.write(' '.join(format(float(x),'.17g') for x in v)+'\n')
        for index,(edge,seeds) in enumerate(sorted(edges.items())):
            if len(seeds)!=2:raise ValueError('Non-manifold spherical Voronoi edge')
            minus,plus=sorted(seeds);n=directions[plus]-directions[minus];n/=np.linalg.norm(n)
            a,b=[voronoi.vertices[j].copy() for j in edge]
            a-=n*np.dot(n,a);a/=np.linalg.norm(a)
            b-=n*np.dot(n,b);b/=np.linalg.norm(b)
            if np.dot(n,np.cross(a,b))<0:a,b=b,a
            tangent=np.cross(n,a)
            theta=math.atan2(np.dot(b,tangent),np.dot(b,a))
            if not 0<theta<math.pi:raise ValueError('Invalid crack-sector angle')
            f.write(f'{index} {minus} {plus} 1 6\n')
            for v in [inner*a,outer*a,outer*b,inner*b]:point(v)
            for v in [np.zeros(3),outer*a,outer*b]:point(v)
            area=0
            for gauss in [-1/math.sqrt(3),1/math.sqrt(3)]:
                angle=.5*theta*(1+gauss);direction=math.cos(angle)*a+math.sin(angle)*tangent
                for radial,w in zip(xi,weight):
                    radius=inner+.5*thickness*(1+radial)
                    measure=.5*theta*.5*thickness*w*radius
                    f.write(' '.join(format(float(x),'.17g') for x in [*(radius*direction),measure])+'\n')
                    area+=measure
            expected=.5*(outer**2-inner**2)*theta
            if abs(area/expected-1)>1e-12:raise ValueError('Physical crack-area quadrature mismatch')
            areas.append(area)
    receipt={'schema':'PD_COHESIVE_BANK_V1','particle_columns':len(directions),'candidate_facets':len(edges),'inner_radius_m':inner,'thickness_m':thickness,'pd_spacing_m':spacing,'area_min_m2':min(areas),'area_max_m2':max(areas),'area_definition':'annular radial section, counted once','barrier_extension_area_charged':False,'material_parameters_used':False,'sha256':hashlib.sha256(output.read_bytes()).hexdigest()}
    output.with_suffix(output.suffix+'.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--inner-radius',type=float,required=True);parser.add_argument('--thickness',type=float,required=True);parser.add_argument('--spacing',type=float,required=True);parser.add_argument('--output',required=True);a=parser.parse_args()
    generate(a.inner_radius,a.thickness,a.spacing,a.output)
