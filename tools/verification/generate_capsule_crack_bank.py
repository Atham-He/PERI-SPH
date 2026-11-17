"""Resolution-compatible radial crack bank; no material/fracture fitting.
Material points are unchanged. Macro Voronoi regions have width 2*horizon,
so fully separated regions retain live NOSB support. Curved radial faces are
subdivided only for geometric quadrature/contact; subdivision does not create
new crack planes inside a material region. Exact annular area pays Gc once.
"""
import argparse, json, math, hashlib
from pathlib import Path
import numpy as np
from scipy.spatial import SphericalVoronoi

def directions(radius,spacing,minimum_half):
    half=max(minimum_half,math.ceil(2*math.pi*radius**2/spacing**2));result=[]
    for i in range(half):
        z=1-i/half;phi=i*math.pi*(3-math.sqrt(5));r=math.sqrt(max(0,1-z*z));n=np.array([r*math.cos(phi),r*math.sin(phi),z]);result.extend([n,-n])
    return np.array(result)

def generate(inner,thickness,spacing,output,horizon=None):
    outer=inner+thickness;horizon=3.015*spacing if horizon is None else horizon
    material=directions(outer,spacing,32);seeds=directions(outer,2*horizon,6)
    labels=np.argmax(material@seeds.T,axis=1);counts=np.bincount(labels,minlength=len(seeds))
    if counts.min()<8:raise ValueError('PD discretization does not resolve the macro crack regions')
    # A frozen geometric support audit is not another physical loading case.
    min_neighbors=1000000;max_affine_error=0
    for radius in [inner+.5*thickness*(1-math.sqrt(3/5)),inner+.5*thickness,inner+.5*thickness*(1+math.sqrt(3/5))]:
        X=radius*material
        for region in range(len(seeds)):
            ids=np.flatnonzero(labels==region);x=X[ids];dist=np.linalg.norm(x[:,None,:]-x[None,:,:],axis=2)
            for local,node in enumerate(ids):
                neighbors=np.flatnonzero((dist[local]>0)&(dist[local]<horizon));min_neighbors=min(min_neighbors,len(neighbors))
                B=(x[neighbors]-x[local])/spacing;n=material[node];target=np.eye(3)-np.outer(n,n)
                fit=np.linalg.lstsq(B.T,target,rcond=1e-12)[0];error=np.linalg.norm(B.T@fit-target)/np.linalg.norm(target);max_affine_error=max(max_affine_error,error)
    if max_affine_error>1e-8:raise ValueError('Fully split macro regions cannot reproduce affine NOSB kinematics')
    voronoi=SphericalVoronoi(seeds,radius=1,center=np.zeros(3));voronoi.sort_vertices_of_regions();edges={}
    for region,vertices in enumerate(voronoi.regions):
        for a,b in zip(vertices,vertices[1:]+vertices[:1]):edges.setdefault(tuple(sorted((a,b))),[]).append(region)
    xi=np.array([-math.sqrt(3/5),0,math.sqrt(3/5)]);weight=np.array([5/9,8/9,5/9]);mid=inner+.5*thickness
    # Keep each analytic quadrature point inside its chordal contact patch.
    sag_limit=.25*thickness*(1-math.sqrt(3/5));theta_limit=2*math.acos(1-sag_limit/outer)
    bank=[];areas=[]
    for parent,(edge,owners) in enumerate(sorted(edges.items())):
        if len(owners)!=2:raise ValueError('Non-manifold Voronoi edge')
        minus,plus=sorted(owners);n=seeds[plus]-seeds[minus];n/=np.linalg.norm(n);a,b=[voronoi.vertices[j].copy() for j in edge]
        a-=n*np.dot(n,a);a/=np.linalg.norm(a);b-=n*np.dot(n,b);b/=np.linalg.norm(b)
        if np.dot(n,np.cross(a,b))<0:a,b=b,a
        tangent=np.cross(n,a);theta=math.atan2(np.dot(b,tangent),np.dot(b,a))
        if not 0<theta<math.pi:raise ValueError('Invalid radial crack-sector angle')
        pieces=max(1,math.ceil(theta/theta_limit));side=material@n
        for piece in range(pieces):
            lo=theta*piece/pieces;hi=theta*(piece+1)/pieces;u=math.cos(lo)*a+math.sin(lo)*tangent;v=math.cos(hi)*a+math.sin(hi)*tangent
            points=[];pairs=[];area=0
            for gauss in [-1/math.sqrt(3),1/math.sqrt(3)]:
                angle=.5*(lo+hi)+.5*(hi-lo)*gauss;direction=math.cos(angle)*a+math.sin(angle)*tangent;scores=material@direction
                im=int(np.argmax(np.where(side< -1e-13,scores,-np.inf)));ip=int(np.argmax(np.where(side>1e-13,scores,-np.inf)))
                for layer,(radial,w) in enumerate(zip(xi,weight)):
                    radius=mid+.5*thickness*radial;measure=.5*(hi-lo)*.5*thickness*w*radius;points.append([*(radius*direction),measure]);pairs.append((layer*len(material)+im,layer*len(material)+ip));area+=measure
            expected=.5*(outer**2-inner**2)*(hi-lo)
            if abs(area/expected-1)>1e-12:raise ValueError('Physical fracture area quadrature mismatch')
            bank.append((parent,[inner*u,outer*u,outer*v,inner*v],[np.zeros(3),outer*u,outer*v],points,pairs));areas.append(area)
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('w') as f:
        f.write(f'PD_COHESIVE_BANK_V2 {len(bank)} {len(material)}\n')
        def point(v):f.write(' '.join(format(float(x),'.17g') for x in v)+'\n')
        for index,(parent,vertices,barrier,points,pairs) in enumerate(bank):
            f.write(f'{index} 1 {len(points)} {len(pairs)}\n')
            for q in vertices:point(q)
            for q in barrier:point(q)
            for q in points:point(q)
            for i,j in pairs:f.write(f'{i} {j}\n')
    receipt={'schema':'PD_COHESIVE_BANK_V2','particle_columns':len(material),'macro_regions':len(seeds),'macro_faces':len(edges),'candidate_subfacets':len(bank),'horizon_m':horizon,'macro_width_m':2*horizon,'minimum_material_columns_per_region':int(counts.min()),'minimum_same_layer_neighbors_if_all_regions_split':min_neighbors,'maximum_affine_functional_error_if_all_regions_split':max_affine_error,'geometric_sag_limit_m':sag_limit,'inner_radius_m':inner,'thickness_m':thickness,'pd_spacing_m':spacing,'area_min_m2':min(areas),'area_max_m2':max(areas),'area_definition':'exact annular section, once','barrier_extension_area_charged':False,'material_parameters_used':False,'sha256':hashlib.sha256(output.read_bytes()).hexdigest()}
    output.with_suffix(output.suffix+'.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inner-radius',type=float,required=True);p.add_argument('--thickness',type=float,required=True);p.add_argument('--spacing',type=float,required=True);p.add_argument('--horizon',type=float);p.add_argument('--output',required=True);a=p.parse_args();generate(a.inner_radius,a.thickness,a.spacing,a.output,a.horizon)
