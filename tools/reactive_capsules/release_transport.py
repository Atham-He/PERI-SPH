"""Conservative diffusion/dissolution on the conforming solid reference cells.

First implementation: fast imbibition after hydraulic access, local solubility
equilibrium, dilute external sink, and one-way mechanics -> transport -> bath
equilibrium. Every one of these reductions is recorded in the case manifest.
No SS inventory is deleted when shell damage appears.
"""
from pathlib import Path
import hashlib
import os
import struct
import sys
if os.environ.get('PERI_REACTIVE_PYTHON_DEPS'):
    sys.path.insert(0,os.environ['PERI_REACTIVE_PYTHON_DEPS'])
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import splu
from scipy.spatial import cKDTree


def load_reference_cells(path):
    """The exact binary format consumed by SnsniReference, not a new mesh."""
    with Path(path).open('rb') as f:
        nc, nf = struct.unpack('<QQ', f.read(16))
        center=np.empty((nc,3));volume=np.empty(nc);owners=[[] for _ in range(nf)]
        for i in range(nc):
            a=np.frombuffer(f.read(13*8),dtype='<f8')
            center[i]=a[:3];volume[i]=a[3]
            n,=struct.unpack('<I',f.read(4))
            ids=np.frombuffer(f.read(n*8),dtype='<i8')
            for fid in ids:owners[abs(fid)-1].append(i)
        faces=[];boundary=[];face_normals=[]
        for k in range(nf):
            normal=np.frombuffer(f.read(24),dtype='<f8').copy()
            is_boundary,n=struct.unpack('<II',f.read(8))
            q=np.frombuffer(f.read(n*32),dtype='<f8').reshape(n,4)
            area=q[:,3].sum();mid=(q[:,:3]*q[:,3,None]).sum(axis=0)/area
            if is_boundary:
                if len(owners[k])!=1:raise ValueError('Boundary ownership mismatch')
                boundary.append((owners[k][0],area,mid,normal))
            else:
                if len(owners[k])!=2:raise ValueError('Internal face ownership mismatch')
                i,j=owners[k];dist=abs(np.dot(center[j]-center[i],normal))
                if dist<=0:raise ValueError('Invalid cell distance')
                faces.append((i,j,area,dist));face_normals.append(normal)
        if f.read(1):raise ValueError('Unexpected geometry trailer')
    edge=np.array(faces)
    return dict(center=center,volume=volume,i=edge[:,0].astype(int),j=edge[:,1].astype(int),
                geometric_conductance=edge[:,2]/edge[:,3],
                face_normal=np.array(face_normals),
                boundary_cell=np.array([b[0] for b in boundary],dtype=int),
                boundary_area=np.array([b[1] for b in boundary]),
                boundary_normal=np.array([b[3] for b in boundary]),
                boundary_point=np.array([b[2] for b in boundary]))


class ConservativeRelease:
    def __init__(self, mesh, rho0, ss_fraction, modulus, porosity,
                 diffusivity, free_diffusivity, solubility_molSi_m3):
        if not (0<=ss_fraction<=1 and 0<porosity<1 and min(rho0,modulus,diffusivity,free_diffusivity,solubility_molSi_m3)>0):
            raise ValueError('Invalid release parameters')
        self.mesh=mesh;self.modulus=modulus
        # kg/mol of Na2O.(SiO2)_modulus. No Na2SiO3 stoichiometry substitution.
        self.formula_kg_mol=.06197853856+modulus*.0600843
        self.initial_ss_kg=rho0*ss_fraction*mesh['volume'].sum()
        self.initial=mesh['volume']*rho0*ss_fraction/self.formula_kg_mol*modulus
        self.solid=self.initial.copy();self.aqueous=np.zeros_like(self.solid)
        self.released_Si=0.;self.wetted=False;self.physical_time=0.
        self.porosity=porosity;self.D=diffusivity;self.Dfree=free_diffusivity
        self.Csat=solubility_molSi_m3;self._cache=None

    def set_barrier(self, barrier, J=None, F=None, relative_density=None):
        directions=np.column_stack([barrier[x] for x in ('nx','ny','nz')])
        query=self.mesh['boundary_point'];query=query/np.linalg.norm(query,axis=1)[:,None]
        _,nearest=cKDTree(directions).query(query)
        fraction=np.atleast_1d(barrier['open_fraction'])[nearest]
        length=np.atleast_1d(barrier['transport_length_m'])[nearest]
        if np.any(fraction<0) or np.any(fraction>1) or np.any(length<=0):raise ValueError('Invalid barrier')
        if J is None:J=np.ones_like(self.solid)
        if np.any(J<=0) or np.any(~np.isfinite(J)):raise ValueError('Invalid transport Jacobian')
        self.conductance=self.D*self.mesh['geometric_conductance']
        area=self.mesh['boundary_area']
        if F is not None:
            if not np.allclose(np.linalg.det(F),J,rtol=1e-9,atol=1e-12):raise ValueError('F/J mismatch')
            inverse_transpose=np.linalg.inv(F).transpose(0,2,1)
            i,j=self.mesh['i'],self.mesh['j'];n=self.mesh['face_normal']
            # Piola pullback of Fick flux: J F^-1 D F^-T. The positive
            # two-point normal projection is used on each shared cell face.
            ni=np.einsum('nij,nj->ni',inverse_transpose[i],n)
            nj=np.einsum('nij,nj->ni',inverse_transpose[j],n)
            self.conductance=self.conductance*.5*(J[i]*(ni*ni).sum(axis=1)+J[j]*(nj*nj).sum(axis=1))
            b=self.mesh['boundary_cell'];nb=self.mesh['boundary_normal']
            area=area*J[b]*np.linalg.norm(np.einsum('nij,nj->ni',inverse_transpose[b],nb),axis=1)
        self.boundary=np.bincount(self.mesh['boundary_cell'],
            weights=self.Dfree*area*fraction/length,minlength=len(self.solid))
        self.exposed_area=float(np.dot(area,fraction))
        porosity=self.porosity if relative_density is None else 1-relative_density
        if np.any(np.asarray(porosity)<=0) or np.any(np.asarray(porosity)>=1):
            raise ValueError('Assumed DPC porosity mapping is outside (0,1); no clipping')
        self.capacity=porosity*self.mesh['volume']*J
        self._cache=None

    def step(self, dt):
        if dt<=0:raise ValueError('Physical time step must be positive')
        if np.any(self.boundary>0):self.wetted=True
        self.physical_time+=dt
        if not self.wetted or self.initial_ss_kg==0:return self.summary()
        # Fast-imbibition limit, followed by local dissolution equilibrium.
        # Recrystallization into the same SS reservoir conserves inventory.
        saturated=self.Csat*self.capacity
        reprecip=np.maximum(self.aqueous-saturated,0)
        self.aqueous-=reprecip;self.solid+=reprecip
        dissolve=np.minimum(self.solid,np.maximum(saturated-self.aqueous,0))
        self.solid-=dissolve;self.aqueous+=dissolve
        if self._cache is None or self._cache[0]!=dt:
            i,j=self.mesh['i'],self.mesh['j'];g=self.conductance
            diagonal=self.capacity/dt+self.boundary+np.bincount(np.r_[i,j],weights=np.r_[g,g],minlength=len(self.solid))
            matrix=coo_matrix((np.r_[diagonal,-g,-g],
                (np.r_[np.arange(len(diagonal)),i,j],np.r_[np.arange(len(diagonal)),j,i])),shape=(len(diagonal),len(diagonal))).tocsc()
            self._cache=(dt,splu(matrix))
        concentration=self._cache[1].solve(self.aqueous/dt)
        if np.any(~np.isfinite(concentration)) or np.min(concentration)<-1e-12:
            raise RuntimeError('Nonpositive/nonfinite diffusion solution; no inventory clipping')
        self.aqueous=self.capacity*concentration
        self.released_Si+=dt*float(np.dot(self.boundary,concentration))
        v=self.summary()
        if abs(v['Si_inventory_residual_mol'])>max(1e-20,1e-9*self.initial.sum()):
            raise RuntimeError('Release violated Si conservation')
        return v

    def summary(self):
        remaining=float(self.solid.sum());aq=float(self.aqueous.sum());total=float(self.initial.sum())
        return dict(physical_time_s=self.physical_time,initial_SS_kg=self.initial_ss_kg,
            solid_Si_mol=remaining,core_aqueous_Si_mol=aq,released_Si_mol=self.released_Si,
            released_SS_formula_mol=self.released_Si/self.modulus,
            released_SS_equivalent_kg=self.released_Si/self.modulus*self.formula_kg_mol,
            released_fraction=self.released_Si/total if total else 0,
            hydraulic_exposed_area_m2=getattr(self,'exposed_area',0.),
            Si_inventory_residual_mol=remaining+aq+self.released_Si-total,
            fast_imbibition_assumed=True,external_dilute_sink_assumed=True,
            mechanical_dissolution_feedback=False)

    def save(self,path):
        p=Path(path);temp=Path(str(p)+'.tmp')
        with temp.open('wb') as f:np.savez_compressed(f,solid=self.solid,aqueous=self.aqueous,
            released_Si=self.released_Si,wetted=self.wetted,physical_time=self.physical_time)
        temp.replace(p)

    def restore(self,path):
        a=np.load(path);self.solid=a['solid'];self.aqueous=a['aqueous']
        self.released_Si=float(a['released_Si']);self.wetted=bool(a['wetted']);self.physical_time=float(a['physical_time'])
