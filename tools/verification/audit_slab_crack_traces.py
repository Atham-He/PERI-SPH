"""Independent frozen-field audit of the formal slab's GMLS crack traces.
No forward solve, parameter modification or damage reassignment.
Polynomial span is invariant under nonsingular whitening; solve the same
weighted constraints in raw coordinates with direct NumPy pseudoinverses.
"""
import sys,json,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
A=Path(sys.argv[1]);dp=25e-6;th=5e-6;Y=16*dp
xi=np.array([-np.sqrt(3/5),0,np.sqrt(3/5)]);wl=np.array([5/9,8/9,5/9]);X=np.array([[(i+.5-8)*dp,(j+.5-8)*dp,.5*th*z] for z in xi for i in range(16) for j in range(16)]);vol=np.repeat(wl,256)
xml=ET.parse(A/'benchmark/growth_frames/final.vtp');Q=np.fromstring(xml.find('.//Points/DataArray').text,sep=' ').reshape(-1,3)
def vertices(k):
 lo=-Y/2+k*Y/4;hi=lo+Y/4;return np.array([[0,lo,-th/2],[0,hi,-th/2],[0,hi,th/2],[0,lo,th/2]])
def visible(a,b,cut_hi):
 if a[0]*b[0]>=0:return True
 t=a[0]/(a[0]-b[0]);y=(a+t*(b-a))[1]
 return y>cut_hi+1e-16
def fit(p,side,cut_hi):
 cand=np.flatnonzero(np.linalg.norm(X-p,axis=1)<3.5*dp)
 if side:cand=cand[side*X[cand,0]>0]
 anchor=cand[np.argmin(np.sum((X[cand]-p)**2,axis=1))];seen={int(anchor)};queue=[int(anchor)]
 for i in queue:
  for j in cand:
   if j not in seen and np.linalg.norm(X[j]-X[i])<3.015*dp and visible(X[i],X[j],cut_hi):seen.add(int(j));queue.append(int(j))
 ids=np.array([j for j in cand if j in seen and visible(X[anchor],X[j],cut_hi)],dtype=int);z=(X[ids]-p)/dp;W=np.exp(-np.sum(z*z,axis=1))*vol[ids];W/=W.sum();B=np.column_stack([np.ones(len(ids)),z,z*z,z[:,0]*z[:,1],z[:,0]*z[:,2],z[:,1]*z[:,2]])
 M=B.T@(W[:,None]*B);inv=np.linalg.pinv(M,rcond=1e-12);w=W*(B@inv[:,0]);assert np.linalg.norm(B.T@w-np.eye(10)[0])<1e-8
 return ids,B,W,inv,w
def transfer(f,target,Qbirth):
 ids,B,W,inv,w=f;w=w.copy();F=Qbirth[ids];R=F-B@(inv@(B.T@(W[:,None]*F)))
 for k in range(2):R-=B@(inv@(B.T@(W[:,None]*R)))
 rms=np.sqrt(np.sum(W[:,None]*R*R,axis=0));scale=np.where(rms>1024*np.finfo(float).eps*np.linalg.norm(F,axis=0)/np.sqrt(len(ids)),1/np.maximum(rms,1e-300),0);weighted=np.sqrt(W)[:,None]*R*scale;right=np.sqrt(W)[:,None]*np.linalg.pinv(weighted,rcond=1e-12).T*scale
 for k in range(5):
  error=target-F.T@w
  if np.linalg.norm(error)<=2e-15:break
  dw=right@error;dw-=W*(B@(inv@(B.T@dw)));w+=dw
 assert np.linalg.norm(F.T@w-target)<2e-15
 return ids,w
base={k:{s:[fit(p,s,1) for p in vertices(k)] for s in [-1,1]} for k in range(4)}
maps={0:{s:[(f[0],f[4]) for f in base[0][s]] for s in [-1,1]}}
for event in [1,2,3]:
 birth=np.loadtxt(A/('benchmark/preinsert_'+str(event)+'.csv'),delimiter=',',skiprows=1)[:,1:4]
 target={k:{s:np.array([birth[ids].T@w for ids,w in maps[k][s]]) for s in [-1,1]} for k in maps}
 v=vertices(event);common=[]
 for point in v:
  f=fit(point,0,v[0,1]);common.append(birth[f[0]].T@f[4])
 target[event]={s:np.array(common) for s in [-1,1]}
 maps={k:{s:[transfer(f,target[k][s][q],birth) for q,f in enumerate(base[k][s])] for s in [-1,1]} for k in target}
result=[]
for k in range(4):
 corners={s:np.array([Q[ids].T@w for ids,w in maps[k][s]]) for s in [-1,1]}
 for a in [-1/np.sqrt(3),1/np.sqrt(3)]:
  s=(a+1)/2
  for z in xi:
   t=(z+1)/2;N=np.array([(1-s)*(1-t),s*(1-t),s*t,(1-s)*t]);Nu=np.array([-(1-t),1-t,t,-t]);Nv=np.array([-(1-s),-s,s,1-s]);m=corners[-1];p=corners[1];n=np.cross(Nu@m,Nv@m);n/=np.linalg.norm(n);jump=N@(p-m);gn=jump@n;effective=jump-min(gn,0)*n
   result.append({'candidate':k,'angular_gauss':float(a),'thickness_gauss':float(z),'normal':n.tolist(),'raw_jump_m':jump.tolist(),'normal_gap_m':float(gn),'effective_opening_m':float(np.linalg.norm(effective))})
receipt={'method':'Independent reconstruction of position-constrained GMLS traces from all three saved events and final material point field','no_forward_solve':True,'no_state_or_parameter_changes':True,'points':result}
(A/'all_final_crack_trace_audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
for r in result:print(r['candidate'],round(r['angular_gauss'],4),round(r['thickness_gauss'],4),'n',np.array(r['normal']),'gn_um',r['normal_gap_m']*1e6,'effective_um',r['effective_opening_m']*1e6)
