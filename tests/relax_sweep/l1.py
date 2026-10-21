import sys; sys.argv=['x','none']; sys.path.insert(0,'.'); sys.path.insert(0,'/home/user/diffusionIMEX/tests')
import ana as A, numpy as np, exact_relaxed as er, glob, os
PEL=[(1.187,7028e5,-1177788.0,0.0,3610.0,6.61e-4),(1.467,0.0,2077616.0,14317.0,955.0,0.0)]
F=er.Fan(PEL,1e5,353.0,1e-2); ex={}
for chem in (False,True):
    (x,p,u,Y,a,T),I=F.full(2.0,3e-3,chem); o=np.argsort(x); ex[chem]=(x[o],p[o],Y[o],T[o])
def err(n,chem):
    d=A.L(A.f(n)); x,p,Y,T=ex[chem]
    pe=np.interp(d[:,0],x,p); Ye=np.interp(d[:,0],x,Y); Te=np.interp(d[:,0],x,T)
    return np.mean(np.abs(d[:,2]-pe)), np.mean(np.abs(d[:,7]-Ye)), np.mean(np.abs(d[:,12]-Te)), np.mean(np.abs(d[:,13]-Te))
print('theta: L1 error vs exact 4-eq pT fan: p, T1, T2 (BE | EXP)')
for th in A.THS:
    r=[err(A.thname(th,fo),False) for fo in (0,1)]; print(th,' | '.join('%.6g %.6g %.6g'%(e[0],e[2],e[3]) for e in r))
print('nu: L1 error vs exact pT->pTg fan: p, Y, T (BE | EXP)')
for nu in A.NUS:
    r=[err(A.nuname(nu,fo),True) for fo in (0,1)]; print(nu,' | '.join('%.6g %.6g %.6g'%(e[0],e[1],e[2]) for e in r))
print('inf proj1', '%.6g %.6g %.6g'%err('E2_pTg_f1_up1',True)[:3])
print('floor: other pTg-limit runs (centre p, Y, alpha2 ; L1 p,Y vs exact)')
for n in ['/home/user/diffusionIMEX/run/pc/pel821_pTg.dat']+sorted(glob.glob(A.R+'*pTg.dat')):
    if 'pTgf' in n: continue
    c=A.ctr(n); e=err(n,True); print(os.path.basename(n),'%.7g %.6g %.6g ; %.6g %.6g'%(c['p'],c['Y2'],c['a2'],e[0],e[1]), open(n[:-4]+'.log').readline()[:0])
