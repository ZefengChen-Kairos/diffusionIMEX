import sys, glob, numpy as np
sys.path.insert(0,'/home/user/diffusionIMEX/blastwave_debug')
from plot_blast import read_vtk

def vtks(t):
    g=glob.glob(f'e2/{t}/plot/*/plt_*.vtk')
    if not g: return None
    g.sort(key=lambda s:int(s.rsplit('_',1)[1].split('.vtk')[0]))
    return g

def F(path):
    d,nx,ny,tt=read_vtk(path)
    K={'rho':'rho','p':'p','u':'velX','v':'velY','Bx':'B_X','By':'B_Y','Bz':'B_Z','eM':'eps_M','eB':'eps_B'}
    return {a:np.asarray(d[b],float).reshape(ny,nx) for a,b in K.items()},tt,nx,ny

LAB={'base':'base (default: eps_B adaptive, HLLEM delta_k on)',
     'noad':'noad  (hllemoff=1: all delta_k anti-diffusion off)',
     'epsB0':'epsB0 (epsBmax=0: eps_B pinned to 0 = single-parameter proxy)'}
R={}
for t in ('base','noad','epsB0'):
    g=vtks(t)
    if not g: print(f"[{t}] 无输出"); continue
    f0,t0,nx,ny=F(g[0]); fN,tN,_,_=F(g[-1])
    R[t]=(f0,t0,fN,tN,nx,ny)

print("="*78)
print("field loop, 128x64, A0=1, CFL 0.2, senv=2, one full period (t=1)")
print("exact solution at t=1 == initial condition (periodic, u=(2,1), domain 2 x 1)")
print("="*78)
for t,(f0,t0,fN,tN,nx,ny) in R.items():
    dx,dy=2.0/nx,1.0/ny
    B0=np.hypot(f0['Bx'],f0['By']); BN=np.hypot(fN['Bx'],fN['By'])
    e0=0.5*(f0['Bx']**2+f0['By']**2).sum()*dx*dy
    eN=0.5*(fN['Bx']**2+fN['By']**2).sum()*dx*dy
    print(f"\n--- {LAB[t]} ---   t0={t0}  tN={tN}")
    print(f"    eps_M: min={fN['eM'].min():.4e} max={fN['eM'].max():.4e}   eps_B: min={fN['eB'].min():.4e} max={fN['eB'].max():.4e}")
    print(f"    E_mag(t)/E_mag(0) = {eN/e0:.6f}      (1 = no numerical dissipation)")
    print(f"    max|B_perp|:  {BN.max():.6e} / {B0.max():.6e} = {BN.max()/B0.max():.6f}")
    for k in ('Bx','By'):
        d=fN[k]-f0[k]; sc=np.abs(f0[k]).max()
        print(f"    err vs exact {k}:  L1={np.abs(d).mean()/sc:.4e}  L2={np.sqrt((d**2).mean())/sc:.4e}  Linf={np.abs(d).max()/sc:.4e}")

if 'base' in R:
    f0b,_,fb,_,nx,ny=R['base']
    for t in ('noad','epsB0'):
        if t not in R: continue
        _,_,fo,_,_,_=R[t]
        print("\n"+"="*78); print(f"PAIRWISE  {t}  vs  base      (roundoff would be ~1e-14)"); print("="*78)
        for k in ('Bx','By','Bz','u','v','rho','p'):
            d=np.abs(fo[k]-fb[k]); sc=max(np.abs(fb[k]).max(),1e-300)
            print(f"    {k:3s}  max|diff|={d.max():.4e}   rel={d.max()/sc:.4e}   L2rel={np.linalg.norm(d)/max(np.linalg.norm(fb[k]),1e-300):.4e}")
