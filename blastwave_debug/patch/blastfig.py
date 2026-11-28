import sys,glob,os,numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,'/home/user/diffusionIMEX/blastwave_debug')
from plot_blast import read_vtk
def last(d):
    g=glob.glob(f'{d}/plot/*/plt_*.vtk')
    return sorted(g,key=lambda s:int(s.rsplit('_',1)[1].split('.vtk')[0]))[-1] if g else None
def load(d):
    f=last(d); assert f, d
    dd,nx,ny,t=read_vtk(f)
    F={k:np.asarray(v,float).reshape(ny,nx) for k,v in dd.items() if np.asarray(v).size==nx*ny}
    F['_t']=t; F['_n']=nx
    F['umag']=np.hypot(F['velX'],F['velY'])
    F['pmag']=0.5*(F['B_X']**2+F['B_Y']**2+F.get('B_Z',np.zeros_like(F['rho']))**2)
    return F
# ---------- 图 1:论文的四面板 ----------
D='cfl08'; F=load(D); n=F['_n']; ext=[-0.5,0.5,-0.5,0.5]
panels=[('rho',r'density $\rho$','viridis'),('umag',r'$|\mathbf{u}|$','inferno'),
        ('p',r'gas pressure $p$','plasma'),('pmag',r'magnetic pressure $\|\mathbf{B}\|^2/2$','cividis')]
fig,ax=plt.subplots(2,2,figsize=(10.4,9.4))
for a,(k,lab,cm) in zip(ax.ravel(),panels):
    v=F[k]
    im=a.imshow(v,origin='lower',extent=ext,cmap=cm,aspect='equal')
    cs=a.contour(np.linspace(-0.5,0.5,n),np.linspace(-0.5,0.5,n),v,levels=12,colors='w',linewidths=0.35,alpha=.55)
    a.set_title(f"{lab}    [{v.min():.3g}, {v.max():.3g}]",fontsize=10)
    a.set_xlabel('x'); a.set_ylabel('y')
    plt.colorbar(im,ax=a,fraction=0.046,pad=0.03)
plt.suptitle(f"MHD blast wave, {n}$\\times${n}, CFL 0.8, $t={F['_t']}$   ($B_0=100/\\sqrt{{4\\pi}}$, $p_{{in}}=1000$)",fontsize=11)
plt.tight_layout(); plt.savefig('blast128_panels.png',dpi=135); print('wrote blast128_panels.png')
# ---------- 图 2:切线对比 cfl 0.2 / 0.8 ----------
A=load('n4g_new'); B=load('cfl08'); x=np.linspace(-0.5,0.5,n); mid=n//2
fig,ax=plt.subplots(2,2,figsize=(12,7.4))
for a,(k,lab) in zip(ax.ravel(),[('rho',r'$\rho$'),('p','$p$'),('umag',r'$|\mathbf{u}|$'),('pmag',r'$\|\mathbf{B}\|^2/2$')]):
    a.plot(x,A[k][mid,:],'-',color='tab:blue',lw=1.2,label='CFL 0.2, cut along $y=0$')
    a.plot(x,B[k][mid,:],'--',color='tab:red',lw=1.2,label='CFL 0.8, cut along $y=0$')
    a.plot(x,A[k][:,mid],'-',color='tab:green',lw=1.0,alpha=.8,label='CFL 0.2, cut along $x=0$')
    a.plot(x,B[k][:,mid],':',color='tab:orange',lw=1.4,label='CFL 0.8, cut along $x=0$')
    a.set_title(lab); a.set_xlabel('position'); a.grid(alpha=.3); a.legend(fontsize=7)
    if k=='p': a.set_yscale('log')
plt.suptitle(f'Blast wave {n}$\\times${n}: centre cuts, CFL 0.2 vs 0.8',fontsize=11)
plt.tight_layout(); plt.savefig('blast128_cuts.png',dpi=135); print('wrote blast128_cuts.png')
for lab,F2 in (('CFL 0.2',A),('CFL 0.8',B)):
    print(f"  {lab}: rho [{F2['rho'].min():.4f}, {F2['rho'].max():.4f}]  p [{F2['p'].min():.3e}, {F2['p'].max():.1f}]  |u| max {F2['umag'].max():.3f}")
