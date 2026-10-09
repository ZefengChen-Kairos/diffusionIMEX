import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, exact5 as e
W=([4.4,1.4],[6e8,0.]); ES=([2.43,1.62],[5.3e9,141e9])
def ex(L,R,x0,t,n=2001):
    x=np.linspace(0,1,n); o,ps,us=e.solve(L,R,x0,t,x); o['x']=x
    Gm,PI=L.G,L.PI; a1=o['a1']; rho=o['rho']; p=o['p']
    with np.errstate(all='ignore'): o['cf']=np.sqrt((a1*o['r1']*Gm[0]*(p+PI[0])/o['r1']+(1-a1)*o['r2']*Gm[1]*(p+PI[1])/o['r2'])/rho)
    with np.errstate(all='ignore'): o['cw']=np.sqrt(1/(rho*(a1/(Gm[0]*(p+PI[0]))+(1-a1)/(Gm[1]*(p+PI[1]))))); return o
def Q(d,o,key):   # (numerical, exact) for a named quantity
    M={'p':(d[:,1]/1e5,o['p']/1e5),'u':(d[:,2],o['u']),'rho':(d[:,3],o['rho']),'a1':(1-d[:,4],o['a1']),'a2':(d[:,4],1-o['a1']),'r2':(d[:,7],o['r2']),
       'Mf':(abs(d[:,2])/d[:,8],abs(o['u'])/o['cf']),'Mw':(abs(d[:,2])/d[:,9],abs(o['u'])/o['cw'])}
    return M[key]
def fig(files,o,panels,title,out,ncol=2,cols=None,ylims=None):
    n=len(panels); nr=(n+ncol-1)//ncol; f,ax=plt.subplots(nr,ncol,figsize=(5.2*ncol,3.5*nr)); ax=np.atleast_1d(ax).ravel()
    for j,(key,lab) in enumerate(panels):
        d=np.loadtxt(files[j] if isinstance(files,list) else files); num,exa=Q(d,o,key)
        ax[j].plot(o['x'],exa,'k-',lw=0.9,label='exact'); ax[j].plot(d[:,0],num,'o',ms=2.6,mfc='none',mec='C0',mew=0.7,label='IMEX')
        ax[j].set_title(lab if cols is None else lab+' — '+cols[j%ncol],fontsize=10); ax[j].set_xlabel('x (m)'); ax[j].set_xlim(0,1)
        if j==0: ax[j].legend(fontsize=8)
        if ylims and key in ylims: ax[j].set_ylim(*ylims[key]); ax[j].ticklabel_format(useOffset=False)
    for j in range(n,len(ax)): ax[j].axis('off')
    f.suptitle(title,fontsize=9.5); f.tight_layout(); f.savefig(out,dpi=110); plt.close(f)
T='6+1 model, IMEX 2nd order (v22), CFL 0.2 — '
# 4.1
L=e.Side([1-1e-8,1e-8],[1000,10],100,1e5,*W); R=e.Side([1e-8,1-1e-8],[1000,10],100,1e5,*W); o=ex(L,R,0.5,2.79e-3)
fig('sp/c41.dat',o,[('p','Pressure (Bar)'),('u','Normal Velocity (m/s)'),('rho','Mixture density (kg/m3)'),('a1','Volume Fraction of water'),('Mf','Mach Number with frozen Speed of Sound'),('Mw','Mach Number with Wood Speed of Sound')],T+'SPB 4.1 interface advection, 200 cells, t = 2.79 ms (layout of Fig. 5)','spb_fig05.png',ylims={'p':(0.99,1.01),'u':(99,101)})
# 4.3.1
L=e.Side([1-1e-6,1e-6],[1000,1],0,1e9,*W); R=e.Side([1e-6,1-1e-6],[1000,1],0,1e5,*W); o=ex(L,R,0.75,240e-6); P4=[('p','Pressure (Bar)'),('u','Normal Velocity (m/s)'),('rho','Mixture density (kg/m3)'),('a1','Volume Fraction of water')]
fig('sp/c431_1000.dat',o,P4,T+'SPB 4.3.1 water–air shock tube, 1000 cells, t = 240 μs (layout of Fig. 8)','spb_fig08.png')
fig('sp/c431_100.dat',o,P4,T+'SPB 4.3.1 water–air shock tube, 100 cells, t = 240 μs (layout of Fig. 9)','spb_fig09.png')
P3=[('p','Pressure (Bar)'),('u','Normal Velocity (m/s)'),('rho','Mixture density (kg/m3)')]
fig(['sp/c431_1000.dat','sp/c44.dat']*3,o,[P3[0],P3[0],P3[1],P3[1],P3[2],P3[2]],T+'SPB 4.4 interface-pressure estimate in the stage-end relaxation, 1000 cells (layout of Fig. 11)','spb_fig11.png',cols=[r'$\hat p_I=p_I^0$',r'$\hat p_I=p$'])
# 4.3.2
L=e.Side([1-1e-6,1e-6],[1000,10],0,1e12,*W); R=e.Side([1e-6,1-1e-6],[1000,10],0,1e5,*W); o=ex(L,R,0.75,8.3e-6)
fig('sp/c432.dat',o,[('p','Pressure (Bar)'),('u','Normal Velocity (m/s)'),('Mf','Mach Number'),('rho','Mixture density (kg/m3)'),('a1','Volume Fraction of water')],T+'SPB 4.3.2 water–air shock tube, extreme conditions, 1000 cells, t = 8.3 μs (layout of Fig. 10)','spb_fig10.png')
# 4.5
L=e.Side([0.99,0.01],[1000,1],-100,1e5,*W); R=e.Side([0.99,0.01],[1000,1],100,1e5,*W); o=ex(L,R,0.5,1.85e-3)
fig('sp/c45.dat',o,[('p','Pressure (Bar)'),('u','Normal Velocity (m/s)'),('rho','Mixture density (kg/m3)'),('a2','Volume Fraction of air')],T+'SPB 4.5 cavitation test, 1000 cells, t = 1.85 ms (layout of Fig. 12)','spb_fig12.png')
# 5.1, 5.2
P6=[('p','Pressure (Bar)'),('u','Normal Velocity (m/s)'),('Mf','Mach Number'),('rho','Mixture density (kg/m3)'),('a1','Volume Fraction of epoxy'),('r2','Density of spinel (kg/m3)')]
a=[0.5954,0.4046]; r=[1185,3622]
for pl,t,fn,out,lab in [(1e10,80e-6,'sp/c51.dat','spb_fig19.png','SPB 5.1 epoxy–spinel shock tube, pressure ratio 1e5, 500 cells, t = 80 μs (layout of Fig. 19)'),(2e11,34e-6,'sp/c52.dat','spb_fig20.png','SPB 5.2 epoxy–spinel shock tube, pressure ratio 2e6, 500 cells, t = 34 μs (layout of Fig. 20; no artificial heat exchange)')]:
    L=e.Side(a,r,0,pl,*ES); R=e.Side(a,r,0,1e5,*ES); o=ex(L,R,0.6,t); fig(fn,o,P6,T+lab,out)
    d=np.loadtxt(fn); i=np.argmin(abs(d[:,0]-0.85)); j=np.argmin(abs(o['x']-0.85)); print(lab[:8],'shocked state at x=0.85: p %.4e (exact %.4e) u %.1f (%.1f) rho %.1f (%.1f) a_epoxy %.4f (%.4f) rho_spinel %.1f (%.1f)'%(d[i,1],o['p'][j],d[i,2],o['u'][j],d[i,3],o['rho'][j],1-d[i,4],o['a1'][j],d[i,7],o['r2'][j]))
    de=np.loadtxt(fn.replace('/c','/e')); print('         explicit:                p %.4e u %.1f rho %.1f a_epoxy %.4f rho_spinel %.1f'%(de[i,1],de[i,2],de[i,3],1-de[i,4],de[i,7]))
