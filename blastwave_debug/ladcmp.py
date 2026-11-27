import sys,glob,os,numpy as np
sys.path.insert(0,'/home/user/diffusionIMEX/blastwave_debug')
from plot_blast import read_vtk
LAB={'ref':'ref   flat1=0 (现状: 3 传感器+cap)','L1':'L1    flat1=1 (统一, env 挂 chi, cap 开)',
     'L2':'L2    flat1=1 cap=0','L3':'L3    flat1=2 cap=0 (无压力传感器)',
     'L4':'L4    flat1=3 cap=0 (env 只挂跳变)','L5':'L5    flat1=4 cap=0 (最小)',
     'L6':'L6    flat1=3 (env 只挂跳变, cap 开)','none':'none  完全无 flattener'}
def noise(a):
    c=a[1:-1,1:-1]
    up=(c>a[:-2,1:-1])&(c>a[2:,1:-1])&(c>a[1:-1,:-2])&(c>a[1:-1,2:])
    dn=(c<a[:-2,1:-1])&(c<a[2:,1:-1])&(c<a[1:-1,:-2])&(c<a[1:-1,2:])
    return 100.0*(up|dn).mean()
print(f"{'run':40s} {'t':>6s} {'rho max':>8s} {'grad rho':>9s} {'shk r':>7s} {'nz rho%':>8s} {'nz p%':>7s} {'rms divu':>9s} {'<=8格%':>7s}")
for t in ('ref','L1','L2','L3','L4','L5','L6','none'):
    g=sorted(glob.glob(f'lad_{t}/plot/*/plt_*.vtk'),key=lambda s:int(s.rsplit('_',1)[1].split('.vtk')[0]))
    if not g or not os.path.exists(f'lad_{t}/DONE'): print(f"{LAB.get(t,t):40s}   —"); continue
    d,nx,ny,tt=read_vtk(g[-1])
    F={k:np.asarray(v,float).reshape(ny,nx) for k,v in d.items() if np.asarray(v).size==nx*ny}
    rho,u,v=F['rho'],F['velX'],F['velY']
    x=np.linspace(-0.5,0.5,nx); X,Y=np.meshgrid(x,x); rr=np.hypot(X,Y); h=x[1]-x[0]
    gr=np.hypot(*np.gradient(rho)); i=np.unravel_index(gr.argmax(),gr.shape)
    du=np.zeros_like(u); du[1:-1,1:-1]=(u[1:-1,2:]-u[1:-1,:-2])/(2*h)+(v[2:,1:-1]-v[:-2,1:-1])/(2*h)
    ix=np.abs(x)<0.1; iy=np.abs(x)<0.05
    blk=du[np.ix_(iy,ix)]; blk=blk-blk.mean(axis=1,keepdims=True)
    sp=np.mean(np.abs(np.fft.rfft(blk,axis=1))**2,axis=0); k=np.fft.rfftfreq(blk.shape[1],d=1.0)
    sp[0]=0; tot=sp[1:].sum(); hi=sp[1:][k[1:]>=1/8.0].sum()
    print(f"{LAB.get(t,t):40s} {tt:6.4f} {rho.max():8.4f} {gr.max():9.4f} {rr[i]:7.4f} {noise(rho):8.3f} {noise(F['p']):7.3f} {np.sqrt((du**2).mean()):9.3f} {100*hi/max(tot,1e-300):7.1f}")
