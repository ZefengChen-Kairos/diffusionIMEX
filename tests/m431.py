import numpy as np, sys, exact_sg as e
# columns: 0 x 1 y 2 p 3 u 4 v 5 rho 6 alpha2 7 Y2 ...
for f in sys.argv[1:]:
    d=np.loadtxt(f); x=d[:,0]
    r,u,p,ps,us=e.solve(1000,0,1e9,4.4,6e8,1,0,1e5,1.4,0,0.75,240e-6,x)
    g=x>0.75+us*240e-6+0.01
    xs=x[np.where(d[:,3]>0.5*us)[0].max()]
    print('%-28s L1 u %6.2f rho %6.2f | gas max u %5.0f (exact %.0f) shock %.3f (0.935) gas L1 p/p* %.3f | u* %.1f p* %.4e'%(f.split('/')[-1],np.mean(abs(d[:,3]-u)),np.mean(abs(d[:,5]-r)),d[g,3].max(),us,xs,np.mean(abs(d[:,2]-p)[g])/ps,us,ps))
