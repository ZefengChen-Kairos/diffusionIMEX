import numpy as np, sys, re, exact_sg as e, exactfan as xf, exact5 as e5
import os; S=os.environ.get('RUN', os.getcwd()+'/run')+'/'; RF=os.environ.get('REF', os.getcwd()+'/ref')+'/'
def st(n,T):
    t=open(S+'su/%s_%s.out'%(n,T)).read(); g=lambda p:re.search(p,t).group(1)
    return 'fb %s nf %s pf %s'%(g(r'RK1 fallback (\d+)'),g(r'Newton fail (\d+)'),g(r'p-floor hits (\d+)'))
for T in sys.argv[1:]:
    print('==',T); L=lambda n:np.loadtxt(S+'su/%s_%s.dat'%(n,T))
    for n in ['wa1000','wa100']:
        d=L(n); x=d[:,0]; r,u,p,ps,us=e.solve(1000,0,1e9,4.4,6e8,1,0,1e5,1.4,0,0.75,240e-6,x); g=x>0.75+us*240e-6+0.01; xs=x[np.where(d[:,2]>0.5*us)[0].max()]
        print('  4.3.1 %-6s: L1 u %6.2f rho %6.2f | gas max u %4.0f (492) shock %.3f (0.935) gas L1 p/p* %.3f | %s'%(n[2:],np.mean(abs(d[:,2]-u)),np.mean(abs(d[:,3]-r)),d[g,2].max(),xs,np.mean(abs(d[:,1]-p)[g])/ps,st(n,T)))
    d=L('w2'); x=d[:,0]; r,u,p,ps,us=e.solve(1000,0,1e12,4.4,6e8,10,0,1e5,1.4,0,0.75,8.3e-6,x); print('  4.3.2: L1 u %.1f rho %.3f max u/u* %.4f | %s'%(np.mean(abs(d[:,2]-u)),np.mean(abs(d[:,3]-r)),d[:,2].max()/us,st('w2',T)))
    d=L('adv'); print('  4.1: max|p/p0-1| %.1e max|u-100| %.1e'%(abs(d[:,1]/1e5-1).max(),abs(d[:,2]-100).max()))
    for rl in [0,2]:
        R=np.loadtxt(RF+'ref3_%d.dat'%rl).reshape(-1,4,6).mean(1); d=L('d%d'%rl)
        print('  dodecane RELAX %d: L1 u %.3f p/1e5 %.3f Yv %.2e | min p on [0.1,0.4] %.3e (ref %.3e) | %s'%(rl,np.mean(abs(d[:,2]-R[:,2])),np.mean(abs(d[:,1]-R[:,1]))/1e5,np.mean(abs(d[:,5]-R[:,5])),d[100:400,1].min(),R[100:400,1].min(),st('d%d'%rl,T)))
    R=np.loadtxt(RF+'ref0_0.dat').reshape(-1,4,6).mean(1); d=L('cav'); print('  4.5: L1 u %.3f alpha %.2e | %s'%(np.mean(abs(d[:,2]-R[:,2])),np.mean(abs(d[:,4]-R[:,4])),st('cav',T)))
    for C,U,t in [(1,2,3.2e-3),(2,500,0.58e-3)]:
        X,P,Uu,Y,A=xf.full(U,t); d=L('c%d'%C); ip=lambda q:np.interp(d[:,0],X,q)
        print('  cav %d pTg: L1 p/1e5 %.2e u/U %.2e Yv/Ymax %.2e | centre p %.1f a %.4f'%(U,np.mean(abs(d[:,1]-ip(P)))/1e5,np.mean(abs(d[:,2]-ip(Uu)))/U,np.mean(abs(d[:,5]-ip(Y)))/Y.max(),d[500,1],d[500,4]))
    a=[0.5954,0.4046]; rr=[1185,3622]; ES=([2.43,1.62],[5.3e9,141e9])
    for nm,pl,t in [('e1',1e10,80e-6),('e2',2e11,34e-6)]:
        d=L(nm); o,ps,us=e5.solve(e5.Side(a,rr,0,pl,*ES),e5.Side(a,rr,0,1e5,*ES),0.6,t,d[:,0]); print('  5.%s: L1 u %.2f rho %.2f p/p* %.4f max p/p* %.3f | %s'%(nm[1],np.mean(abs(d[:,2]-o['u'])),np.mean(abs(d[:,3]-o['rho'])),np.mean(abs(d[:,1]-o['p']))/ps,d[:,1][d[:,0]>0.65].max()/ps,st(nm,T)))
    ref=np.loadtxt(S+'su/sm0.003125_%s.dat'%T); pr=None; out=[]
    for c in ['0.2','0.1','0.05']:
        d=np.loadtxt(S+'su/sm%s_%s.dat'%(c,T)); er=np.mean(abs(d[:,3]-ref[:,3]))/np.mean(ref[:,3]); out.append('%.2e'%er+('' if pr is None else ' (%.2f)'%np.log2(pr/er))); pr=er
    print('  smooth, time refinement:',' | '.join(out))
