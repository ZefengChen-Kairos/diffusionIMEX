import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, exactfan as xf
S='/tmp/claude-0/-home-claude/31e0ebec-a735-5ef0-b582-99f4c98f4569/scratchpad/ex2/'
X,P,U,Y,A=xf.full(500,0.58e-3)
f,ax=plt.subplots(1,2,figsize=(13,5.4))
a=ax[0]
for N,c in [(500,'#f4b183'),(1000,'#e07b39'),(2000,'#a84300')]:
    d=np.loadtxt(S+'e_%d.dat'%N); a.plot(d[:,0],d[:,1],color=c,lw=1.3,label='explicit HLLC (paper-type method), N = %d'%N)
for N,c in [(500,'#9ecae1'),(2000,'#08519c')]:
    d=np.loadtxt(S+'i_%d.dat'%N); a.plot(d[:,0],d[:,1],color=c,lw=1.3,label='IMEX, N = %d'%N)
a.plot(X,P,'k-',lw=1.0,label='exact (equilibrium model)')
a.plot([0.5],[1.0e4],'s',color='k',mfc='none',ms=8,label='P-S Fig. 2 centre value, N = 5000 (read off, ≈ 1×10⁴ Pa)')
a.set_yscale('log'); a.set_xlim(0,1); a.set_ylim(1e2,3e7); a.set_xlabel('x [m]'); a.set_ylabel('p [Pa]'); a.legend(fontsize=7.5,loc='upper center',ncol=2); a.set_title('Pressure, p–pT–pTg, log scale',fontsize=10)
b=ax[1]; Ns=np.array([500,1000,2000]); pe=np.array([np.loadtxt(S+'e_%d.dat'%N)[N//2,1] for N in Ns]); pi=np.array([np.loadtxt(S+'i_%d.dat'%N)[N//2,1] for N in Ns])
k,c0=np.polyfit(np.log(Ns),np.log(pe),1); NN=np.array([500,5000]); b.plot(NN,np.exp(c0)*NN**k,'--',color='#e07b39',lw=1,label='power-law fit to explicit, slope %.2f'%k)
b.plot(Ns,pe,'o-',color='#a84300',label='explicit, centre pressure'); b.plot(Ns,pi,'o-',color='#08519c',label='IMEX, centre pressure')
b.plot([5000],[1.0e4],'s',color='k',mfc='none',ms=8,label='P-S Fig. 2 (read off)'); b.axhline(P[len(P)//2],color='k',lw=1,label='exact, 664 Pa')
b.set_xscale('log'); b.set_yscale('log'); b.set_xlabel('number of cells N'); b.set_ylabel('p at x = 0.5 [Pa]'); b.legend(fontsize=8); b.set_title('Centre pressure against resolution',fontsize=10); b.set_xticks([500,1000,2000,5000]); b.set_xticklabels(['500','1000','2000','5000'])
f.suptitle('Water cavitation tube, |u| = 500 m/s, t = 0.58 ms: why the pressure differs from P-S Fig. 2',fontsize=10.5); f.tight_layout(); f.savefig('ps_fig2_pressure.png',dpi=115)
