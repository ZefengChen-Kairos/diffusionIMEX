import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, exactfan as e
ST=[(0,'p relaxation',dict(color='m',ls='-',lw=1.2)),(1,'p–pT',dict(color='r',ls='-.',lw=1.2)),(2,'p–pT–pTg',dict(color='b',ls='none',marker='.',ms=3.5)),(3,'p–pTg',dict(color='g',ls='--',lw=1.4))]
def fig(C,U,t,title,lims,out,ysc):
    X=e.full(U,t); ex={1:X[1],2:X[2],4:X[4],5:X[3]}
    f,ax=plt.subplots(2,2,figsize=(10,7.6));ax=ax.ravel()
    for j,(c,lab) in enumerate([(1,'p [Pa]'),(2,'u [m/s]'),(4,r'$\alpha_v$'),(5,r'$Y_v$')]):
        for R,l,kw in ST:
            d=np.loadtxt('f%d_%d.dat'%(C,R));ax[j].plot(d[:,0],d[:,c]/ysc[j],label=l,**kw)
        ax[j].plot(X[0],ex[c]/ysc[j],'k-',lw=0.9,label='exact, p–pT–pTg equilibrium',zorder=5)
        ax[j].set_xlim(0,1);ax[j].set_ylim(*lims[j]);ax[j].set_xlabel('x [m]');ax[j].set_ylabel(lab);ax[j].legend(fontsize=7.5)
        if ysc[j]!=1:ax[j].text(0,1.01,r'$\times 10^{%d}$'%np.log10(ysc[j]),transform=ax[j].transAxes,fontsize=9)
    f.suptitle(title,fontsize=10);f.tight_layout();f.savefig(out,dpi=120)
fig(1,2,3.2e-3,'Water cavitation tube, |u| = 2 m/s, t = 3.2 ms — 6+1 model, IMEX 2nd order, N = 1000, CFL = 0.2; black: exact equilibrium solution',[(0,10),(-2,2),(0,0.8),(0,8)],'ps_fig1_imex_exact.png',[1e4,1,1,1e-4])
fig(2,500,0.58e-3,'Water cavitation tube, |u| = 500 m/s, t = 0.58 ms — 6+1 model, IMEX 2nd order, N = 1000, CFL = 0.2; black: exact equilibrium solution',[(0,10),(-500,500),(0,1),(0,0.16)],'ps_fig2_imex_exact.png',[1e4,1,1,1])
