import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ST=[(0,'p relaxation',dict(color='m',ls='-',lw=1.2)),(1,'p–pT',dict(color='r',ls='-.',lw=1.2)),(2,'p–pT–pTg',dict(color='b',ls='none',marker='.',ms=3.5)),(3,'p–pTg',dict(color='g',ls='--',lw=1.4))]
def fig(C,title,lims,out,ysc):
    f,ax=plt.subplots(2,2,figsize=(10,7.6));ax=ax.ravel()
    for j,(c,lab) in enumerate([(1,'p [Pa]'),(2,'u [m/s]'),(4,r'$\alpha_v$'),(5,r'$Y_v$')]):
        for R,l,kw in ST:
            d=np.loadtxt('f%d_%d.dat'%(C,R));ax[j].plot(d[:,0],d[:,c]/ysc[j],label=l,**kw)
        ax[j].set_xlim(0,1);ax[j].set_ylim(*lims[j]);ax[j].set_xlabel('x [m]');ax[j].set_ylabel(lab);ax[j].legend(fontsize=8,loc=lims[j][2] if len(lims[j])>2 else 'best')
        if ysc[j]!=1:ax[j].text(0,1.01,r'$\times 10^{%d}$'%np.log10(ysc[j]),transform=ax[j].transAxes,fontsize=9)
    f.suptitle(title,fontsize=10);f.tight_layout();f.savefig(out,dpi=120)
fig(1,'Water cavitation tube, |u| = 2 m/s, t = 3.2 ms — 6+1 model, IMEX 2nd order, N = 1000, CFL = 0.2 (32 steps)',[(0,10),(-2,2),(0,0.35),(0,1.4)],'ps_fig1_imex.png',[1e4,1,1,1e-4])
fig(1,'Same, own axis range (the P-S range clips the centre)',[(0,10),(-2,2),(0,0.8),(0,8)],'ps_fig1_imex_full.png',[1e4,1,1,1e-4])
fig(2,'Water cavitation tube, |u| = 500 m/s, t = 0.58 ms — 6+1 model, IMEX 2nd order, N = 1000, CFL = 0.2 (1450 steps)',[(0,10),(-500,500),(0,1),(0,0.16)],'ps_fig2_imex.png',[1e4,1,1,1])
