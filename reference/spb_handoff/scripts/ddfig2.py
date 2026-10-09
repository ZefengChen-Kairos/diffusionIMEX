import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ST=[(0,'p relaxation',dict(color='m',ls='-',lw=1.3)),(2,'p–pT–pTg',dict(color='b',ls='none',marker='.',ms=3.5)),(3,'p–pTg',dict(color='g',ls='--',lw=1.4))]
G=[2.35,1.025];PI=[4e8,0.];CV=[1077.7,1956.45]
f,ax=plt.subplots(3,2,figsize=(10,11.5));ax=ax.ravel()
P=[(3,r'$\rho$ [kg/m$^3$]',1,True,(1,1e3)),(2,'u [m/s]',1,False,(0,350)),(4,r'$\alpha_v$',1,False,(0,1)),(5,r'$Y_v$',1,False,(0,1)),(1,'p [bar]',1e5,True,(1,1e3))]
for j,(c,lab,sc,lg,yl) in enumerate(P):
    for R,l,kw in ST:
        d=np.loadtxt('d_i2_%d.dat'%R);ax[j].plot(d[:,0],d[:,c]/sc,label=l,**kw)
    if lg:ax[j].set_yscale('log')
    ax[j].set_xlim(0,1);ax[j].set_ylim(*yl);ax[j].set_xlabel('x [m]');ax[j].set_ylabel(lab);ax[j].legend(fontsize=8)
a=ax[5]
for R,l,kw in ST:
    d=np.loadtxt('d_i2_%d.dat'%R);p=d[:,1];rv=d[:,5]*d[:,3]/d[:,4];rl=(1-d[:,5])*d[:,3]/(1-d[:,4])
    a.plot(d[:,0],(p+PI[1])/((G[1]-1)*CV[1]*rv)-(p+PI[0])/((G[0]-1)*CV[0]*rl),label=l,**kw)
a.set_xlim(0,1);a.set_xlabel('x [m]');a.set_ylabel(r'$T_v-T_l$ [K]');a.legend(fontsize=8)
f.suptitle('Dodecane liquid–vapor shock tube, t = 473 μs — 6+1 model, IMEX 2nd order, N = 1000, CFL = 0.2',fontsize=10);f.tight_layout();f.savefig('ps_fig3_imex.png',dpi=115)
