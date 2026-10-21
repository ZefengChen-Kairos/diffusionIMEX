import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.argv=[sys.argv[0],'none']
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import ana as A
plt.rcParams.update({'font.size':9,'axes.titlesize':9,'legend.fontsize':7})
fig,AX=plt.subplots(2,2,figsize=(10,7.2),gridspec_kw=dict(height_ratios=[1.3,1]))
ax=AX[0]
import l1 as LL
def xs(lst,infx):
    return [infx if t=='inf' else float(t) for t in lst]
# 左：theta
TH=[t for t in A.THS if t!='0']; x=xs(TH,1e12)
for fo,sty,lab in ((0,dict(color='r',ls='-',marker='s',ms=4),'backward Euler (relaxform=0)'),(1,dict(color='k',ls='none',marker='o',mfc='none',ms=7,mew=1),'exponential (relaxform=1)')):
    ax[0].plot(x,[A.ctr(A.thname(t,fo))['p'] for t in TH],label=lab,**sty)
E=A.EX['pel']
ax[0].axhline(E['pT']['p'],color='k',ls='--',lw=0.9,label='exact 4-eq. (pT) %.5g Pa'%E['pT']['p'])
p0=A.ctr(A.PC+'pel821_p.dat')['p']
ax[0].axhline(p0,color='darkviolet',ls=':',lw=1.1,label='p-relax (theta=0), IMEX %.5g Pa'%p0)
ax[0].axhline(E['p5']['p'],color='darkviolet',ls='-.',lw=0.8,label='exact 5-eq. (p) %.5g Pa'%E['p5']['p'])
ax[0].set_xscale('log'); ax[0].set_xlabel(r'$\theta$ (inf plotted at 1e12, unify=1)'); ax[0].set_ylabel('centre pressure (Pa)')
ax[0].set_title('relax=1, pT(f): centre p vs heat-exchange rate'); ax[0].legend(loc='center right')
NU=A.NUS; x=xs(NU,1e6)
for fo,sty,lab in ((0,dict(color='r',ls='-',marker='s',ms=4),'backward Euler (relaxform=0)'),(1,dict(color='k',ls='none',marker='o',mfc='none',ms=7,mew=1),'exponential (relaxform=1)')):
    ax[1].plot(x,[A.ctr(A.nuname(t,fo))['p'] for t in NU],label=lab,**sty)
ax[1].axhline(E['pTg']['p'],color='k',ls='--',lw=0.9,label='exact pT->pTg %.5g Pa'%E['pTg']['p'])
ax[1].axhline(E['pT']['p'],color='k',ls=':',lw=0.9,label='exact 4-eq. (pT) %.5g Pa'%E['pT']['p'])
ax[1].set_xscale('log'); ax[1].set_xlabel(r'$\nu$ (inf plotted at 1e6, unify=1)'); ax[1].set_ylabel('centre pressure (Pa)')
ax[1].set_title('relax=2, pTg(f): centre p vs mass-transfer rate'); ax[1].legend(loc='center right')
fig.suptitle('Pelanti 2022 8.2.1 (NASG water, |u| = 2 m/s, t = 3 ms, 5000 cells): one finite-rate formula from frozen to equilibrium',fontsize=9)
BE=dict(color='r',ls='-',marker='s',ms=4,label='backward Euler'); EXP=dict(color='k',ls='none',marker='o',mfc='none',ms=7,mew=1,label='exponential')
for fo,sty in ((0,BE),(1,EXP)):
    AX[1,0].loglog(xs(TH,1e12),[LL.err(A.thname(t,fo),False)[0] for t in TH],**sty)
    AX[1,1].loglog(xs(NU,1e6),[LL.err(A.nuname(t,fo),True)[0] for t in NU],**sty)
AX[1,0].set_title('L1 error of p vs exact 4-eq. (pT) fan'); AX[1,1].set_title('L1 error of p vs exact pT->pTg fan')
for a,l in ((AX[1,0],r'$\theta$'),(AX[1,1],r'$\nu$')): a.set_xlabel(l+' (inf at right end)'); a.set_ylabel('mean |p - p_exact| (Pa)'); a.legend(); a.grid(True,which='major',lw=0.3)
fig.tight_layout(); out='/home/user/diffusionIMEX/docs/figures/expo_relax.png'; fig.savefig(out,dpi=120); print('wrote',out)
