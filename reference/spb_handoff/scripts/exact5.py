# exact Riemann solver of the 5-equation (mechanical equilibrium) model for a two-SG mixture with uniform mass fractions on each side,
# shock relations of Saurel et al. (II.5): each phase follows its own Hugoniot at the common pressure; rarefactions: each phase isentropic.
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
class Side:
    def __init__(s,a,r,u,p,G,PI):   # a, r: phase volume fractions and densities
        s.a=np.array(a,float); s.r=np.array(r,float); s.u=u; s.p=p; s.G=np.array(G,float); s.PI=np.array(PI,float)
        s.rho=(s.a*s.r).sum(); s.Y=s.a*s.r/s.rho; s.v0=1/s.r
    def vk(s,p,shock):
        P=p+s.PI; P0=s.p+s.PI; g=s.G
        if shock: return s.v0*((g-1)*P+(g+1)*P0)/((g+1)*P+(g-1)*P0)
        return s.v0*(P0/P)**(1/g)
    def v(s,p,shock): return (s.Y*s.vk(p,shock)).sum()
    def cw_rho(s,p):  # rho*c (Wood) along the isentrope
        vk=s.vk(p,False); dv=-(s.Y*vk/(s.G*(p+s.PI))).sum(); return np.sqrt(-1/dv)
    def f(s,p):       # |u* - u| with sign: f>0 for compression
        if p>s.p: return np.sqrt((p-s.p)*(1/s.rho-s.v(p,True)))
        return -quad(lambda z:np.exp(z)/s.cw_rho(np.exp(z)),np.log(p),np.log(s.p),epsabs=0,epsrel=1e-10,limit=200)[0]   # integrate in ln p (p may span many decades)
def solve(L,R,x0,t,x):
    F=lambda p:L.f(p)+R.f(p)+R.u-L.u
    lo=max(max(-L.PI.min(),-R.PI.min())*(1-1e-9),0.0)+1e-10; vac=F(lo)>0   # vac: the two rarefactions do not meet (escape velocity exceeded) -> vacuum between the two tails
    if vac: ps=lo; usL=L.u-L.f(lo); usR=R.u+R.f(lo); us=0.5*(usL+usR)
    else: ps=brentq(F,lo,10*max(L.p,R.p),xtol=1e-9,rtol=1e-13); us=0.5*(L.u+R.u)+0.5*(R.f(ps)-L.f(ps)); usL=usR=us
    xi=(x-x0)/t; out={k:np.zeros_like(x) for k in ['p','u','rho','a1','r1','r2']}
    def fill(m,S,p,u):
        sh=p>S.p*(1+1e-14); vk=S.vk(p,sh) if p!=S.p else S.v0; v=(S.Y*vk).sum(); out['p'][m]=p; out['u'][m]=u; out['rho'][m]=1/v; out['a1'][m]=S.Y[0]*vk[0]/v; out['r1'][m]=1/vk[0]; out['r2'][m]=1/vk[1]
    for S,sg,ust in ((L,-1,usL),(R,1,usR)):
        side=(xi<us) if sg<0 else (xi>=us)
        if ps>S.p:
            W=S.u+sg*np.sqrt((ps-S.p)/(1/S.rho-S.v(ps,True)))/S.rho; fill(side&(sg*(xi-W)>0),S,S.p,S.u); fill(side&(sg*(xi-W)<=0),S,ps,us)
        else:
            c0=S.cw_rho(S.p)/S.rho; vs=S.v(ps,False); cs=S.cw_rho(ps)*vs; H=S.u+sg*c0; T=ust+sg*cs
            fill(side&(sg*(xi-H)>0),S,S.p,S.u); fill(side&(sg*(xi-T)<=0),S,ps,ust)
            fan=side&(sg*(xi-H)<=0)&(sg*(xi-T)>0)
            if fan.any():
                P=np.geomspace(S.p,ps,1500); U=np.array([S.u+sg*S.f(q) for q in P]); C=np.array([S.cw_rho(q)*S.v(q,False) for q in P]); X=U+sg*C
                o=np.argsort(X); pf=np.interp(xi[fan],X[o],P[o])
                for i,q in zip(np.where(fan)[0],pf): fill(np.array([i]),S,q,np.interp(q,P[::-1],U[::-1]))
    if vac:
        m=(xi>usL)&(xi<usR); out['u'][m]=np.nan; out['p'][m]=0.0; out['rho'][m]=0.0; out['a1'][m]=0.0
    return out,ps,us
