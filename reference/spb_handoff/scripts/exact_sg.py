# exact Riemann solver, two stiffened gases separated by a contact (left fluid gl,pil; right gr,pir)
import numpy as np
from scipy.optimize import brentq
def solve(rl,ul,pl,gl,pil,rr,ur,pr,gr,pir,x0,t,x):
    def f(p,r,pk,g,pi):
        P=p+pi;Pk=pk+pi;c=np.sqrt(g*Pk/r)
        if p>pk: A=2/((g+1)*r);B=(g-1)/(g+1)*Pk;return (p-pk)*np.sqrt(A/(P+B))
        return 2*c/(g-1)*((P/Pk)**((g-1)/(2*g))-1)
    F=lambda p:f(p,rl,pl,gl,pil)+f(p,rr,pr,gr,pir)+ur-ul
    ps=brentq(F,max(-pil,-pir)*(1-1e-12)+1e-6,10*max(pl,pr),xtol=1e-12,rtol=1e-14)
    us=0.5*(ul+ur)+0.5*(f(ps,rr,pr,gr,pir)-f(ps,rl,pl,gl,pil))
    def side(r,u,pk,g,pi,sgn,xi):   # returns rho,u,p on that side of contact for xi=(x-x0)/t ; sgn=-1 left, +1 right
        Pk=pk+pi;P=ps+pi;c=np.sqrt(g*Pk/r)
        if ps>pk:
            rs=r*((P/Pk+(g-1)/(g+1))/((g-1)/(g+1)*P/Pk+1));S=u+sgn*c*np.sqrt((g+1)/(2*g)*P/Pk+(g-1)/(2*g))
            return np.where(sgn*(xi-S)>0,r,rs),np.where(sgn*(xi-S)>0,u,us),np.where(sgn*(xi-S)>0,pk,ps)
        rs=r*(P/Pk)**(1/g);cs=c*(P/Pk)**((g-1)/(2*g));H=u+sgn*c;Tl=us+sgn*cs
        cf=2/(g+1)*(c-sgn*(g-1)/2*(u-xi)) ; cf=np.maximum(cf,1e-30); rf=r*(cf/c)**(2/(g-1));uf=2/(g+1)*(-sgn*c+(g-1)/2*u+xi);pf=Pk*(cf/c)**(2*g/(g-1))-pi
        out=sgn*(xi-H)>0;ins=sgn*(xi-Tl)<0
        return np.where(out,r,np.where(ins,rs,rf)),np.where(out,u,np.where(ins,us,uf)),np.where(out,pk,np.where(ins,ps,pf))
    xi=(x-x0)/t;L=side(rl,ul,pl,gl,pil,-1,xi);R=side(rr,ur,pr,gr,pir,+1,xi);m=xi<us
    return np.where(m,L[0],R[0]),np.where(m,L[1],R[1]),np.where(m,L[2],R[2]),ps,us
