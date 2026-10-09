import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
G=[2.35,1.43];PI=[1e9,0.];Q=[-1167e3,2030e3];CV=[1816.,1040.];QP=[0.,-23.4e3]
v=lambda k,p,T:(G[k]-1)*CV[k]*T/(p+PI[k]); e=lambda k,p,T:CV[k]*T*(p+G[k]*PI[k])/(p+PI[k])+Q[k]
g=lambda k,p,T:(G[k]*CV[k]-QP[k])*T-CV[k]*T*np.log(T**G[k]/(p+PI[k])**(G[k]-1))+Q[k]
Ts=lambda p:brentq(lambda T:g(0,p,T)-g(1,p,T),250,700,xtol=1e-12)
p0,T0,av=1e5,354.728,1e-2
m=[(1-av)/v(0,p0,T0),av/v(1,p0,T0)];rho0=sum(m);Y=m[1]/rho0
mix=lambda p,T,y:((1-y)*v(0,p,T)+y*v(1,p,T),(1-y)*e(0,p,T)+y*e(1,p,T))
# stage 1: pT isentrope, state (T,u) vs p (decreasing)
def r1(p,s):
    T=s[0];h=1e-5
    vp=(mix(p*(1+h),T,Y)[0]-mix(p*(1-h),T,Y)[0])/(2*h*p);ep=(mix(p*(1+h),T,Y)[1]-mix(p*(1-h),T,Y)[1])/(2*h*p)
    vT=(mix(p,T*(1+h),Y)[0]-mix(p,T*(1-h),Y)[0])/(2*h*T);eT=(mix(p,T*(1+h),Y)[1]-mix(p,T*(1-h),Y)[1])/(2*h*T)
    dT=-(ep+p*vp)/(eT+p*vT);vv=mix(p,T,Y)[0];dv=vp+vT*dT;c=np.sqrt(-vv*vv/dv)
    return [dT,1/(c/vv)*(-1)*(-1)*(-1)]   # du/dp = -1/(rho c) for left wave as p decreases u increases
ev=lambda p,s:s[0]-Ts(p);ev.terminal=True
c_at=lambda p,T:(lambda vp,vT,dT,vv:np.sqrt(-vv*vv/(vp+vT*dT)))(*(lambda h:((mix(p*(1+h),T,Y)[0]-mix(p*(1-h),T,Y)[0])/(2*h*p),(mix(p,T*(1+h),Y)[0]-mix(p,T*(1-h),Y)[0])/(2*h*T),r1(p,[T,0])[0],mix(p,T,Y)[0]))(1e-5))
s=solve_ivp(r1,[p0,1e3],[T0,-500.0],events=ev,rtol=1e-10,atol=1e-12)
ps,Tsn,us=s.t[-1],s.y[0,-1],s.y[1,-1]
# stage 2: pTg isentrope, state (Y,u) vs p
def st(p,y): T=Ts(p);return mix(p,T,y)
def r2(p,s):
    y=s[0];h=1e-5
    a=st(p*(1+h),y);b=st(p*(1-h),y);vp=(a[0]-b[0])/(2*h*p);ep=(a[1]-b[1])/(2*h*p)
    hy=1e-7;a=st(p,y+hy);b=st(p,y-hy);vy=(a[0]-b[0])/(2*hy);ey=(a[1]-b[1])/(2*hy)
    dy=-(ep+p*vp)/(ey+p*vy);vv=st(p,y)[0];c=np.sqrt(-vv*vv/(vp+vy*dy));return [dy,-vv/c]
r2.c=lambda p,y:(-st(p,y)[0]/r2(p,[y,0])[1])
ev2=lambda p,s:s[1];ev2.terminal=True

def profile(U0,t,n=600):
    """exact equilibrium (pT -> pTg) solution of the symmetric expansion tube, left half; returns x,p,u,Yv,alpha_v sorted in x (0..0.5)"""
    s=solve_ivp(r1,[p0,1e3],[T0,-U0],events=ev,rtol=1e-10,atol=1e-12,dense_output=True)
    ps,Tsn,us=s.t[-1],s.y[0,-1],s.y[1,-1]
    P1=np.geomspace(p0,ps,n); S1=s.sol(P1); x1=np.array([0.5+(S1[1,i]-c_at(P1[i],S1[0,i]))*t for i in range(n)])
    a1=np.array([Y*v(1,P1[i],S1[0,i])/mix(P1[i],S1[0,i],Y)[0] for i in range(n)])
    s2=solve_ivp(r2,[ps,1.0],[Y,us],events=ev2,rtol=1e-10,atol=1e-14,dense_output=True); pc=s2.t[-1]
    P2=np.geomspace(ps,pc,n); S2=s2.sol(P2); x2=np.array([0.5+(S2[1,i]-r2.c(P2[i],S2[0,i]))*t for i in range(n)])
    a2=np.array([S2[0,i]*v(1,P2[i],Ts(P2[i]))/st(P2[i],S2[0,i])[0] for i in range(n)])
    x=np.concatenate([[0.0],x1,x2,[0.5]]); p=np.concatenate([[p0],P1,P2,[pc]]); u=np.concatenate([[-U0],S1[1],S2[1],[0.0]])
    y=np.concatenate([[Y],np.full(n,Y),S2[0],[S2[0,-1]]]); a=np.concatenate([[av],a1,a2,[a2[-1]]])
    return x,p,u,y,a
def full(U0,t):
    x,p,u,y,a=profile(U0,t); return np.concatenate([x,1-x[::-1]]),np.concatenate([p,p[::-1]]),np.concatenate([u,-u[::-1]]),np.concatenate([y,y[::-1]]),np.concatenate([a,a[::-1]])
