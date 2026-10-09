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
s=solve_ivp(r1,[p0,1e3],[T0,-2.0],events=ev,rtol=1e-10,atol=1e-12)
ps,Tsn,us=s.t[-1],s.y[0,-1],s.y[1,-1]
print('pT c at IC %.2f m/s; saturation onset: p %.5e T %.3f u %.4f  c_pT %.2f'%(c_at(p0,T0),ps,Tsn,us,c_at(ps,Tsn)))
# stage 2: pTg isentrope, state (Y,u) vs p
def st(p,y): T=Ts(p);return mix(p,T,y)
def r2(p,s):
    y=s[0];h=1e-5
    a=st(p*(1+h),y);b=st(p*(1-h),y);vp=(a[0]-b[0])/(2*h*p);ep=(a[1]-b[1])/(2*h*p)
    hy=1e-7;a=st(p,y+hy);b=st(p,y-hy);vy=(a[0]-b[0])/(2*hy);ey=(a[1]-b[1])/(2*hy)
    dy=-(ep+p*vp)/(ey+p*vy);vv=st(p,y)[0];c=np.sqrt(-vv*vv/(vp+vy*dy));return [dy,-vv/c]
r2.c=lambda p,y:(-st(p,y)[0]/r2(p,[y,0])[1])
ev2=lambda p,s:s[1];ev2.terminal=True
s2=solve_ivp(r2,[ps,1e3],[Y,us],events=ev2,rtol=1e-10,atol=1e-14)
pc,yc=s2.t[-1],s2.y[0,-1];vc=st(pc,yc)[0];T=Ts(pc)
print('centre: p %.5e Yv %.4e av %.4f T %.3f rho %.2f'%(pc,yc,yc*v(1,pc,T)/vc,T,1/vc))
print('c_HEM at onset %.3f, at centre %.3f m/s; fan head x=0.5-/+ (c-u)t: head %.4f tail %.4f (t=3.2ms)'%(r2.c(ps,Y),r2.c(pc,yc),0.5-(r2.c(ps,Y)-us)*3.2e-3,0.5-r2.c(pc,yc)*3.2e-3))
