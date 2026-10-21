import sys, numpy as np, json
sys.path.insert(0,'/home/user/diffusionIMEX/tests')
import exact5 as e5, exact_relaxed as er
def ic_rho(P, p, T): return [1 / er.Fan(P, p, T, 0.5).v(k, p, T) for k in (0, 1)]
def ex5min(P, L, R, t):
    x=np.linspace(0,1,4001); G=[P[0][0],P[1][0]]; PI=[P[0][1],P[1][1]]; B=[P[0][5],P[1][5]]
    side=lambda a0,r0,r1,u,p: e5.Side([a0,1-a0],[r0,r1],u,p,G,PI,B)
    o,ps,us=e5.solve(side(*L),side(*R),0.5,t,x); m=len(x)//2
    r=[o['r1'][m],o['r2'][m]]; p=o['p'][m]
    T=[(1/r[k]-P[k][5])*(p+P[k][1])/((P[k][0]-1)*P[k][4]) for k in (0,1)]
    return dict(p=float(p),a2=float(1-o['a1'][m]),T1=float(T[0]),T2=float(T[1]))
out={}
PEL=[(1.187,7028e5,-1177788.0,0.0,3610.0,6.61e-4),(1.467,0.0,2077616.0,14317.0,955.0,0.0)]
PS=[(2.35,1e9,-1167e3,0.0,1816.0,0.0),(1.43,0.0,2030e3,-23.4e3,1040.0,0.0)]
for nm,P,T0,U0,t in (('pel',PEL,353.0,2.,3e-3),('ps1',PS,354.728,2.,3.2e-3),('ps2',PS,354.728,500.,0.58e-3)):
    r0=ic_rho(P,1e5,T0); d={}
    try: d['p5']=ex5min(P,(0.99,r0[0],r0[1],-U0,1e5),(0.99,r0[0],r0[1],U0,1e5),t)
    except Exception as ex: d['p5']=str(ex)
    F=er.Fan(P,1e5,T0,1e-2)
    for chem in (False,True):
        try:
            (x,p,u,Y,a,T),I=F.full(U0,t,chem); m=len(x)//2-1
            d['pTg' if chem else 'pT']=dict(p=float(p[m]),Y=float(Y[m]),a2=float(a[m]),T=float(T[m]),centre=bool(I['centre']))
        except Exception as ex: d['pTg' if chem else 'pT']=str(ex)
    out[nm]=d
print(json.dumps(out,indent=1)); json.dump(out,open(sys.argv[1],'w'),indent=1)
