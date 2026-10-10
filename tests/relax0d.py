"""0D relaxation tests for src/SPB6_TV_IMEX.cpp (case 40: uniform non-equilibrium state, u = 0, periodic).
With no gradients the IMEX stage reduces to the local closure, so the runs test the relaxation part alone:
  (1) frozen limit  : rate = 0  -> state unchanged (bit-for-bit up to round-off)
  (2) AP limit      : rate -> inf at fixed dt -> converges to the instantaneous closure, error O(1/(tau rate))
  (3) finite rate   : time convergence of the IMEX-RK against a Radau solution of the relaxation ODE
Reference ODEs (SG water/vapour of Pelanti-Shyue 2014, Table 1):
  finite mu      : d alpha1/dt = mu (p1-p2),  d eps1/dt = -pI mu (p1-p2),  pI = (Z2 p1 + Z1 p2)/(Z1+Z2)
  finite theta   : mu = inf; each phase isentropic except the heat Q = theta (T2-T1):  m_k T_k ds_k/dt = +-Q
  finite nu      : mu = theta = inf; dY2/dt = -(nu/rho)(g2 - g1) on the (p,T)-equilibrium manifold
usage: python3 tests/relax0d.py [BIN]"""
import numpy as np, subprocess, sys, os, tempfile
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

G = [2.35, 1.43]; PI = [1e9, 0.0]; Q = [-1167e3, 2030e3]; QP = [0.0, -23.4e3]; CV = [1816.0, 1040.0]
def v_pT(k, p, T): return (G[k] - 1) * CV[k] * T / (p + PI[k])
def e_pT(k, p, T): return CV[k] * T * (p + G[k] * PI[k]) / (p + PI[k]) + Q[k]
def s_pT(k, p, T): return CV[k] * (G[k] * np.log(T) - (G[k] - 1) * np.log(p + PI[k])) + QP[k]
def g_pT(k, p, T): return G[k] * CV[k] * T + Q[k] - T * s_pT(k, p, T)
def p_ame(k, a, m, e): return (G[k] - 1) * (e - m * Q[k]) / a - G[k] * PI[k]
def T_amp(k, a, m, p): return (a / m) * (p + PI[k]) / ((G[k] - 1) * CV[k])
# initial state of case 40
a2 = 0.2; a1 = 1 - a2; T1, T2, p1, p2 = 350.0, 380.0, 1.2e5, 0.8e5
m1 = a1 / v_pT(0, p1, T1); m2 = a2 / v_pT(1, p2, T2); e1 = m1 * e_pT(0, p1, T1); e2 = m2 * e_pT(1, p2, T2)
rho = m1 + m2; E = e1 + e2

def ref_mu(mu, t):
    def f(_, y):
        a, e = y; pa = p_ame(0, a, m1, e); pb = p_ame(1, 1 - a, m2, E - e)
        Z1 = np.sqrt(G[0] * (pa + PI[0]) * m1 / a); Z2 = np.sqrt(G[1] * (pb + PI[1]) * m2 / (1 - a))
        pI = (Z2 * pa + Z1 * pb) / (Z1 + Z2); R = mu * (pa - pb); return [R, -pI * R]
    s = solve_ivp(f, [0, t], [a1, e1], method='Radau', rtol=1e-12, atol=[1e-15, 1e-6])
    a, e = s.y[:, -1]; return dict(p1=p_ame(0, a, m1, e), p2=p_ame(1, 1 - a, m2, E - e), a2=1 - a)

def p_equil_from_s(s1, s2):   # p such that m1 v1(p,s1) + m2 v2(p,s2) = 1 (SG isentropes)
    def T_ps(k, p, s): return np.exp((s - QP[k] + CV[k] * (G[k] - 1) * np.log(p + PI[k])) / (G[k] * CV[k]))
    f = lambda lp: m1 * v_pT(0, np.exp(lp), T_ps(0, np.exp(lp), s1)) + m2 * v_pT(1, np.exp(lp), T_ps(1, np.exp(lp), s2)) - 1
    p = np.exp(brentq(f, np.log(1.0), np.log(1e9), xtol=1e-14)); return p, T_ps(0, p, s1), T_ps(1, p, s2)
def ref_theta(theta, t, s0):
    def f(_, y):
        p, Ta, Tb = p_equil_from_s(*y); Qh = theta * (Tb - Ta); return [Qh / (m1 * Ta), -Qh / (m2 * Tb)]
    s = solve_ivp(f, [0, t], list(s0), method='Radau', rtol=1e-12, atol=1e-12)
    p, Ta, Tb = p_equil_from_s(*s.y[:, -1]); return dict(p=p, T1=Ta, T2=Tb)

def pT_state(Y):              # (p,T) equilibrium at fixed rho, rho e, Y
    ee = E / rho
    def T_of_p(p): return (ee - (1 - Y) * Q[0] - Y * Q[1]) / ((1 - Y) * CV[0] * (p + G[0] * PI[0]) / (p + PI[0]) + Y * CV[1] * (p + G[1] * PI[1]) / (p + PI[1]))
    f = lambda lp: rho * ((1 - Y) * v_pT(0, np.exp(lp), T_of_p(np.exp(lp))) + Y * v_pT(1, np.exp(lp), T_of_p(np.exp(lp)))) - 1
    p = np.exp(brentq(f, np.log(1.0), np.log(1e10), xtol=1e-14)); return p, T_of_p(p)
def ref_nu(nu, t, Y0):
    def f(_, y):
        p, T = pT_state(y[0]); return [-(nu / rho) * (g_pT(1, p, T) - g_pT(0, p, T))]
    s = solve_ivp(f, [0, t], [Y0], method='Radau', rtol=1e-12, atol=1e-14)
    p, T = pT_state(s.y[0, -1]); return dict(p=p, T=T, Y=s.y[0, -1])

BIN = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src', 'spb6')
def run(args):
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, 'o.dat'); r = subprocess.run([BIN, 'case=40', 'N=4', 'mode=1', 'out=' + f] + args, capture_output=True, text=True)
        if r.returncode != 0: return None
        x = np.loadtxt(f); return dict(p=x[0, 2], a2=x[0, 6], Y=x[0, 7], T1=x[0, 12], T2=x[0, 13], p1=x[0, 14], p2=x[0, 15])

tend = 1e-4
print('initial: p1 %.6e p2 %.6e T1 %.3f T2 %.3f alpha2 %.4f Y2 %.6e' % (p1, p2, T1, T2, a2, m2 / rho))
# ---------------- (1) frozen limit
o = run(['relax=0', 'mu=0', 'rk=2', 'dt=1e-5', 'tend=%g' % tend])
print('\n(1) frozen, mu = 0 (RK2): p1 %.10e p2 %.10e alpha2 %.10f  | change p1 %.2e p2 %.2e' % (o['p1'], o['p2'], o['a2'], o['p1'] - p1, o['p2'] - p2))
o = run(['relax=1', 'theta=0', 'rk=2', 'dt=1e-5', 'tend=%g' % tend])
print('    theta = 0 (mu = inf): T1 %.6f T2 %.6f (pressure-equilibrated, temperatures stay apart)' % (o['T1'], o['T2']))
# ---------------- (2) AP limit (mu)
print('\n(2) AP limit at fixed dt = 1e-5 (one RK2 step sequence), mu -> inf: |p1-p2| and distance to the mu = inf result')
oinf = run(['relax=0', 'mu=inf', 'rk=2', 'dt=1e-5', 'tend=%g' % tend])
for mu in [1e-6, 1e-4, 1e-2, 1, 1e2, 1e4]:
    o = run(['relax=0', 'mu=%g' % mu, 'rk=2', 'dt=1e-5', 'tend=%g' % tend])
    if o is None: print('   mu %g failed' % mu); continue
    print('   mu %8.0e: |p1-p2| %.3e  |p - p(inf)| %.3e  |alpha2 - alpha2(inf)| %.3e' % (mu, abs(o['p1'] - o['p2']), abs(o['p'] - oinf['p']), abs(o['a2'] - oinf['a2'])))
print('   AP limit for theta and nu:')
oT = run(['relax=1', 'theta=inf', 'rk=2', 'dt=1e-5', 'tend=%g' % tend]); print('     theta = inf: T1 %.6f T2 %.6f p %.6e' % (oT['T1'], oT['T2'], oT['p']))
for th in [1e4, 1e6, 1e8, 1e10]:
    o = run(['relax=1', 'theta=%g' % th, 'rk=2', 'dt=1e-5', 'tend=%g' % tend])
    if o: print('     theta %8.0e: |T1-T2| %.3e  |p - p(inf)| %.3e' % (th, abs(o['T1'] - o['T2']), abs(o['p'] - oT['p'])))
# ---------------- (3) finite-rate convergence in dt
def conv(label, args, refval, keys):
    print('\n(3) %s: error vs dt (RK1 / RK2 / RK3)' % label)
    for rk in [1, 2, 3]:
        errs = []
        for dt in [2e-5, 1e-5, 5e-6, 2.5e-6, 1.25e-6]:
            o = run(args + ['rk=%d' % rk, 'dt=%g' % dt, 'tend=%g' % tend])
            if o is None: errs.append(np.nan); continue
            errs.append(max(abs(o[k] - refval[k]) / abs(refval[k]) for k in keys))
        orders = [np.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)]
        print('   rk%d: ' % rk + ' '.join('%.2e' % e for e in errs) + '   orders ' + ' '.join('%.2f' % q for q in orders))
mu = 1e-5
r = ref_mu(mu, tend); print('\nreference finite mu = %g: p1 %.8e p2 %.8e alpha2 %.8f' % (mu, r['p1'], r['p2'], r['a2']))
conv('finite mu = %g (theta = nu = 0)' % mu, ['relax=0', 'mu=%g' % mu], r, ['p1', 'p2', 'a2'])
# entropies of the initial state after instantaneous mechanical relaxation are what the theta-ODE starts from: take them from the scheme's
# own mu = inf state at t ~ 0 (one tiny step) so that only the thermal part is compared
o0 = run(['relax=1', 'theta=0', 'rk=1', 'dt=1e-9', 'tend=1e-9'])
s0 = (s_pT(0, o0['p'], o0['T1']), s_pT(1, o0['p'], o0['T2']))
th = 5e6
r = ref_theta(th, tend - 1e-9, s0); print('\nreference finite theta = %g (from the p-relaxed state): p %.8e T1 %.6f T2 %.6f' % (th, r['p'], r['T1'], r['T2']))
conv('finite theta = %g (mu = inf)' % th, ['relax=1', 'theta=%g' % th], r, ['p', 'T1', 'T2'])

# ---------------- finite nu (case 41: pT equilibrium, superheated, T > Tsat): mu = theta = inf
T0, p0 = 380.0, 0.8e5
m1 = 0.8 / v_pT(0, p0, T0); m2 = 0.2 / v_pT(1, p0, T0); rho = m1 + m2; E = m1 * e_pT(0, p0, T0) + m2 * e_pT(1, p0, T0)
def run41(args):
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, 'o.dat'); r = subprocess.run([BIN, 'case=41', 'N=4', 'mode=1', 'out=' + f] + args, capture_output=True, text=True)
        if r.returncode != 0: return None
        x = np.loadtxt(f); return dict(p=x[0, 2], Y=x[0, 7], T=x[0, 12])
nu = 1e-3
r = ref_nu(nu, tend, m2 / rho); print('\nreference finite nu = %g: p %.8e T %.6f Y2 %.8e' % (nu, r['p'], r['T'], r['Y']))
print('(3) finite nu = %g (mu = theta = inf): error vs dt (RK1 / RK2 / RK3)' % nu)
for rk in [1, 2, 3]:
    errs = []
    for dt in [2e-5, 1e-5, 5e-6, 2.5e-6, 1.25e-6]:
        o = run41(['relax=2', 'nu=%g' % nu, 'rk=%d' % rk, 'dt=%g' % dt, 'tend=%g' % tend])
        errs.append(np.nan if o is None else max(abs(o['p'] - r['p']) / r['p'], abs(o['Y'] - r['Y']) / r['Y']))
    print('   rk%d: ' % rk + ' '.join('%.2e' % e for e in errs) + '   orders ' + ' '.join('%.2f' % np.log2(errs[i] / errs[i + 1]) for i in range(len(errs) - 1)))
