"""Finite-rate closure forms in the implicit stage, 0D (case 40, u = 0): backward Euler (relaxform=0) against the frozen-coefficient
exponential (relaxform=1). Left: distance to the infinite-rate state at fixed dt = 1e-5, t = 1e-4, against the rate. Right: LSDIRK2 error
against a Radau solution of the relaxation ODE, against dt (mu = 1e-5; theta = 5e6, from the p-relaxed state).
usage: python3 tests/plot_relax_forms.py BIN OUT.png"""
import sys, os, subprocess, tempfile
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
BIN, OUT = sys.argv[1], sys.argv[2]
G = [2.35, 1.43]; PI = [1e9, 0.0]; Q = [-1167e3, 2030e3]; QP = [0.0, -23.4e3]; CV = [1816.0, 1040.0]
def v_pT(k, p, T): return (G[k] - 1) * CV[k] * T / (p + PI[k])
def e_pT(k, p, T): return CV[k] * T * (p + G[k] * PI[k]) / (p + PI[k]) + Q[k]
def s_pT(k, p, T): return CV[k] * (G[k] * np.log(T) - (G[k] - 1) * np.log(p + PI[k])) + QP[k]
def p_ame(k, a, m, e): return (G[k] - 1) * (e - m * Q[k]) / a - G[k] * PI[k]
a2 = 0.2; a1 = 1 - a2; T1, T2, p1, p2 = 350.0, 380.0, 1.2e5, 0.8e5
m1 = a1 / v_pT(0, p1, T1); m2 = a2 / v_pT(1, p2, T2); e1 = m1 * e_pT(0, p1, T1); e2 = m2 * e_pT(1, p2, T2); E = e1 + e2
def ref_mu(mu, t):
    def f(_, y):
        a, e = y; pa = p_ame(0, a, m1, e); pb = p_ame(1, 1 - a, m2, E - e)
        Z1 = np.sqrt(G[0] * (pa + PI[0]) * m1 / a); Z2 = np.sqrt(G[1] * (pb + PI[1]) * m2 / (1 - a))
        pI = (Z2 * pa + Z1 * pb) / (Z1 + Z2); R = mu * (pa - pb); return [R, -pI * R]
    s = solve_ivp(f, [0, t], [a1, e1], method='Radau', rtol=1e-12, atol=[1e-15, 1e-6])
    a, e = s.y[:, -1]; return dict(p1=p_ame(0, a, m1, e), p2=p_ame(1, 1 - a, m2, E - e), a2=1 - a)
def p_equil_from_s(s1, s2):
    def T_ps(k, p, s): return np.exp((s - QP[k] + CV[k] * (G[k] - 1) * np.log(p + PI[k])) / (G[k] * CV[k]))
    f = lambda lp: m1 * v_pT(0, np.exp(lp), T_ps(0, np.exp(lp), s1)) + m2 * v_pT(1, np.exp(lp), T_ps(1, np.exp(lp), s2)) - 1
    p = np.exp(brentq(f, np.log(1.0), np.log(1e9), xtol=1e-14)); return p, T_ps(0, p, s1), T_ps(1, p, s2)
def ref_theta(theta, t, s0):
    def f(_, y):
        p, Ta, Tb = p_equil_from_s(*y); Qh = theta * (Tb - Ta); return [Qh / (m1 * Ta), -Qh / (m2 * Tb)]
    s = solve_ivp(f, [0, t], list(s0), method='Radau', rtol=1e-12, atol=1e-12)
    p, Ta, Tb = p_equil_from_s(*s.y[:, -1]); return dict(p=p, T1=Ta, T2=Tb)
def run(args):
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, 'o.dat'); r = subprocess.run([BIN, 'case=40', 'N=4', 'out=' + f] + args, capture_output=True, text=True)
        if r.returncode != 0: return None
        x = np.loadtxt(f); return dict(p=x[0, 2], a2=x[0, 6], T1=x[0, 12], T2=x[0, 13], p1=x[0, 14], p2=x[0, 15])
tend = 1e-4; FORMS = [(0, 'backward Euler'), (1, 'exponential')]
# limit at fixed dt
mus = np.logspace(-7, 3, 21); ths = np.logspace(2, 12, 21)
oinf = run(['relax=0', 'mu=inf', 'rk=2', 'dt=1e-5', 'tend=%g' % tend]); tinf = run(['relax=1', 'theta=inf', 'rk=2', 'dt=1e-5', 'tend=%g' % tend])
lim = {}
for f, _ in FORMS:
    lim[('mu', f)] = [abs(run(['relax=0', 'mu=%g' % m, 'rk=2', 'dt=1e-5', 'tend=%g' % tend, 'relaxform=%d' % f])['a2'] - oinf['a2']) for m in mus]
    lim[('th', f)] = [(lambda o: max(abs(o['T1'] - tinf['T1']), abs(o['T2'] - tinf['T2'])) / tinf['T1'])(run(['relax=1', 'theta=%g' % t, 'rk=2', 'dt=1e-5', 'tend=%g' % tend, 'relaxform=%d' % f])) for t in ths]
# convergence in dt
dts = np.array([2e-5, 1e-5, 5e-6, 2.5e-6, 1.25e-6, 6.25e-7])
rm = ref_mu(1e-5, tend)
o0 = run(['relax=1', 'theta=0', 'rk=1', 'dt=1e-9', 'tend=1e-9']); s0 = (s_pT(0, o0['p'], o0['T1']), s_pT(1, o0['p'], o0['T2'])); rt = ref_theta(5e6, tend - 1e-9, s0)
cv = {}
for f, _ in FORMS:
    cv[('mu', f)] = [max(abs(o[k] - rm[k]) / abs(rm[k]) for k in ('p1', 'p2', 'a2')) for o in (run(['relax=0', 'mu=1e-5', 'rk=2', 'dt=%g' % d, 'tend=%g' % tend, 'relaxform=%d' % f]) for d in dts)]
    cv[('th', f)] = [max(abs(o[k] - rt[k]) / abs(rt[k]) for k in ('p', 'T1', 'T2')) for o in (run(['relax=1', 'theta=5e6', 'rk=2', 'dt=%g' % d, 'tend=%g' % tend, 'relaxform=%d' % f]) for d in dts)]
for key, e in cv.items(): print(key, ' '.join('%.2e' % x for x in e), 'orders', ' '.join('%.2f' % np.log2(e[i] / e[i + 1]) for i in range(len(e) - 1)))
for key, e in lim.items(): print('limit', key, ' '.join('%.1e' % x for x in e))
C = {0: dict(color='#2a78d6', marker='o', ls='-'), 1: dict(color='#eb6834', marker='s', ls='--')}
plt.rcParams.update({'font.size': 9, 'legend.fontsize': 8, 'axes.edgecolor': '#888', 'axes.labelcolor': '#222', 'xtick.color': '#444', 'ytick.color': '#444'})
fig, ax = plt.subplots(2, 2, figsize=(10, 7.2))
floor = 1e-16
for r, (rel, xs, xl, yl) in enumerate([('mu', mus, r'$\mu$ (Pa$^{-1}$ s$^{-1}$)', r'$|\alpha_2-\alpha_2^{\mu=\infty}|$'),
                                        ('th', ths, r'$\theta$ (W m$^{-3}$ K$^{-1}$)', r'$\max_k|T_k-T_k^{\theta=\infty}|/T$')]):
    a = ax[r, 0]
    for f, lab in FORMS:
        a.loglog(xs, np.maximum(lim[(rel, f)], floor), lw=1.4, ms=4.5, mfc='none', mew=1.0, label=lab, **C[f])
    a.set_xlabel(xl); a.set_ylabel(yl); a.grid(alpha=0.25, which='major'); a.legend(loc='lower left', frameon=False)
    a.set_title('%s: distance to the infinite-rate closure, $\\Delta t=10^{-5}$' % ('pressure relaxation' if rel == 'mu' else 'thermal relaxation'), fontsize=9)
    b = ax[r, 1]
    for f, lab in FORMS:
        e = np.array(cv[(rel, f)]); q = np.polyfit(np.log(dts), np.log(e), 1)[0]   # least-squares slope over all dt
        b.loglog(dts, e, lw=1.4, ms=4.5, mfc='none', mew=1.0, label='%s (order %.2f)' % (lab, q), **C[f])
    for p_, y0 in ((1, cv[(rel, 1)][0]), (2, cv[(rel, 0)][0])):
        b.loglog(dts, y0 * 0.5 * (dts / dts[0]) ** p_, color='#999', lw=0.8, ls=':')
        b.text(dts[3], y0 * 0.5 * (dts[3] / dts[0]) ** p_ * 0.8, 'slope %d' % p_, fontsize=7.5, color='#666', ha='left', va='top')
    b.set_xlabel(r'$\Delta t$ (s)'); b.set_ylabel('relative error vs Radau'); b.grid(alpha=0.25); b.legend(loc='upper left', frameon=False)
    b.set_title('LSDIRK2, %s' % (r'$\mu=10^{-5}$' if rel == 'mu' else r'$\theta=5\times10^{6}$'), fontsize=9)
fig.suptitle('0D relaxation (case 40): backward-Euler and exponential closures in the implicit stage', fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(OUT, dpi=140); print('wrote', OUT)
