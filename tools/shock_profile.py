"""Shock-profile ODE for the 5-eq (Kapila) model with a local energy-defect split lambda_k.
Per unit mass of phase k (SG EOS): de_k = -p dv_k + dq_k, Y_k dq_k = lam_k dq, dq = (p - P(v)) dv,
P(v) = p0 + j^2 (v0 - v).  Integrate in v from v0 until p = P(v) again -> post-shock state.
Compare the end partition with the (II.5) closure (each phase on its own Hugoniot at p1)."""
import numpy as np, sys
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

class Mix:
    def __init__(self, ph, p0=1e5):
        self.g = np.array([a[0] for a in ph]); self.pi = np.array([a[1] for a in ph])
        rho = np.array([a[2] for a in ph]); al = np.array([a[3] for a in ph])
        self.m0 = al * rho; self.rho0 = self.m0.sum(); self.Y = self.m0 / self.rho0
        self.v0 = 1 / self.rho0; self.vk0 = 1 / rho; self.al0 = al; self.p0 = p0
    def K(self, p): return self.g * (p + self.pi)
    def G(self): return (self.g + 1) / 2
    def e(self, p, vk): return (p + self.g * self.pi) * vk / (self.g - 1)
    def T(self, p, vk):   # cv from reference: T = (p+pi) v /((g-1) cv); use cv=1 -> relative T only
        return (p + self.pi) * vk / (self.g - 1)

def lam_rule(rule, mx, p, vk):
    Y = mx.Y; al = Y * vk / (Y * vk).sum()
    if rule == 'mass': w = Y
    elif rule == 'spb': w = al / (mx.g - 1)
    elif rule == 'weak': w = al * mx.G() / mx.K(p) ** 2
    elif rule == 'aK': w = al / mx.K(p)
    elif rule == 'weak0': w = mx.al0 * mx.G() / mx.K(mx.p0) ** 2     # frozen at pre-shock state (non-local)
    else: raise ValueError(rule)
    return w / w.sum()

def profile(mx, j2, rule):
    Y = mx.Y; g = mx.g; pi = mx.pi; p0 = mx.p0; v0 = mx.v0
    def rhs(v, y):
        p, v1 = y; v2 = (v - Y[0] * v1) / Y[1]; vk = np.array([v1, v2])
        P = p0 + j2 * (v0 - v); lam = lam_rule(rule, mx, p, vk)
        # vk p' + g(p+pi) vk' = (g-1)(lam/Y)(p-P)   (k=1,2),   Y1 v1' + Y2 v2' = 1
        A = np.array([[vk[0], g[0] * (p + pi[0]), 0.0], [vk[1], 0.0, g[1] * (p + pi[1])], [0.0, Y[0], Y[1]]])
        b = np.array([(g[0] - 1) * lam[0] / Y[0] * (p - P), (g[1] - 1) * lam[1] / Y[1] * (p - P), 1.0])
        s = np.linalg.solve(A, b); return [s[0], s[1]]
    def ev(v, y): return y[0] - (p0 + j2 * (v0 - v))
    ev.terminal = True; ev.direction = 1
    # start a hair inside the profile along the isentrope (p - P = 0 at v0 exactly)
    dv = -1e-9 * v0; y0 = np.array([p0, mx.vk0[0]]); y0 = y0 + dv * np.array(rhs(v0, y0))
    sol = solve_ivp(rhs, [v0 + dv, 0.2 * v0], y0, events=ev, rtol=1e-10, atol=1e-14, max_step=1e-3 * v0)
    v1 = sol.t_events[0][0]; p1, v11 = sol.y_events[0][0]; v21 = (v1 - Y[0] * v11) / Y[1]
    return p1, v1, np.array([v11, v21])

def hug_vk(mx, p1):      # phase Hugoniot centred at (p0, vk0): e - e0 = -(p0+p1)/2 (v - v0)
    g, pi, p0 = mx.g, mx.pi, mx.p0
    return mx.vk0 * ((g - 1) * (p1 + pi) + (g + 1) * (p0 + pi)) / ((g + 1) * (p1 + pi) + (g - 1) * (p0 + pi))

def ii5_from_p1(mx, p1):
    vk = hug_vk(mx, p1); v1 = (mx.Y * vk).sum(); j2 = (p1 - mx.p0) / (mx.v0 - v1); return j2, v1, vk

def end_state(mx, p1, vk):
    al = mx.Y * vk / (mx.Y * vk).sum(); T = mx.T(p1, vk) / mx.T(mx.p0, mx.vk0)
    return al, T

if __name__ == '__main__':
    cases = {'epoxy-spinel': [(2.43, 5.3e9, 1185, 0.5954), (1.62, 141e9, 3622, 0.4046)],
             'water+air a=0.5': [(4.4, 6e8, 1000, 0.5), (1.4, 0, 1, 0.5)],
             'water+air a=0.1': [(4.4, 6e8, 1000, 0.9), (1.4, 0, 1, 0.1)],
             'water+1e-3 air': [(4.4, 6e8, 1000, 1 - 1e-3), (1.4, 0, 1, 1e-3)]}
    for nm, ph in cases.items():
        mx = Mix(ph)
        for p1 in [1e6, 1e8, 1e9, 1e10, 2e11]:
            j2, v1, vk = ii5_from_p1(mx, p1); al, T = end_state(mx, p1, vk); us = np.sqrt(j2) * mx.v0
            print('%-16s p1 %.0e  (II.5): us %9.1f  rho1/rho0 %.4f  a1 %.5f  T1/T0 %9.2f T2/T0 %9.2f' % (nm, p1, us, mx.v0 / v1, al[0], T[0], T[1]))
            for rule in ['weak', 'weak0', 'mass', 'spb', 'aK']:
                try:
                    pp, vv, vkk = profile(mx, j2, rule); a, TT = end_state(mx, pp, vkk)
                    print('   %-6s same us: p1 %.4e (%+.2e rel) rho1/rho0 %.4f a1 %.5f (%+.1e) T1/T0 %9.2f T2/T0 %9.2f' % (rule, pp, pp / p1 - 1, mx.v0 / vv, a[0], a[0] - al[0], TT[0], TT[1]))
                except Exception as ex: print('   %-6s failed: %s' % (rule, ex))

# ---- (II.5)-exact constant share c along the profile, vs closed forms ----
def profile_c(mx, j2, c):
    Y = mx.Y; g = mx.g; pi = mx.pi; p0 = mx.p0; v0 = mx.v0; lam = np.array([c, 1 - c])
    def rhs(v, y):
        p, v1 = y; v2 = (v - Y[0] * v1) / Y[1]; vk = np.array([v1, v2]); P = p0 + j2 * (v0 - v)
        A = np.array([[vk[0], g[0] * (p + pi[0]), 0.0], [vk[1], 0.0, g[1] * (p + pi[1])], [0.0, Y[0], Y[1]]])
        b = np.array([(g[0] - 1) * lam[0] / Y[0] * (p - P), (g[1] - 1) * lam[1] / Y[1] * (p - P), 1.0])
        s = np.linalg.solve(A, b); return [s[0], s[1]]
    def ev(v, y): return y[0] - (p0 + j2 * (v0 - v))
    ev.terminal = True; ev.direction = 1
    dv = -1e-9 * v0; y0 = np.array([p0, mx.vk0[0]]); y0 = y0 + dv * np.array(rhs(v0, y0))
    sol = solve_ivp(rhs, [v0 + dv, 0.2 * v0], y0, events=ev, rtol=1e-10, atol=1e-14, max_step=1e-3 * v0)
    v1 = sol.t_events[0][0]; p1, v11 = sol.y_events[0][0]; return p1, v11

def gap_share(mx, p1):   # alpha_k rho_k (e_H - e_S)(p0 -> p1) per mixture volume
    g, pi = mx.g, mx.pi; vH = hug_vk(mx, p1); vS = mx.vk0 * ((mx.p0 + pi) / (p1 + pi)) ** (1 / g)
    q = mx.m0 * (p1 + g * pi) * (vH - vS) / (g - 1); return q / q.sum()

def Tds_share(mx, p1, pref):   # T_k(pref on isentrope) * Delta s_k^H   (cv = 1 units)
    g, pi = mx.g, mx.pi; vH = hug_vk(mx, p1)
    ds = np.log((p1 + pi) * vH ** g / ((mx.p0 + pi) * mx.vk0 ** g))          # cv * ln(...)
    vS = mx.vk0 * ((mx.p0 + pi) / (pref + pi)) ** (1 / g); T = (pref + pi) * vS / (g - 1)
    q = mx.m0 * T * ds; return q / q.sum()

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'c':
    cases = {'epoxy-spinel': [(2.43, 5.3e9, 1185, 0.5954), (1.62, 141e9, 3622, 0.4046)],
             'water+air a=0.5': [(4.4, 6e8, 1000, 0.5), (1.4, 0, 1, 0.5)],
             'water+air a=0.1': [(4.4, 6e8, 1000, 0.9), (1.4, 0, 1, 0.1)],
             'water+1e-3 air': [(4.4, 6e8, 1000, 1 - 1e-3), (1.4, 0, 1, 1e-3)],
             'water+1e-6 air': [(4.4, 6e8, 1000, 1 - 1e-6), (1.4, 0, 1, 1e-6)]}
    print('phase-1 share:  c_exact | weak(p0)  gapH(p1)  TdS(pbar)  TdS(p0)  TdS(p1)')
    for nm, ph in cases.items():
        mx = Mix(ph)
        for p1 in [3e5, 1e6, 1e7, 1e8, 1e9, 1e10, 2e11]:
            j2, v1, vk = ii5_from_p1(mx, p1)
            f = lambda c: profile_c(mx, j2, c)[1] - vk[0]
            try: c = brentq(f, 1e-12, 1 - 1e-12, xtol=1e-12)
            except Exception as ex: c = float('nan')
            w = lam_rule('weak0', mx, mx.p0, mx.vk0)[0]; gh = gap_share(mx, p1)[0]
            print('%-16s p1 %.0e  %.5f | %.5f  %.5f  %.5f  %.5f  %.5f' % (nm, p1, c, w, gh, Tds_share(mx, p1, 0.5 * (mx.p0 + p1))[0], Tds_share(mx, p1, mx.p0)[0], Tds_share(mx, p1, p1)[0]))
