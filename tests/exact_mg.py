"""Exact Riemann solution of the Euler equations with two Mie-Grueneisen-type EOS (NASG, Cochran-Chan, JWL),
for the MG shock tubes (V12): SPB 2009 4.2 (nitromethane, Cochran-Chan) and Pelanti 2022 8.1.1 (JWL gas / NASG water).
Rarefactions: integrate de/drho = p/rho^2 (isentrope) and du = -/+ c/rho drho; shocks: Hugoniot
e - e0 + (p + p0)(1/rho - 1/rho0)/2 = 0.  usage: python3 tests/exact_mg.py CASE x0 x1 N [out.dat]  (columns x rho u p)"""
import sys, numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


class MG:
    """p = G rho (e - er(rho)) + pr(rho), G constant, er' = pr/rho^2 (A2.1)."""
    def __init__(self, kind, G, rr, c1, c2, e1, e2, e0=0.0):
        self.kind, self.G, self.rr, self.k1, self.k2, self.e1, self.e2, self.e0 = kind, G, rr, c1, c2, e1, e2, e0

    def ref(self, r):
        if self.kind == 'CC':
            x = r / self.rr; t1, t2 = self.k1 * x ** self.e1, self.k2 * x ** self.e2
            return t1 - t2, (self.e1 * t1 - self.e2 * t2) / r, (t1 / (self.e1 - 1) - t2 / (self.e2 - 1)) / r
        t1, t2 = self.k1 * np.exp(-self.e1 * self.rr / r), self.k2 * np.exp(-self.e2 * self.rr / r)
        return t1 + t2, (self.e1 * t1 + self.e2 * t2) * self.rr / r ** 2, t1 / (self.e1 * self.rr) + t2 / (self.e2 * self.rr) - self.e0

    def p(self, r, e): pr, _, er = self.ref(r); return self.G * r * (e - er) + pr
    def e(self, r, p): pr, _, er = self.ref(r); return er + (p - pr) / (self.G * r)
    def c2(self, r, p): pr, dpr, _ = self.ref(r); return (self.G + 1) * (p - pr) / r + dpr


class NASG:
    def __init__(self, g, pi, q=0.0, b=0.0): self.g, self.pi, self.q, self.b = g, pi, q, b
    def p(self, r, e): return (self.g - 1) * r * (e - self.q) / (1 - r * self.b) - self.g * self.pi
    def e(self, r, p): return (p + self.g * self.pi) * (1 - r * self.b) / ((self.g - 1) * r) + self.q
    def c2(self, r, p): return self.g * (p + self.pi) / (r * (1 - r * self.b))


class Side:
    """one side of the Riemann problem; sgn = -1 left (waves u - c), +1 right"""
    def __init__(self, eos, r, u, p, sgn):
        self.E, self.r, self.u, self.p, self.s = eos, r, u, p, sgn
        self.e = eos.e(r, p); self.c = np.sqrt(eos.c2(r, p))
        self._rare = None

    def rare(self, rmin):
        # isentrope from (r, e): y = [e, u]; d e/d r = p/r^2; d u/d r = sgn c / r  (left: u decreases as r decreases ... with sgn)
        E, s = self.E, self.s
        def f(r, y):
            p = E.p(r, y[0]); c = np.sqrt(max(E.c2(r, p), 0.0))
            return [p / r ** 2, s * c / r]
        sol = solve_ivp(f, [self.r, rmin], [self.e, self.u], method='DOP853', rtol=1e-12, atol=[1e-6 * abs(self.e) + 1e-6, 1e-9], dense_output=True)
        self._rare = sol; return sol

    def star(self, ps):
        """velocity on the contact side for pressure ps, and the star density"""
        E = self.E
        if ps >= self.p:   # shock
            f = lambda r: E.e(r, ps) - self.e + 0.5 * (ps + self.p) * (1 / r - 1 / self.r)
            rs = brentq(f, self.r * (1 + 1e-14), self.r * 50, xtol=1e-14 * self.r, maxiter=500)
            du = np.sqrt(max((ps - self.p) * (1 / self.r - 1 / rs), 0.0))
            return self.u + self.s * du, rs
        sol = self._rare
        g = lambda r: E.p(r, sol.sol(r)[0]) - ps
        rs = brentq(g, sol.t[-1], self.r, xtol=1e-14 * self.r, maxiter=500)
        return sol.sol(rs)[1], rs


def solve(L, R):
    L.rare(L.r * 0.05); R.rare(R.r * 0.05)
    pmin = max(L.E.p(L._rare.t[-1], L._rare.y[0][-1]), R.E.p(R._rare.t[-1], R._rare.y[0][-1])) + 1.0
    F = lambda p: L.star(p)[0] - R.star(p)[0]
    lo, hi = pmin, max(L.p, R.p) * 2
    while F(hi) > 0: hi *= 2
    ps = brentq(F, lo, hi, xtol=1e-12 * hi, rtol=1e-15, maxiter=500)
    us, rL = L.star(ps); _, rR = R.star(ps)
    return ps, us, rL, rR


def sample(L, R, ps, us, rLs, rRs, xi):
    out = np.zeros((len(xi), 3))
    for i, x in enumerate(xi):
        if x < us: S, rs = L, rLs
        else: S, rs = R, rRs
        s = S.s
        if ps >= S.p:   # shock
            m = np.sqrt((ps - S.p) / (1 / S.r - 1 / rs)); W = S.u + s * m / S.r
            inside = (x < W) if s < 0 else (x > W)
            out[i] = [S.r, S.u, S.p] if inside else [rs, us, ps]
        else:
            cs = np.sqrt(S.E.c2(rs, ps)); head = S.u + s * S.c; tail = us + s * cs
            if (s < 0 and x <= head) or (s > 0 and x >= head): out[i] = [S.r, S.u, S.p]
            elif (s < 0 and x >= tail) or (s > 0 and x <= tail): out[i] = [rs, us, ps]
            else:
                sol = S._rare
                def h(r):
                    e, u = sol.sol(r); p = S.E.p(r, e); return u + s * np.sqrt(S.E.c2(r, p)) - x
                r = brentq(h, rs, S.r, xtol=1e-13 * S.r)
                e, u = sol.sol(r); out[i] = [r, u, S.E.p(r, e)]
    return out


CASES = {   # case: (EOS_L, state L, EOS_R, state R, x_interface, t)
    60: (lambda: MG('CC', 1.19, 1134.0, 0.819181e9, 1.50835e9, 4.52969, 1.42144), (1134.0, 0.0, 20e9),
         lambda: MG('CC', 1.19, 1134.0, 0.819181e9, 1.50835e9, 4.52969, 1.42144), (1200.0, 0.0, 0.2e6), 0.5, 67e-6),
    61: (lambda: MG('JWL', 0.25, 1840.0, 854.5e9, 20.5e9, 4.6, 1.35), (2000.0, 0.0, 4.6406e10),
         lambda: NASG(1.187, 7028e5, -1177788.0, 6.61e-4), (1044.0, 0.0, 1e5), 0.0, 0.9e-3),
}

if __name__ == '__main__':
    case, x0, x1, N = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
    eL, sL, eR, sR, xd, t = CASES[case]
    L = Side(eL(), *sL, -1); R = Side(eR(), *sR, +1)
    ps, us, rLs, rRs = solve(L, R)
    print('p* = %.8e Pa, u* = %.8e m/s, rho*L = %.6f, rho*R = %.6f' % (ps, us, rLs, rRs))
    x = x0 + (np.arange(N) + 0.5) * (x1 - x0) / N
    prof = sample(L, R, ps, us, rLs, rRs, (x - xd) / t)
    if len(sys.argv) > 5: np.savetxt(sys.argv[5], np.column_stack([x, prof]), fmt='%.10e', header='x rho u p')
