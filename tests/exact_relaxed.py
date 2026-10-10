"""Exact solutions of relaxed limit models for the symmetric expansion tube (u = -U0 | +U0 at x = 0.5, uniform p0, T0, alpha_v):
  pT  : pT-relaxed (4-equation) model, frozen Y: one pT-isentropic fan down to u = 0 (Pelanti 2022 Fig. 4 "exact 4-eq.");
  pTg : activation-limited model (P-S 2014 p-pT-pTg, Pelanti 2022 pTg): pT isentrope while T < Tsat(p), then, after a constant
        state, the pTg (HEM) isentrope with T = Tsat(p) down to u = 0 (the characteristic speed drops from c_pT to c_pTg).
Phases: NASG (g, pi, q, qp, cv, b) as in src/spb6_physics.h (b = 0: stiffened gas); phase 0 liquid, phase 1 vapour.
Generalises tests/exactfan.py (P-S SG water) to NASG and to the 4-equation model; derivatives are analytic.
usage: F = Fan(PH, p0, T0, av); x, p, u, Y, a, T = F.full(U0, t, chem)"""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp


class Fan:
    def __init__(s, PH, p0, T0, av):
        s.PH = [tuple(map(float, P)) for P in PH]; s.p0, s.T0, s.av = p0, T0, av
        m = [(1 - av) / s.v(0, p0, T0), av / s.v(1, p0, T0)]; s.rho0 = sum(m); s.Y0 = m[1] / s.rho0

    # phase EOS (as src/spb6_physics.h): v(p,T), e(p,T), s(p,T), g = h - T s, and partial derivatives
    def v(s, k, p, T): g, pi, q, qp, cv, b = s.PH[k]; return (g - 1) * cv * T / (p + pi) + b
    def e(s, k, p, T): g, pi, q, qp, cv, b = s.PH[k]; return cv * T * (p + g * pi) / (p + pi) + q
    def s_(s, k, p, T): g, pi, q, qp, cv, b = s.PH[k]; return cv * (g * np.log(T) - (g - 1) * np.log(p + pi)) + qp
    def g(s, k, p, T): G, pi, q, qp, cv, b = s.PH[k]; return G * cv * T + b * p + q - T * s.s_(k, p, T)
    def dv(s, k, p, T): g, pi, q, qp, cv, b = s.PH[k]; return -(g - 1) * cv * T / (p + pi) ** 2, (g - 1) * cv / (p + pi)          # (v_p, v_T)
    def de(s, k, p, T): g, pi, q, qp, cv, b = s.PH[k]; return cv * T * (1 - g) * pi / (p + pi) ** 2, cv * (p + g * pi) / (p + pi)  # (e_p, e_T)
    def Ts(s, p): return brentq(lambda T: s.g(0, p, T) - s.g(1, p, T), 150.0, 1500.0, xtol=1e-13, rtol=1e-15)
    def dTs(s, p):   # Clausius-Clapeyron: dTs/dp = (v_v - v_l) / (s_v - s_l)
        T = s.Ts(p); return T, (s.v(1, p, T) - s.v(0, p, T)) / (s.s_(1, p, T) - s.s_(0, p, T))

    # stage 1: pT isentrope at frozen Y; state (T, u) as a function of p
    def _pT(s, p, T, Y):
        w = (1 - Y, Y); vp = vT = ep = eT = 0.0
        for k in (0, 1):
            a, b = s.dv(k, p, T); c, d = s.de(k, p, T); vp += w[k] * a; vT += w[k] * b; ep += w[k] * c; eT += w[k] * d
        dT = -(ep + p * vp) / (eT + p * vT); v = w[0] * s.v(0, p, T) + w[1] * s.v(1, p, T)
        c = np.sqrt(-v * v / (vp + vT * dT)); return dT, v, c

    # stage 2: pTg isentrope (T = Ts(p)); state (Y, u) as a function of p
    def _pTg(s, p, Y):
        T, Tp = s.dTs(p); w = (1 - Y, Y); vp = ep = 0.0
        for k in (0, 1):
            a, b = s.dv(k, p, T); c, d = s.de(k, p, T); vp += w[k] * (a + b * Tp); ep += w[k] * (c + d * Tp)
        vy = s.v(1, p, T) - s.v(0, p, T); ey = s.e(1, p, T) - s.e(0, p, T)
        dY = -(ep + p * vp) / (ey + p * vy); v = w[0] * s.v(0, p, T) + w[1] * s.v(1, p, T)
        c = np.sqrt(-v * v / (vp + vy * dY)); return dY, v, c, T

    def left(s, U0, t, chem=True, n=800, pmin=1.0):
        """left half (x <= 0.5): arrays x, p, u, Y, alpha_v, T sorted in x; also returns a dict of key states"""
        def r1(p, y): dT, v, c = s._pT(p, y[0], s.Y0); return [dT, -v / c]
        ev_c = lambda p, y: y[1]; ev_c.terminal = True                              # centre reached (u = 0)
        ev_a = lambda p, y: y[0] - s.Ts(p); ev_a.terminal = True                     # activation T = Tsat(p)
        S1 = solve_ivp(r1, [s.p0, pmin], [s.T0, -U0], events=[ev_c, ev_a] if chem else [ev_c], rtol=1e-11, atol=1e-12, dense_output=True)
        p1 = S1.t[-1]; P = np.geomspace(s.p0, p1, n); Z = S1.sol(P)
        C = np.array([s._pT(P[i], Z[0, i], s.Y0)[2] for i in range(n)]); V = np.array([s._pT(P[i], Z[0, i], s.Y0)[1] for i in range(n)])
        A = np.array([s.Y0 * s.v(1, P[i], Z[0, i]) for i in range(n)]) / V
        xs = [np.array([0.0]), 0.5 + (Z[1] - C) * t]; ps = [[s.p0], P]; us = [[-U0], Z[1]]; Ys = [[s.Y0], np.full(n, s.Y0)]; As = [[s.av], A]; Ts = [[s.T0], Z[0]]
        info = {'p_end1': p1, 'u_end1': Z[1, -1], 'T_end1': Z[0, -1]}
        activated = chem and len(S1.t_events) > 1 and len(S1.t_events[1]) > 0
        if activated:
            def r2(p, y): dY, v, c, T = s._pTg(p, y[0]); return [dY, -v / c]
            ev_y = lambda p, y: y[0] - 1.0; ev_y.terminal = True                     # complete evaporation (not handled further)
            S2 = solve_ivp(r2, [p1, pmin], [s.Y0, Z[1, -1]], events=[ev_c, ev_y], rtol=1e-11, atol=1e-14, dense_output=True)
            p2 = S2.t[-1]; P2 = np.geomspace(p1, p2, n); Z2 = S2.sol(P2); R = [s._pTg(P2[i], Z2[0, i]) for i in range(n)]
            C2 = np.array([r[2] for r in R]); V2 = np.array([r[1] for r in R]); T2 = np.array([r[3] for r in R])
            A2 = Z2[0] * np.array([s.v(1, P2[i], T2[i]) for i in range(n)]) / V2
            xs.append(0.5 + (Z2[1] - C2) * t); ps.append(P2); us.append(Z2[1]); Ys.append(Z2[0]); As.append(A2); Ts.append(T2)
            info.update(p_end2=p2, u_end2=Z2[1, -1], Y_end2=Z2[0, -1], evaporated=len(S2.t_events[1]) > 0)
        # constant state up to the centre if u = 0 was reached; otherwise the fan tail borders a vacuum (pmin reached) or a
        # completely evaporated state (not computed): the profile then ends at the tail
        info['centre'] = abs(us[-1][-1]) < 1e-9 * U0
        if info['centre']:
            xs.append([0.5]); ps.append([ps[-1][-1]]); us.append([us[-1][-1]]); Ys.append([Ys[-1][-1]]); As.append([As[-1][-1]]); Ts.append([Ts[-1][-1]])
        out = [np.concatenate(z) for z in (xs, ps, us, Ys, As, Ts)]
        info['activated'] = activated; return out, info

    def full(s, U0, t, chem=True):
        (x, p, u, Y, a, T), info = s.left(U0, t, chem)
        r = lambda z: z[::-1]; gap = [] if info['centre'] else [np.nan]   # NaN: no line across a vacuum
        return (np.concatenate([x, gap, 1 - r(x)]), np.concatenate([p, gap, r(p)]), np.concatenate([u, gap, -r(u)]),
                np.concatenate([Y, gap, r(Y)]), np.concatenate([a, gap, r(a)]), np.concatenate([T, gap, r(T)])), info


if __name__ == '__main__':   # checks against tests/exactfan.py (P-S SG water, +-2 m/s, t = 3.2 ms) and Pelanti 2022 8.2.1 (NASG)
    PS = [(2.35, 1e9, -1167e3, 0.0, 1816.0, 0.0), (1.43, 0.0, 2030e3, -23.4e3, 1040.0, 0.0)]
    F = Fan(PS, 1e5, 354.728, 1e-2)
    for chem in (False, True):
        (x, p, u, Y, a, T), I = F.full(2.0, 3.2e-3, chem); m = len(x) // 2 - 1
        print('P-S +-2 m/s %-4s centre p %.5e Y %.4e alpha %.4f T %.3f | activation p %.5e' % ('pTg' if chem else 'pT', p[m], Y[m], a[m], T[m], I['p_end1']))
    PEL = [(1.187, 7028e5, -1177788.0, 0.0, 3610.0, 6.61e-4), (1.467, 0.0, 2077616.0, 14317.0, 955.0, 0.0)]
    F = Fan(PEL, 1e5, 353.0, 1e-2); print('Pelanti Tsat(1e5) = %.2f K' % F.Ts(1e5))
    for chem in (False, True):
        (x, p, u, Y, a, T), I = F.full(2.0, 3e-3, chem); m = len(x) // 2 - 1
        print('Pelanti 8.2.1 %-4s centre p %.5e Y %.4e alpha %.4f T %.3f | %s' % ('pTg' if chem else 'pT', p[m], Y[m], a[m], T[m], I))
