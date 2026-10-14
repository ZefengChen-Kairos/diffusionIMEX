"""Exact HEM (p-T-g relaxed, instantaneous heat and mass transfer) solution of the Bartak blowdown (case 32, relax=2) at the probe,
before waves reflected from the closed end come back.  Initial state: NASG water (Pelanti 2022 Table 5), 12.5 MPa, 563.15 K, at rest.
The opening sends a right-facing rarefaction into the tube: first an isentropic liquid wave down to psat (speed c_l ~ 1590 m/s), then
an equilibrium two-phase fan (HEM sound speed c_e(p) ~ 30-230 m/s) that ends at the sonic (choked) state at the open end
(x = 0, u + c = 0 with u < 0).  In the fan x/t = u + c_e, |u|(p) = |u_l| + int_p^psat v / c_e dp.  The trace vapour (alpha 1e-3) is
neglected for the liquid wave.  The liquid wave reflects at the closed end (x = 1.7 m, t ~ 1.07 ms) into liquid that is already at psat:
the reflected rarefaction flashes at once and travels at c_e <= 30 m/s, so it stays near the wall for the whole run (it would need
~60 ms to reach the probe).  The simple wave is therefore exact at x = 48 mm for 0 <= t <= 17.5 ms.
usage: python3 tests/exact_hem_bartak.py [x=0.048]"""
import math, sys
import numpy as np
from scipy.optimize import brentq
L = (1.387, 8899e5, -1244191.0, 0.0, 3202.0, 4.78e-4); V = (1.954, 0.0, 2287484.0, 6417.0, 462.0, 0.0)
def v(P, p, T): g, pi, q, qp, cv, b = P; return (g - 1) * cv * T / (p + pi) + b
def s(P, p, T): g, pi, q, qp, cv, b = P; return cv * math.log(T ** g / (p + pi) ** (g - 1)) + qp
def gfun(P, p, T): g, pi, q, qp, cv, b = P; return g * cv * T + b * p + q - T * s(P, p, T)
def Tsat(p): return brentq(lambda T: gfun(L, p, T) - gfun(V, p, T), 300, 700)
p0, T0, x = 12.5e6, 563.15, float(sys.argv[1]) if len(sys.argv) > 1 else 0.048
s0 = s(L, p0, T0)
# liquid isentrope (NASG: (p + pi) v^g... use T from s = const) down to the saturation point
def T_liq(p): return brentq(lambda T: s(L, p, T) - s0, 300, 700)
ps = brentq(lambda p: T_liq(p) - Tsat(p), 1e6, p0)
def c_liq(p): g, pi, q, qp, cv, b = L; vv = v(L, p, T_liq(p)); return math.sqrt(g * (p + pi) * vv * vv / (vv - b))
ul = 0.0; pp = np.linspace(p0, ps, 2001)
for a, b in zip(pp[:-1], pp[1:]): m = 0.5 * (a + b); ul += (a - b) * v(L, m, T_liq(m)) / c_liq(m)
def mix(p):
    T = Tsat(p); Y = (s0 - s(L, p, T)) / (s(V, p, T) - s(L, p, T)); return (1 - Y) * v(L, p, T) + Y * v(V, p, T), Y
def ce(p): dp = p * 1e-6; v1, _ = mix(p); v2, _ = mix(p - dp); return v1 * math.sqrt(dp / (v2 - v1))
# tabulate the fan: p, |u|, c_e, x/t = c_e - |u|
tab = []; U = ul; p = ps * (1 - 1e-9); dp = ps / 20000
while p > 1.2e5:
    c = ce(p); tab.append((p, U, c, c - U, mix(p)[1]))
    if c - U < 0: break
    m = p - 0.5 * dp; U += dp * mix(m)[0] / ce(m); p -= dp
tab = np.array(tab)
print('psat(T0) on the liquid isentrope %.4f MPa (T %.2f K), |u| after the liquid wave %.2f m/s, liquid speed %.0f m/s' % (ps / 1e6, T_liq(ps), ul, c_liq(ps)))
k = np.argmax(tab[:, 3] < 0); print('sonic (choked) state at the open end: p %.4f MPa, Y %.4f, |u| %.1f m/s' % (tab[k, 0] / 1e6, tab[k, 4], tab[k, 1]))
print('fan head speed %.2f m/s' % tab[0, 3])
print('# t (ms)  p (MPa)  Y  at x = %.3f m' % x)
for t in [0.5e-3, 1e-3, 1.5e-3, 2e-3, 2.5e-3, 3e-3, 5e-3, 10e-3, 17.5e-3]:
    xi = x / t
    if xi >= tab[0, 3]: print('%6.2f %.4f %.5f   (liquid plateau, before the two-phase fan)' % (t * 1e3, ps / 1e6, 0.0)); continue
    j = np.argmax(tab[:, 3] <= xi); print('%6.2f %.4f %.5f%s' % (t * 1e3, tab[j, 0] / 1e6, tab[j, 4], ''))
