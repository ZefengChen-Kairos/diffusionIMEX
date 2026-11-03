#!/usr/bin/env python3
"""Is the dual-parameter explicit split system hyperbolic for eps_M != eps_B?
   Magnetosonic quartic of the split system (paper 2 / code SplitMagnetosonicCoeffs):
     P(l) = (nuM - eM^2 a^2)(nuB - eB^2 can^2) - eB^2 bp^2 nuM,  nuM = (l-u)(l-eM u), nuB = (l-u)(l-eB u)
   For eM = eB the roots are real (closed form).  Sample random states and report complex roots."""
import numpy as np
rng = np.random.default_rng(1)
def roots(u, a2, can2, bp2, eM, eB):
    s1 = (1+eM)*u; a1 = (1+eB)*u; s0 = eM*u*u - eM*eM*a2; a0 = eB*u*u - eB*eB*can2; q0 = eM*u*u; e = eB*eB*bp2
    c3 = -(s1+a1); c2 = s1*a1 + s0 + a0 - e; c1 = -(s0*a1 + s1*a0) + e*s1; c0 = s0*a0 - e*q0
    return np.roots([1.0, c3, c2, c1, c0])
N = 200000; bad = 0; worst = (0, None)
for _ in range(N):
    u = rng.uniform(-3, 3); a2 = rng.uniform(0.01, 9); can2 = rng.uniform(0, 9); bp2 = rng.uniform(0, 9)
    eM = rng.uniform(0, 1); eB = rng.uniform(0, 1)
    r = roots(u, a2, can2, bp2, eM, eB); scale = abs(u) + np.sqrt(a2 + can2 + bp2)
    im = np.max(np.abs(r.imag))/scale
    if im > 1e-8:
        bad += 1
        if im > worst[0]: worst = (im, (u, a2, can2, bp2, eM, eB, r))
print('random states: %d of %d have complex magnetosonic roots (rel. imag > 1e-8)' % (bad, N))
if worst[1] is not None:
    u, a2, can2, bp2, eM, eB, r = worst[1]
    print('worst: u=%.3f a=%.3f c_an=%.3f b_perp=%.3f eps_M=%.3f eps_B=%.3f  roots=%s  rel.imag=%.3g' % (u, np.sqrt(a2), np.sqrt(can2), np.sqrt(bp2), eM, eB, np.round(r, 4), worst[0]))
# a structured scan: u = 1, a = 1, c_an = 1, b_perp = 1, eps grid
print('\nscan u=1, a=c_an=b_perp=1: max rel. imag part over the (eps_M, eps_B) grid')
g = np.linspace(0, 1, 101); M = np.zeros((101, 101))
for i, eM in enumerate(g):
    for j, eB in enumerate(g):
        r = roots(1.0, 1.0, 1.0, 1.0, eM, eB); M[i, j] = np.max(np.abs(r.imag))/3.0
print('max = %.3g at eps_M=%.2f eps_B=%.2f' % (M.max(), g[np.unravel_index(M.argmax(), M.shape)[0]], g[np.unravel_index(M.argmax(), M.shape)[1]]))
# the blast-wave-like state: rho=1, p=100 (a=11.8), Bx=19.7 (c_an=19.7), By small, u~5, eps_M=1, eps_B=0.2 and vice versa
for (u, a2, can2, bp2, eM, eB, lab) in [(5.0, 140.0, 390.0, 1.0, 1.0, 0.2, 'blast ring eM=1,eB=0.2'), (5.0, 140.0, 390.0, 1.0, 0.1, 1.0, 'eM=0.1,eB=1'), (5.0, 140.0, 1.0, 390.0, 1.0, 0.2, 'transverse field eM=1,eB=0.2')]:
    r = roots(u, a2, can2, bp2, eM, eB); print('%-32s roots = %s' % (lab, np.round(r, 3)))
