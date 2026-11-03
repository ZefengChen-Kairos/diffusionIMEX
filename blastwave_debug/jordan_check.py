#!/usr/bin/env python3
"""Exact (rational) Jordan-structure check of the three sub-system Jacobians of the three-way split
   (paper 1 eq. 13-15 / paper 2 eq. 12), 1D x-direction, conservative variables U=(rho,mx,my,mz,E,Bx,By,Bz).
   For each sub-system: eigenvalues, algebraic and geometric multiplicities (8x8 with the Bx column kept,
   and 7x7 with the Bx column removed, i.e. using d_x Bx = 0 as in paper 2, Thm 2.1)."""
import sympy as sp
rho, mx, my, mz, E, Bx, By, Bz = sp.symbols('rho mx my mz E Bx By Bz')
g = sp.Rational(5, 3)
U = sp.Matrix([rho, mx, my, mz, E, Bx, By, Bz])
u, v, w = mx/rho, my/rho, mz/rho
k = (mx**2 + my**2 + mz**2)/(2*rho); B2 = Bx**2 + By**2 + Bz**2; p = (g-1)*(E - k - B2/2)
h = (p/(g-1) + p)/rho   # specific enthalpy e + p/rho
Fc = sp.Matrix([rho*u, rho*u*u, rho*u*v, rho*u*w, k*u, 0, 0, 0])
Fp = sp.Matrix([0, p, 0, 0, rho*h*u, 0, 0, 0])
FB = sp.Matrix([0, B2/2 - Bx**2, -Bx*By, -Bx*Bz, B2*u - (u*Bx + v*By + w*Bz)*Bx, 0, u*By - v*Bx, u*Bz - w*Bx])
state = {rho: sp.Rational(7, 5), mx: sp.Rational(3, 7), my: sp.Rational(-2, 9), mz: sp.Rational(1, 4), E: sp.Rational(23, 3), Bx: sp.Rational(5, 6), By: sp.Rational(-4, 5), Bz: sp.Rational(2, 7)}
def analyse(name, F, drop_Bx):
    A = F.jacobian(U).subs(state)
    if drop_Bx:
        keep = [0, 1, 2, 3, 4, 6, 7]; A = A.extract(keep, keep)
    n = A.shape[0]; lam = sp.symbols('lam')
    cp = sp.factor((A - lam*sp.eye(n)).det())
    print('%s (%dx%d%s): char poly = %s' % (name, n, n, ', Bx column removed' if drop_Bx else '', cp))
    for ev, mult in sp.roots(sp.Poly((A - lam*sp.eye(n)).det(), lam)).items():
        if ev.is_real:
            geo = n - (A - ev*sp.eye(n)).rank()
            print('   eigenvalue %-10s alg. mult = %d  geom. mult = %d  %s' % (sp.nsimplify(ev) if ev.is_rational else sp.N(ev, 6), mult, geo, 'DEFECTIVE' if geo < mult else ''))
print('state: u = %s, v = %s, w = %s' % (sp.N(state[mx]/state[rho], 5), sp.N(state[my]/state[rho], 5), sp.N(state[mz]/state[rho], 5)))
for drop in [False, True]:
    analyse('advective F^c', Fc, drop); analyse('pressure  F^p', Fp, drop); analyse('magnetic  F^B', FB, drop); print()
