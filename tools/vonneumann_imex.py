"""Linear von Neumann analysis of one backward-Euler stage of the TV-split IMEX (B3 + B4.1) for single-phase
acoustics around (rho0, u0 > 0, p0), first order in space.  Non-dimensional: p~ = p'/Z, Z = rho0 c.
  C  = u0 tau/dx (material CFL), M = u0/c, nu = c tau/dx = C/M.
Explicit Rusanov with a = u0 is pure upwind:  X* = g_c X^n,  g_c = 1 - C (1 - e^{-i th}).
Implicit stage (P implicit, u* explicit):  p~^{n+1} = p~* - nu delta[u_f],  u^{n+1} = u* - nu delta[p~_f]
  u_f  = avg(u*) - (nu + beta/2) D p~^{n+1}         (D = forward difference)
  p~_f = avg(p~^{n+1}) - (beta/2) D u*
(increment-flux terms are O(u0 * du) and do not enter at linear order beyond the convection already in g_c)
beta: 0 = old scheme;  (1-nu)^+ = Godunov top-up (B4.1).
Reference: explicit first-order Godunov for the full system (waves u0-c, u0, u0+c), stable iff (1+M) nu <= 1."""
import numpy as np
def G(C, nu, th, beta):
    S = np.sin(th / 2) ** 2; s = np.sin(th); gc = 1 - C * (1 - np.exp(-1j * th))
    D = 1 + 4 * S * (nu * nu + beta * nu / 2)
    # p^{n+1} = (gc p^n - i nu s gc u^n)/D ;  u^{n+1} = (1 - 2 beta nu S) gc u^n - i nu s p^{n+1}
    A = np.array([[gc / D, -1j * nu * s * gc / D],
                  [-1j * nu * s * gc / D, (1 - 2 * beta * nu * S) * gc - nu * nu * s * s * gc / D]])
    return A
def godunov(C, nu, th):
    lam = []
    for cfl in (C - nu, C + nu):           # acoustic waves u0 -/+ c (upwind for each sign)
        if cfl >= 0: lam.append(abs(1 - cfl * (1 - np.exp(-1j * th))))
        else: lam.append(abs(1 - cfl * (np.exp(1j * th) - 1)))
    return max(lam)
th = np.linspace(1e-3, np.pi, 1500)
print('  C     M     nu   | old: max|lam|  |lam(pi)|  d_eff  | top-up: max|lam| |lam(pi)| d_eff | Godunov (if stable): |lam(pi)| d_eff')
for C, M in [(0.2, 0.01), (0.2, 0.09), (0.2, 0.2), (0.2, 0.5), (0.2, 1.3), (0.05, 0.09), (0.05, 0.5), (0.5, 0.09), (0.2, 3.0)]:
    nu = C / M
    out = []
    for beta in (0.0, max(0.0, 1 - nu)):
        lam = np.array([max(abs(np.linalg.eigvals(G(C, nu, t, beta)))) for t in th])
        k = th < 0.2; d = np.polyfit(th[k] ** 2, 1 - lam[k], 1)[0]   # 1-|lam| ~ d th^2 (d in dx^2/tau units)
        out.append((lam.max(), lam[-1], d))
    gs = ''
    if (1 + M) * nu <= 1.0 + 1e-12:
        lg = np.array([godunov(C, nu, t) for t in th]); k = th < 0.2; dg = np.polyfit(th[k] ** 2, 1 - lg[k], 1)[0]
        gs = f'{lg[-1]:.3f}  {dg:.4f}'
    else: gs = 'unstable'
    print(f'{C:5.2f} {M:5.2f} {nu:6.2f} | {out[0][0]:.6f} {out[0][1]:.4f} {out[0][2]:.4f} | {out[1][0]:.6f} {out[1][1]:.4f} {out[1][2]:.4f} | {gs}')
