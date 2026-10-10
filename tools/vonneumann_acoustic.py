"""Von Neumann analysis of the implicit acoustic step of the TV-split IMEX (B4.1), linear acoustics, u0 = 0
(so the explicit Rusanov part vanishes: a_Omega = 0).  One backward-Euler stage tau, r = tau/dx, nu = c tau/dx.
  u_f   = (u_i+u_{i+1})/2 - (r/rho) dP - beta dP/(2Z)          (P implicit, u explicit = u^n)
  p_f   = (P_i+P_{i+1})/2 - beta (Z/2) du                        (rho uniform -> density weights 1/2)
  P^{n+1} = P^n - rho c^2 r delta[u_f],  u^{n+1} = u^n - (r/rho) delta[p_f]
beta = 0 : old scheme;  beta = (1-nu)^+ : Godunov top-up.  Non-dimensional p~ = P/Z.
Compare with first-order acoustic Godunov |g|^2 = 1 - 4 nu (1-nu) sin^2(theta/2) (explicit, nu <= 1)."""
import numpy as np
def G(nu, th, beta):
    S = np.sin(th / 2) ** 2; s = np.sin(th); D = 1 + 4 * S * (nu * nu + beta * nu / 2)
    return np.array([[1 / D, -1j * nu * s / D], [-1j * nu * s / D, (1 - 2 * beta * nu * S) - nu * nu * s * s / D]])
th = np.linspace(1e-3, np.pi, 2000)
print(' nu   | max|lam| old (theta)      | top-up max|lam|  | Godunov max|g| | damping per step at theta=pi: old / top-up / Godunov')
for nu in [0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 3.0, 10.0]:
    lo = [max(abs(np.linalg.eigvals(G(nu, t, 0.0)))) for t in th]
    beta = max(0.0, 1 - nu)
    lt = [max(abs(np.linalg.eigvals(G(nu, t, beta)))) for t in th]
    gd = np.sqrt(np.maximum(0, 1 - 4 * nu * (1 - nu) * np.sin(th / 2) ** 2)) if nu <= 1 else np.full_like(th, np.nan)
    # compare damping mode by mode: is top-up at least as damping as Godunov?  ratio of (1-|lam|) to (1-|g|) at mid/high frequencies
    k = th > 0.5
    ok = np.all(np.array(lt)[k] <= gd[k] + 1e-12) if nu <= 1 else True
    print(f"{nu:5.2f} | {max(lo):.6f} ({th[int(np.argmax(lo))]:.2f})          | {max(lt):.6f}        | {np.nanmax(gd) if nu<=1 else float('nan'):.6f}       | "
          f"{max(abs(np.linalg.eigvals(G(nu, np.pi, 0.0)))):.4f} / {max(abs(np.linalg.eigvals(G(nu, np.pi, beta)))):.4f} / {np.sqrt(max(0,1-4*nu*(1-nu))) if nu<=1 else float('nan'):.4f}   top-up<=Godunov(theta>0.5): {ok}")
