"""Low-wavenumber acoustic viscosity of ONE implicit stage (backward Euler, step tau) of the IMEX pressure system (B4.1),
linear acoustics rho = c = h = 1, no advection, for the complement weight phi.  D_eff = -ln|g| / (k^2 tau), reported in units
of the upwind (Godunov) viscosity c h / 2.  Result: phi = (1 - nu)^+ gives exactly max(1, nu) = max(ch/2, c^2 tau/2)
(the stage reaches the upwind viscosity and adds nothing once backward Euler alone exceeds it); phi = 1 gives 1 + nu.
usage: python3 tools/stage_viscosity.py"""
import numpy as np
def stage(theta, nu, phiP, phiU, Z=1.0, K=1.0):
    tau = nu; e = np.exp(1j * theta); A = (1 + e) / 2; B = e - 1; Dm = 1 - 1 / e
    den = 1 - K * tau * (tau + phiP / (2 * Z)) * B * Dm
    a11 = 1 / den; a12 = -K * tau * A * Dm / den
    a21 = -tau * Dm * A * a11; a22 = 1 - tau * Dm * (A * a12 - phiU * Z / 2 * B)
    return np.array([[a11, a12], [a21, a22]])
if __name__ == '__main__':
    th = 0.01
    print(' nu   phi=(1-nu)+  phi=1   phi=0    |g(pi)| with (1-nu)+')
    for nu in (0.1, 0.2, 0.5, 0.8, 1.0, 2.0, 5.0):
        r = []
        for phi in (max(0.0, 1 - nu), 1.0, 0.0):
            g = np.linalg.eigvals(stage(th, nu, phi, phi)); r.append((-np.log(np.abs(g)) / (th ** 2 * nu)).mean() / 0.5)
        gp = np.abs(np.linalg.eigvals(stage(np.pi, nu, max(0, 1 - nu), max(0, 1 - nu))))
        print('%4.1f   %8.3f   %6.3f  %6.3f    %s' % (nu, *r, np.round(np.sort(gp), 3)))
