"""Von Neumann analysis of one implicit stage of the sequential magnetic -> pressure split
for the linear perpendicular fast-magnetosonic wave (B = B0 z, k along x, uniform rho, p0).

  magnetic step : b1 = b - h B0 dx u_m ,  u_m = u - h/rho (dx p^e + B0 dx b1)      (paper Eq. 27)
  momentum      : u** = u - h/rho B0 dx b1
  pressure step : p1 = p - h g p0 dx u1 ,  u1 = u** - h/rho dx p1                    (Eq. 37-38)

Scaled with a = c_A h k, s = c_s h k (spectral, dx -> i k).  'pe' selects the pressure used
in the magnetic step: 'n' (explicit, as in the paper/code) or 'none' (dropped).
A magnetic Laplacian nu k^2 h (nu = kappa_B) can be added to the b equation (Eq. 83).
"""
import numpy as np

def amp(a, s, pe='n', d=0.0):
    # state (U, Bn, P), Bn = b/sqrt(rho), P = p/sqrt(g p0 rho)
    M = np.zeros((3, 3), complex)
    for col in range(3):
        U, Bn, P = np.eye(3)[col].astype(complex)
        cp = 1.0 if pe == 'n' else 0.0
        # (1 + a^2 + d) Bn1 = Bn - i a (U - i s cp P)
        Bn1 = (Bn - 1j*a*(U - 1j*s*cp*P)) / (1 + a*a + d)
        Uss = U - 1j*a*Bn1
        P1 = (P - 1j*s*Uss) / (1 + s*s)
        U1 = Uss - 1j*s*P1
        M[:, col] = [U1, Bn1, P1]
    return np.max(np.abs(np.linalg.eigvals(M)))

if __name__ == '__main__':
    print("max |lambda| of the stage amplification matrix (backward-Euler stage)")
    print(" a=cA h k   s=cs h k   explicit p^e   no p^e    p^e + nu=cA dx (d=a*kdx)")
    for a in [0.1, 0.3, 1, 3, 10]:
        for s in [0.1, 0.3, 1, 3]:
            print(f"{a:8.2f} {s:8.2f}   {amp(a,s,'n'):10.4f}   {amp(a,s,'none'):8.4f}   {amp(a,s,'n',d=a*np.pi):10.4f}")
