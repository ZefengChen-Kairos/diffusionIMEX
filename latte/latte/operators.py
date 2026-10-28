"""Metric-consistent diffusion-type operator on the sigma mesh.

For a scalar phi and coefficient c (1/rho for the pressure projection, mu for viscosity)
the discrete energy

  Q(phi) = sum_faces T_f d_f^2 + sum_edges W_e [ -2 z_e a_e b_e + z_e^2 b_e^2 ]

approximates  int c |grad phi|^2 dV  written in sigma coordinates,
  |grad phi|^2 = phi_x^2 + phi_y^2 + phi_z^2 - 2 z_x phi_x phi_z - 2 z_y phi_y phi_z
                 + (z_x^2 + z_y^2) phi_z^2,
with phi_x the derivative along the layer, phi_z the true vertical derivative and z_x the
slope of the layer surfaces.  d_f is the difference of phi across face f, a_e an average of
the lateral differences at the two faces that meet edge e (volume weighted), b_e the average
vertical derivative of the two interfaces meeting at e.  Each lateral face gives half of its
volume to the edge above and below it, so by Jensen  Q >= sum_e W_e (a_e - z_e b_e)^2 >= 0:
the operator is symmetric positive semi-definite by construction (positive definite with a
Dirichlet face).  Q = d^T B d, the face flux of c grad phi is F = B d, and with G the
face-difference matrix the cell operator is L = G^T B G.
"""
import numpy as np
from scipy import sparse


class MetricOperator:
    def __init__(self, mesh, geom, cx, cy, cz, surface="dirichlet", walls=False, bottom=False, c_wall=None,
                 surface_mask=None):
        m, g = mesh, geom
        # surface_mask: columns whose free-surface face is a Dirichlet face (default all); the others
        # (inlet patch columns) carry a prescribed flux and are not degrees of freedom
        smask = np.ones(m.ncols, dtype=bool) if surface_mask is None else np.asarray(surface_mask, bool)
        scols = np.nonzero(smask)[0] if surface == "dirichlet" else np.zeros(0, dtype=np.int64)
        spos = -np.ones(m.ncols, dtype=np.int64)
        spos[scols] = np.arange(scols.size)
        self.surface_cols = scols
        N, nc, npx, npy = m.N, m.ncols, m.npx, m.npy
        nwx, nwy = m.wx_col.size, m.wy_col.size
        self.surface, self.walls, self.bottom = surface, walls, bottom
        names = ["x", "y", "z"]
        sizes = [N * npx, N * npy, (N - 1) * nc]
        if surface == "dirichlet":
            names.append("s"); sizes.append(scols.size)
        if walls:
            names += ["wx", "wy"]; sizes += [N * nwx, N * nwy]
        if bottom:
            names.append("b"); sizes.append(nc)
        offs = np.concatenate([[0], np.cumsum(sizes)]).astype(np.int64)
        self.nf = int(offs[-1])
        slot = {nm: int(offs[i]) for i, nm in enumerate(names)}
        self.slices = {nm: slice(int(offs[i]), int(offs[i + 1])) for i, nm in enumerate(names)}
        o_x, o_y, o_z = slot["x"], slot["y"], slot["z"]
        o_s = slot.get("s")
        K = np.arange(N)[:, None]
        rows, cols, vals = [], [], []

        def G_add(r, c, v):
            rows.append(np.ravel(r)); cols.append(np.ravel(c)); vals.append(np.broadcast_to(v, np.shape(np.ravel(r))).astype(float))

        # ---- face-difference matrix G (faces x cells), ghost matrix Gg (faces x ncols)
        rx = o_x + K * npx + np.arange(npx)[None, :]
        G_add(rx, m.cell(K, m.xp[None, :]), 1.0)
        G_add(rx, m.cell(K, m.xm[None, :]), -1.0)
        ry = o_y + K * npy + np.arange(npy)[None, :]
        G_add(ry, m.cell(K, m.yp[None, :]), 1.0)
        G_add(ry, m.cell(K, m.ym[None, :]), -1.0)
        Ki = np.arange(N - 1)[:, None]
        cl = np.arange(nc)[None, :]
        rz = o_z + Ki * nc + cl
        G_add(rz, m.cell(Ki + 1, cl), 1.0)
        G_add(rz, m.cell(Ki, cl), -1.0)
        grow, gcol = [], []
        if surface == "dirichlet":
            rs = o_s + np.arange(scols.size)
            G_add(rs, m.cell(N - 1, scols), -1.0)
            grow, gcol = rs, scols
        if walls:
            rwx = slot["wx"] + K * nwx + np.arange(nwx)[None, :]
            G_add(rwx, m.cell(K, m.wx_col[None, :]), -1.0)
            rwy = slot["wy"] + K * nwy + np.arange(nwy)[None, :]
            G_add(rwy, m.cell(K, m.wy_col[None, :]), -1.0)
        if bottom:
            rb = slot["b"] + np.arange(nc)
            G_add(rb, m.cell(0, np.arange(nc)), -1.0)
        self.G = sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                                   shape=(self.nf, m.nc))
        self.Gg = sparse.csr_matrix((np.ones(len(grow)), (np.asarray(grow, dtype=np.int64), np.asarray(gcol, dtype=np.int64))),
                                    shape=(self.nf, nc))

        # ---- symmetric face matrix B
        br, bc, bv = [], [], []

        def B_add(r, c, v):
            r = np.ravel(r); c = np.ravel(c); v = np.ravel(v)
            br.append(r); bc.append(c); bv.append(v)

        A = m.A
        B_add(rx, rx, cx * g.afx / m.dx)
        B_add(ry, ry, cy * g.afy / m.dy)
        B_add(rz, rz, cz[:-1] * A / g.dz_int[:-1])
        if surface == "dirichlet":
            B_add(rs, rs, cz[-1][scols] * A / g.dz_int[-1][scols])
        if walls:
            cw = c_wall
            B_add(rwx, rwx, cw * g.ax_wall / (0.5 * m.dx))
            B_add(rwy, rwy, cw * g.ay_wall / (0.5 * m.dy))
        if bottom:
            B_add(rb, rb, c_wall * A / (0.5 * g.h[0]))

        def edges(axis):
            if axis == "x":
                npair, mm, pp, Vf, cf, zed, dh, roff = npx, m.xm, m.xp, g.Vfx, cx, g.zx_edge, m.dx, o_x
            else:
                npair, mm, pp, Vf, cf, zed, dh, roff = npy, m.ym, m.yp, g.Vfy, cy, g.zy_edge, m.dy, o_y
            if npair == 0:
                return
            half = 0.5 * Vf * cf                                   # (N, npair) half-volume weights
            P = np.arange(npair)[None, :]
            # interior levels k+1/2, k = 0..N-2
            if N > 1:
                k = np.arange(N - 1)[:, None]
                f0 = roff + k * npair + P
                f1 = roff + (k + 1) * npair + P
                gm = o_z + k * nc + mm[None, :]
                gp = o_z + k * nc + pp[None, :]
                z = zed[:-1]
                We = half[:-1] + half[1:]
                bm = 0.5 / g.dz_int[k, mm[None, :]]
                bp = 0.5 / g.dz_int[k, pp[None, :]]
                for f, hw in ((f0, half[:-1]), (f1, half[1:])):
                    for gg, bb in ((gm, bm), (gp, bp)):
                        v = -z * (hw / dh) * bb
                        B_add(f, gg, v)
                        B_add(gg, f, v)
                for g1, b1 in ((gm, bm), (gp, bp)):
                    for g2, b2 in ((gm, bm), (gp, bp)):
                        B_add(g1, g2, We * z ** 2 * b1 * b2)
            if surface == "dirichlet":
                ok = (spos[mm] >= 0) & (spos[pp] >= 0)            # surface edges between two Dirichlet columns
                f0 = roff + (N - 1) * npair + P[0][ok]
                gm = o_s + spos[mm[ok]]
                gp = o_s + spos[pp[ok]]
                z = zed[-1][ok]
                We = half[-1][ok]
                half_s = half[-1][ok]
                bm = 0.5 / g.dz_int[-1, mm[ok]]
                bp = 0.5 / g.dz_int[-1, pp[ok]]
                for gg, bb in ((gm, bm), (gp, bp)):
                    v = -z * (half_s / dh) * bb
                    B_add(f0, gg, v)
                    B_add(gg, f0, v)
                for g1, b1 in ((gm, bm), (gp, bp)):
                    for g2, b2 in ((gm, bm), (gp, bp)):
                        B_add(g1, g2, We * z ** 2 * b1 * b2)

        edges("x")
        edges("y")
        self.B = sparse.csr_matrix((np.concatenate(bv), (np.concatenate(br), np.concatenate(bc))),
                                   shape=(self.nf, self.nf))
        self.B.sum_duplicates()
        self._L = None

    @property
    def L(self):
        if self._L is None:
            self._L = (self.G.T @ (self.B @ self.G)).tocsr()
        return self._L

    def ghost_rhs(self, ghost):
        """G^T B Gg ghost  (contribution of Dirichlet ghost values to the cell equations)."""
        return self.G.T @ (self.B @ (self.Gg @ ghost))

    def flux(self, phi, ghost=None):
        d = self.G @ phi
        if ghost is not None:
            d = d + self.Gg @ ghost
        return self.B @ d
