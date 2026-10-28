"""Mesh, topology and geometry for the layered (sigma-coordinate) cup solver.

Horizontal: uniform Cartesian columns (dx, dy); the cup cross-section is a set of
fluid columns (staircase mask).  Vertical: N terrain-following layers over a flat
bottom z = 0 that follow the free surface, z_{k+1/2} = S_{k+1/2} D(x, y, t), where
D is the local liquid depth and S the cumulative layer fractions.

Index conventions
-----------------
cell (k, col)   flat id  k * ncols + col,   k = 0 is the bottom layer
x-face (k, p)   between columns xm[p] (minus side) and xp[p] (plus side), +x oriented
y-face (k, q)   between columns ym[q] and yp[q], +y oriented
interface (k, col)  top face of cell (k, col), +z oriented; k = N-1 is the free surface
wall faces      lateral faces of fluid columns that touch a solid column / domain wall
"""
import numpy as np


class Mesh:
    def __init__(self, nx, ny, dx, dy, fractions, x0=0.0, y0=0.0,
                 periodic_x=False, periodic_y=False, mask=None):
        self.nx, self.ny = int(nx), int(ny)
        self.dx, self.dy = float(dx), float(dy)
        self.x0, self.y0 = float(x0), float(y0)
        self.periodic_x, self.periodic_y = bool(periodic_x), bool(periodic_y)
        s = np.asarray(fractions, dtype=float)
        if s.ndim != 1 or np.any(s <= 0) or abs(s.sum() - 1.0) > 1e-14:
            raise ValueError("layer fractions must be positive and sum to one")
        self.s = s
        self.N = s.size
        self.S = np.concatenate([[0.0], np.cumsum(s)])
        self.S[-1] = 1.0
        self.A = self.dx * self.dy

        fluid = np.ones((self.ny, self.nx), dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
        if fluid.shape != (self.ny, self.nx):
            raise ValueError("mask shape must be (ny, nx)")
        self.fluid = fluid
        self.cj, self.ci = np.nonzero(fluid)
        self.ncols = self.ci.size
        col_of = -np.ones((self.ny, self.nx), dtype=np.int64)
        col_of[self.cj, self.ci] = np.arange(self.ncols)
        self.col_of = col_of
        self.xc = self.x0 + (self.ci + 0.5) * self.dx
        self.yc = self.y0 + (self.cj + 0.5) * self.dy
        self.nc = self.N * self.ncols

        # neighbour columns: >= 0 fluid column, -1 wall, -2 periodic self (no face, no wall)
        self.right = self._neighbour(+1, 0)
        self.left = self._neighbour(-1, 0)
        self.up = self._neighbour(0, +1)
        self.down = self._neighbour(0, -1)

        cols = np.arange(self.ncols)
        hx = self.right >= 0
        self.xm, self.xp = cols[hx], self.right[hx]
        hy = self.up >= 0
        self.ym, self.yp = cols[hy], self.up[hy]
        self.npx, self.npy = self.xm.size, self.ym.size
        # wall faces: (column, outward sign); outward normal is sign * e_axis
        self.wx_col = np.concatenate([cols[self.right == -1], cols[self.left == -1]])
        self.wx_sgn = np.concatenate([np.ones((self.right == -1).sum()), -np.ones((self.left == -1).sum())])
        self.wy_col = np.concatenate([cols[self.up == -1], cols[self.down == -1]])
        self.wy_sgn = np.concatenate([np.ones((self.up == -1).sum()), -np.ones((self.down == -1).sum())])
        # pair index of the x-face on the left/right of each column (-1 if none)
        self.xpair_right = -np.ones(self.ncols, dtype=np.int64)
        self.xpair_left = -np.ones(self.ncols, dtype=np.int64)
        self.xpair_right[self.xm] = np.arange(self.npx)
        self.xpair_left[self.xp] = np.arange(self.npx)
        self.ypair_up = -np.ones(self.ncols, dtype=np.int64)
        self.ypair_down = -np.ones(self.ncols, dtype=np.int64)
        self.ypair_up[self.ym] = np.arange(self.npy)
        self.ypair_down[self.yp] = np.arange(self.npy)
        self._build_transport_topology()

    # ------------------------------------------------------------------ topology
    def _neighbour(self, di, dj):
        out = np.empty(self.ncols, dtype=np.int64)
        i2 = self.ci + di
        j2 = self.cj + dj
        if di != 0:
            if self.periodic_x:
                i2 = i2 % self.nx
            valid = (i2 >= 0) & (i2 < self.nx)
        else:
            if self.periodic_y:
                j2 = j2 % self.ny
            valid = (j2 >= 0) & (j2 < self.ny)
        i2c = np.clip(i2, 0, self.nx - 1)
        j2c = np.clip(j2, 0, self.ny - 1)
        nb = np.where(valid, self.col_of[j2c, i2c], -1)
        out[:] = nb
        out[out == np.arange(self.ncols)] = -2
        return out

    def cell(self, k, col):
        return np.asarray(k) * self.ncols + np.asarray(col)

    def _build_transport_topology(self):
        """Face lists [x faces, y faces, interior interfaces] with the cells on both
        sides and the next cells upstream of each side (for limited reconstruction)."""
        N, nc = self.N, self.ncols
        K = np.arange(N)[:, None]
        cm_x = self.cell(K, self.xm[None, :]).ravel()
        cp_x = self.cell(K, self.xp[None, :]).ravel()
        lm = self.left[self.xm]
        rp = self.right[self.xp]
        um_x = np.where(lm[None, :] >= 0, self.cell(K, np.maximum(lm, 0)[None, :]), -1).ravel()
        up_x = np.where(rp[None, :] >= 0, self.cell(K, np.maximum(rp, 0)[None, :]), -1).ravel()
        cm_y = self.cell(K, self.ym[None, :]).ravel()
        cp_y = self.cell(K, self.yp[None, :]).ravel()
        dm = self.down[self.ym]
        upp = self.up[self.yp]
        um_y = np.where(dm[None, :] >= 0, self.cell(K, np.maximum(dm, 0)[None, :]), -1).ravel()
        up_y = np.where(upp[None, :] >= 0, self.cell(K, np.maximum(upp, 0)[None, :]), -1).ravel()
        Ki = np.arange(N - 1)[:, None]
        colsr = np.arange(nc)[None, :]
        cm_z = self.cell(Ki, colsr).ravel()
        cp_z = self.cell(Ki + 1, colsr).ravel()
        um_z = np.where(Ki >= 1, self.cell(np.maximum(Ki - 1, 0), colsr), -1).ravel()
        up_z = np.where(Ki + 2 <= N - 1, self.cell(np.minimum(Ki + 2, N - 1), colsr), -1).ravel()
        self.tf_cm = np.concatenate([cm_x, cm_y, cm_z])
        self.tf_cp = np.concatenate([cp_x, cp_y, cp_z])
        self.tf_um = np.concatenate([um_x, um_y, um_z])
        self.tf_up = np.concatenate([up_x, up_y, up_z])
        self.n_tf = self.tf_cm.size

    # ------------------------------------------------------------------ helpers
    def column_field(self, f2d):
        """Sample a (ny, nx) array on fluid columns."""
        return np.asarray(f2d)[self.cj, self.ci]

    def to_grid(self, colvals, fill=np.nan):
        g = np.full((self.ny, self.nx), fill, dtype=float)
        g[self.cj, self.ci] = colvals
        return g

    def column_at(self, x, y):
        """Column id containing horizontal point(s) (x, y); -1 outside fluid."""
        x = np.atleast_1d(np.asarray(x, dtype=float))
        y = np.atleast_1d(np.asarray(y, dtype=float))
        fi = (x - self.x0) / self.dx
        fj = (y - self.y0) / self.dy
        if self.periodic_x:
            fi = np.mod(fi, self.nx)
        if self.periodic_y:
            fj = np.mod(fj, self.ny)
        i = np.floor(fi).astype(np.int64)
        j = np.floor(fj).astype(np.int64)
        inside = (i >= 0) & (i < self.nx) & (j >= 0) & (j < self.ny)
        out = -np.ones(i.shape, dtype=np.int64)
        out[inside] = self.col_of[j[inside], i[inside]]
        return out

    def area_total(self):
        return self.ncols * self.A


class Geometry:
    """Geometric quantities of the layered mesh for a given depth field D (per column)."""

    def __init__(self, mesh, D):
        m = mesh
        D = np.asarray(D, dtype=float)
        if D.shape != (m.ncols,):
            raise ValueError("D must have one value per fluid column")
        self.D = D
        s, S = m.s[:, None], m.S
        self.h = s * D[None, :]                                   # (N, ncols)
        self.zc = (S[:-1, None] + 0.5 * s) * D[None, :]
        self.zi = S[1:, None] * D[None, :]                        # top interface heights
        self.V = self.h * m.A
        Dfx = 0.5 * (D[m.xm] + D[m.xp])
        Dfy = 0.5 * (D[m.ym] + D[m.yp])
        self.hfx = s * Dfx[None, :]
        self.hfy = s * Dfy[None, :]
        self.afx = self.hfx * m.dy
        self.afy = self.hfy * m.dx
        self.Vfx = self.afx * m.dx
        self.Vfy = self.afy * m.dy
        dz = np.empty_like(self.h)
        dz[:-1] = 0.5 * (self.h[:-1] + self.h[1:])
        dz[-1] = 0.5 * self.h[-1]
        self.dz_int = dz                                          # centre distance across each top face
        # free-surface slope at x/y faces and metric slope at edges (level k+1/2, k=0..N-1)
        self.etax_face = (D[m.xp] - D[m.xm]) / m.dx
        self.etay_face = (D[m.yp] - D[m.ym]) / m.dy
        self.zx_edge = S[1:, None] * self.etax_face[None, :]
        self.zy_edge = S[1:, None] * self.etay_face[None, :]
        # column-centred surface slope: average of both sides, a wall side counts as zero slope
        gx = np.zeros(m.ncols)
        np.add.at(gx, m.xm, 0.5 * self.etax_face)
        np.add.at(gx, m.xp, 0.5 * self.etax_face)
        gy = np.zeros(m.ncols)
        np.add.at(gy, m.ym, 0.5 * self.etay_face)
        np.add.at(gy, m.yp, 0.5 * self.etay_face)
        self.etax_col, self.etay_col = gx, gy
        self.zx_int = S[1:, None] * gx[None, :]                   # slope of top interface of each cell
        self.zy_int = S[1:, None] * gy[None, :]
        self.ax_wall = (s * D[m.wx_col][None, :]) * m.dy          # wall face areas (N, nwx)
        self.ay_wall = (s * D[m.wy_col][None, :]) * m.dx
        self.max_slope = float(np.max(np.abs(np.concatenate([self.etax_face, self.etay_face, [0.0]]))))


def layer_fractions(N, kind="uniform", ratio=1.0):
    """Layer fractions from bottom to top.  'geometric': each layer is `ratio` times the one
    below it (ratio < 1 refines towards the surface)."""
    if kind == "uniform":
        return np.full(N, 1.0 / N)
    if kind == "geometric":
        w = ratio ** np.arange(N)
        return w / w.sum()
    raise ValueError(kind)


def cylinder_mesh(n, radius, fractions, center=(0.0, 0.0)):
    """Straight-walled cylindrical cup on an n x n grid covering [-R, R]^2 (staircase)."""
    dx = 2.0 * radius / n
    x0, y0 = center[0] - radius, center[1] - radius
    xc = x0 + (np.arange(n) + 0.5) * dx
    X, Y = np.meshgrid(xc, y0 + (np.arange(n) + 0.5) * dx)
    mask = (X - center[0]) ** 2 + (Y - center[1]) ** 2 < radius ** 2
    return Mesh(n, n, dx, dx, fractions, x0=x0, y0=y0, mask=mask)
