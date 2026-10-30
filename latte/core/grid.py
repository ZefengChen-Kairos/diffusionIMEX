"""Cartesian grids with a disc wet mask, kernels on local windows, coarse/fine transfer."""
import numpy as np
from scipy import ndimage


class Grid:
    def __init__(self, L, N, R_cup):
        self.L, self.N, self.R = L, N, R_cup
        self.dx = L / N
        self.dA = self.dx * self.dx
        # exactly antisymmetric cell centres (mirror tests are bitwise on the geometry)
        self.xc = (np.arange(N) - 0.5 * (N - 1)) * self.dx
        self.X, self.Y = np.meshgrid(self.xc, self.xc)          # [j, i] = (y, x)
        self.wet = (self.X ** 2 + self.Y ** 2) < R_cup ** 2
        w = self.wet
        # open interior faces
        self.ox = w[:, :-1] & w[:, 1:]                             # (N, N-1), face i+1/2
        self.oy = w[:-1, :] & w[1:, :]                             # (N-1, N)
        # closed faces bordering a wet cell: right/left (x), top/bottom (y)
        oxp = np.zeros((N, N + 1), bool)
        oxp[:, 1:-1] = self.ox
        oyp = np.zeros((N + 1, N), bool)
        oyp[1:-1, :] = self.oy
        self.oxp, self.oyp = oxp, oyp
        self.wall_xr = w & ~oxp[:, 1:]
        self.wall_xl = w & ~oxp[:, :-1]
        self.wall_yt = w & ~oyp[1:, :]
        self.wall_yb = w & ~oyp[:-1, :]
        # nearest wet cell for every cell (extrapolation outside the mask)
        _, idx = ndimage.distance_transform_edt(~w, return_indices=True)
        self.near_j, self.near_i = idx[0], idx[1]
        self.dry = ~w

    def index_of(self, x):
        """Continuous index coordinates (j, i) of physical points x[..., (x, y)]."""
        i = (x[..., 0] + 0.5 * self.L) / self.dx - 0.5
        j = (x[..., 1] + 0.5 * self.L) / self.dx - 0.5
        return j, i

    def window(self, x0, rad):
        """Index slices of the cells within a square of half-width rad around x0."""
        N = self.N
        i0 = max(int(np.floor((x0[0] - rad + 0.5 * self.L) / self.dx)), 0)
        i1 = min(int(np.ceil((x0[0] + rad + 0.5 * self.L) / self.dx)) + 1, N)
        j0 = max(int(np.floor((x0[1] - rad + 0.5 * self.L) / self.dx)), 0)
        j1 = min(int(np.ceil((x0[1] + rad + 0.5 * self.L) / self.dx)) + 1, N)
        return slice(j0, max(j1, j0)), slice(i0, max(i1, i0))

    def kernel(self, x0, kind, a, b=0.0):
        """Kernel normalised on the wet cells, integral K dA = 1.

        kind 'ring': exp(-((rho - b)/a)^2) with ring radius b and width a
        kind 'gauss': exp(-rho^2 / (2 a^2))
        Returns (window slices, values on the window).
        """
        if kind == "ring":
            rad = b + 3.5 * a
        else:
            rad = 4.0 * a
        # box with a one-cell margin; truncation by distance keeps the kernel mirror symmetric
        sj, si = self.window(x0, rad + self.dx)
        X = self.X[sj, si] - x0[0]
        Y = self.Y[sj, si] - x0[1]
        rho = np.sqrt(X * X + Y * Y)
        if kind == "ring":
            K = np.exp(-((rho - b) / a) ** 2)
        else:
            K = np.exp(-0.5 * (rho / a) ** 2)
        K *= self.wet[sj, si] & (rho <= rad)
        s = K.sum() * self.dA
        if s <= 0:
            return sj, si, K
        return sj, si, K / s

    def extrapolate(self, f):
        """Fill dry cells with the value of the nearest wet cell (in place, last two axes)."""
        f[..., self.dry] = f[..., self.near_j[self.dry], self.near_i[self.dry]]
        return f


def bilinear(field, j, i, order=1):
    """Sample a cell-centred 2D field at continuous indices (edge clamped)."""
    return ndimage.map_coordinates(field, [j, i], order=order, mode="nearest")


class Bilin:
    """Bilinear sampling weights at fixed continuous indices, shared by many fields."""

    def __init__(self, j, i, n):
        j = np.clip(j, 0.0, n - 1.0)
        i = np.clip(i, 0.0, n - 1.0)
        j0 = np.minimum(np.floor(j).astype(np.int64), n - 2)
        i0 = np.minimum(np.floor(i).astype(np.int64), n - 2)
        self.wj = j - j0
        self.wi = i - i0
        self.k00 = (j0 * n + i0).ravel()
        self.n = n
        self.shape = j.shape

    def __call__(self, f):
        """f: (..., n, n) -> (..., *shape)."""
        lead = f.shape[:-2]
        ff = f.reshape(lead + (-1,))
        k = self.k00
        n = self.n
        a = ff[..., k]
        b = ff[..., k + 1]
        c = ff[..., k + n]
        d = ff[..., k + n + 1]
        wj = self.wj.ravel()
        wi = self.wi.ravel()
        out = (1 - wj) * ((1 - wi) * a + wi * b) + wj * ((1 - wi) * c + wi * d)
        return out.reshape(lead + self.shape)

    def minmax(self, f):
        lead = f.shape[:-2]
        ff = f.reshape(lead + (-1,))
        k = self.k00
        n = self.n
        a, b, c, d = ff[..., k], ff[..., k + 1], ff[..., k + n], ff[..., k + n + 1]
        lo = np.minimum(np.minimum(a, b), np.minimum(c, d)).reshape(lead + self.shape)
        hi = np.maximum(np.maximum(a, b), np.maximum(c, d)).reshape(lead + self.shape)
        return lo, hi


def _interp_matrix(n_f, n_c):
    r = n_f / n_c
    idx = np.clip((np.arange(n_f) + 0.5) / r - 0.5, 0.0, n_c - 1.0)
    i0 = np.minimum(np.floor(idx).astype(int), n_c - 2)
    w = idx - i0
    W = np.zeros((n_f, n_c))
    W[np.arange(n_f), i0] = 1 - w
    W[np.arange(n_f), i0 + 1] += w
    return W


class FineMap:
    """Coarse -> fine sampling helpers (fine cell centres in coarse index space)."""

    def __init__(self, coarse: Grid, fine: Grid):
        self.c, self.f = coarse, fine
        r = fine.N / coarse.N
        idx = (np.arange(fine.N) + 0.5) / r - 0.5
        self.J, self.I = np.meshgrid(idx, idx, indexing="ij")
        self.r = int(round(r))
        self.W = _interp_matrix(fine.N, coarse.N)      # separable bilinear (edge clamped)

    def to_fine(self, field):
        return self.W @ field @ self.W.T

    def to_fine_window(self, field, sj, si):
        return self.W[sj] @ field @ self.W[si].T

    def restrict(self, ff):
        r = self.r
        N = self.c.N
        return ff.reshape(N, r, N, r).mean(axis=(1, 3))
