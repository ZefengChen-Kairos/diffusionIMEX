"""Shared helpers for the stage-1 suite."""
import copy
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from ..mesh import Geometry
from ..solver import Model, Numerics
from ..transport import transport


def merge(base, over):
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = merge(out[k], v)
        else:
            out[k] = v
    return out


def make_model(cfg, **over):
    d = merge(cfg["model"], over)
    d["rho"] = tuple(d["rho"])
    return Model(**d)


def make_numerics(cfg, **over):
    d = merge(cfg["numerics"], over)
    return Numerics(**d)


class Prescribed:
    """Transport-only driver with prescribed face fluxes (G2 tests)."""

    def __init__(self, mesh, D, C, rho, v=None, scheme="fct", recon="minmod", ho_limiter="superbee"):
        self.m = mesh
        self.g = Geometry(mesh, np.broadcast_to(np.asarray(D, float), (mesh.ncols,)).copy())
        C = np.asarray(C, float).reshape(3, mesh.nc)
        self.rho = np.asarray(rho, float)
        self.Vtot = self.g.V.ravel().copy()
        self.Vc = C[:2] * self.Vtot[None, :]
        self.M = (self.rho @ C) * self.Vtot
        v = np.zeros((3, mesh.nc)) if v is None else np.asarray(v, float).reshape(3, mesh.nc)
        self.P = self.M[None, :] * v
        self.scheme, self.recon, self.ho = scheme, recon, ho_limiter
        self.t = 0.0
        self.steps = 0
        self.inflow_tot = np.zeros(3)
        self.outflow_tot = np.zeros(3)
        self.max_cfl = 0.0

    def C(self):
        VG = self.Vtot - self.Vc[0] - self.Vc[1]
        return np.stack([self.Vc[0], self.Vc[1], VG]) / self.Vtot[None, :]

    def step(self, dV, dt, inflow=None, outflow=None):
        r = transport(self.m, self.Vc, self.Vtot, self.P, self.M, dV, self.rho, inflow, outflow,
                      self.scheme, self.recon, 1e-12, self.ho)
        self.Vc, self.Vtot, self.M, self.P = r.Vc, r.Vtot, r.M, r.P
        self.inflow_tot += r.in_comp.sum(1)
        self.outflow_tot += r.out_comp.sum(1)
        self.t += dt
        self.steps += 1
        self.max_cfl = max(self.max_cfl, r.cfl)
        return r


def l1(a, b, w):
    return float(np.sum(np.abs(a - b) * w) / np.sum(w))


def rate(errs, grids):
    e = np.asarray(errs, float); n = np.asarray(grids, float)
    return [float(np.log(e[i] / e[i + 1]) / np.log(n[i + 1] / n[i])) for i in range(len(e) - 1)]


def savefig(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def disc(x, y, xc, yc, r):
    return ((x - xc) ** 2 + (y - yc) ** 2 < r ** 2).astype(float)


def subcell_disc(mesh, xc, yc, r, sub=8):
    """Cell-average of the indicator of a disc (sub x sub sampling per cell)."""
    o = (np.arange(sub) + 0.5) / sub - 0.5
    acc = np.zeros(mesh.ncols)
    for a in o:
        for b in o:
            acc += disc(mesh.xc + a * mesh.dx, mesh.yc + b * mesh.dy, xc, yc, r)
    return acc / sub ** 2
