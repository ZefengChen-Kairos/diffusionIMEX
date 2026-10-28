"""Handoff section Gamma_J and the inlet input contract (V3 Sec. 2.2, 5.1).

Internally an arrival is a set of quadrature samples {x_q, dA_q, j_q, v_q, C_q}: the
parametric ellipse below produces them, and a future MPM adapter can produce the same
samples from particles crossing the handoff section.  Volume, component volumes, mass
and advective momentum all come from the same samples and the same quadrature.

Sign convention: n is the outward normal of the cup-liquid domain on Gamma_J;
the inflow flux density is j = -(v_J - w_J) . n >= 0.
"""
from dataclasses import dataclass, field
import numpy as np

from .failures import SimFailure

FLOW, NO_FLOW, UNKNOWN = "flow", "no_flow", "unknown"


@dataclass
class InletSection:
    center: np.ndarray            # (3,) world position of the section centre
    normal: np.ndarray            # (3,) unit outward normal of the liquid domain
    axis1: np.ndarray             # (3,) unit in-plane axis of semi-axis a
    a: float
    b: float
    velocity: np.ndarray          # (3,) material arrival velocity (peak value for 'parabolic')
    composition: np.ndarray       # (3,) volume fractions (C_B, C_M, C_G)
    section_velocity: np.ndarray = field(default_factory=lambda: np.zeros(3))  # w_J
    profile: str = "uniform"      # 'uniform' or 'parabolic' (v = velocity * (1 - rho^2))

    def validate(self):
        n = np.asarray(self.normal, float)
        e1 = np.asarray(self.axis1, float)
        if not (np.all(np.isfinite(n)) and np.all(np.isfinite(e1)) and np.all(np.isfinite(self.center))):
            raise SimFailure("GEOMETRY_INVALID", "non-finite inlet geometry")
        if abs(np.linalg.norm(n) - 1) > 1e-12 or abs(np.linalg.norm(e1) - 1) > 1e-12 or abs(n @ e1) > 1e-12:
            raise SimFailure("GEOMETRY_INVALID", "inlet normal/axis must be orthonormal")
        if not (self.a > 0 and self.b > 0):
            raise SimFailure("GEOMETRY_INVALID", "inlet semi-axes must be positive")
        C = np.asarray(self.composition, float)
        if C.shape != (3,) or np.any(C < 0) or abs(C.sum() - 1) > 1e-12:
            raise SimFailure("INPUT_INVALID", "inlet composition must be a point of the simplex")
        if C[0] + C[1] <= 0:
            raise SimFailure("INPUT_INVALID", "inlet has no liquid fraction")
        if self.profile not in ("uniform", "parabolic"):
            raise SimFailure("INPUT_INVALID", f"unknown velocity profile {self.profile}")
        if not np.all(np.isfinite(self.velocity)) or not np.all(np.isfinite(self.section_velocity)):
            raise SimFailure("INPUT_INVALID", "non-finite inlet velocity")

    def samples(self, n_r=8, n_theta=32):
        """Polar Gauss quadrature on the ellipse.  Returns points (n,3), area weights (n,),
        flux density j (n,), material velocity (n,3)."""
        self.validate()
        n = np.asarray(self.normal, float)
        e1 = np.asarray(self.axis1, float)
        e2 = np.cross(n, e1)
        xg, wg = np.polynomial.legendre.leggauss(n_r)
        r = 0.5 * (xg + 1.0)
        wr = 0.5 * wg * r                       # weight of r dr on [0, 1]
        th = (np.arange(n_theta) + 0.5) * 2 * np.pi / n_theta
        wt = 2 * np.pi / n_theta
        R, TH = np.meshgrid(r, th, indexing="ij")
        W = (wr[:, None] * np.full(n_theta, wt)[None, :]) * self.a * self.b
        pts = (np.asarray(self.center, float)[None, :]
               + (self.a * R * np.cos(TH)).reshape(-1, 1) * e1[None, :]
               + (self.b * R * np.sin(TH)).reshape(-1, 1) * e2[None, :])
        shape = np.ones_like(R) if self.profile == "uniform" else (1.0 - R ** 2)
        v = shape.reshape(-1, 1) * np.asarray(self.velocity, float)[None, :]
        j = -((v - np.asarray(self.section_velocity, float)[None, :]) @ n)
        return pts, W.ravel(), j, v


@dataclass
class InletState:
    status: str
    section: InletSection = None
    scan_speed: float = 0.0       # horizontal speed of the section centre (for the time-step limit)


class NoInlet:
    """Provider for runs without any arrival (confirmed zero flow at all times)."""

    def sample(self, t):
        return InletState(NO_FLOW)

    def next_event(self, t):
        return np.inf


class PiecewiseInlet:
    """Inlet nodes in physical seconds.

    Interval [t_i, t_{i+1}) takes the status of node i.  Section parameters are linearly
    interpolated only when nodes i and i+1 are both 'flow'; a 'flow' node followed by a
    status change is held constant up to the change.  'unknown' intervals and times
    outside [t_0, t_last] raise INPUT_INVALID (unknown is never treated as zero).
    The solver must not step across node times (see next_event).
    """

    def __init__(self, nodes):
        if len(nodes) < 1:
            raise SimFailure("INPUT_INVALID", "inlet needs at least one node")
        self.t = np.array([nd["t"] for nd in nodes], dtype=float)
        if np.any(np.diff(self.t) <= 0):
            raise SimFailure("INPUT_INVALID", "inlet node times must increase strictly")
        self.nodes = nodes
        for nd in nodes:
            if nd["status"] not in (FLOW, NO_FLOW, UNKNOWN):
                raise SimFailure("INPUT_INVALID", f"bad inlet status {nd['status']}")
            if nd["status"] == FLOW:
                nd["section"].validate()

    def _index(self, t):
        if not np.isfinite(t) or t < self.t[0] or t > self.t[-1]:
            raise SimFailure("INPUT_INVALID", f"t={t} outside the inlet record [{self.t[0]}, {self.t[-1]}]")
        i = int(np.searchsorted(self.t, t, side="right") - 1)
        return min(i, len(self.t) - 1)

    def sample(self, t):
        i = self._index(t)
        nd = self.nodes[i]
        if nd["status"] == UNKNOWN:
            raise SimFailure("INPUT_INVALID", f"inlet state unknown at t={t}")
        if nd["status"] == NO_FLOW:
            return InletState(NO_FLOW)
        sec = nd["section"]
        if i + 1 < len(self.t) and self.nodes[i + 1]["status"] == FLOW:
            s2 = self.nodes[i + 1]["section"]
            w = (t - self.t[i]) / (self.t[i + 1] - self.t[i])
            lerp = lambda a, b: (1 - w) * np.asarray(a, float) + w * np.asarray(b, float)
            n = lerp(sec.normal, s2.normal)
            n /= np.linalg.norm(n)
            e1 = lerp(sec.axis1, s2.axis1)
            e1 -= (e1 @ n) * n
            e1 /= np.linalg.norm(e1)
            new = InletSection(center=lerp(sec.center, s2.center), normal=n, axis1=e1,
                               a=float(lerp(sec.a, s2.a)), b=float(lerp(sec.b, s2.b)),
                               velocity=lerp(sec.velocity, s2.velocity),
                               composition=lerp(sec.composition, s2.composition),
                               section_velocity=(np.asarray(s2.center, float) - np.asarray(sec.center, float))
                               / (self.t[i + 1] - self.t[i]),
                               profile=sec.profile)
            scan = float(np.linalg.norm(new.section_velocity[:2]))
            return InletState(FLOW, new, scan)
        return InletState(FLOW, sec, float(np.linalg.norm(np.asarray(sec.section_velocity)[:2])))

    def next_event(self, t):
        later = self.t[self.t > t]
        return float(later[0]) if later.size else np.inf


@dataclass
class SurfaceSource:
    """Per-column arrival rates through the free surface (all zero if no flow)."""
    vol_rate: np.ndarray          # (ncols,) volume per second
    comp_rate: np.ndarray         # (3, ncols)
    mass_rate: np.ndarray         # (ncols,)
    mom_rate: np.ndarray          # (3, ncols) advective momentum per second
    Q: float = 0.0
    abs_rate: np.ndarray = None   # (ncols,) absolute volume flux of material through the patch, sum v_J.n dA (outward +)
    cover: np.ndarray = None      # (ncols,) horizontally projected covered area / column area


def zero_source(ncols):
    return SurfaceSource(np.zeros(ncols), np.zeros((3, ncols)), np.zeros(ncols), np.zeros((3, ncols)), 0.0,
                         np.zeros(ncols), np.zeros(ncols))


def surface_source(state, mesh, D, rho_components, n_r=8, n_theta=32, connect_tol=None, max_spacing=None):
    """Map an arrival onto surface columns by the horizontal position of each quadrature
    sample.  Positive flux outside the fluid columns, or a section that is not connected
    to the liquid surface, is an error; nothing is renormalised back into the cup."""
    if state.status == NO_FLOW:
        return zero_source(mesh.ncols)
    if state.status != FLOW:
        raise SimFailure("INPUT_INVALID", f"inlet status {state.status}")
    if max_spacing is not None:
        # refine the quadrature so that samples are at most `max_spacing` apart (column shares then
        # converge to the geometric overlap instead of aliasing onto the columns that happen to hold a sample)
        R = max(state.section.a, state.section.b)
        n_r = max(n_r, int(np.ceil(R / max_spacing)))
        n_theta = max(n_theta, int(np.ceil(2 * np.pi * R / max_spacing)))
    pts, dA, j, v = state.section.samples(n_r, n_theta)
    if np.any(j < -1e-14 * max(1.0, np.abs(j).max())):
        raise SimFailure("INPUT_INVALID", "inlet has outflow (j < 0) on the handoff section")
    j = np.maximum(j, 0.0)
    w = j * dA
    col = mesh.column_at(pts[:, 0], pts[:, 1])
    if np.any((col < 0) & (w > 0)):
        raise SimFailure("OUT_OF_SCOPE", "inlet footprint leaves the fluid columns (no renormalisation)")
    if connect_tol is not None:
        gap = np.abs(pts[:, 2] - D[np.maximum(col, 0)])
        if np.any((gap > connect_tol) & (w > 0)):
            raise SimFailure("GEOMETRY_INVALID", "handoff section is not connected to the liquid surface",
                             max_gap=float(gap[w > 0].max()))
    C = np.asarray(state.section.composition, float)
    rho_in = float(np.asarray(rho_components, float) @ C)
    nc = mesh.ncols
    vol = np.bincount(col[w > 0], weights=w[w > 0], minlength=nc)
    comp = C[:, None] * vol[None, :]
    mom = np.stack([np.bincount(col[w > 0], weights=(rho_in * w * v[:, d])[w > 0], minlength=nc) for d in range(3)])
    n = np.asarray(state.section.normal, float)
    vn = v @ n                                            # absolute normal velocity of the arriving material
    pos = w > 0
    absr = np.bincount(col[pos], weights=(vn * dA)[pos], minlength=nc)
    cover = np.bincount(col[pos], weights=(dA * abs(n[2]))[pos], minlength=nc) / mesh.A
    out = SurfaceSource(vol, comp, rho_in * vol, mom, float(w.sum()), absr, cover)
    out.quadrature = (n_r, n_theta)
    return out


def surface_source_subgrid(state, mesh, D, rho_components, sub=8, connect_tol=None, n_r=16, n_theta=64):
    """Mesh-aligned distribution of an arrival: every column in the footprint bounding box is sampled
    with sub x sub points, each projected vertically onto the section plane and tested against the
    ellipse (no aliasing between sample spacing and columns).  Totals come from the exact polar Gauss
    integral; one global factor (1 + O(sub^-2) boundary quadrature error) reconciles the two.  Samples
    with positive flux outside the fluid columns are an error (nothing is moved back into the cup)."""
    if state.status == NO_FLOW:
        return zero_source(mesh.ncols)
    if state.status != FLOW:
        raise SimFailure("INPUT_INVALID", f"inlet status {state.status}")
    sec = state.section
    sec.validate()
    n = np.asarray(sec.normal, float); e1 = np.asarray(sec.axis1, float); e2 = np.cross(n, e1)
    c = np.asarray(sec.center, float)
    if abs(n[2]) < 0.2:
        raise SimFailure("GEOMETRY_INVALID", "handoff section too steep to map onto the free surface (|n_z| < 0.2)")
    hx = np.sqrt((sec.a * e1[0]) ** 2 + (sec.b * e2[0]) ** 2)
    hy = np.sqrt((sec.a * e1[1]) ** 2 + (sec.b * e2[1]) ** 2)
    i0 = int(np.floor((c[0] - hx - mesh.x0) / mesh.dx)) - 1; i1 = int(np.floor((c[0] + hx - mesh.x0) / mesh.dx)) + 1
    j0 = int(np.floor((c[1] - hy - mesh.y0) / mesh.dy)) - 1; j1 = int(np.floor((c[1] + hy - mesh.y0) / mesh.dy)) + 1
    o = (np.arange(sub) + 0.5) / sub
    xs = mesh.x0 + (np.arange(i0, i1 + 1)[:, None] + o[None, :]).ravel() * mesh.dx
    ys = mesh.y0 + (np.arange(j0, j1 + 1)[:, None] + o[None, :]).ravel() * mesh.dy
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    X = X.ravel(); Y = Y.ravel()
    Z = c[2] - (n[0] * (X - c[0]) + n[1] * (Y - c[1])) / n[2]
    P = np.stack([X - c[0], Y - c[1], Z - c[2]], axis=1)
    xi = (P @ e1) / sec.a; et = (P @ e2) / sec.b
    rr = xi ** 2 + et ** 2
    inside = rr < 1.0
    X, Y, Z, rr = X[inside], Y[inside], Z[inside], rr[inside]
    shape = np.ones_like(rr) if sec.profile == "uniform" else (1.0 - rr)
    v = shape[:, None] * np.asarray(sec.velocity, float)[None, :]
    j = -((v - np.asarray(sec.section_velocity, float)[None, :]) @ n)
    if np.any(j < -1e-14 * max(1.0, np.abs(j).max(initial=0.0))):
        raise SimFailure("INPUT_INVALID", "inlet has outflow (j < 0) on the handoff section")
    j = np.maximum(j, 0.0)
    dAs = (mesh.dx / sub) * (mesh.dy / sub) / abs(n[2])              # section area per sample
    w = j * dAs
    col = mesh.column_at(X, Y)
    if np.any((col < 0) & (w > 0)):
        raise SimFailure("OUT_OF_SCOPE", "inlet footprint leaves the fluid columns (no renormalisation)")
    if connect_tol is not None:
        gap = np.abs(Z - D[np.maximum(col, 0)])
        if np.any((gap > connect_tol) & (w > 0)):
            raise SimFailure("GEOMETRY_INVALID", "handoff section is not connected to the liquid surface",
                             max_gap=float(gap[w > 0].max()))
    exact = section_totals(sec, rho_components, n_r, n_theta)
    Qs = float(w.sum())
    f = exact["Q"] / Qs if Qs > 0 else 0.0
    C = np.asarray(sec.composition, float)
    rho_in = float(np.asarray(rho_components, float) @ C)
    nc = mesh.ncols
    pos = w > 0
    vol = f * np.bincount(col[pos], weights=w[pos], minlength=nc)
    mom = f * np.stack([np.bincount(col[pos], weights=(rho_in * w * v[:, d])[pos], minlength=nc) for d in range(3)])
    absr = f * np.bincount(col[pos], weights=((v @ n) * dAs)[pos], minlength=nc)
    cover = np.bincount(col[pos], weights=np.full(pos.sum(), dAs * abs(n[2])), minlength=nc) / mesh.A
    out = SurfaceSource(vol, C[:, None] * vol[None, :], rho_in * vol, mom, float(vol.sum()), absr, cover)
    out.quadrature = ("subgrid", sub, f)
    return out


def section_totals(section, rho_components, n_r=8, n_theta=32):
    """Exact-quadrature totals over the section (independent of any mesh)."""
    pts, dA, j, v = section.samples(n_r, n_theta)
    C = np.asarray(section.composition, float)
    rho_in = float(np.asarray(rho_components, float) @ C)
    Q = float((j * dA).sum())
    return {"Q": Q, "comp": Q * C, "mass": rho_in * Q, "mom": rho_in * ((j * dA)[:, None] * v).sum(0)}
