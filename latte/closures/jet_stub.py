"""Analytic jet stub: pour sample -> impact event (section 1.4).

The solver only sees an ImpactEvent; phase two replaces this module by MPM boundary
crossing events with the same fields.  Height enters only through v_j (CF2).
"""
from dataclasses import dataclass
import math
import numpy as np


@dataclass
class ImpactEvent:
    x: np.ndarray        # impact point [m]
    mdot: float          # incoming mass rate [kg/s]
    Q: float             # volume rate [m^3/s]
    v_j: float           # impact speed [m/s]
    D_j: float           # jet diameter at impact [m]
    r: float             # non-penetrating (retained) share; 1 - r penetrates (drives E, tau_j)
    phi_f: float         # foam yield of the retained share: D = phi_f r mdot
    u_d: np.ndarray      # horizontal velocity of the retained milk [m/s]
    v_h: np.ndarray      # horizontal velocity of the penetrating jet [m/s]
    eta: float           # fast-move narrowing factor
    tag: int


def jet_speed(Q, h, P):
    v0 = Q / P.A_spout
    return math.sqrt(v0 * v0 + 2.0 * P.g * h)


def jet_diameter(Q, v_j):
    return math.sqrt(4.0 * Q / (math.pi * v_j)) if Q > 0 else 0.0


def retention(v_j, D_j, P, rho_bar=None):
    if P.retention == "sigmoid":
        z = (v_j - P.v_e) / P.dv
        return 1.0 - 1.0 / (1.0 + math.exp(-z))
    rho_bar = P.rho_c if rho_bar is None else rho_bar
    g_dep = P.g * (P.rho_c - P.rho_m) / P.rho_c
    Pi = v_j ** 2 / (g_dep * D_j + P.sigma / (rho_bar * D_j))
    return 1.0 / (1.0 + Pi / P.Pi0)


def impact_event(sample, P, v_j_override=None):
    """sample: actions.scripts.Sample.  v_j_override is used by CF2 only."""
    if sample.Q <= 0.0:
        return None
    v_j = jet_speed(sample.Q, sample.h, P) if v_j_override is None else v_j_override
    D_j = jet_diameter(sample.Q, v_j)
    r = retention(v_j, D_j, P)
    e = np.asarray(sample.e, float)
    v_p = np.asarray(sample.v, float)
    v_h = P.k_h * v_j * e + v_p
    u_d = P.chi * v_h
    eta = 1.0 / (1.0 + float(np.linalg.norm(v_p)) / P.v_c)
    return ImpactEvent(np.asarray(sample.x, float), P.rho_milk * sample.Q, sample.Q, v_j, D_j,
                       r, P.phi_f, u_d, v_h, eta, sample.tag)


def footprint(ev, P):
    """Kernel geometry (all lengths in m)."""
    r0 = P.r0_fac * ev.D_j * ev.r * ev.eta
    sd = max(P.sig_d_min, P.sig_d_fac * ev.D_j * ev.eta)
    sk = max(ev.D_j, P.sig_k_min)
    sr = P.sig_r_fac * sk
    return r0, sd, sk, sr
