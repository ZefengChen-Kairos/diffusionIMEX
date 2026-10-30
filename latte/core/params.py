"""Shared parameter sets.

PhysParams is the single physics parameter set used by all six patterns and all gates.
Gates and pattern scripts construct it with no arguments; only counterfactual / sensitivity
runs override fields, and they record the override.
"""
from dataclasses import dataclass, asdict, replace
import hashlib
import json
import os


@dataclass(frozen=True)
class PhysParams:
    # fixed material properties
    rho_c: float = 970.0        # coffee bulk density [kg/m^3]
    rho_m: float = 550.0        # microfoam effective density (with air) [kg/m^3]
    rho_milk: float = 1000.0    # liquid milk density of the incoming jet [kg/m^3]
    sigma: float = 0.04         # surface tension [N/m]
    g: float = 9.81
    # geometry / inputs
    R_cup: float = 0.042
    h0: float = 0.004           # initial active layer thickness (pure coffee)
    A_spout: float = 2e-5
    V_bulk0: float = 180e-6     # base liquid in the cup when the pattern starts [m^3]
    bulk_return: str = "white"  # R returns the bulk's white (whipped) material share ("white")
                                # or its total milk share ("all_milk", task-book form)
    # calibrated physical closures
    gamma: float = 7.0          # drag on the bulk [1/s]
    tau_y: float = 0.0          # foam yield stress on the base [Pa] (0 = spec model)
    tau_ex: float = 0.0         # layer-thickness exchange time with the bulk [s] (0 = off)
    chi: float = 0.35           # momentum retention of the deposited milk
    k_h: float = 0.30           # horizontal fraction of the jet speed
    nu_s: float = 2e-5          # kinematic surface viscosity [m^2/s]; T_s = 2 nu_s M D_h
    # empirical closures
    retention: str = "sigmoid"  # "pi" or "sigmoid"
    Pi0: float = 207.0          # only for retention == "pi"
    phi_f: float = 0.14         # foam yield: share of the non-penetrating pour that surfaces as D
                                # (rest drains quietly to the bulk); calibrated from white-area growth
    v_e: float = 0.80
    dv: float = 0.10
    C_E: float = 0.3            # task-book start 0.5 (calibrated on the pull-through, REPORT F5)
    C_t: float = 0.3
    v_c: float = 0.025
    sig_d_min: float = 1.1e-3
    sig_d_fac: float = 0.75
    r0_fac: float = 0.6
    sig_k_min: float = 2.0e-3   # task-book start 3.5 mm (a thin high-pour stream, D_j ~ 2 mm)
    sig_r_fac: float = 2.5
    # optics
    Phi_star: float = 1.2       # white scale [kg/m^2]
    pressure: str = "elastic"   # "buoyant" (spec main, eq. D) or "elastic" (spec fallback)
    K_s: float = 970.0 * 9.81   # hydrostatic stiffness rho_c g (pressure == "elastic")

    def to_dict(self):
        return asdict(self)

    def digest(self):
        return hashlib.sha1(json.dumps(self.to_dict(), sort_keys=True).encode()).hexdigest()[:12]

    def with_(self, **kw):
        return replace(self, **kw)


@dataclass(frozen=True)
class NumParams:
    L: float = 0.096            # square domain side
    N: int = int(os.environ.get("LATTE_N", 128))   # M/P grid (LATTE_N only for smoke tests)
    ratio: int = 3              # fine material grid = ratio*N
    cfl: float = 0.35
    dt_max: float = 2e-3
    mat_every: int = 4          # material (fine grid) step every k coarse steps
    kappa_fac: float = 0.25     # interface dissipation kappa = kappa_fac * c' * dx
    limiter: str = "mc"         # MUSCL limiter for the coarse fluxes
    n_tags: int = 6
    check_every: int = 1        # conservation assertion cadence (steps)
    tol: float = 1e-6

    def to_dict(self):
        return asdict(self)

    def with_(self, **kw):
        return replace(self, **kw)
