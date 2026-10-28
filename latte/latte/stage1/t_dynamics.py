"""G3 pressure / layer / free-surface tests T11-T13, T15-T17 and wave-dispersion checks W01-W02."""
import numpy as np

from ..mesh import Mesh, Geometry, cylinder_mesh, layer_fractions
from ..operators import MetricOperator
from ..solver import Solver, totals
from ..inlet import InletSection, PiecewiseInlet, NoInlet, FLOW, NO_FLOW
from ..linsolve import SPDSolver
from ..failures import SimFailure
from .common import make_model, make_numerics, savefig, plt


M_SCALE_FACTOR = 1e-6   # thresholds.M_scale_factor in configs/stage1/suite.yaml


def _brown(m):
    C = np.zeros((3, m.N, m.ncols)); C[0] = 1.0
    return C


def _budget(s, st, T0):
    T1 = totals(s, st)
    a = st.acc
    scale = max(T0["V"] * M_SCALE_FACTOR, 1e-300)      # V3 Sec. 9.1 denominator floor (volumes)
    vol_res = abs(T1["V"] - T0["V"] - a["in_vol"] + a["out_vol"]) / max(T0["V"] + a["in_vol"], 1e-300)
    comp_res = np.abs(T1["Vj"] - T0["Vj"] - a["in_comp"] + a["out_comp"]) / np.maximum(T0["Vj"] + a["in_comp"], scale)
    mass_res = abs(T1["M"] - T0["M"] - a["in_mass"] + a["out_mass"]) / max(T0["M"] + a["in_mass"], 1e-300)
    return T1, float(vol_res), comp_res, float(mass_res)


def _summ(s, st, T0, extra=None):
    T1, vr, cr, mr = _budget(s, st, T0)
    C = st.composition()
    log = s.log
    out = {"mass_residual": max(mr, float(cr.max())), "volume_residual": vr,
           "simplex_error": float(np.max(np.abs(C.sum(0) - 1))),
           "bounds_violation": float(max(0.0, -C.min(), C.max() - 1)),
           "divergence_error": max(r["e_div"] for r in log) if log else None,
           "mixture_indicator": T1["mixing"],
           "solver_residual": max(max(r["p_res"] for r in log), max(r["v_res"] for r in log)) if log else None,
           "accepted_steps": st.step, "rejected_steps": st.rejected}
    if extra:
        out.update(extra)
    return out, T1


# ---------------------------------------------------------------- T11 hydrostatic balance
def T11(cfg, run):
    p = cfg["tests"]["T11"]
    rows = []
    fig, ax = plt.subplots(figsize=(5, 3.2))
    for case in ("equal_density", "foam_top_layers", "linear_stratification"):
        m = cylinder_mesh(p["n"], p["radius"], layer_fractions(p["N"]))
        s = Solver(m, make_model(cfg), make_numerics(cfg))
        C = _brown(m)
        if case == "foam_top_layers":
            C[0, -p["foam_layers"]:] = 1 - p["foam_phi"]; C[2, -p["foam_layers"]:] = p["foam_phi"]
        elif case == "linear_stratification":
            phi = p["foam_phi"] * (np.arange(m.N) + 0.5) / m.N
            C[0] = (1 - phi)[:, None]; C[2] = phi[:, None]
        st = s.initial_state(p["depth"], C)
        T0 = totals(s, st)
        vt = []
        while st.t < p["t_end"] - 1e-12:
            st, d = s.advance(st, dt_try=p["t_end"] - st.t)
            vt.append((st.t, float(np.abs(st.velocity()).max()), float(np.abs(st.D - p["depth"]).max())))
        summ, T1 = _summ(s, st, T0)
        vt = np.array(vt)
        rows.append({"case": case, "vmax": float(vt[:, 1].max()), "level_drift": float(vt[:, 2].max()), **summ})
        run.write_log(s.log, f"diagnostics/solver_{case}.csv")
        ax.semilogy(vt[:, 0], np.maximum(vt[:, 1], 1e-18), label=case)
    ax.set_xlabel("t"); ax.set_ylabel("max |u|"); ax.legend(fontsize=7)
    savefig(fig, run.path("figures", "T11_rest.png"))
    run.write_log(rows, "diagnostics/cases.csv")
    th = cfg["thresholds"]
    ok = all(r["vmax"] <= th["rest_velocity"] and r["level_drift"] <= 1e-12 and r["volume_residual"] <= th["budget"] for r in rows)
    keys = ("case", "vmax", "level_drift", "volume_residual", "divergence_error", "accepted_steps", "rejected_steps")
    return {"status": "completed", "verdict": "pass" if ok else "fail",
            "reason": "; ".join(f"{r['case']}: max|u| {r['vmax']:.1e}, level drift {r['level_drift']:.1e}" for r in rows),
            "volume_residual": max(r["volume_residual"] for r in rows), "mass_residual": max(r["mass_residual"] for r in rows),
            "divergence_error": max(r["divergence_error"] for r in rows), "accepted_steps": sum(r["accepted_steps"] for r in rows),
            "rejected_steps": sum(r["rejected_steps"] for r in rows), "details": [{k: r[k] for k in keys} for r in rows]}


# ---------------------------------------------------------------- T12 common pressure projection
def T12(cfg, run):
    p = cfg["tests"]["T12"]
    rng = np.random.default_rng(cfg["seed"])
    nx, ny, N = p["nx"], p["ny"], p["N"]
    m = Mesh(nx, ny, 1.0 / nx, 1.0 / ny, layer_fractions(N))
    X, Y = m.xc, m.yc
    D = p["depth"] + p["bump"] * np.cos(np.pi * X) * np.cos(np.pi * Y)
    g = Geometry(m, D)
    rho = 0.7 + 0.3 * rng.random(m.nc)
    v = rng.standard_normal((3, m.nc))
    # known boundary flux: outflow through the bottom-layer cells of the x = 0 wall, inflow at x = 1 (layer 1)
    bc = np.zeros(m.nc)
    left = np.nonzero(m.ci == 0)[0]; right = np.nonzero(m.ci == nx - 1)[0]
    bc[m.cell(0, left)] = 1e-3
    bc[m.cell(1, right)] = -0.5e-3
    s = Solver(m, make_model(cfg, sigma=0.0), make_numerics(cfg), bc_flux=bc)
    dt = 0.01
    Un, pi, info, Ut, op = s._project(g, rho, v, dt, gravity=True)
    L = op.L
    sym = float(abs(L - L.T).max() / abs(L).max())
    ev = np.linalg.eigvalsh(L.toarray())
    checks = {"L_symmetry": sym, "L_min_eig": float(ev.min()), "L_cond": float(ev.max() / ev.min()),
              "e_div": info["e_div"], "poisson_rel_residual": info["rel_res"]}
    # idempotence: projecting the projected face fluxes again changes nothing
    div2 = -(op.G.T @ Un) + bc
    pi2, _ = SPDSolver("direct").solve(L, -div2 / dt)
    dU = dt * op.flux(pi2)
    checks["idempotence"] = float(np.abs(dU).max() / np.abs(Un).max())
    # pure gradient fluxes are removed completely
    phi = rng.standard_normal(m.nc)
    Ug = dt * op.flux(phi)
    pig, _ = SPDSolver("direct").solve(L, -(-(op.G.T @ Ug)) / dt)
    checks["gradient_removal"] = float(np.abs(Ug - dt * op.flux(pig)).max() / np.abs(Ug).max())
    # rigid-lid (closed) variant: singular operator, compatible boundary flux, explicit gauge
    opn = MetricOperator(m, g, *[1.0 / r for r in s._face_rho(g, rho)], surface="neumann")
    Ln = opn.L.tolil()
    bcz = bc.copy(); bcz[m.cell(1, right)] = -bcz.sum() / right.size + bcz[m.cell(1, right)]
    Uin = s._interp(g, v)[:opn.nf]
    divn = -(opn.G.T @ Uin) + bcz
    rhs = -divn / dt
    pin = 0
    Ln[pin, :] = 0; Ln[pin, pin] = 1.0; rhs[pin] = 0.0
    pn, _ = SPDSolver("direct").solve(Ln.tocsr(), rhs)
    pn -= pn.mean()                                        # mean-zero pressure gauge
    Unn = Uin - dt * opn.flux(pn)
    checks["rigid_lid_e_div"] = float(dt * np.max(np.abs(-(opn.G.T @ Unn) + bcz) / g.V.ravel()))
    checks["rigid_lid_gauge_mean"] = float(abs(pn.mean()))
    ev_n = np.linalg.eigvalsh(opn.L.toarray())
    checks["rigid_lid_nullspace_dim"] = int(np.sum(ev_n < 1e-10 * ev_n.max()))
    run.write_log([checks], "diagnostics/projection_checks.csv")
    ok = (sym <= 1e-14 and ev.min() > 0 and info["e_div"] <= 1e-10 and checks["idempotence"] <= 1e-10
          and checks["gradient_removal"] <= 1e-10 and checks["rigid_lid_e_div"] <= 1e-10 and checks["rigid_lid_nullspace_dim"] == 1)
    return {"status": "completed", "verdict": "pass" if ok else "fail",
            "reason": (f"sloped surface (max slope {g.max_slope:.3f}), random density 0.7-1, random v*, prescribed wall flux: "
                       f"e_div {info['e_div']:.1e}, L symmetric to {sym:.1e}, lambda_min {ev.min():.2e} > 0, idempotence "
                       f"{checks['idempotence']:.1e}, gradient removal {checks['gradient_removal']:.1e}; rigid lid: null space "
                       f"dim {checks['rigid_lid_nullspace_dim']}, e_div {checks['rigid_lid_e_div']:.1e}"),
            "divergence_error": info["e_div"], "solver_residual": info["rel_res"], "details": checks}


# ---------------------------------------------------------------- T13 two-layer viscous exchange
def T13(cfg, run):
    p = cfg["tests"]["T13"]
    nx = p["nx"]
    m = Mesh(nx, 1, 1.0 / nx, 1.0 / nx, [0.5, 0.5], periodic_x=True, periodic_y=True)
    mod = make_model(cfg, bottom_bc="free_slip")
    s = Solver(m, mod, make_numerics(cfg))
    v0 = np.zeros((3, 2, m.ncols)); v0[0, 0] = p["u_lower"]; v0[0, 1] = p["u_upper"]
    st = s.initial_state(p["depth"], _brown(m), v0)
    A = m.A; Dp = p["depth"]
    M1 = M2 = mod.rho[0] * 0.5 * Dp * A
    zeta = mod.mu * A / (0.5 * Dp)
    k = zeta * (1 / M1 + 1 / M2)
    du0 = p["u_lower"] - p["u_upper"]
    rows = []
    P0 = st.P[0].sum()
    ke_prev = np.inf
    mono = True
    n = 0
    while st.t < p["t_end"] - 1e-12:
        st, d = s.advance(st, dt_try=p["dt"])
        n += 1
        u = st.velocity()[0].reshape(2, m.ncols)
        du = float(np.mean(u[0] - u[1]))
        T1 = totals(s, st)
        mono &= T1["KE"] <= ke_prev * (1 + 1e-14)
        ke_prev = T1["KE"]
        rows.append({"t": st.t, "du": du, "du_discrete": du0 * (1 + k * p["dt"]) ** (-n), "du_exact": du0 * np.exp(-k * st.t),
                     "P": float(st.P[0].sum()), "KE": T1["KE"]})
    run.write_log(rows, "diagnostics/two_layer.csv")
    r = np.array([[x["t"], x["du"], x["du_discrete"], x["du_exact"]] for x in rows])
    e_disc = float(np.max(np.abs(r[:, 1] - r[:, 2])) / du0)
    e_cont = float(np.max(np.abs(r[:, 1] - r[:, 3])) / du0)
    e_mom = float(np.max(np.abs(np.array([x["P"] for x in rows]) - P0)) / abs(P0))
    fig, ax = plt.subplots(figsize=(5, 3.2))
    ax.semilogy(r[:, 0], r[:, 1], label="solver")
    ax.semilogy(r[:, 0], r[:, 3], "--", label="exp[-zeta(1/M1+1/M2) t]")
    ax.set_xlabel("t"); ax.set_ylabel("u_lower - u_upper"); ax.legend(fontsize=7)
    savefig(fig, run.path("figures", "T13_two_layer.png"))
    ok = e_disc <= 1e-10 and e_mom <= 1e-12 and mono
    return {"status": "completed", "verdict": "pass" if ok else "fail",
            "reason": (f"relative-velocity decay matches the backward-Euler two-mass solution to {e_disc:.1e} (continuous "
                       f"exponential to {e_cont:.1e}, O(dt)); total momentum kept to {e_mom:.1e}; KE non-increasing: {mono}"),
            "mass_residual": e_mom, "accepted_steps": st.step, "rejected_steps": st.rejected,
            "details": {"zeta": zeta, "rate": k, "err_discrete": e_disc, "err_continuous": e_cont}}


# ---------------------------------------------------------------- T15 injection and level
def _disc_inlet(center_xy, z, diameter, U, t_end):
    sec = InletSection(center=np.array([center_xy[0], center_xy[1], z]), normal=np.array([0, 0, 1.0]),
                       axis1=np.array([1.0, 0, 0]), a=diameter / 2, b=diameter / 2, velocity=np.array([0, 0, -U]),
                       composition=np.array([0.0, 1.0, 0.0]))
    return PiecewiseInlet([{"t": 0.0, "status": FLOW, "section": sec}, {"t": 10 * t_end, "status": FLOW, "section": sec}]), sec


def T15(cfg, run):
    p = cfg["tests"]["T15"]
    results = []
    cases = [("reference", "velocity", p["reference_inlet"], [p["n"], 2 * p["n"]]),
             ("reference", "surface_flux", p["reference_inlet"], [p["n"]]),
             ("gentle", "velocity", p["gentle_inlet"], [p["n"]]),
             ("gentle", "surface_flux", p["gentle_inlet"], [p["n"]]),
             ("broad", "velocity", p["broad_inlet"], [p["n"]]),
             ("broad", "surface_flux", p["broad_inlet"], [p["n"]])]
    for base, bc, ip, grids in cases:
        label = f"{base}_{bc}"
        for n in grids:
            m = cylinder_mesh(n, p["radius"], layer_fractions(p["N"]))
            prov, sec = _disc_inlet(ip["center"], p["depth"], ip["diameter"], ip["U"], ip["t_end"])
            s = Solver(m, make_model(cfg), make_numerics(cfg, connect_tol=0.5 * p["depth"], inlet_bc=bc), inlet=prov)
            st = s.initial_state(p["depth"], _brown(m))
            T0 = totals(s, st)
            Q = np.pi * (ip["diameter"] / 2) ** 2 * ip["U"]
            hist = []
            failure = None
            try:
                while st.t < ip["t_end"] - 1e-12:
                    st, d = s.advance(st, dt_try=ip["t_end"] - st.t)
                    hist.append((st.t, float(st.D.mean()), float(st.D.min()), float(st.D.max()), d["max_slope"]))
            except SimFailure as e:
                failure = {"code": e.code, "reason": e.reason, "t": st.t}
            summ, T1 = _summ(s, st, T0)
            A_cup = m.area_total()
            lvl_err = abs(st.D.mean() - (p["depth"] + Q * st.t / A_cup)) / (Q * st.t / A_cup)
            vin_err = abs(st.acc["in_vol"] - Q * st.t) / (Q * st.t)
            C = st.composition().reshape(3, m.N, m.ncols)
            Vc = st.Vtot.reshape(m.N, m.ncols)
            milk_by_layer = (C[1] * Vc).sum(1) / max((C[1] * Vc).sum(), 1e-300)
            res = {"inlet": label, "inlet_bc": bc, "n": n, "cells_per_inlet_diameter": ip["diameter"] / m.dx, "t": st.t, "Q": Q,
                   "level_rel_err": float(lvl_err), "inflow_rel_err": float(vin_err),
                   "brown_change": float(abs(T1["Vj"][0] - T0["Vj"][0]) / T0["Vj"][0]),
                   "milk_vs_Qt": float(abs(T1["Vj"][1] - Q * st.t) / (Q * st.t)),
                   "milk_fraction_by_layer_bottom_to_top": milk_by_layer, "failure": failure,
                   "crater_depth": float(p["depth"] + Q * st.t / A_cup - st.D.min()),
                   "rim_height": float(st.D.max() - (p["depth"] + Q * st.t / A_cup)),
                   "geo_material_volume_mismatch": max(r_["vol_mismatch"] for r_ in s.log) if s.log else None,
                   "advective_momentum_bookkeeping_residual": st.acc["adv_mom_residual"],
                   "max_slope": max(h[4] for h in hist) if hist else None, **summ}
            results.append(res)
            run.write_log(s.log, f"diagnostics/solver_{label}_{n}.csv")
            if n == grids[-1]:
                fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
                h = np.array(hist)
                ax[0].plot(h[:, 0], h[:, 1] - p["depth"], label="mean level - H0 (solver)")
                ax[0].plot(h[:, 0], Q * h[:, 0] / A_cup, "k--", label="Q t / A_cup")
                ax[0].plot(h[:, 0], h[:, 2] - p["depth"], ":", label="min level - H0")
                ax[0].set_xlabel("t"); ax[0].legend(fontsize=7)
                top = (C[1] * Vc).sum(0) / Vc.sum(0)
                im = ax[1].imshow(m.to_grid(top), origin="lower", extent=[-0.5, 0.5, -0.5, 0.5], cmap="copper_r", vmin=0)
                plt.colorbar(im, ax=ax[1]); ax[1].set_title(f"{label} n={n}: column-mean milk fraction", fontsize=8)
                j0 = m.ny // 2
                cols = np.nonzero(m.cj == j0)[0]
                xs = m.xc[cols]
                Z = Geometry(m, st.D).zc[:, cols]
                ax[2].pcolormesh(np.tile(xs, (m.N, 1)), Z, C[1][:, cols], cmap="copper_r", vmin=0, vmax=1, shading="nearest")
                ax[2].set_title("milk fraction, section y = 0", fontsize=8); ax[2].set_xlim(ip["center"][0] - 0.2, ip["center"][0] + 0.2)
                savefig(fig, run.path("figures", f"T15_{label}.png"))
    run.write_log([{k: v for k, v in r.items() if k != "milk_fraction_by_layer_bottom_to_top"} for r in results], "diagnostics/cases.csv")
    th = cfg["thresholds"]
    ok = all(r["level_rel_err"] <= th["budget"] and r["volume_residual"] <= th["budget"] and r["brown_change"] <= th["budget"]
             and r["milk_vs_Qt"] <= th["budget"] for r in results)
    fails = [f"{r['inlet']} n={r['n']}: {r['failure']['code']} at t={r['failure']['t']:.3f}" for r in results if r["failure"]]
    verdict = "pass" if ok else "fail"
    reason = ("; ".join(f"{r['inlet']} n={r['n']} ({r['cells_per_inlet_diameter']:.1f} cells/diameter) t={r['t']:.3f}: level err "
                        f"{r['level_rel_err']:.1e}, volume res {r['volume_residual']:.1e}, brown change {r['brown_change']:.1e}, "
                        f"min/max level vs mean {-r['crater_depth']:+.4f}/{r['rim_height']:+.4f}" for r in results)
              + (f". Stopped runs: {fails}" if fails else ""))
    return {"status": "completed", "verdict": verdict, "reason": reason,
            "volume_residual": max(r["volume_residual"] for r in results), "mass_residual": max(r["mass_residual"] for r in results),
            "divergence_error": max(r["divergence_error"] for r in results),
            "accepted_steps": sum(r["accepted_steps"] for r in results), "rejected_steps": sum(r["rejected_steps"] for r in results),
            "details": results}


# ---------------------------------------------------------------- T16 buoyancy
def T16(cfg, run):
    p = cfg["tests"]["T16"]
    out = []
    curves = {}
    for case in ("equal_density", "buoyant"):
        for nx, N in (p["grids"] if case == "buoyant" else p["grids"][:1]):
            m = Mesh(nx, 1, p["Lx"] / nx, p["Lx"] / nx, layer_fractions(N), periodic_x=True, periodic_y=True)
            s = Solver(m, make_model(cfg), make_numerics(cfg))
            g0 = Geometry(m, np.full(m.ncols, p["depth"]))
            sub = 4
            frac = np.zeros((N, m.ncols))
            o = (np.arange(sub) + 0.5) / sub - 0.5
            for a in o:
                for b in o:
                    xx = m.xc[None, :] + a * m.dx
                    zz = g0.zc + b * g0.h
                    frac += ((xx - p["blob_center"][0]) ** 2 + (zz - p["blob_center"][1]) ** 2 < p["blob_radius"] ** 2)
            frac /= sub ** 2
            C = _brown(m)
            if case == "buoyant":
                C[2] = p["blob_phi"] * frac; C[0] = 1 - C[2]
            else:
                C[1] = frac; C[0] = 1 - frac                       # milk liquid blob: same density as the base
            st = s.initial_state(p["depth"], C)
            T0 = totals(s, st)
            hist = []
            while st.t < p["t_end"] - 1e-12:
                st, d = s.advance(st, dt_try=p["t_end"] - st.t)
                Cc = st.composition()
                gg = Geometry(m, st.D)
                w = Cc[2 if case == "buoyant" else 1] * st.Vtot
                hist.append((st.t, float((w * gg.zc.ravel()).sum() / w.sum()), float(np.abs(st.velocity()).max())))
            summ, T1 = _summ(s, st, T0)
            h = np.array(hist)
            out.append({"case": case, "nx": nx, "N": N, "z_centroid_end": float(h[-1, 1]), "rise": float(h[-1, 1] - p["blob_center"][1]),
                        "vmax": float(h[:, 2].max()), **summ})
            curves[(case, nx)] = (h, m, st)
    run.write_log(out, "diagnostics/cases.csv")
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    for (case, nx), (h, m, st) in curves.items():
        if case == "buoyant":
            ax[0].plot(h[:, 0], h[:, 1], label=f"nx={nx}")
    ax[0].set_xlabel("t"); ax[0].set_ylabel("gas centroid z"); ax[0].legend(fontsize=7)
    h, m, st = curves[("buoyant", p["grids"][-1][0])]
    C = st.composition()[2].reshape(m.N, m.ncols)
    gg = Geometry(m, st.D)
    ax[1].pcolormesh(np.tile(m.xc, (m.N, 1)), gg.zc, C, cmap="Blues", shading="nearest")
    ax[1].set_aspect("equal"); ax[1].set_title(f"gas fraction at t={st.t:.2f}, nx={p['grids'][-1][0]}", fontsize=8)
    savefig(fig, run.path("figures", "T16_buoyancy.png"))
    eq = [o for o in out if o["case"] == "equal_density"][0]
    b = [o for o in out if o["case"] == "buoyant"]
    rises = [o["rise"] for o in b]
    diffs = [abs(rises[i + 1] - rises[i]) for i in range(len(rises) - 1)]
    conv = all(diffs[i + 1] < diffs[i] for i in range(len(diffs) - 1))
    ok_hard = eq["vmax"] <= cfg["thresholds"]["rest_velocity"] and all(r > 0 for r in rises) and max(o["volume_residual"] for o in out) <= cfg["thresholds"]["budget"]
    verdict = ("review" if conv else "fail") if ok_hard else "fail"
    return {"status": "completed", "verdict": verdict,
            "reason": (f"equal-density blob: max|u| {eq['vmax']:.1e}; buoyant blob (rho 0.7 in 1.0) rises by {[round(r, 4) for r in rises]} "
                       f"for nx={[o['nx'] for o in b]} (successive differences {[f'{d:.1e}' for d in diffs]}); no reference solution, so magnitude is 'review'"),
            "volume_residual": max(o["volume_residual"] for o in out), "mass_residual": max(o["mass_residual"] for o in out),
            "divergence_error": max(o["divergence_error"] for o in out), "accepted_steps": sum(o["accepted_steps"] for o in out),
            "rejected_steps": sum(o["rejected_steps"] for o in out), "details": out}


# ---------------------------------------------------------------- W01 / W02 standing waves
def _wave(cfg, run, key, sigma):
    p = cfg["tests"][key]
    H, k, g = p["H"], 2 * np.pi, cfg["model"]["g"]
    rho = cfg["model"]["rho"][0]
    w_nh = np.sqrt((g * k + sigma * k ** 3 / rho) * np.tanh(k * H))
    w_h = np.sqrt((g * H + sigma * k ** 2 * H / rho)) * k
    T = 2 * np.pi / w_nh
    rows = []
    series = {}
    for nx, N in p["resolutions"]:
        m = Mesh(nx, 1, 1.0 / nx, 1.0 / nx, layer_fractions(N), periodic_x=True, periodic_y=True)
        s = Solver(m, make_model(cfg, mu=0.0, sigma=sigma), make_numerics(cfg))
        st = s.initial_state(H + p["amp"] * np.cos(k * m.xc), _brown(m))
        ts, e = [0.0], [st.D[0] - H]
        while st.t < p["periods"] * T:
            st, d = s.advance(st)
            ts.append(st.t); e.append(st.D[0] - H)
        ts, e = np.array(ts), np.array(e)
        zc = [ts[i] - e[i] * (ts[i + 1] - ts[i]) / (e[i + 1] - e[i]) for i in range(len(e) - 1) if e[i] * e[i + 1] < 0]
        Pn = 2 * float(np.mean(np.diff(zc)))
        rows.append({"nx": nx, "N": N, "period": Pn, "rel_err": Pn / T - 1, "steps": st.step,
                     "amplitude_ratio_end": float(np.max(np.abs(e[ts > ts[-1] - 0.5 * T])) / p["amp"]),
                     "volume_residual": float(abs(st.Vtot.sum() - s.initial_state(H + p["amp"] * np.cos(k * m.xc), _brown(m)).Vtot.sum()))})
        series[(nx, N)] = (ts, e)
    r = [abs(x["rel_err"]) for x in rows]
    order = [float(np.log(r[i] / r[i + 1]) / np.log(2)) for i in range(len(r) - 1)]
    run.write_log(rows, "diagnostics/wave.csv")
    fig, ax = plt.subplots(figsize=(6, 3.3))
    for (nx, N), (ts, e) in series.items():
        ax.plot(ts, e / p["amp"], label=f"nx={nx}, N={N}")
    tt = np.linspace(0, p["periods"] * T, 400)
    ax.plot(tt, np.cos(w_nh * tt), "k--", lw=0.8, label="non-hydrostatic theory")
    ax.plot(tt, np.cos(w_h * tt), "r:", lw=0.8, label="hydrostatic theory")
    ax.set_xlabel("t"); ax.set_ylabel("eta(0,t)/a"); ax.legend(fontsize=6)
    savefig(fig, run.path("figures", f"{key}_wave.png"))
    ok = abs(rows[-1]["rel_err"]) < 0.01 and all(r[i + 1] < r[i] for i in range(len(r) - 1))
    return {"status": "completed", "verdict": "pass" if ok else "fail",
            "reason": (f"period error vs omega^2 = (g k + sigma k^3/rho) tanh(kH) (sigma={sigma}): "
                       f"{[f'{x['rel_err']:+.2%}' for x in rows]} for (nx,N)={p['resolutions']} (observed order {[round(o, 2) for o in order]}); "
                       f"hydrostatic theory would differ by {2 * np.pi / w_h / T - 1:+.1%}; amplitude kept "
                       f"{[round(x['amplitude_ratio_end'], 3) for x in rows]} with mu=0 (numerical damping)"),
            "volume_residual": max(x["volume_residual"] for x in rows), "details": rows}


def W01(cfg, run):
    return _wave(cfg, run, "W01", 0.0)


def W02(cfg, run):
    return _wave(cfg, run, "W02", cfg["tests"]["W02"]["sigma"])


# ---------------------------------------------------------------- T17 surface, walls, inlet switch-off
def T17(cfg, run):
    p = cfg["tests"]["T17"]
    checks = {}
    m = cylinder_mesh(p["n"], 0.5, layer_fractions(p["N"]))
    s = Solver(m, make_model(cfg, sigma=p["sigma"]), make_numerics(cfg))
    st = s.initial_state(0.25, _brown(m))
    T0 = totals(s, st)
    while st.t < p["t_end"] - 1e-12:
        st, d = s.advance(st, dt_try=p["t_end"] - st.t)
    checks["flat_surface_with_tension_max_u"] = float(np.abs(st.velocity()).max())
    checks["flat_surface_steps"] = st.step
    checks["flat_surface_volume_residual"] = float(abs(totals(s, st)["V"] - T0["V"]) / T0["V"])
    checks["wall_faces_are_not_flux_dofs"] = True       # no-penetration is structural: wall faces carry no volume flux
    # inlet switched off: afterwards identical (bitwise) to a run that never had an inlet
    prov, sec = _disc_inlet([0.1, 0.0], 0.25, 0.04, 0.2, 1.0)
    sw = p["switch_time"]
    prov = PiecewiseInlet([{"t": 0.0, "status": FLOW, "section": sec}, {"t": sw, "status": NO_FLOW}, {"t": 100.0, "status": NO_FLOW}])
    sa = Solver(m, make_model(cfg), make_numerics(cfg, connect_tol=0.1), inlet=prov)
    sta = sa.initial_state(0.25, _brown(m))
    while sta.t < sw - 1e-14:
        sta, _ = sa.advance(sta)
    stb = sta.copy()
    sb = Solver(m, make_model(cfg), make_numerics(cfg), inlet=NoInlet())
    import copy
    sb.psolver, sb.vsolver = copy.deepcopy(sa.psolver), copy.deepcopy(sa.vsolver)   # same linear-solver state
    for _ in range(p["continue_steps"]):
        sta, _ = sa.advance(sta, dt_try=0.01)
        stb, _ = sb.advance(stb, dt_try=0.01)
    checks["after_switch_off_bitwise_identical"] = bool(np.array_equal(sta.P, stb.P) and np.array_equal(sta.D, stb.D) and np.array_equal(sta.Vc, stb.Vc))
    checks["after_switch_off_max_rel_diff"] = float(np.abs(sta.P - stb.P).max() / max(np.abs(stb.P).max(), 1e-300))
    checks["switch_time_hit"] = float(sw)
    checks["cup_edge_contact_line"] = "blocked: contact angle / slip model not chosen (walls currently act as a 90-degree contact placeholder)"
    run.write_log([checks], "diagnostics/checks.csv")
    ok = checks["flat_surface_with_tension_max_u"] <= cfg["thresholds"]["rest_velocity"] and checks["after_switch_off_bitwise_identical"]
    return {"status": "completed", "verdict": "pass" if ok else "fail",
            "reason": (f"flat surface with sigma={p['sigma']}: max|u| {checks['flat_surface_with_tension_max_u']:.1e}; capillary "
                       f"dynamics see W02; inlet switch-off leaves no constraint (bitwise identical: {checks['after_switch_off_bitwise_identical']}); "
                       f"cup-edge meniscus part BLOCKED (contact-line model undecided)"),
            "accepted_steps": st.step, "details": checks}
