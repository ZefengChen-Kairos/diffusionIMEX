"""G2 transport tests T05-T10 and the interlayer exchange test T14."""
import numpy as np

from ..mesh import Mesh, Geometry, layer_fractions
from ..solver import Solver, totals
from ..remap import remap_column
from .common import Prescribed, make_model, make_numerics, l1, rate, savefig, plt, subcell_disc


def _mixing(C, V):
    liq = C[0] + C[1]
    xM = np.where(liq > 0, C[1] / liq, 0.0)
    return float((liq * V * xM * (1 - xM)).sum())


# ---------------------------------------------------------------- T05 static stripes (full solver)
def T05(cfg, run):
    p = cfg["tests"]["T05"]
    n, N = p["n"], p["N"]
    m = Mesh(n, n, 1.0 / n, 1.0 / n, layer_fractions(N), periodic_x=True, periodic_y=True)
    s = Solver(m, make_model(cfg), make_numerics(cfg))
    band = (np.floor(m.xc * 2 * p["stripes"]) % 2 == 0).astype(float)      # sharp stripes of fixed width
    C = np.zeros((3, N, m.ncols)); C[1] = band[None, :]; C[0] = 1 - C[1]
    st = s.initial_state(p["depth"], C)
    C0 = st.composition().copy()
    T0 = totals(s, st)
    while st.t < p["t_end"] - 1e-12:
        st, d = s.advance(st, dt_try=p["t_end"] - st.t)
    dC = float(np.max(np.abs(st.composition() - C0)))
    T1 = totals(s, st)
    run.write_log(s.log)
    inv = float(np.max(np.abs(T1["Vj"] - T0["Vj"])) / T0["V"])
    ok = dC <= cfg["thresholds"]["local_update"] and inv <= cfg["thresholds"]["budget"]
    return {"status": "completed", "verdict": "pass" if ok else "fail",
            "reason": f"max|C(t)-C(0)| = {dC:.2e}, max|u| = {T1['vmax']:.2e} after {st.step} steps",
            "volume_residual": inv, "bounds_violation": 0.0, "divergence_error": max(r["e_div"] for r in s.log),
            "mixture_indicator": [T0["mixing"], T1["mixing"]], "solver_residual": max(r["p_res"] for r in s.log),
            "accepted_steps": st.step, "rejected_steps": st.rejected}


# ---------------------------------------------------------------- prescribed flux helpers
def _uniform_flux(m, g, u, v):
    return np.concatenate([(g.afx * u).ravel(), (g.afy * v).ravel(), np.zeros((m.N - 1) * m.ncols)])


def _psi_flux(m, g, psi):
    """Discretely divergence-free fluxes from a streamfunction psi(x, y) at cell corners (one layer set)."""
    i = m.ci[m.xp].astype(float); j = m.cj[m.xp].astype(float)
    xf = m.x0 + i * m.dx
    if m.periodic_x:
        xf = np.where(m.ci[m.xp] == 0, m.x0, xf)   # wrap face at the domain edge
    yb, yt = m.y0 + j * m.dy, m.y0 + (j + 1) * m.dy
    Fx = psi(xf, yt) - psi(xf, yb)
    i2 = m.ci[m.yp].astype(float); j2 = m.cj[m.yp].astype(float)
    yf = m.y0 + j2 * m.dy
    xl, xr = m.x0 + i2 * m.dx, m.x0 + (i2 + 1) * m.dx
    Fy = -(psi(xr, yf) - psi(xl, yf))
    hx = g.hfx; hy = g.hfy
    return np.concatenate([(hx * Fx[None, :]).ravel(), (hy * Fy[None, :]).ravel(), np.zeros((m.N - 1) * m.ncols)])


# ---------------------------------------------------------------- T06 periodic translation
def T06(cfg, run):
    p = cfg["tests"]["T06"]
    rho = cfg["model"]["rho"]
    u, v = p["velocity"]
    res = {"disc": [], "bell": []}
    figs = {}
    for n in p["grids"]:
        m = Mesh(n, n, 1.0 / n, 1.0 / n, [1.0], periodic_x=True, periodic_y=True)
        g = Geometry(m, np.ones(m.ncols))
        for shape in ("disc", "bell"):
            if shape == "disc":
                f = subcell_disc(m, 0.5, 0.5, p["radius"])
            else:
                r = np.hypot(m.xc - 0.5, m.yc - 0.5)
                f = np.where(r < 0.25, 0.5 * (1 + np.cos(np.pi * r / 0.25)), 0.0)
            # foamy milk (liquid 0.7, gas 0.3) in a brown base: density jump 0.7 : 1
            C = np.stack([1 - f, 0.7 * f, 0.3 * f])
            vel = np.zeros((3, m.nc)); vel[0], vel[1] = u, v
            pr = Prescribed(m, 1.0, C, rho, v=vel)
            dt0 = p["cfl"] / (abs(u) / m.dx + abs(v) / m.dy)
            steps = int(np.ceil(p["t_end"] / dt0)); dt = p["t_end"] / steps
            F = _uniform_flux(m, g, u, v)
            M0 = pr.M.sum(); V0 = pr.Vc.sum(1)
            for _ in range(steps):
                pr.step(dt * F, dt)
            Cn = pr.C()
            vel_err = float(np.max(np.abs(pr.P / pr.M[None, :] - vel)))
            res[shape].append({"n": n, "L1": l1(Cn[1], C[1], np.ones(m.nc)), "steps": steps,
                               "mass_err": abs(pr.M.sum() - M0) / M0, "vol_err": float(np.max(np.abs(pr.Vc.sum(1) - V0)) / m.nc / m.A),
                               "velocity_uniformity_err": vel_err, "cmin": float(Cn.min()), "cmax": float(Cn.max())})
            if n == p["grids"][-1]:
                figs[shape] = (m, C[1], Cn[1])
    for shape in res:
        r = rate([x["L1"] for x in res[shape]], p["grids"])
        for i, x in enumerate(res[shape][1:]):
            x["rate"] = r[i]
    run.write_log(res["disc"] + res["bell"], "diagnostics/convergence.csv")
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.4))
    for k, shape in enumerate(("disc", "bell")):
        mm, c0, cn = figs[shape]
        ax[k].contour(mm.to_grid(c0), levels=[0.05, 0.5, 0.95], colors="k", linewidths=0.6)
        ax[k].contour(mm.to_grid(cn), levels=[0.05, 0.5, 0.95], colors="r", linewidths=0.6)
        ax[k].set_title(f"{shape}: initial (k) vs after one period (r), n={p['grids'][-1]}", fontsize=8)
        ax[k].set_aspect("equal")
    for shape in res:
        ax[2].loglog(p["grids"], [x["L1"] for x in res[shape]], "o-", label=shape)
    ax[2].set_xlabel("n"); ax[2].set_ylabel("L1 error of C_M"); ax[2].legend()
    savefig(fig, run.path("figures", "T06_translation.png"))
    dec = all(res[s][i + 1]["L1"] < res[s][i]["L1"] for s in res for i in range(len(p["grids"]) - 1))
    worst_mass = max(x["mass_err"] for s in res for x in res[s])
    worst_vel = max(x["velocity_uniformity_err"] for s in res for x in res[s])
    bnd = max(max(-x["cmin"], x["cmax"] - 1, 0) for s in res for x in res[s])
    ok = dec and worst_mass <= 1e-12 and worst_vel <= 1e-12 and bnd <= cfg["thresholds"]["bounds"]
    return {"status": "completed", "verdict": "pass" if ok else "fail",
            "reason": (f"L1 decreases with n: disc {[round(x['L1'], 5) for x in res['disc']]} (rates "
                       f"{[round(x.get('rate', 0), 2) for x in res['disc'][1:]]}), bell "
                       f"{[round(x['L1'], 6) for x in res['bell']]} (rates {[round(x.get('rate', 0), 2) for x in res['bell'][1:]]}); "
                       f"uniform velocity across the 0.7:1 density jump kept to {worst_vel:.1e}"),
            "mass_residual": worst_mass, "bounds_violation": bnd, "details": res}


# ---------------------------------------------------------------- T07 reversed vortex
def T07(cfg, run):
    p = cfg["tests"]["T07"]
    rho = cfg["model"]["rho"]
    T = p["period"]
    out = []
    snaps = {}
    for n in p["grids"]:
        m = Mesh(n, n, 1.0 / n, 1.0 / n, [1.0])
        g = Geometry(m, np.ones(m.ncols))
        f = subcell_disc(m, p["center"][0], p["center"][1], p["radius"])
        C = np.stack([1 - f, f, 0 * f])
        pr = Prescribed(m, 1.0, C, rho)
        base = _psi_flux(m, g, lambda x, y: np.sin(np.pi * x) ** 2 * np.sin(np.pi * y) ** 2 / np.pi)
        umax = 1.0
        dt0 = p["cfl"] * m.dx / (2 * umax)
        steps = int(np.ceil(T / dt0)); dt = T / steps
        mix = []
        half = None
        for k in range(steps):
            tm = (k + 0.5) * dt
            pr.step(dt * np.cos(np.pi * tm / T) * base, dt)
            if k % max(1, steps // 100) == 0:
                mix.append((pr.t, _mixing(pr.C(), pr.Vtot)))
            if half is None and pr.t >= T / 2 - 1e-12:
                half = pr.C()[1].copy()
        Cn = pr.C()
        # mean filament width at T/2 from area / perimeter of the C_M = 0.5 contour
        import matplotlib.pyplot as _plt
        fig_ = _plt.figure(); cs = _plt.contour((np.arange(n) + 0.5) / n, (np.arange(n) + 0.5) / n, m.to_grid(half), levels=[0.5])
        per = 0.0
        for path in cs.get_paths():
            v_ = path.vertices
            per += float(np.sum(np.hypot(np.diff(v_[:, 0]), np.diff(v_[:, 1]))))
        _plt.close(fig_)
        area = float((half >= 0.5).sum() * m.A)
        out.append({"n": n, "L1": l1(Cn[1], C[1], np.ones(m.nc)), "steps": steps,
                    "mass_err": float(abs(pr.Vc[1].sum() - C[1].sum() * m.A) / (C[1].sum() * m.A)),
                    "mixing_final": _mixing(Cn, pr.Vtot), "mixing_max": max(x[1] for x in mix),
                    "mean_width_T2": 2 * area / per if per > 0 else None, "cmin": float(Cn.min()), "cmax": float(Cn.max())})
        snaps[n] = (m, C[1], half, Cn[1], mix)
    r = rate([x["L1"] for x in out], p["grids"])
    for i, x in enumerate(out[1:]):
        x["rate"] = r[i]
    run.write_log(out, "diagnostics/convergence.csv")
    fig, ax = plt.subplots(1, 4, figsize=(14, 3.5))
    mm, c0, ch, cn, mix = snaps[p["grids"][-1]]
    for a_, f_, ttl in ((ax[0], c0, "t=0"), (ax[1], ch, "t=T/2 (max stretch)"), (ax[2], cn, "t=T (returned)")):
        a_.imshow(mm.to_grid(f_), origin="lower", extent=[0, 1, 0, 1], vmin=0, vmax=1, cmap="copper_r")
        a_.contour((np.arange(mm.nx) + 0.5) / mm.nx, (np.arange(mm.ny) + 0.5) / mm.ny, mm.to_grid(c0), levels=[0.5], colors="c", linewidths=0.6)
        a_.set_title(f"{ttl}, n={p['grids'][-1]}", fontsize=8)
    for n in p["grids"]:
        mx = np.array(snaps[n][4])
        ax[3].plot(mx[:, 0], mx[:, 1], label=f"n={n}")
    ax[3].set_xlabel("t"); ax[3].set_ylabel("mixing indicator M"); ax[3].legend(fontsize=7)
    savefig(fig, run.path("figures", "T07_reversed_vortex.png"))
    dec = all(out[i + 1]["L1"] < out[i]["L1"] for i in range(len(out) - 1))
    bnd = max(max(-x["cmin"], x["cmax"] - 1, 0) for x in out)
    ok = dec and max(x["mass_err"] for x in out) <= 1e-12 and bnd <= cfg["thresholds"]["bounds"]
    return {"status": "completed", "verdict": "pass" if ok else "fail",
            "reason": (f"return error L1 {[round(x['L1'], 5) for x in out]} for n={p['grids']} (rates "
                       f"{[round(x['rate'], 2) for x in out[1:]]}); mean filament width at T/2 "
                       f"{[round(x['mean_width_T2'], 4) for x in out]}; final M {[round(x['mixing_final'], 5) for x in out]}"),
            "mass_residual": max(x["mass_err"] for x in out), "bounds_violation": bnd,
            "mixture_indicator": [x["mixing_final"] for x in out], "details": out}


# ---------------------------------------------------------------- T08 piston displacement
def T08(cfg, run):
    p = cfg["tests"]["T08"]
    rho = np.asarray(cfg["model"]["rho"], float)
    U, tend = p["U"], p["t_end"]
    out, prof = [], {}
    for n in p["grids"]:
        m = Mesh(n, 1, 1.0 / n, 1.0 / n, [1.0], periodic_x=False, periodic_y=True)
        g = Geometry(m, np.ones(m.ncols))
        C = np.stack([np.ones(m.nc), np.zeros(m.nc), np.zeros(m.nc)])
        vel = np.zeros((3, m.nc)); vel[0] = U
        pr = Prescribed(m, 1.0, C, rho, v=vel)
        a = m.dy * 1.0
        dt0 = p["cfl"] * m.dx / U
        steps = int(np.ceil(tend / dt0)); dt = tend / steps
        F = _uniform_flux(m, g, U, 0.0)
        left = np.array([m.col_of[0, 0]]); right = np.array([m.col_of[0, n - 1]])
        qv = U * a * dt
        rho_m = float(rho[1])
        for _ in range(steps):
            inflow = {"cell": left, "vol": np.array([qv]), "comp": np.array([[0.0], [qv], [0.0]]),
                      "mass": np.array([rho_m * qv]), "mom": np.array([[rho_m * qv * U], [0.0], [0.0]])}
            outflow = {"cell": right, "vol": np.array([qv])}
            pr.step(dt * F, dt, inflow, outflow)
        Cn = pr.C()
        xM = Cn[1] / (Cn[0] + Cn[1])
        x = m.xc
        exact = np.clip((U * tend - (x - 0.5 * m.dx)) / m.dx, 0, 1)
        def cross(level):
            idx = np.nonzero((xM[:-1] - level) * (xM[1:] - level) <= 0)[0]
            if idx.size == 0:
                return np.nan
            i = idx[0]
            return float(x[i] + (level - xM[i]) * (x[i + 1] - x[i]) / (xM[i + 1] - xM[i]))
        xf = cross(0.5)
        down = xM[x > U * tend + 0.1]
        out.append({"n": n, "front_x": xf, "front_err": abs(xf - U * tend), "width_10_90": cross(0.1) - cross(0.9),
                    "L1": float(np.abs(xM - exact).mean()), "max_xM_beyond_front_plus_0.1": float(down.max()),
                    "min_xM_behind_front_minus_0.1": float(xM[x < U * tend - 0.1].min()),
                    "milk_err": abs(pr.Vc[1].sum() - (pr.inflow_tot[1] - pr.outflow_tot[1])) / pr.inflow_tot[1]})
        prof[n] = (x, xM, exact)
    r = rate([o["L1"] for o in out], p["grids"])
    for i, o in enumerate(out[1:]):
        o["rate"] = r[i]
    run.write_log(out, "diagnostics/piston.csv")
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    for n in p["grids"]:
        x, xm_, ex = prof[n]
        ax[0].plot(x, xm_, label=f"n={n}")
    ax[0].plot(prof[p["grids"][-1]][0], prof[p["grids"][-1]][2], "k--", lw=0.8, label="exact x_f = U t")
    ax[0].set_xlim(0.3, 0.7); ax[0].set_xlabel("x"); ax[0].set_ylabel("milk fraction of liquid"); ax[0].legend(fontsize=7)
    ax[1].loglog(p["grids"], [o["width_10_90"] for o in out], "o-", label="10-90 width")
    ax[1].loglog(p["grids"], [o["front_err"] + 1e-16 for o in out], "s-", label="|x_f - Ut|")
    ax[1].loglog(p["grids"], [o["L1"] for o in out], "^-", label="L1")
    ax[1].set_xlabel("n"); ax[1].legend(fontsize=7)
    savefig(fig, run.path("figures", "T08_piston.png"))
    ok = (all(o["front_err"] <= 2.0 / o["n"] for o in out) and all(out[i + 1]["width_10_90"] < out[i]["width_10_90"] for i in range(len(out) - 1))
          and max(o["max_xM_beyond_front_plus_0.1"] for o in out) <= 1e-6 and max(o["milk_err"] for o in out) <= 1e-12)
    return {"status": "completed", "verdict": "pass" if ok else "fail",
            "reason": (f"front x_f vs Ut=0.5: errors {[f'{o['front_err']:.1e}' for o in out]}; 10-90 width "
                       f"{[round(o['width_10_90'], 4) for o in out]} (shrinks with n); downstream milk beyond x_f+0.1 <= "
                       f"{max(o['max_xM_beyond_front_plus_0.1'] for o in out):.1e} (no whole-pipe whitening)"),
            "mass_residual": max(o["milk_err"] for o in out), "details": out}


# ---------------------------------------------------------------- T09 simplex
def T09(cfg, run):
    p = cfg["tests"]["T09"]
    rho = cfg["model"]["rho"]
    rng = np.random.default_rng(cfg["seed"])
    n = p["n"]
    m = Mesh(n, n, 1.0 / n, 1.0 / n, [1.0], periodic_x=True, periodic_y=True)
    g = Geometry(m, np.ones(m.ncols))
    C = rng.dirichlet([0.5, 0.5, 0.5], size=m.nc).T
    C[:, rng.random(m.nc) < 0.2] = np.array([[0.0], [1.0], [0.0]])     # pure states as hard bounds
    pr = Prescribed(m, 1.0, C, rho)
    base = _psi_flux(m, g, lambda x, y: np.sin(np.pi * x) ** 2 * np.sin(np.pi * y) ** 2 / np.pi)
    dt = p["cfl"] * m.dx / 2
    V0 = pr.Vc.sum(1); VG0 = (pr.Vtot - pr.Vc[0] - pr.Vc[1]).sum()
    lo0, hi0 = C.min(1), C.max(1)
    worst_out = 0.0
    for _ in range(p["steps"]):
        pr.step(dt * base, dt)
        Cn = pr.C()
        worst_out = max(worst_out, float(np.max(lo0[:, None] - Cn)), float(np.max(Cn - hi0[:, None])))
    Cn = pr.C()
    cons = float(max(np.max(np.abs(pr.Vc.sum(1) - V0)), abs((pr.Vtot - pr.Vc[0] - pr.Vc[1]).sum() - VG0)) / pr.Vtot.sum())
    simplex = float(np.max(np.abs(Cn.sum(0) - 1)))
    ok = cons <= 1e-12 and max(worst_out, 0) <= cfg["thresholds"]["bounds"] and simplex <= cfg["thresholds"]["simplex"]
    return {"status": "completed", "verdict": "pass" if ok else "fail",
            "reason": (f"{p['steps']} FCT steps on random simplex data: component conservation {cons:.1e}, "
                       f"max excursion beyond initial per-component range {max(worst_out, 0):.1e}, "
                       f"min fraction {Cn.min():.2e}; sum-to-one {simplex:.1e} (structural: third component derived)"),
            "mass_residual": cons, "simplex_error": simplex, "bounds_violation": max(worst_out, 0.0)}


# ---------------------------------------------------------------- T10 remap
def T10(cfg, run):
    p = cfg["tests"]["T10"]
    rng = np.random.default_rng(cfg["seed"])

    def prof(z):
        cm = 0.5 + 0.4 * np.tanh((z - 0.6) / 0.06)
        cg = 0.1 * z
        return np.stack([1 - cm - cg, cm, cg])

    def cell_avg(zf, f, k=64):
        out = []
        for a, b in zip(zf[:-1], zf[1:]):
            xg, wg = np.polynomial.legendre.leggauss(k)
            zz = 0.5 * (b - a) * xg + 0.5 * (a + b)
            out.append((0.5 * (b - a) * wg[None, :] * f(zz)).sum(1) / (b - a))
        return np.array(out).T

    rows = []
    for order in (1, 2):
        for N in p["N_list"]:
            zu = np.linspace(0, 1, N + 1)
            fr = layer_fractions(N, "geometric", 0.8 ** (8.0 / N)); zg = np.concatenate([[0], np.cumsum(fr)]); zg[-1] = 1
            rnd = rng.random(N) + 0.5; zr = np.concatenate([[0], np.cumsum(rnd / rnd.sum())]); zr[-1] = 1
            C0 = cell_avg(zu, prof)
            h = np.diff(zu)
            vol = C0 * h[None, :]
            mom = (np.sin(2 * np.pi * 0.5 * (zu[1:] + zu[:-1])) * h)[None, :]
            tot0, m0 = vol.sum(1), mom.sum()
            z = zu
            worst_b = 0.0
            for _ in range(p["cycles"]):
                for znew in (zg, zr, zu):
                    vol, mom = remap_column(z, znew, vol, mom, order=order)
                    z = znew
                    Cn = vol / np.diff(z)[None, :]
                    worst_b = max(worst_b, float(-Cn.min()), float(Cn.max() - 1))
            Cn = vol / np.diff(zu)[None, :]
            rows.append({"order": order, "N": N, "profile_L1": float(np.abs(Cn - C0).sum(0).mean()),
                         "comp_conservation": float(np.max(np.abs(vol.sum(1) - tot0))),
                         "mom_conservation": float(abs(mom.sum() - m0)), "bounds_excursion": max(worst_b, 0.0),
                         "simplex": float(np.max(np.abs(Cn.sum(0) - 1)))})
    run.write_log(rows, "diagnostics/remap.csv")
    fig, ax = plt.subplots(figsize=(5, 3.4))
    for order in (1, 2):
        rr = [r for r in rows if r["order"] == order]
        ax.loglog([r["N"] for r in rr], [r["profile_L1"] for r in rr], "o-", label=f"order {order}")
    ax.set_xlabel("layers N"); ax.set_ylabel(f"profile L1 after {p['cycles']} x 3 remaps"); ax.legend()
    savefig(fig, run.path("figures", "T10_remap.png"))
    dec = all(rows[i + 1]["profile_L1"] < rows[i]["profile_L1"] for i in range(len(rows) - 1) if rows[i]["order"] == rows[i + 1]["order"])
    cons = max(max(r["comp_conservation"], r["mom_conservation"]) for r in rows)
    bnd = max(r["bounds_excursion"] for r in rows)
    ok = dec and cons <= 1e-13 and bnd <= cfg["thresholds"]["bounds"]
    r2 = [r for r in rows if r["order"] == 2]
    return {"status": "completed", "verdict": "pass" if ok else "fail",
            "reason": (f"30 remaps (uniform->geometric->random->uniform x10): conservation {cons:.1e}, bounds {bnd:.1e}; "
                       f"profile L1 (order 2) {[f'{r['profile_L1']:.2e}' for r in r2]} for N={p['N_list']}"),
            "mass_residual": cons, "bounds_violation": bnd, "simplex_error": max(r["simplex"] for r in rows), "details": rows}


# ---------------------------------------------------------------- T14 interlayer exchange
def T14(cfg, run):
    p = cfg["tests"]["T14"]
    rho = np.asarray(cfg["model"]["rho"], float)
    m = Mesh(1, 1, 1.0, 1.0, [0.5, 0.5], periodic_x=True, periodic_y=True)
    # lower layer: brown at rest; upper layer: foamy milk moving in x
    C = np.array([[1.0, 0.0], [0.0, 0.7], [0.0, 0.3]])
    vel = np.zeros((3, 2)); vel[0, 1] = 1.0
    checks = []
    for direction in (-1, +1):
        pr = Prescribed(m, 1.0, C, rho, v=vel)
        V3 = lambda: np.stack([pr.Vc[0], pr.Vc[1], pr.Vtot - pr.Vc[0] - pr.Vc[1]])
        before = V3().copy(); P0 = pr.P.copy(); M0 = pr.M.copy()
        dV = np.array([direction * p["dV_fraction"] * pr.Vtot[0]])
        donor = 1 if direction < 0 else 0
        Cd = (V3()[:, donor] / pr.Vtot[donor]).copy(); vd = (pr.P[:, donor] / pr.M[donor]).copy()
        pr.step(dV, 1.0)
        after = V3()
        lost = before[:, donor] - after[:, donor]
        gained = after[:, 1 - donor] - before[:, 1 - donor]
        checks.append({"direction": "down" if direction < 0 else "up",
                       "loss_equals_gain": float(np.max(np.abs(lost - gained))),
                       "carried_composition_err": float(np.max(np.abs(lost / abs(dV[0]) - Cd))),
                       "momentum_transfer_err": float(np.max(np.abs((P0[:, donor] - pr.P[:, donor]) - (rho @ Cd) * abs(dV[0]) * vd))),
                       "total_momentum_err": float(np.max(np.abs(pr.P.sum(1) - P0.sum(1)))),
                       "total_mass_err": float(abs(pr.M.sum() - M0.sum()))})
    run.write_log(checks, "diagnostics/exchange.csv")
    worst = max(max(c.values()) if isinstance(c, dict) else 0 for c in [{k: v for k, v in c.items() if k != "direction"} for c in checks])
    ok = worst <= 1e-14
    return {"status": "completed", "verdict": "pass" if ok else "fail",
            "reason": f"one exchange step each way: loss=gain, donor composition and donor momentum carried, worst error {worst:.1e}",
            "mass_residual": max(c["total_mass_err"] for c in checks), "details": checks}
