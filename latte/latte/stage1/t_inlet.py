"""G1 inlet tests T02 (flux integral) and T03 (input contract)."""
import numpy as np

from ..inlet import (InletSection, PiecewiseInlet, NoInlet, surface_source, section_totals,
                     FLOW, NO_FLOW, UNKNOWN)
from ..mesh import cylinder_mesh, layer_fractions
from ..solver import Solver
from ..failures import SimFailure
from .common import make_model, make_numerics, savefig, plt


def _tilted(c, tilt_deg):
    th = np.radians(tilt_deg)
    n = np.array([np.sin(th), 0.0, np.cos(th)])
    e1 = np.array([np.cos(th), 0.0, -np.sin(th)])
    return n, e1


def T02(cfg, run):
    p = cfg["tests"]["T02"]
    rho = np.asarray(cfg["model"]["rho"], float)
    n, e1 = _tilted(p["center"], p["tilt_deg"])
    e2 = np.cross(n, e1)
    a, b, U = p["a"], p["b"], p["U"]
    area = np.pi * a * b
    C = np.array([0.0, 0.8, 0.2])                       # milk with effective gas fraction 0.2
    rhoJ = float(rho @ C)
    cases = {
        "uniform_normal": dict(v=-U * n, prof="uniform", w=np.zeros(3), Q=area * U, mom=rhoJ * area * U * (-U * n)),
        "uniform_oblique": dict(v=-U * n + p["tangential"] * U * e1, prof="uniform", w=np.zeros(3), Q=area * U,
                                mom=rhoJ * area * U * (-U * n + p["tangential"] * U * e1)),
        "parabolic_oblique": dict(v=-U * n + p["tangential"] * U * e1, prof="parabolic", w=np.zeros(3), Q=area * U / 2,
                                  mom=rhoJ * U * area / 3 * (-U * n + p["tangential"] * U * e1)),
        "moving_inplane": dict(v=-U * n, prof="uniform", w=p["w_inplane"] * (e1 + e2) / np.sqrt(2), Q=area * U,
                               mom=rhoJ * area * U * (-U * n)),
        "moving_normal": dict(v=-U * n, prof="uniform", w=-p["w_normal"] * n, Q=area * (U - p["w_normal"]),
                              mom=rhoJ * area * (U - p["w_normal"]) * (-U * n)),
    }
    rows = []
    worst = 0.0
    for name, cs in cases.items():
        for nr, nt in p["quadratures"]:
            sec = InletSection(center=np.array(p["center"]), normal=n, axis1=e1, a=a, b=b, velocity=cs["v"],
                               composition=C, section_velocity=cs["w"], profile=cs["prof"])
            tot = section_totals(sec, rho, nr, nt)
            eQ = abs(tot["Q"] - cs["Q"]) / cs["Q"]
            eM = float(np.linalg.norm(tot["mom"] - cs["mom"]) / np.linalg.norm(cs["mom"]))
            eC = float(np.max(np.abs(tot["comp"] - cs["Q"] * C)) / cs["Q"])
            rows.append({"case": name, "n_r": nr, "n_theta": nt, "err_Q": eQ, "err_mom": eM, "err_comp": eC})
        worst = max(worst, rows[-1]["err_Q"], rows[-1]["err_mom"], rows[-1]["err_comp"])
    run.write_log(rows, "diagnostics/quadrature.csv")

    # mapping onto surface columns: exact totals and convergence of the per-column shares
    ng = p["map_grid"]
    mesh = cylinder_mesh(ng, 0.5, layer_fractions(4))
    D = np.full(mesh.ncols, 0.25)
    sec = InletSection(center=np.array(p["center"]), normal=n, axis1=e1, a=a, b=b, velocity=-U * n, composition=C)
    from ..inlet import InletState
    ref = surface_source(InletState(FLOW, sec), mesh, D, rho, *p["map_reference"])
    map_rows = []
    for nr, nt in p["map_quadratures"]:
        s = surface_source(InletState(FLOW, sec), mesh, D, rho, nr, nt)
        map_rows.append({"n_r": nr, "n_theta": nt, "sum_err": abs(s.vol_rate.sum() - area * U) / (area * U),
                         "share_L1_vs_ref": float(np.abs(s.vol_rate - ref.vol_rate).sum() / ref.vol_rate.sum()),
                         "columns_hit": int((s.vol_rate > 0).sum()),
                         "mom_sum_err": float(np.linalg.norm(s.mom_rate.sum(1) - rhoJ * area * U * (-U * n)) / (rhoJ * area * U * U))})
    run.write_log(map_rows, "diagnostics/mapping.csv")
    share = [r["share_L1_vs_ref"] for r in map_rows]
    mono = all(share[i + 1] <= share[i] * 1.05 for i in range(len(share) - 1))
    # mesh-aligned sub-grid mapping (solver default)
    from ..inlet import surface_source_subgrid
    refs = surface_source_subgrid(InletState(FLOW, sec), mesh, D, rho, sub=p["subgrid_reference"])
    sub_rows = []
    for sb in p["subgrid_levels"]:
        ss = surface_source_subgrid(InletState(FLOW, sec), mesh, D, rho, sub=sb)
        sub_rows.append({"sub": sb, "sum_err": abs(ss.vol_rate.sum() - area * U) / (area * U),
                         "share_L1_vs_ref": float(np.abs(ss.vol_rate - refs.vol_rate).sum() / refs.vol_rate.sum()),
                         "share_max_vs_ref": float(np.abs(ss.vol_rate - refs.vol_rate).max() / refs.vol_rate.max()),
                         "mom_sum_err": float(np.linalg.norm(ss.mom_rate.sum(1) - rhoJ * area * U * (-U * n)) / (rhoJ * area * U * U))})
    run.write_log(sub_rows, "diagnostics/mapping_subgrid.csv")
    sshare = [r["share_L1_vs_ref"] for r in sub_rows]
    mono = mono and all(sshare[i + 1] < sshare[i] for i in range(len(sshare) - 1))
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    for name in cases:
        rr = [r for r in rows if r["case"] == name]
        ax[0].semilogy([r["n_r"] * r["n_theta"] for r in rr], [max(r["err_Q"], 1e-17) for r in rr], "o-", label=name)
    ax[0].set_xlabel("quadrature points"); ax[0].set_ylabel("relative error of Q"); ax[0].legend(fontsize=7)
    ax[1].loglog([r["n_r"] * r["n_theta"] for r in map_rows], share, "s-", label="polar samples")
    ax[1].loglog([(r["sub"] ** 2) * (refs.vol_rate > 0).sum() for r in sub_rows], sshare, "o-", label="mesh-aligned sub-grid (default)")
    ax[1].set_xlabel("samples"); ax[1].set_ylabel("L1 of column shares vs reference"); ax[1].legend(fontsize=7)
    savefig(fig, run.path("figures", "T02_quadrature.png"))
    ok = (worst <= cfg["thresholds"]["local_update"] and max(r["sum_err"] for r in map_rows + sub_rows) <= 1e-12 and mono
          and max(r["mom_sum_err"] for r in sub_rows) <= 1e-12)
    return {"status": "completed", "verdict": "pass" if ok else "fail",
            "reason": (f"worst finest-quadrature relative error {worst:.2e}; column sum exact; "
                       f"column-share L1 vs reference: polar {share[0]:.2e} -> {share[-1]:.2e}, "
                       f"sub-grid {sshare[0]:.2e} -> {sshare[-1]:.2e} (sub {p['subgrid_levels']})"),
            "volume_residual": max(r["sum_err"] for r in map_rows),
            "details": {"quadrature": rows, "mapping": map_rows, "mapping_subgrid": sub_rows}}


def _expect(fn, code):
    try:
        fn()
    except SimFailure as e:
        return e.code == code, e.code, e.reason
    return False, None, "no failure raised"


def T03(cfg, run):
    p = cfg["tests"]["T03"]
    rho = np.asarray(cfg["model"]["rho"], float)
    mesh = cylinder_mesh(p["grid"], 0.5, layer_fractions(4))
    D = np.full(mesh.ncols, p["depth"])
    zs = p["depth"]
    good = lambda **kw: InletSection(**{**dict(center=np.array([0.0, 0.0, zs]), normal=np.array([0, 0, 1.0]),
                                             axis1=np.array([1.0, 0, 0]), a=0.02, b=0.02,
                                             velocity=np.array([0, 0, -1.0]), composition=np.array([0, 1.0, 0])), **kw})
    from ..inlet import InletState
    checks = []

    def add(name, ok, got, why):
        checks.append({"case": name, "ok": bool(ok), "got": got, "detail": why})

    s0 = surface_source(InletState(NO_FLOW), mesh, D, rho)
    add("zero_flow_returns_zero_source", s0.vol_rate.sum() == 0 and s0.mom_rate.sum() == 0, "zero", "")
    prov = PiecewiseInlet([{"t": 0.0, "status": FLOW, "section": good()}, {"t": 0.5, "status": UNKNOWN},
                           {"t": 1.0, "status": NO_FLOW}])
    add("unknown_interval_is_error", *_expect(lambda: prov.sample(0.7), "INPUT_INVALID"))
    add("outside_record_is_error", *_expect(lambda: prov.sample(1.5), "INPUT_INVALID"))
    add("bad_semi_axis", *_expect(lambda: good(a=-0.01).validate(), "GEOMETRY_INVALID"))
    add("non_orthonormal_axes", *_expect(lambda: good(axis1=np.array([1.0, 0, 0.1])).validate(), "GEOMETRY_INVALID"))
    add("bad_composition", *_expect(lambda: good(composition=np.array([0.2, 0.9, 0.0])).validate(), "INPUT_INVALID"))
    add("outflow_through_section", *_expect(
        lambda: surface_source(InletState(FLOW, good(velocity=np.array([0, 0, 1.0]))), mesh, D, rho), "INPUT_INVALID"))
    add("footprint_leaves_cup", *_expect(
        lambda: surface_source(InletState(FLOW, good(center=np.array([0.49, 0.0, zs]))), mesh, D, rho), "OUT_OF_SCOPE"))
    add("no_liquid_connection", *_expect(
        lambda: surface_source(InletState(FLOW, good(center=np.array([0.0, 0.0, zs + 0.2]))), mesh, D, rho,
                               connect_tol=p["connect_tol"]), "GEOMETRY_INVALID"))
    # interpolation within flow nodes, hold before a status change, no interpolation across it
    s1, s2 = good(), good(center=np.array([0.1, 0.0, zs]))
    pv = PiecewiseInlet([{"t": 0.0, "status": FLOW, "section": s1}, {"t": 1.0, "status": FLOW, "section": s2},
                         {"t": 2.0, "status": NO_FLOW}, {"t": 3.0, "status": NO_FLOW}])
    mid = pv.sample(0.5).section.center
    add("linear_interpolation_between_flow_nodes", np.allclose(mid, [0.05, 0, zs]), mid.tolist(), "")
    held = pv.sample(1.5).section.center
    add("hold_before_status_change", np.allclose(held, [0.1, 0, zs]), held.tolist(), "")
    add("no_flow_after_change", pv.sample(2.5).status == NO_FLOW, pv.sample(2.5).status, "")
    add("next_event", pv.next_event(0.3) == 1.0 and pv.next_event(1.0) == 2.0, [pv.next_event(0.3), pv.next_event(1.0)], "")
    add("section_velocity_from_geometry", np.allclose(pv.sample(0.5).section.section_velocity, [0.1, 0, 0]),
        pv.sample(0.5).section.section_velocity.tolist(), "w_J = d(center)/dt; translation is not counted as inflow")

    # solver never steps across an event and the delivered volume equals Q * t_switch exactly
    T = p["switch_time"]
    sec = good()
    prov = PiecewiseInlet([{"t": 0.0, "status": FLOW, "section": sec}, {"t": T, "status": NO_FLOW},
                           {"t": 10.0, "status": NO_FLOW}])
    solver = Solver(mesh, make_model(cfg), make_numerics(cfg, connect_tol=p["connect_tol"]), inlet=prov)
    Cc = np.zeros((3, mesh.N, mesh.ncols)); Cc[0] = 1
    st = solver.initial_state(p["depth"], Cc)
    times = []
    while st.t < 2 * T:
        st, d = solver.advance(st)
        times.append(st.t)
    hit = min(abs(np.array(times) - T))
    Q = np.pi * 0.02 ** 2
    err = abs(st.acc["in_vol"] - Q * T) / (Q * T)
    add("step_lands_on_event", hit < 1e-14, float(hit), "")
    add("delivered_volume_equals_QT", err < 1e-12, float(err), "")
    # a no-flow provider gives bit-identical steps to having no inlet at all
    st_a = solver.initial_state(p["depth"], Cc)
    sa = Solver(mesh, make_model(cfg), make_numerics(cfg), inlet=PiecewiseInlet(
        [{"t": 0.0, "status": NO_FLOW}, {"t": 10.0, "status": NO_FLOW}]))
    sb = Solver(mesh, make_model(cfg), make_numerics(cfg), inlet=NoInlet())
    xa, xb = st_a.copy(), st_a.copy()
    xa.P[0] += 1e-3 * np.sin(np.arange(mesh.nc)); xb.P[0] += 1e-3 * np.sin(np.arange(mesh.nc))
    for _ in range(5):
        xa, _ = sa.advance(xa, dt_try=0.01)
        xb, _ = sb.advance(xb, dt_try=0.01)
    add("no_flow_identical_to_no_inlet", np.array_equal(xa.P, xb.P) and np.array_equal(xa.D, xb.D), "bitwise", "")
    run.write_log(checks, "diagnostics/contract_checks.csv")
    ok = all(c["ok"] for c in checks)
    failed = [c["case"] for c in checks if not c["ok"]]
    return {"status": "completed", "verdict": "pass" if ok else "fail",
            "reason": f"{sum(c['ok'] for c in checks)}/{len(checks)} contract checks pass" + (f"; failed: {failed}" if failed else ""),
            "volume_residual": float(err), "accepted_steps": len(times), "details": checks}
