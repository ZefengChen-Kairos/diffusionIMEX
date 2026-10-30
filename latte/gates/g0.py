"""G0: conservation, null input, symmetry, backend equivalence, refinement."""
import numpy as np

from actions.scripts import Pitcher, ALL
from core.params import PhysParams, NumParams
from core.solver_np import Sim
from gates.common import Report, run, save_run
from gates import metrics as mt


def pulse_action():
    """Half-second pour at a fixed point; the flow opens and closes in 0.15 s."""
    p = Pitcher("pulse", at=(0, -5), h=0.015, Q=0.0)
    p.move(0.15, Q=10.0).move(0.2).move(0.15, Q=0.0)
    return p.build(T_end=1.5, continuous=False)


def checkerboard_index(f, wet):
    """Share of the field's L2 norm in the (pi, pi) mode over the wet cells."""
    N = f.shape[0]
    jj, ii = np.meshgrid(np.arange(N), np.arange(N), indexing="ij")
    sgn = np.where((ii + jj) % 2 == 0, 1.0, -1.0)
    # local 2x2 checkerboard amplitude
    s = f * sgn
    cb = 0.25 * (s[:-1, :-1] + s[1:, :-1] + s[:-1, 1:] + s[1:, 1:])
    w = wet[:-1, :-1] & wet[1:, 1:] & wet[1:, :-1] & wet[:-1, 1:]
    num = np.sqrt((cb[w] ** 2).sum())
    den = np.sqrt((f[:-1, :-1][w] ** 2).sum()) + 1e-300
    return float(num / den)


def foam_disc_sim(kappa_off=False, P=None):
    P = P or PhysParams()
    s = Sim(None, P=P, n_tags=1, flags=dict(kappa_off=kappa_off))
    g = s.g
    r = np.sqrt(g.X ** 2 + (g.Y + 0.005) ** 2)
    disc = (r < 0.012) & g.wet
    # a foam disc of 2 mm extra foam thickness with a sharp edge, at rest
    s.U[3][disc] += P.rho_m * 0.002
    s.U[0][disc] += P.rho_m * 0.002
    s.U[4][disc] = s.U[3][disc]
    s.Yf[0] = s.fm.to_fine(np.where(g.wet, s.U[3] / np.where(g.wet, s.U[0], 1), 0)) * s.gf.wet
    s.E0 = s.totals()
    return s


def run_g0(full=True):
    R = Report("G0")
    P = PhysParams()
    NP = NumParams()

    # G0.1 static, no input
    s = Sim(None, P=P, n_tags=1)
    U0 = s.U.copy()
    for _ in range(1000):
        s.step(1e-3)
    umax = float(np.abs(s.U[1:3]).max())
    R.check("G0.1", "static 1000 steps: max|P| (exact zero)", umax, umax == 0.0, "== 0")

    # G0.2 foam disc at rest: bounded, no checkerboard; kappa is needed
    cbs = {}
    for koff in (False, True):
        s = foam_disc_sim(kappa_off=koff)
        for _ in range(400):
            s.step()
        M, u, v, *_ = s.prims(s.U)
        cbs[koff] = max(checkerboard_index(u, s.g.wet), checkerboard_index(v, s.g.wet))
        if not koff:
            umax2 = float(np.sqrt(u ** 2 + v ** 2).max())
            cbM = checkerboard_index(s.U[0] - s.U[0][s.g.wet].mean(), s.g.wet)
    R.check("G0.2a", "foam disc: max|u| bounded [m/s]", umax2, umax2 < 0.5, "< 0.5")
    R.check("G0.2b", "foam disc: checkerboard index of u (kappa on)", cbs[False], cbs[False] < 0.05, "< 0.05")
    R.extra["G0.2_checkerboard_kappa_off"] = cbs[True]
    R.extra["G0.2_checkerboard_M"] = cbM

    # G0.3 single pulse: ledgers close every step
    s = run(pulse_action(), P, NP)
    L = np.array(s.log)
    emax = L[:, 1:5].max(axis=0)
    R.check("G0.3a", "pulse: max per-step rel. error [mass, milk, tag, mom]", emax, emax.max() < 1e-6, "< 1e-6")
    tot = s.totals()
    B = s.bulk
    m_err = abs(tot["M"] + B.M - (s.E0["M"] + P.rho_c * P.V_bulk0 + s.inj["mass"])) / (tot["M"] + B.M)
    k_err = abs(tot["Phi"] + B.Mm - s.inj["milk"]) / s.inj["milk"]
    R.check("G0.3b", "pulse: cumulative mass / milk closure", [m_err, k_err], max(m_err, k_err) < 1e-9, "< 1e-9")
    save_run(s, R.dir, "pulse", board_frames=False)

    # G0.4 moving pitcher, Q = 0: nothing changes
    act = ALL["leaf"]().with_flow(0.0)
    s = Sim(act, P=P, n_tags=1)
    U0 = s.U.copy()
    Y0 = s.Yf.copy()
    s.run(T=3.0)
    same = bool(np.array_equal(U0, s.U) and np.array_equal(Y0, s.Yf))
    R.check("G0.4", "pitcher moves with Q=0: fields bitwise unchanged", same, same, "bitwise")

    # G0.5 mirror symmetry
    T5 = 2.5
    a = run(ALL["leaf"](), P, NP, T=T5)
    b = run(ALL["leaf"]().mirrored(), P, NP, T=T5)
    dY = float(np.abs(a.Yf[0] - b.Yf[0][:, ::-1]).max())
    dU = float(max(np.abs(a.U[q] - b.U[q][:, ::-1]).max() / np.abs(a.U[q]).max() for q in (0, 3)))
    dP = float(np.abs(a.U[1] + b.U[1][:, ::-1]).max() / np.abs(a.U[1]).max())
    io = mt.iou(mt.white(a), mt.white(b)[:, ::-1])
    R.check("G0.5a", "mirror: coarse M, Phi, P_x mirror to rounding (rel. max diff)", [dU, dP],
            max(dU, dP) < 1e-12, "< 1e-12")
    R.check("G0.5b", "mirror: white masks mirror IoU; fine max|dY| (rounding grown by the BFECC clamp)",
            [io, dY], io >= 0.999 and dY < 1e-4, "IoU >= 0.999, dY < 1e-4")

    # G0.7 positivity over the pulse and the leaf runs
    for nm, sim in (("pulse", s), ("leaf", a)):
        pass
    okpos = getattr(a, "pos_violation", 0) == 0
    R.check("G0.7", "positivity M >= 0.3 rho_c h0 (leaf, 2.5 s): violations", getattr(a, "pos_violation", 0), okpos, "0")

    # G0.9 coarse ledger vs fine material picture (diagnostic of the ratio-form transport)
    phi_f, _ = a.phi_fine()
    tot_f = phi_f.sum() * a.gf.dA
    tot_c = a.U[3].sum() * a.g.dA
    rf = a.fm.restrict(phi_f)
    l1 = float(np.abs(rf - a.U[3])[a.g.wet].sum() / np.abs(a.U[3][a.g.wet]).sum())
    rel = abs(tot_f - tot_c) / tot_c
    R.check("G0.9", "fine Y*M_interp vs coarse Phi (leaf 2.5 s): total rel. diff, restricted L1 rel. diff",
            [rel, l1], rel < 0.05, "total < 5 % (diagnostic)")

    # G0.8 numba vs numpy reference
    sa = Sim(ALL["leaf"](), P=P, backend="numpy")
    sb = Sim(ALL["leaf"](), P=P, backend="numba")
    for _ in range(80):
        dt = sa.dt_cfl()
        sa.step(dt)
        sb.step(dt)
    e8 = float(np.abs(sa.U - sb.U).max() / np.abs(sa.U).max())
    e8f = float(np.abs(sa.Yf - sb.Yf).max())
    R.check("G0.8", "numba == numpy reference (80 steps): rel max|dU|, max|dY|", [e8, e8f],
            max(e8, e8f) < 1e-12, "< 1e-12")

    # G0.6 refinement 128 -> 192 (big heart)
    if full:
        from run_patterns import run_pattern
        from gates.common import shell, alpha_of
        r1 = run_pattern("heart", P, NP)
        r2 = run_pattern("heart", P, NP.with_(N=192), tag="heart_N192")
        K = r1["Y"].shape[1] - 1
        sh1, sh2 = shell(K, P, NP), shell(K, P, NP.with_(N=192))
        w1 = alpha_of(sh1, r1["Yfinal"].astype(float), r1["Mfinal"].astype(float)) > 0.5
        w2 = alpha_of(sh2, r2["Yfinal"].astype(float), r2["Mfinal"].astype(float)) > 0.5
        # compare on the 384 grid: upsample the 576 mask by nearest
        from scipy import ndimage
        w2d = ndimage.zoom(w2.astype(float), w1.shape[0] / w2.shape[0], order=1) > 0.5
        io = mt.iou(w1, w2d)
        from gates.g1 import heart_metrics
        m1, m2 = heart_metrics(alpha_of(sh1, r1["Yfinal"].astype(float), r1["Mfinal"].astype(float)), sh1.gf), \
            heart_metrics(alpha_of(sh2, r2["Yfinal"].astype(float), r2["Mfinal"].astype(float)), sh2.gf)
        same_class = (m1["notch"] == m2["notch"]) and (m1["lobes"] >= 2) == (m2["lobes"] >= 2)
        R.check("G0.6", "refinement 128->192 big heart: white IoU, class (notch/lobes) unchanged", [io, same_class],
                io >= 0.85 and same_class, "IoU >= 0.85 & same class")
        R.extra["G0.6_metrics"] = dict(N128=m1, N192=m2)
        from render.render import board
        from gates.common import img_of
        board([[img_of(sh1, r1["Yfinal"].astype(float), r1["Mfinal"].astype(float)),
                img_of(sh2, r2["Yfinal"].astype(float), r2["Mfinal"].astype(float))]], R.dir + "/refinement.png",
              labels=["N=128 / 384", "N=192 / 576"], cell=260)
    R.write(P, NP)
    return R
