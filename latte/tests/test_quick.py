"""Fast regression checks (a few seconds).  The full verification is latte.stage1.suite."""
import numpy as np
import pytest

from latte.mesh import Mesh, Geometry, cylinder_mesh, layer_fractions
from latte.operators import MetricOperator
from latte.transport import transport
from latte.inlet import InletSection, section_totals, surface_source, InletState, FLOW
from latte.solver import Solver, Model, Numerics
from latte.failures import SimFailure


def sloped_mesh(n=10, N=5):
    m = Mesh(n, n - 2, 1.0 / n, 1.0 / (n - 2), layer_fractions(N), periodic_x=True)
    D = 0.4 + 0.03 * np.sin(2 * np.pi * m.xc) + 0.02 * np.cos(np.pi * m.yc)
    return m, Geometry(m, D)


def test_operator_symmetric_positive_and_hydrostatic_exact():
    m, g = sloped_mesh()
    N = m.N
    op = MetricOperator(m, g, np.ones((N, m.npx)), np.ones((N, m.npy)), np.ones((N, m.ncols)))
    L = op.L.toarray()
    assert np.abs(L - L.T).max() < 1e-14 * np.abs(L).max()
    assert np.linalg.eigvalsh(L).min() > 0
    F = op.flux((-1.3 * g.zc).ravel(), -1.3 * g.D)
    Fx = F[op.slices["x"]]
    assert np.abs(Fx).max() < 1e-13           # no spurious lateral flux of a pure vertical gradient


def test_transport_conserves_and_stays_bounded():
    rng = np.random.default_rng(0)
    m = Mesh(16, 16, 1 / 16, 1 / 16, [0.5, 0.5], periodic_x=True, periodic_y=True)
    g = Geometry(m, np.ones(m.ncols))
    C = rng.dirichlet([0.4, 0.4, 0.4], size=m.nc).T
    V = g.V.ravel().copy()
    Vc = C[:2] * V
    rho = np.array([1.0, 1.0, 1e-3])
    M = (rho @ C) * V
    P = np.zeros((3, m.nc)); P[0] = M
    dV = np.concatenate([0.10 * g.Vfx.ravel(), 0.05 * g.Vfy.ravel(), np.zeros(m.ncols)])   # Courant 0.1 / 0.05
    tot0 = Vc.sum(1)
    for _ in range(20):
        r = transport(m, Vc, V, P, M, dV, rho)
        Vc, V, M, P = r.Vc, r.Vtot, r.M, r.P
    Cn = np.stack([Vc[0], Vc[1], V - Vc[0] - Vc[1]]) / V
    assert np.abs(Vc.sum(1) - tot0).max() < 1e-13
    assert Cn.min() > -1e-12 and Cn.max() < 1 + 1e-12
    assert np.abs(P[0] / M - 1.0).max() < 1e-12      # uniform velocity kept across density jumps


def test_inlet_totals_exact():
    n = np.array([0.0, 0.0, 1.0])
    sec = InletSection(center=np.array([0.0, 0.0, 0.25]), normal=n, axis1=np.array([1.0, 0, 0]), a=0.03, b=0.02,
                       velocity=-n, composition=np.array([0.0, 0.8, 0.2]), profile="parabolic")
    t = section_totals(sec, [1.0, 1.0, 1e-3])
    assert abs(t["Q"] - np.pi * 0.03 * 0.02 / 2) < 1e-15
    m = cylinder_mesh(32, 0.5, layer_fractions(4))
    s = surface_source(InletState(FLOW, sec), m, np.full(m.ncols, 0.25), [1.0, 1.0, 1e-3])
    assert abs(s.vol_rate.sum() - t["Q"]) < 1e-15
    off = InletSection(center=np.array([0.495, 0.0, 0.25]), normal=n, axis1=np.array([1.0, 0, 0]), a=0.03, b=0.02,
                       velocity=-n, composition=np.array([0.0, 1.0, 0.0]))
    with pytest.raises(SimFailure) as e:
        surface_source(InletState(FLOW, off), m, np.full(m.ncols, 0.25), [1.0, 1.0, 1e-3])
    assert e.value.code == "OUT_OF_SCOPE"


def test_rest_state_and_wave_period():
    m = cylinder_mesh(12, 0.5, layer_fractions(4))
    s = Solver(m, Model(), Numerics())
    C = np.zeros((3, 4, m.ncols)); C[0] = 1; C[0, -1] = 0.7; C[2, -1] = 0.3
    st = s.initial_state(0.25, C)
    for _ in range(5):
        st, _ = s.advance(st)
    assert np.abs(st.velocity()).max() < 1e-12
    # coarse non-hydrostatic standing wave: period within 1.5 % of theory at 32 x 8
    m = Mesh(32, 1, 1 / 32, 1 / 32, layer_fractions(8), periodic_x=True, periodic_y=True)
    s = Solver(m, Model(mu=0.0), Numerics())
    C = np.zeros((3, 8, m.ncols)); C[0] = 1
    H, k = 0.25, 2 * np.pi
    st = s.initial_state(H + 2.5e-4 * np.cos(k * m.xc), C)
    T = 2 * np.pi / np.sqrt(k * np.tanh(k * H))
    ts, e = [0.0], [st.D[0] - H]
    while st.t < 1.6 * T:
        st, _ = s.advance(st)
        ts.append(st.t); e.append(st.D[0] - H)
    zc = [ts[i] - e[i] * (ts[i + 1] - ts[i]) / (e[i + 1] - e[i]) for i in range(len(e) - 1) if e[i] * e[i + 1] < 0]
    assert abs(2 * np.mean(np.diff(zc)) / T - 1) < 0.015
