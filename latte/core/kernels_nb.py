"""Numba kernels mirroring the numpy reference operators one-to-one.

Each function is a plain loop nest over cells / faces, i.e. the same structure the phase-two
Warp kernels will have.  Equivalence with the numpy reference is gate G0.8.
"""
import os
import numpy as np
import numba as nb

# the CI / sandbox container gives ~1 core; prange loops run serially unless LATTE_PARALLEL=1
PAR = os.environ.get("LATTE_PARALLEL", "0") == "1"


@nb.njit(cache=True, inline="always")
def _mc(dl, dr):
    if dl * dr <= 0.0:
        return 0.0
    a = min(2.0 * abs(dl), 2.0 * abs(dr), 0.5 * abs(dl + dr))
    return a if dl > 0 else -a


@nb.njit(cache=True, inline="always")
def _minmod(dl, dr):
    if dl * dr <= 0.0:
        return 0.0
    a = min(abs(dl), abs(dr))
    return a if dl > 0 else -a


@nb.njit(cache=True, parallel=PAR)
def rhs_kernel(U, W, h, rb, Ps, c, ox, oy, wxr, wxl, wyt, wyb, wet, dx, kap, nu, use_mc, dU, wall):
    """W: primitive stack [M, u, v, Y, Yk...]; fills dU and wall[2] (force on the layer)."""
    nv, N, _ = W.shape
    nc = U.shape[0]
    for q in range(nc):
        for j in range(N):
            for i in range(N):
                dU[q, j, i] = 0.0
    # ---- x faces (rows are independent)
    for j in nb.prange(N):
        WL = np.empty(nv)
        WR = np.empty(nv)
        fl = np.empty(nc)
        for i in range(N - 1):
            if not ox[j, i]:
                continue
            for q in range(nv):
                d0 = W[q, j, i + 1] - W[q, j, i]
                dlm = (W[q, j, i] - W[q, j, i - 1]) if (i >= 1 and ox[j, i - 1]) else 0.0
                drp = (W[q, j, i + 2] - W[q, j, i + 1]) if (i + 2 < N and ox[j, i + 1]) else 0.0
                if use_mc:
                    sl = _mc(dlm, d0)
                    sr = _mc(d0, drp)
                else:
                    sl = _minmod(dlm, d0)
                    sr = _minmod(d0, drp)
                WL[q] = W[q, j, i] + 0.5 * sl
                WR[q] = W[q, j, i + 1] - 0.5 * sr
            un = 0.5 * (WL[1] + WR[1])
            F = un * (WL[0] if un > 0 else WR[0])
            cf = max(c[j, i], c[j, i + 1])
            F -= kap * cf * 0.5 * (rb[j, i] + rb[j, i + 1]) * (h[j, i + 1] - h[j, i]) / dx
            up = F > 0
            Mf = 0.5 * (W[0, j, i] + W[0, j, i + 1])
            fl[0] = F
            fl[1] = F * (WL[1] if up else WR[1]) + 0.5 * (Ps[j, i] + Ps[j, i + 1]) \
                - nu * Mf * (W[1, j, i + 1] - W[1, j, i]) / dx
            fl[2] = F * (WL[2] if up else WR[2]) - nu * Mf * (W[2, j, i + 1] - W[2, j, i]) / dx
            for q in range(3, nc):
                fl[q] = F * (WL[q] if up else WR[q])
            for q in range(nc):
                dU[q, j, i] -= fl[q] / dx
                dU[q, j, i + 1] += fl[q] / dx
    # ---- y faces (columns are independent)
    for i in nb.prange(N):
        WL = np.empty(nv)
        WR = np.empty(nv)
        fl = np.empty(nc)
        for j in range(N - 1):
            if not oy[j, i]:
                continue
            for q in range(nv):
                d0 = W[q, j + 1, i] - W[q, j, i]
                dlm = (W[q, j, i] - W[q, j - 1, i]) if (j >= 1 and oy[j - 1, i]) else 0.0
                drp = (W[q, j + 2, i] - W[q, j + 1, i]) if (j + 2 < N and oy[j + 1, i]) else 0.0
                if use_mc:
                    sl = _mc(dlm, d0)
                    sr = _mc(d0, drp)
                else:
                    sl = _minmod(dlm, d0)
                    sr = _minmod(d0, drp)
                WL[q] = W[q, j, i] + 0.5 * sl
                WR[q] = W[q, j + 1, i] - 0.5 * sr
            un = 0.5 * (WL[2] + WR[2])
            F = un * (WL[0] if un > 0 else WR[0])
            cf = max(c[j, i], c[j + 1, i])
            F -= kap * cf * 0.5 * (rb[j, i] + rb[j + 1, i]) * (h[j + 1, i] - h[j, i]) / dx
            up = F > 0
            Mf = 0.5 * (W[0, j, i] + W[0, j + 1, i])
            fl[0] = F
            fl[1] = F * (WL[1] if up else WR[1]) - nu * Mf * (W[1, j + 1, i] - W[1, j, i]) / dx
            fl[2] = F * (WL[2] if up else WR[2]) + 0.5 * (Ps[j, i] + Ps[j + 1, i]) \
                - nu * Mf * (W[2, j + 1, i] - W[2, j, i]) / dx
            for q in range(3, nc):
                fl[q] = F * (WL[q] if up else WR[q])
            for q in range(nc):
                dU[q, j, i] -= fl[q] / dx
                dU[q, j + 1, i] += fl[q] / dx
    wall[0] = 0.0
    wall[1] = 0.0
    for j in range(N):
        for i in range(N):
            if not wet[j, i]:
                for q in range(nc):
                    dU[q, j, i] = 0.0
                continue
            fx = -Ps[j, i] * ((1.0 if wxr[j, i] else 0.0) - (1.0 if wxl[j, i] else 0.0)) / dx
            fy = -Ps[j, i] * ((1.0 if wyt[j, i] else 0.0) - (1.0 if wyb[j, i] else 0.0)) / dx
            dU[1, j, i] += fx
            dU[2, j, i] += fy
            wall[0] += fx
            wall[1] += fy


@nb.njit(cache=True, inline="always")
def _bil(F, q, j, i, n):
    if j < 0.0:
        j = 0.0
    if i < 0.0:
        i = 0.0
    if j > n - 1.0:
        j = n - 1.0
    if i > n - 1.0:
        i = n - 1.0
    j0 = min(int(np.floor(j)), n - 2)
    i0 = min(int(np.floor(i)), n - 2)
    wj = j - j0
    wi = i - i0
    return (1 - wj) * ((1 - wi) * F[q, j0, i0] + wi * F[q, j0, i0 + 1]) + \
        wj * ((1 - wi) * F[q, j0 + 1, i0] + wi * F[q, j0 + 1, i0 + 1])


@nb.njit(cache=True, parallel=PAR)
def trace_kernel(cvel, sign, L, dxf, R, out_j, out_i):
    """RK2 back-trace (sign=+1) or forward trace (sign=-1) in index space, clamped into the cup."""
    n = cvel.shape[1]
    rmax = R - 0.5 * dxf
    for j in nb.prange(n):
        for i in range(n):
            jm = j - 0.5 * sign * cvel[1, j, i]
            im = i - 0.5 * sign * cvel[0, j, i]
            um = _bil(cvel, 0, jm, im, n)
            vm = _bil(cvel, 1, jm, im, n)
            jd = j - sign * vm
            idp = i - sign * um
            xc = (idp + 0.5) * dxf - 0.5 * L
            yc = (jd + 0.5) * dxf - 0.5 * L
            rr = np.sqrt(xc * xc + yc * yc)
            sc = rmax / rr if rr > rmax else 1.0
            out_j[j, i] = (yc * sc + 0.5 * L) / dxf - 0.5
            out_i[j, i] = (xc * sc + 0.5 * L) / dxf - 0.5


@nb.njit(cache=True, parallel=PAR)
def sample_kernel(F, tj, ti, out):
    nq, n, _ = F.shape
    for j in nb.prange(n):
        for i in range(n):
            for q in range(nq):
                out[q, j, i] = _bil(F, q, tj[j, i], ti[j, i], n)


@nb.njit(cache=True)
def extrapolate_kernel(F, dry_j, dry_i, near_j, near_i):
    for k in range(dry_j.shape[0]):
        for q in range(F.shape[0]):
            F[q, dry_j[k], dry_i[k]] = F[q, near_j[k], near_i[k]]


@nb.njit(cache=True, parallel=PAR)
def bfecc_finish_kernel(F, corr, tj, ti, wet, out):
    """out = clamp(sample(corr), local min/max of F at the departure cell); tags <= Y."""
    nq, n, _ = F.shape
    for j in nb.prange(n):
        for i in range(n):
            if not wet[j, i]:
                for q in range(nq):
                    out[q, j, i] = 0.0
                continue
            jj = min(max(tj[j, i], 0.0), n - 1.0)
            ii = min(max(ti[j, i], 0.0), n - 1.0)
            j0 = min(int(np.floor(jj)), n - 2)
            i0 = min(int(np.floor(ii)), n - 2)
            for q in range(nq):
                v = _bil(corr, q, jj, ii, n)
                a = F[q, j0, i0]
                b = F[q, j0, i0 + 1]
                c = F[q, j0 + 1, i0]
                d = F[q, j0 + 1, i0 + 1]
                lo = min(min(a, b), min(c, d))
                hi = max(max(a, b), max(c, d))
                if v < lo:
                    v = lo
                if v > hi:
                    v = hi
                out[q, j, i] = v
            for q in range(1, nq):
                if out[q, j, i] > out[0, j, i]:
                    out[q, j, i] = out[0, j, i]
