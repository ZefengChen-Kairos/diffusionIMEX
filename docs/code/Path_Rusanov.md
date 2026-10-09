# `Path_Rusanov.cpp`: code notes

Source: `e7a02682-Path_Rusanov.cpp` (1884 lines, C++/Eigen). The line numbers below refer to that file.

**Summary.** The file solves a single-velocity **6-equation two-phase model** (Saurel–Petitpas–Berry / Pelanti–Shyue type) with stiffened-gas (SG) EOS. It uses a **Toro–Vázquez-type flux splitting**:

- **Explicit convective subsystem.** α-transport, phasic mass, momentum convection and the *kinetic-energy* part of the phasic energy fluxes. All of it is written in quasi-linear form `A(Q) ∂x Q` and discretised with a **path-conservative Rusanov** scheme. The path is linear in conserved variables, the path integral uses 3-point Gauss, and the dissipation speed is the global max |u|.
- **Implicit pressure subsystem.** Enthalpy flux, pressure gradient, and the non-conservative term Σ. The unknowns are (α1p1, α2p2), and the system is solved matrix-free with GMRES (two Picard passes).
- **Instantaneous exact pressure relaxation.** Pelanti–Shyue Appendix C.1 quadratic.

The grid is written for 2D (x, y). The only test is run in **1D** (`nyCells = 1`): a shock hitting a helium bubble.

---

## 1. Model and equations

### 1.1 State vector (`NUM_VARS = 7`, l.19)

| k | conserved `u_con` (`Q`) | primitive `u_pri` |
|---|---|---|
| 0 | α1 | α1 |
| 1 | α1ρ1 | ρ1 |
| 2 | α2ρ2 | ρ2 |
| 3 | ρu | u |
| 4 | ρv | v |
| 5 | α1ρ1E1 (phasic **total** energy, E_k = e_k + ½\|u\|²) | α1p1 (α-weighted!) |
| 6 | α2ρ2E2 | α2p2 |

- α2 = 1 − α1, and ρ = α1ρ1 + α2ρ2.
- Conversions are in `ini_con2pri` (l.294–322) and `ini_pri2con` (l.324–342).
- The primitive pressure slots hold **α_k p_k**, not p_k. That keeps them well defined as α_k → 0.
- `ini_con2pri` has a vanishing-phase guard (l.308–311). If α_k < 1e-4, it sets α_k p_k = α_k · p_other.

### 1.2 EOS: stiffened gas (`StiffenedGasEOS`, l.187–254)

`p = (γ−1)ρe − γπ`. In α-weighted form (l.207–222):

```cpp
// α_k p_k from α_kρ_kE_k (subtract phasic kinetic energy with the shared velocity)
alphaRhoe = alphaRhoE - 0.5*alphaRho*(u*u+v*v);
return alphaRhoe*(gamma-1) - alpha*gamma*pInf;          // l.209-211
```

- Phasic sound speed: c_k² = γ_k(p_k+π_k)/ρ_k (l.234).
- The magnetosonic `calcf_x/y` (l.238–246) is left over from an MHD code and is unused.
- **Frozen mixture sound speed**, used for dt and the shock detector:
  - `calMixSoundSpeed` (l.757–784) computes c_f² = Σ Y_k c_k².
  - `cellSound` (l.884–896) computes the equivalent ρc_f² = Σ γ_k(α_k p_k + α_k π_k). It divides only by ρ and never by α.

### 1.3 PDE (as implied by the split)

```
∂t α1 + u·∇α1 = μ (p1 − p2)                              (μ → ∞: instantaneous relaxation)
∂t(α_kρ_k) + ∇·(α_kρ_k u) = 0
∂t(ρu) + ∇·(ρu⊗u) + ∇p = 0,           p = α1p1 + α2p2
∂t(α_kρ_kE_k) + ∇·(α_kρ_kE_k u + α_k p_k u) ± Σ = ∓ μ p_I (p1 − p2)
Σ = −u·( Y2 ∇(α1p1) − Y1 ∇(α2p2) )                         (Pelanti–Shyue 2014)
```

The only non-conservative terms are:

- `u·∇α` (explicit, treated with the path);
- Σ (implicit, centred differences inside the pressure system, l.1066–1068).

There is no p_I∇α term in the momentum equation, because the model has a single velocity.

### 1.4 How the terms are split

| Term | Where it is computed | Treatment |
|---|---|---|
| u·∇α, ∇·(α_kρ_k u), ∇·(ρu⊗u), ∇·(α_kρ_k K u) with K = ½\|u\|² | `buildA_advX/Y` (l.374–410), `Path_RusanovC_*` | explicit, path-conservative Rusanov |
| ∇p (mixture) | `PressureSubsystem` momentum update (l.1416–1447) | implicit, central |
| ∇·(α_kρ_k h_k u), with α_kρ_k(e_k + p_k/ρ_k) = Y_k h_k ρu | pressure system RHS + Laplacian (l.1400–1402, 1123–1124) | implicit |
| Σ | `sigma` in `Pressure_Subsystem_AX` (l.1066–1068) | implicit, linearised |
| μ(p1−p2) | `PressureRelaxation_Exact` (l.1525–1583) | instantaneous |

The explicit energy rows carry **only the phasic kinetic-energy flux**, Y_k ρ K u = α_kρ_k K u. The comment at l.538 says the same: "能量行已是动能通量…enthalpy 留给 implicit 步" ("the energy rows are already the kinetic-energy flux; the enthalpy is left to the implicit step").

---

## 2. Path-conservative Rusanov, as implemented

### 2.1 Quasi-linear matrix `A_x(Q)` (`buildA_advX`, l.374–391)

This is the Jacobian of the convective flux `F = (uα1?, m1 u, m2 u, ρu², ρuv, m1 K u, m2 K u)`. The α row is the non-conservative `u ∂xα`. With Y_k = m_k/ρ and K = ½(u²+v²):

```cpp
A[0][0]=u;                                                         // u ∂x α1 (non-conservative)
A[1][1]=u*(1-Y1); A[1][2]=-u*Y1;      A[1][3]=Y1;                  // ∂(m1 u)/∂Q
A[2][1]=-u*Y2;    A[2][2]=u*(1-Y2);   A[2][3]=Y2;
A[3][1]=-u*u;     A[3][2]=-u*u;       A[3][3]=2u;                  // ∂(ρu²)/∂Q
A[4][1]=-u*v;     A[4][2]=-u*v;       A[4][3]=v;  A[4][4]=u;       // ∂(ρuv)/∂Q
A[5][1]=K*u*(1-3Y1); A[5][2]=-3Y1*K*u; A[5][3]=Y1*(u*u+K); A[5][4]=Y1*u*v;  // ∂(m1 K u)/∂Q
A[6][...] analogous with Y2;  A[5][5]=A[6][6]=0
```

I checked these entries by differentiating F by hand, and they are correct. `buildA_advY` (l.394–410) is the same matrix with u and v swapped.

### 2.2 Path integral: linear path in **conserved** variables, 3-point Gauss–Legendre (l.413–441)

```cpp
s = {0.5 - sqrt(15)/10, 0.5, 0.5 + sqrt(15)/10};  w = {5/18, 8/18, 5/18};
Qs = QL + s_g (QR - QL);   Ahat += w_g * A(Qs);           // l.415-424
```

So `Â = ∫₀¹ A(QL + s(QR−QL)) ds ≈ Σ w_g A(Q(s_g))`.

**How ū is chosen.** The code never selects an explicit face velocity. In the α row, ū = Σ_g w_g · (ρu)(s_g)/ρ(s_g), i.e. the Gauss average of u along the linear path in (ρu, ρ). This is **not** ½(uL+uR) unless ρL = ρR. If you want ū = ½(uL+uR) exactly, use a linear path in primitive variables (u linear in s). You can also apply the path only to the α row, where ∫u ds = ū then holds exactly.

### 2.3 Fluctuations (`Path_RusanovC_x`, l.444–460)

```cpp
double S = fmax(fabs(vxL), fabs(vxR));
S = vx_max;                                    // l.453: OVERWRITTEN by global max |u| (from caldt)
Dm[r] = 0.5*(Ahat*dQ)[r] - 0.5*S*dQ[r];        // D^-_{i+1/2}
Dp[r] = 0.5*(Ahat*dQ)[r] + 0.5*S*dQ[r];        // D^+_{i+1/2}
```

`D^∓ = ½ Â ΔQ ∓ ½ S ΔQ`, with S = max over the domain of |u| (or of |v| in y).

- S is a **material-speed** (Lax–Friedrichs-type) dissipation. There is no sound speed in it, because the acoustics are implicit.
- The local value `max(|uL|,|uR|)` is computed and then thrown away.
- `vx_max` is a global value set in `caldt` from `u_pri` at time t^n. It is **not** updated per stage.

### 2.4 Cell update with the cell-interior term for second order (`Path_RusanovC_Update`, l.482–546)

```
L_i = −(1/Δx)[ D^-_{i+1/2} + D^+_{i−1/2} + Â(w^L_i → w^R_i)(w^R_i − w^L_i) ]   (+ same in y)
Q_i = Qs_i + adt · L_i ;   α clamped to [ALPHA_MIN, 1−ALPHA_MIN]
```

- Right face (l.495–498): `QL = w_xR(i)` and `QR = w_xL(i+1)`. The cell receives `−Dm/dx`.
- Left face (l.500–503): `QL = w_xR(i−1)` and `QR = w_xL(i)`. The cell receives `−Dp/dx`.
- **In-cell term** `B` (l.515–535): the path average `Â` from the reconstructed left edge to the right edge of the same cell, times `(wR−wL)`. This is the standard Castro/Parés/Dumbser smooth-part term. It is zero at first order.
- The y faces and y in-cell term (l.506–513, 527–535) contribute nothing in 1D. The ghost rows are transmissive copies, so ΔQ = 0, and `vy_max = 0`.
- The update writes `u_con = Qs + adt·L` (l.539–540). Here `adt = A_I(s,s)·dt`, not an explicit tableau weight (see §3.3).

---

## 3. Reconstruction, limiter, time stepping, CFL

### 3.1 Reconstruction (`TVDReconstruction`, l.549–612)

- The reconstruction is done in **primitive** variables (α, ρ1, ρ2, u, v, α1p1, α2p2). The goal is to keep p and u continuous across material interfaces (l.552–553).
- The limiter is **minmod**: `w^{L/R} = w ∓ ½ minmod(Δ⁻, Δ⁺)` (l.575–586).
- The edge states are mapped back to conserved variables with `ini_pri2con` (l.594–609).
- `Reconstruction recon = First_Order` (l.45) is the active setting. `Second_Order` is a switch.
- The input is `Qe` (the explicit stage state).

### 3.2 CFL / dt (`caldt`, l.786–824)

```cpp
dt = CFL * fmin(dx,dy) / fmax(a_max_x, a_max_y);   // CFL = 0.3 (l.44)
```

- `Fullydt`: a = |u|_max + c_f,max.
- `IMEXdt` (active, l.31): uses |u|+c_f for the **first 5 steps** only, then a = |u|_max. This is a material-wave CFL only.
- The comment at l.44 says second-order MUSCL with IMEX needs CFL ≲ 0.4. Otherwise it overshoots at shocks and contacts.

### 3.3 IMEX Runge–Kutta (main, l.1650–1765)

| `RKs` | Explicit tableau A_E | Implicit tableau A_I | b |
|---|---|---|---|
| **1 (active, l.43)** | 0 | 1 | 1 (forward Euler for the convective part, backward Euler for the pressure part) |
| 2 | Ae(1,0) = 1/(2β), β = 1−1/√2 | [[β,0],[1−β,β]] | (1−β, β) |
| 3 | [[0],[½],[½,½]] | diag(¼,¼); last row (⅓,⅓,⅓) | (⅓,⅓,⅓) |

Each stage runs these steps in order:

1. `Qe = Qⁿ + dt Σ_r Ae(s,r) H_r` and `Qs = Qⁿ + dt Σ_r Ai(s,r) H_r` (l.1729–1738).
2. Reconstruct `Qe`.
3. Explicit update `u = Qs + adt·L(Qe)`.
4. Pressure subsystem (h = adt).
5. Exact relaxation.
6. `H_s = (u − Qs)/adt` (l.1760).

**Caveat.** H_s lumps the explicit, implicit and relaxation increments into one stage derivative, which is then reused with *both* A_E and A_I. `bRK` is never used: the last stage value is taken as Qⁿ⁺¹. This is a consistent IMEX-RK only for `RKs = 1`. For RKs = 2 and 3 it is not the textbook IMEX scheme: the explicit last row ≠ b, and there are no separate explicit and implicit stage derivatives. Do not reuse it as is for higher order.

---

## 4. Implicit pressure subsystem (context for the explicit part)

This is `PressureSubsystem` (l.1328–1517) together with `Pressure_Subsystem_AX` (l.1022–1160).

**Unknowns.** q_k = α_k p_k per cell, stored interleaved (`M_P`, l.975–1020). The BC is mapped by index aliasing.

**Operator** (matrix-free, l.1123–1124):

```
q_k/(γ_k−1) − h² ∇·( Y_k h_k ∇p ) ± h Σ + [add_Evisc] −h/(γ_k−1) ∇·(κ ∇q_k)  =  RHS_k
p = q_1 + q_2 ;  Y_k h_k face-averaged ;  compact 5-point Laplacian
Σ = −Y2(u*·∇q1 − (h/ρ)∇p·∇q1ⁿ) + Y1(u*·∇q2 − (h/ρ)∇p·∇q2ⁿ)   (central 2Δx, l.1053–1068)
```

**RHS** (l.1400–1402):

```
RHS_k = α_kρ_kE_k* − α_kγ_kπ_k/(γ_k−1) − m_k K* − h ∇·( Y_k h_k ρu* )   (central face averages)
```

**Shock-activated dissipation** (`add_Evisc = true`, l.47; implemented at l.1128–1157):

- κ = χ · |u|_max · Δx.
- χ comes from the Balsara-type detector `ShockDetector` (l.899–938).
- Formula: χ̃ = clamp(−(δ∇·u + k1 a_min)/(k1 a_min), 0, 1), with k1 = 1e-3. χ is then the 3×3 max of χ̃.

**Solver.** `GMRES_Solver` (l.1260–1326) is a restarted GMRES(m) with m = 400, at most 2000 outer iterations and an **absolute** tolerance of 1e-8. It has no preconditioner. The initial guess is the q from `Qe`.

**Two Picard passes:**

1. **Pass 1** (l.1413). After the solve, momentum is updated with the central mixture-pressure gradient: `ρu −= h/Δx·(½(p_{i+1}+p_i) − ½(p_i+p_{i−1}))` (l.1416–1447).
2. **Pass 2** (l.1452–1503). The kinetic energy is recomputed from the new momentum, the system is solved again, and the phasic energies are rebuilt from q_k (l.1505–1516).

**Momentum is not updated after pass 2.** The final ρu therefore uses p⁽¹⁾, while the energies use p⁽²⁾. In addition, the energy Laplacian is compact while the momentum gradient uses the wide 2Δx central stencil. These two discrete pressure-work forms do not match, so total energy is not exactly consistent. Check this before you rely on energy conservation.

---

## 5. Relaxation procedures

Only **mechanical (pressure) relaxation** is implemented. There is **no thermal or chemical (T/g) relaxation**.

### 5.1 `PressureRelaxation_Exact` (active, l.1525–1583): instantaneous, Pelanti–Shyue App. C.1

The procedure works cell by cell. It skips near-pure cells (α_k < 10·ALPHA_MIN). Phasic masses, momentum and the mixture internal energy are frozen during the step.

```
C_k = 2γ_kπ_k + (γ_k−1)p_k
A = 1 + γ2α1 + γ1α2
B = C1α2 + C2α1 − (1+γ2)α1p1 − (1+γ1)α2p2
D = −(C2α1p1 + C1α2p2)
p*₀ = (−B + √max(B²−4AD,0)) / (2A)                       // closed form, no iteration
α1* = α1((γ1−1)p*₀ + 2p1 + C1) / ((γ1+1)p*₀ + C1)          // clamped
p*  = (ρe_mix − Σα_k*γ_kπ_k/(γ_k−1)) / Σα_k*/(γ_k−1)       // recomputed so internal energy is conserved exactly
α_kρ_ke_k* = α_k*(p* + γ_kπ_k)/(γ_k−1)
```

The root solver is the analytic quadratic formula. Recomputing p* from the mixture EOS at the end is a good idea: it gives exact energy conservation even though α1* was clamped.

### 5.2 `PressureRelaxation_Boscheri` (l.1171–1251; finite rate μ, **commented out** at l.1751)

The scheme is a semi-implicit linearisation with one division per cell:

- p_I = (Z2p1 + Z1p2)/(Z1+Z2), with Z_k = ρ_k c_k.
- K_k = −(1/α_k)[γ_k(p_k+π_k) + (γ_k−1)(p_I − p_k)].
- Δα = dt·μ(p1−p2)/(1 − dt·μ(K1+K2)). The division is applied only if the denominator is > 0.5.
- Δα is limited by |Δα_eq|, where Δα_eq = (p1−p2)/(−(K1+K2)).
- Energies are updated along `de = −p_I dv`, using the same p_I for both phases so the total is conserved.
- `mu = 1e30` (l.51), and the routine uses the global `dt`.

### 5.3 Vanishing-phase handling

- `reset_vanishing_phase` (l.1590–1623) runs at the start of every step on all cells including ghosts. If α_k < 1e-4, it sets ρ_k := ρ_other and p_k := p_other, then rewrites m_k and α_kρ_kE_k. This **changes the trace-phase mass**, so it is not conservative, though the effect is small.
- After the explicit update, α is clamped to [1e-8, 1−1e-8] (l.543).

---

## 6. Test case and switches

### 6.1 `testshockbubble` (the only test; l.622–643 and l.1628–1635)

This is a 1D air-shock / helium-bubble problem.

- Domain [0,1], `nxCells = 1000`, `nyCells = 1`, 3 ghost cells, `t_end = 0.0014`.
- γ1 = 1.4 (air), γ2 = 1.67 (He), π1 = π2 = 0. Units are SI-like (p ~ 1e5).

| region | α1 | ρ1 | ρ2 | u | p |
|---|---|---|---|---|---|
| x ≤ 0.05 (post-shock air) | 1−1e-6 | 1.3333 | 1e3 (dummy) | 0.3535·√1e5 ≈ 111.8 | 1.5e5 |
| 0.05 < x ≤ 0.4 | 1−1e-6 | 1.0 | 1e3 | 0 | 1e5 |
| 0.4 < x ≤ 0.6 (He) | 1e-6 | 1e3 (dummy) | 0.1379 | 0 | 1e5 |
| x > 0.6 | 1−1e-6 | 1.0 | 1e3 | 0 | 1e5 |

- The post-shock state satisfies Rankine–Hugoniot for p2/p1 = 1.5 with γ = 1.4. I checked this: ρ ratio 4/3, u = √(0.5e5·0.25).
- The dummy trace densities are immediately overwritten by `reset_vanishing_phase`.
- A diagnostic dump goes to stderr every 3e-4 s for x ∈ [0.40, 0.58] (l.1772–1785).
- Output files go to hard-coded `/Users/chenzefeng/Desktop/*.dat` (l.1790–1883).

### 6.2 Switches (l.19–63)

- **Active:**
  - `caltimestep = IMEXdt`
  - `x_bc = y_bc = Transmissive`
  - `RKs = 1`
  - `CFL = 0.3`
  - `recon = First_Order`
  - `add_Evisc = true`
  - `ALPHA_MIN = 1e-8`
  - GMRES: tol 1e-8 (absolute), 2000 iterations, restart 400
- **Declared but unused or dead:**
  - `formula`
  - `PressureSolver` (the code always does two passes)
  - `add_shock_detector`, `shock_kappa1`, `shock_threshold`, `threshold_Left/Right`
  - `add_IEvisc`, `add_restart`
  - `Bx_p0`, `By_0p`, `calcf_*`, the `triplets_p` sparse matrix, `calkinetic`

---

## 7. Reuse and avoid

### Reuse

- The **TV-type split** with only kinetic-energy fluxes in the explicit energy rows. The explicit CFL then depends on |u| only, which is the right structure for the IMEX target.
- `buildA_advX/Y` (l.374–410). The quasi-linear matrix is checked and correct.
- The **fluctuation form** `D^∓ = ½ÂΔQ ∓ ½SΔQ`, together with the in-cell term `Â(wL→wR)(wR−wL)` for second order (l.444–460, 515–535).
- Reconstruction in primitive variables with α_k p_k slots, plus the vanishing-phase guard in `ini_con2pri`.
- `PressureRelaxation_Exact`: closed form, robust, and it re-imposes exact energy conservation.
- `cellSound`: the form without division by α.

### Avoid or fix

1. **Conservation loss.** All rows, including mass and momentum, go through the 3-point Gauss path integral of a *rational* Jacobian. As a result `ÂΔQ ≠ F(QR) − F(QL)`, and mass and momentum are not exactly conserved at strong jumps. In this test ρ jumps by about 7×. Fix: use conservative Rusanov fluxes for rows 1–6 and keep the path only for the α row (and any genuinely non-conservative terms). Equivalently, use `ÂΔQ := F(QR)−F(QL)` for those rows.
2. **ū is not the face velocity.** It is a Gauss average of ρu/ρ along a linear path in conserved variables. If you want ū = ½(uL+uR), use a primitive-variable linear path, or set the α row to ū Δα directly.
3. **Global dissipation `S = vx_max`** (l.453). It is very diffusive away from the maximum and is frozen at tⁿ for all stages. Consider local S = max(|uL|,|uR|), which the code computes and then overwrites.
4. **IMEX-RK for RKs > 1** lumps the stage derivatives, `bRK` is unused, and the explicit residual is scaled by A_I(s,s). Only RKs = 1 is consistent. Write a proper ARS/SSP IMEX tableau with separate H^E and H^I.
5. **Pressure system.**
   - Momentum is not updated after the second Picard pass.
   - The Laplacian (compact) and the momentum gradient (2Δx) are inconsistent.
   - The GMRES tolerance is absolute while the residual is O(1e5).
   - There is no preconditioner, and the restart dimension of 400 is costly in 2D.
6. **Out-of-bounds write.** The reflective BCs write component index 7 (l.686–687, 738–739) although `NUM_VARS = 7`. They are unused now but will corrupt memory if enabled.
7. **Garbage output.** `shock_eta` in `main` is never filled, but it is written to `shock_de.dat` (l.1653, 1790–1798).
8. **Overwritten boundary edges.** `bound_cond_u` on the edge arrays (l.610–611) overwrites the boundary-adjacent ghost edges that were just reconstructed. For example, at second order, `w_xR(gc−1)` is replaced by `w_xR(gc)`. This is harmless at first order or with transmissive BCs, but wrong for reflective BCs.
9. **Performance and portability.**
   - `TVDReconstruction` allocates 5 Grids per call.
   - `Grid` has no copy or move semantics.
   - Output paths are hard-coded.
   - `calkinetic` reads `con[5]` as vz (l.968).
