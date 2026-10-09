# Code note: `6-equation-diffuse-IMEX-WB.cpp` (Pelanti–Shyue 6-eq, IMEX, well-balanced)

- **Primary file (analysed line by line):** `696c4f8d-6-equation-diffuse-IMEX-WB.cpp`, 2597 lines. In this note, `L123` means a line in this file.
- **Variant (compared structurally):** `de47bc68-6-equation-diffuse-IMEX_WB.cpp`, 3501 lines. `O:123` means a line in that file. See §11.
- **Model:** Pelanti & Shyue, *A mixture-energy-consistent six-equation two-phase numerical model…* (JCP 2014). The equation numbers (6a)–(6g) refer to the reading note `docs/papers/paper_6.md`.
- **Numerics borrowed from:** Battisti–Boscheri / Boscheri-style pressure-based semi-implicit IMEX. The code cites "B-B (12d),(18),(32b),(39)/(40)", "Boscheri (18)" and "eq 52–54" for the Balsara-type shock detector and the χ-Laplacian. Those papers are not in `docs/papers`, so the numbering is the code's own.
- **Dependencies:** a single translation unit with Eigen (Dense, SparseLU, IterativeLinearSolvers, unsupported/IterativeSolvers). Output goes to hard-coded `/Users/chenzefeng/Desktop/*.dat` (L2504–2597).

---

## 1. State vector, EOS, model, layout

### 1.1 Model solved
This is the single-velocity 6-equation model, with phasic **total** energies (P-S (6a)–(6g)) and 2D Cartesian geometry:

```
∂t α1 + u·∇α1                                   = μ(p1-p2)
∂t (α_k ρ_k) + ∇·(α_k ρ_k u)                   = 0
∂t (ρu) + ∇·(ρu⊗u) + ∇(α1p1+α2p2)              = 0
∂t (α_k E_k) + ∇·(α_k(E_k+p_k)u) ± Σ           = ∓ μ p_I (p1-p2)
Σ = -u·( Y2 ∇(α1p1) - Y1 ∇(α2p2) )             (+Σ for k=1, -Σ for k=2)
```

Relaxation is mechanical only (§5). There is no heat or mass transfer in this file.

### 1.2 Storage (`NUM_VARS = 7`, L19)

| idx | conservative `u_con` | primitive `u_pri` (`ini_con2pri` L476) |
|---|---|---|
| 0 | α1 | α1 |
| 1 | α1ρ1 | ρ1 |
| 2 | α2ρ2 | ρ2 |
| 3 | ρu | u |
| 4 | ρv | v |
| 5 | α1E1 (phasic total energy, including the phase's kinetic energy α1ρ1·½\|u\|²) | **α1p1** |
| 6 | α2E2 | **α2p2** |

The primitive "pressure" slots hold α_k p_k, not p_k. So the mixture pressure is `pri[5]+pri[6]`, and p_k = pri[5+k]/α_k.

### 1.3 Data structures
- `Grid` (L107) stores a flat `double[(nx+2g)(ny+2g)·7]`, indexed `(i,j,k) → 7·(i + j·(nx+2g)) + k`, so `i` is the fast index. There are 3 ghost cells (L41).
- `Array` (L135) is a scalar field of the same shape. Note that it is **not zero-initialised**: ghost values are garbage until `bound_cond_array` runs. The variant fixes this at O:518.
- `Tensor` (L78) adds an RK-stage index and holds the stage increments `Hs(r,i,j,k)`.
- `M_P(i,j,k)` (L1353) maps a cell to an index in the implicit unknown vector. The layout is **interleaved**, `2·cell + k` with k=0 for phase 1 and k=1 for phase 2. A ghost index is mapped back onto an interior cell according to the BC (transmissive clamps, reflective mirrors, periodic wraps). This is how the implicit stencil handles boundaries.
- 1D runs use `nyCells = 1`. The y-loops still execute, and `dy = y-range` (for example 1.0), so the y-terms are negligible but are not removed.

### 1.4 EOS hierarchy (L163–433)
`EquationOfState` is the abstract interface. `MieGruneisenEOS` implements the generic `p = p_ref(ρ) + Γ(ρ)ρ(e − e_ref(ρ))`. The concrete classes are `StiffenedGasEOS(γ,π)` with `p_ref=-γπ, e_ref=0, Γ=γ-1`, `IdealGasEOS`, `JWLEOS` and `CochranChanEOS`. That last class re-implements the MG pieces directly instead of deriving from `MieGruneisenEOS`.

The implicit system depends on a single hook, `e_linear_coeff(ρ, c, d)` with `e = c·p + d` at fixed ρ (L193). `get_phase_linear_coeffs` (L1400) turns it into the per-volume form

```
α_k ρ_k e_k = A_k (α_k p_k) + α_k B_k ,   A_k = ρc = 1/Γ_k(ρ_k),   B_k = ρd = ρ(e_ref - p_ref/(Γρ))
SG:  A_k = 1/(γ_k-1),  B_k = γ_k π_k/(γ_k-1)
```

The relation is linear in q_k = α_k p_k when ρ_k is frozen. **That is the entire linearisation of the pressure system** (§3.3).

Sound speeds come from three places:
- `calcs` per phase.
- `calMixSoundSpeed` and `cellSound` (L1078, L1208) both use the **Wood** formula. The comment at L1077 says "frozen", but the code computes Wood.
- The SG `calcs` floors `γ(p+π)` and `ρ` at 1e-12 for positivity (L243).

---

## 2. Explicit / implicit split

### 2.1 Explicit part: transport without pressure (path-conservative Rusanov)
The explicit operator is the quasi-linear advection matrix `A_x(Q)` (`buildA_advX` L552, `buildA_advY` L572) of the system **without any pressure terms**:
- row 0: `u` (α advection)
- rows 1, 2: d(m_k u)
- row 3: d(ρu²), i.e. `-u² dρ + 2u d(ρu)`, with no ∇p
- rows 5, 6: d(m_k K u) with K = ½|u|². The energy rows carry **only the phasic kinetic-energy flux**, and enthalpy transport is left to the implicit step (comment L749).

The scheme is path-conservative (DOT-type) with a linear path and a 3-point Gauss–Legendre average of `A` (`pathAvgAx` L591). Fluctuations at each interface are computed as follows (L637–662):

```cpp
Dm[r] = 0.5*Av - 0.5*S*kjmp[r];   Dp[r] = 0.5*Av + 0.5*S*kjmp[r];
// S = vx_max  (global max |u|, NOT u+c);  kjmp = ΔQ except rows 5,6 = Δ(m_k K)
```

There are three points to note here:
- **Dissipation speed:** the dissipation uses the global maximum material speed `vx_max` from `caldt`. The local `max(|uL|,|uR|)` is computed and then overwritten (L645–646). With a quiescent initial condition, the first step therefore has **zero** Rusanov dissipation.
- **Energy rows:** the dissipation acts only on the kinetic-energy jump (L647–655), so internal energy is never diffused explicitly. This is essential for well-balancing (§4).
- **Update:** `Path_RusanovC_Update` (L693) adds the interface fluctuations D⁻ and D⁺ and the in-cell term `∫A dQ` along `w_L(i) → w_R(i)`, then computes `u_con = Qs + adt·L(Qe)` and clamps α to [1e-8, 1-1e-8].

**Reconstruction** is done by `TVDReconstruction` (L760) in primitive variables (α, ρ1, ρ2, u, v, α1p1, α2p2):
- Each component gets a minmod slope.
- A Zhang–Shu-type **positivity scaling θ** (`posTheta`, L777) is applied, with one common θ per direction for all components. It enforces α∈[αmin, 1−αmin], ρ_k>1e-10 and α_k(p_k+π_k)>1e-6.
- The edge values are mapped back to conservative variables.

Because minmod is positively homogeneous and θ is shared, a uniform p gives `slope(α_k p_k) = p·slope(α_k)`. So p_k is preserved exactly at reconstructed faces, which is needed for WB.

### 2.2 Implicit part: pressure, enthalpy flux and Σ (`PressureSubsystem` L1866)
The implicit part treats three things:
- the momentum pressure gradient ∇(α1p1+α2p2);
- the enthalpy part of each phasic energy flux, ∇·(α_k ρ_k h_k u) with ρ_k h_k = ρ_k e_k + p_k;
- the non-conservative Σ.

The volume fraction, the phase masses and the convective momentum are **not** changed in the implicit step. They come only from the explicit step and from relaxation.

Notation:
- `u_con` on entry is the post-explicit state: α^{n+1}, m_k^{n+1}, (ρu)*, (α_kE_k)*.
- `Qe` and `Qs` are the RK explicit and implicit stage bases (§7).
- `h = adt = Ai(s,s)·dt`.

---

## 3. The pressure equation

### 3.1 Unknowns
The unknowns are `x = q_k = α_k p_k` for every cell, so `dim = 2·nx·ny`. The initial guess is `Qe_pri(…,5/6)` (L2031), which only matters for the iterative fallback.

### 3.2 Exact discrete form (`Pressure_Subsystem_AX` L1407–1616, RHS L2023–2025)
Write `P = q1 + q2` (mixture pressure) and s_1 = +1, s_2 = −1. For each interior cell (i,j) and phase k:

```
LHS_k = A_k q_k
        - h² [ (W_k,i+½ (P_{i+1}-P_i) - W_k,i-½ (P_i-P_{i-1}))/dx²  + (same in y)/dy² ]
        + s_k · h · σ
        - h · ∇_h·(κ ∇_h (A_k q_k))                       (shock χ-Laplacian, §6)

RHS_k = (α_kE_k)*  -  ½ |(ρu)*|²/ρ² · m_k  -  h · EnthDiv_k  -  α_k B_k
```

**Coefficients** (all frozen):
- `A_k` and `B_k` are evaluated from `ρ_k = m_k^{n+1}/α_k^{n+1}` (L1536).
- Elliptic weight `C_k = α_k^{n+1} (ρ_k h_k)^{Qe} / ρ^{n+1}` ("Boscheri (18)", L1514–1528), with `(ρ_k h_k)^{Qe} = ((α_kρ_ke_k)+(α_kp_k))/α_k` evaluated at `Qe`. The face value is the **arithmetic mean**, `W_k,i+½ = ½(C_k,i + C_k,i+1)` (L1530).
- The **Σ term** `σ` (L1441–1496) uses the expanded (6g) form, `Σ = -u·[(Y2p1+Y1p2)∇α1 + α1Y2∇p1 − α2Y1∇p2]`, rewritten through α_k∇p_k = ∇q_k − p_k∇α_k:

  ```cpp
  Aadv   = (Qs(i,j,0) - u_con(i,j,0)) / h;        // (u·∇α1)^PC  = explicit α increment
  Pw     = Y2*q1/α1 + Y1*q2/α2;                     // Y2 p1 + Y1 p2  (diagonal only)
  S_impl = u*·(Y2 ∇c q1 − Y1 ∇c q2) + Pw*(Aadv − u*·∇c α1);    // ∇c = wide central diff
  Bx_n   = Y2 ∂x(α1p1)^{Qe} − Y1 ∂x(α2p2)^{Qe}               // frozen bracket
  sigma  = −S_impl + (h/ρ) (∇c P · B_n);            // Picard: u^{n+1} = u* − (h/ρ)∇P
  ```

  Y_k, α_k and u* are taken from `u_con`. The sign convention (+σ in row 1, −σ in row 2) matches P-S, where Σ sits on the LHS of (6e).
- **Enthalpy flux** `EnthDiv_k` (L1927–1978) is in **conservative face-flux** form:

  ```
  EnthDiv_k = [Ĥ_k,i+½ Ψ_k,i+½ − Ĥ_k,i−½ Ψ_k,i−½]/dx + (y)
  Ĥ_k,i+½ = ½(H_k,i + H_k,i+1),   H_k = (ρ_k h_k)^{Qe}
  Ψ_1 = Φ + â1 (ū*_face − û),   Ψ_2 = (û − Φ) + (1−â1)(ū*_face − û)
  Φ   = ½ û (α1L+α1R) − ½ vx_max (α1R−α1L)      // SAME α Rusanov face flux as the explicit step
  û   = path-averaged velocity (pathAvgVel, L622), â1 = ½(α1L+α1R), ū*_face = ½(u*_i+u*_{i+1})
  ```

**Derivation of the h² term.** The implicit enthalpy flux needs u^{n+1} = u* − (h/ρ)∇P. Substituting gives `−h∇·(α_kρ_kh_k u*) + h²∇·(C_k∇P)`; the first part is EnthDiv and the second is moved to the LHS.

**Mixture-energy conservation.** Summing the two rows gives the following. The Laplacians are fluxes of (W1+W2)∇P, the ±σ cancel, EnthDiv is a face-flux difference, and the χ-Laplacian is a flux of ∇(A1q1+A2q2), because A_k sits **inside** the divergence (L1596–1605). So the discrete mixture total energy is conservative.

### 3.3 Linearisation
There is no Newton iteration. The operator is **linear in q** because:
- A_k, B_k, C_k, H_k, Y_k, u*, B_n and χ are all frozen at either `Qe` (the explicit stage state) or `u_con` (the post-explicit α^{n+1}, m^{n+1}, u*);
- `Pw` uses q_k/α_k on the diagonal only.

So the system is the single linear solve "A(frozen)·q = b", repeated twice (§3.5).

### 3.4 Linear solver: `PressureSolve_Direct` (L1787)
Despite the GMRES globals, **the matrix-free GMRES is never called**.
- `GMRES_Solver` (L1713) is a hand-written restarted GMRES(m) with modified Gram–Schmidt, Givens rotations and no preconditioner. Its settings are m = `GMRESrest` = 400, outer = `GMRESiter` = 2000 and tol = `GMREStol`·‖b‖ = 1e-8·‖b‖. It is dead code.

What actually runs:
1. **Assembly by coloured probing** of the matvec `Pressure_Subsystem_AX`. Probing is valid because the operator is linear with A(0)=0.
   - Colour = `p + 2·((i-g)%7) + 14·((j-g)%7)`, giving 98 colours (S=7). In 1D only 14 are non-empty.
   - For each probe, entries within radius D=2 cells are read back into triplets.
   - A one-time self-check compares `‖A·x_rand − matvec(x_rand)‖/‖matvec‖` and prints it to stderr.
2. **Symmetric Jacobi equilibration** `As = D^{-1/2} A D^{-1/2}` with `d_k = 1/sqrt|a_kk|`. This rescales the liquid/gas rows, whose magnitudes differ enormously.
3. `Eigen::SparseLU` on `As`, followed by un-scaling.
4. **Fallback** if LU fails or the result is non-finite: `Eigen::GMRES<SparseMatrix, DiagonalPreconditioner>` with restart 200, maxIter 5000, tol 1e-9, warm-started from the scaled guess.

Periodic BC caveat: `M_P` wraps ghost indices to the opposite boundary, but the probing read-back window is ±2 cells, so wrapped couplings are **not extracted**. With `x_Periodic` the assembled matrix is incomplete. Only the self-check message would reveal this.

### 3.5 Correction iterations ("no / once / all")
The enum `SolvePressure {no_correction, once_correction, all_correction}` and `PressureSolver = all_correction` (L26, L35) are **unused**. The code hard-wires two solves, which amounts to one Picard correction:

1. **Solve 1** (L1981–2036). The RHS kinetic energy uses u* (`rhovx_star`).
2. **Momentum update** (L2039–2070), using face pressure ½(P_i+P_{i±1}), which is the wide central gradient:

   ```cpp
   u_con(i,j,3) -= h/dx * (0.5*(P_i+P_{i+1}) - 0.5*(P_{i-1}+P_i));
   ```

3. **Solve 2** (L2075–2129). The RHS is rebuilt with kinetic energy from the updated momentum, ½|ρu^{(1)}|²/ρ²·m_k. EnthDiv still uses u*, which is correct because the h² term accounts for u*→u^{n+1}.
4. **The velocity is not updated again.** The final u^{n+1} uses P^{(1)} and the energies use P^{(2)}.

**Energy update after the solve** (L2131–2147). The energy is reconstructed from the solved pressure:

```cpp
u_con(i,j,5) = eos1.getAlphaE_fromAlphaP(q1, m1, u, v, α1);   // = A1 q1 + α1 B1 + m1 K
```

This is a "pressure → energy" update. Because the rows are flux-form energy equations, it is still conservative for the mixture total energy, as shown in §3.2.

**Volume fractions** are not touched by the pressure solve. α changes only in the explicit step and in relaxation (§5).

**Inconsistencies worth knowing:**
- In solve 2, `S_impl` reads `vxs = u_con(3)/ρ` (L1457), which is now u^{(1)}, while the Picard term `(h/ρ)∇P·B_n` is still added. The Σ velocity correction is therefore counted roughly twice in the second pass.
- The energy Laplacian is **compact** (P_{i+1}−P_i), but the momentum gradient is the **wide** central (P_{i+1}−P_{i−1})/2dx. The velocity implied by the energy equation differs from the actual u^{n+1}, and the momentum is open to odd–even decoupling. A Rhie–Chow-type issue is also flagged in `jcp_2025_semi_implicit_6eq.md` L292.

---

## 4. The well-balanced (WB) mechanism

**Target state:** uniform p and u, with α (and ρ_k) jumping across a contact. The goal is to keep p^{n+1} = p and u^{n+1} = u exactly. Each piece below contributes.

1. **Reconstruction.** Primitive minmod with a shared θ keeps the face values of u and p_k = (α_kp_k)/α_k exactly uniform (§2.1).
2. **Explicit step.**
   - The momentum row of `A` is `-u²dρ + 2u d(ρu)`. With uniform u, `d(ρu) = u dρ`, so the momentum fluctuation equals u × (mass fluctuation) and **u\* = u**.
   - Energy rows 5 and 6 carry m_kK only, with dissipation on the m_kK jump only, so internal energy is not perturbed.
   - Rusanov dissipation on rows 1, 2 and 3 is consistent for the same reason as the momentum row.
3. **Σ uses the explicit α operator.** The α-gradient part of Σ is not re-discretised. It is **the explicit increment itself**, `Aadv = (α1^{Qs} − α1^{n+1})/h` (L1460), so it inherits the reconstruction, path, Rusanov dissipation and clamp bit for bit.
4. **Enthalpy flux uses the same α face flux.** The conservative face flux uses `Φ` = the α Rusanov flux of the explicit step: same states, same û, same S. When û is uniform, `(Φ_{i+½}−Φ_{i−½})/dx ≡ (α^{Qs}−α^{n+1})/h` (comment L1923–1926). The velocity correction `â(ū*−û)` is zero there because u* = û.
5. **Algebraic cancellation at the contact (SG).**
   - Uniform p means ∇P = 0, so the h² and Picard terms vanish. Also `u*·(Y2∇q1−Y1∇q2) = u*(Y2p+Y1p)∇α1`, which cancels against `−Pw·u*·∇α1`, leaving `σ = −p·Aadv`.
   - With Δα = α^{n+1} − α^{Qs}, H = γ(p+π)/(γ−1) and A q + αB = α(p+γπ)/(γ−1), row k reads

     ```
     α^{n+1}(p+γπ)/(γ−1) + pΔα  =  α^{Qs}(p+γπ)/(γ−1) + γ(p+π)/(γ−1)·Δα
     ⇔ Δα[(p+γπ)/(γ−1) + p − γ(p+π)/(γ−1)] = 0  ✓   (identically)
     ```

     So q_k = α_k^{n+1}p solves the system and p^{n+1} = p.
   - This needs H_k = ρ_kh_k to be uniform across the contact. For SG, ρ_kh_k depends only on p, so this holds. For general MG EOS (JWL, NASG), where ρh depends on ρ_k, WB is **only approximate**.
6. **The χ-Laplacian must be off at contacts.** `∇(A_kq_k) = ∇(α_kp/(γ_k−1)) ≠ 0` at a contact. WB therefore relies on χ = 0 there.
   - At a contact ∇·u = 0, which gives `χ̃ = clamp(−1,0,1) = 0`.
   - The pressure-driven propagation step requires a strict P(i) > P(i±1), which never holds for uniform P.
   - The variant documents that an always-on, unweighted version breaks WB (O:2322–2328).
7. **Relaxation.** p1 = p2 already holds, so `da = 0`.

Two earlier approaches ("Rcorr" and other mismatch corrections) were removed. The comments at L1446–1449 and L1915–1916 explain that WB now follows from the discretisation itself.

---

## 5. Relaxation procedures

Only the **mechanical (pressure)** relaxation is active in this file. **pT and pTg (thermal and chemical) relaxations do not exist here**; they live only in the variant (§11).

### 5.1 `PressureRelaxation_Boscheri` (L1627): the active one
This is called **inside every RK stage**, after the pressure solve (L2480), with `mu = 1e60` (L50). It is cell-local, takes one linearised step and has no root solver.

It skips near-pure cells (α_k < 1e-7). Otherwise:

```
p_I = (Z2 p1 + Z1 p2)/(Z1+Z2),  Z_k = ρ_k c_k
K_k = −(1/α_k)[ ρ_k c_k² + Γ_k (p_I − p_k) ],  Γ_k = 1/(c_lin ρ_k)
denom = 1 − dt·μ(K1+K2);   da = dt·μ(p1−p2);  if denom>0.5: da /= denom
da_eq = (p1−p2)/max(−(K1+K2),1e-12);  if |da|>|da_eq|: da = da_eq     // never overshoot
α1 ← clamp(α1+da);   ρ_k = m_k/α_k;   e_k ← e_k − p_I (1/ρ_k^new − 1/ρ_k^old)
```

- **Step size:** with μ = 1e60 the limiter always binds, so this is **one Newton step toward p1=p2** with frozen slopes, using the full dt rather than adt. Equilibrium is only reached to linearisation accuracy.
- **Conservation:** the mixture internal energy is conserved exactly, because Σ m_k Δ(1/ρ_k) = Δα1 + Δα2 = 0. Masses and momentum are untouched.
- **Interaction with the RK stages:** the relaxation increment is included in `Hs` and so is re-weighted by the Butcher coefficients. The variant's comments (O:3248–3262) say exactly this "relaxation inside Hs" design caused pumping (pmax 1e11 → 1e14) and replaced it with Strang half-steps.

### 5.2 `PressureRelaxation_Exact` (L2211): defined, never called
This is Pelanti–Shyue App. C.1 for SG, an instantaneous p1=p2=p\*.
1. Solve the quadratic `A p*² + B p* + D = 0` with `C_k = 2γ_kπ_k + (γ_k−1)p_k`, `A = 1+γ2α1+γ1α2`, `B = C1α2 + C2α1 − (1+γ2)α1p1 − (1+γ1)α2p2` and `D = −(C2α1p1 + C1α2p2)`. Take the positive root.
2. Compute `α1* = α1((γ1−1)p* + 2p1 + C1)/((γ1+1)p* + C1)`.
3. Recompute p\* from the mixture EOS `Eint = p*Σα_k/(γ_k−1) + Σα_kγ_kπ_k/(γ_k−1)`, so energy is exactly consistent, then redistribute the phasic energies.

It is SG-only and skips cells with α_k < 1e-4.

### 5.3 Other per-stage fix-ups (main loop L2442–2482)
- `clampQ` on `Qe` and `Qs`: α is clamped to [1e-8, 1−1e-8], and m_k is set to `≥ 1e-8·ρ` if non-finite or too small. The variant calls this mass floor wrong (O:3300–3307), because it can push ρ_k past the EOS limits.
- **The same clamp on `u_con` after the explicit step** (L2462).
- `reset_vanishing_phase` (L2152): where α_k < 1e-4, the minor-phase pressure is set equal to the dominant-phase pressure while conserving mixture internal energy. Its formula is **SG-only**, written with γ and π.
- `clamp_pressure_floor` (L2190): enforces p_k+π_k ≥ max(1e2, 1e-6π_k). It is **not** conservative.

---

## 6. Shock detector and implicit artificial viscosity (χ-Laplacian, "eq 52–54")

### 6.1 `ShockDetector(chi, Qe)` (L1236)
It is called once per stage inside `PressureSubsystem` (L1882) on the explicit-stage state `Qe`.

1. Compute U, V and the Wood sound speed `cs` on all cells including ghosts.
2. Compute the raw indicator `χ̃ = clamp( −(δ ∇·u + k1 a_min)/(k1 a_min), 0, 1 )`, where:
   - δ = min(dx,dy);
   - ∇·u is a central difference;
   - a_min is the minimum `cs` over the von Neumann neighbourhood;
   - **k1 = 1e-3 is hard-coded** (L1240).

   χ̃ becomes positive when δ∇·u < −k1·a_min and saturates at 1 when δ∇·u ≤ −2k1·a_min. This is a very sensitive compression trigger.
3. Take the 3×3 max of χ̃.
4. **Propagate** χ toward lower mixture pressure `P = α1p1+α2p2`, with 2 sweeps per direction. A cell with χ=0 inherits its neighbour's χ if P_nb > P_cell.
   - Because the sweep runs i-ascending and writes in place, rightward/upward propagation can chain arbitrarily far in one sweep along a monotonically falling P.
   - Leftward/downward propagation moves at most one cell per sweep.
   - So the spreading is **asymmetric**.
5. Apply `bound_cond_array(chi)`.

### 6.2 Use in the pressure operator (L1542–1613)

```cpp
kx = chi(i,j) * vx_max * dx;    ky = chi(i,j) * vy_max * dy;      // κ, then face average ½(κ_i+κ_{i+1})
lap_k = Σ_faces κ_f (A_k q_k)_nb − (A_k q_k)_i  / dx²  (+y)
Ax(k) += −h * lap_k;           // i.e. LHS gets −h ∇·(κ ∇(A_k α_k p_k))
```

For SG, `A_k q_k = α_kp_k/(γ_k−1)` = α_kρ_ke_k − α_kB_k. So the operator is an **implicit diffusion of each phase's "pressure energy"**, not literally of the mixture pressure. The sum over phases is ∇·(κ∇Σα_kp_k/(γ_k−1)), which is conservative.

κ uses the **global** max material speed and **no ½ factor**. The variant instead uses an always-on `0.5·vx_max·dx` without χ.

### 6.3 Unused knobs
None of these is referenced anywhere:
- `add_Evisc = true` (L47; the comment says "eq 53-54");
- `add_shock_detector = No_Second`;
- `shock_kappa1 = 0.001`;
- `shock_threshold = 0.9`;
- `threshold_Left` and `threshold_Right`;
- the global `g_chi` (L74).

The χ-term is **always active**. The detector's k1 is hard-coded and happens to equal `shock_kappa1`. Some commented-out κ choices remain at L1557–1584: zero, always-on `vmax·dx`, and local `|u|dx`.

---

## 7. IMEX-RK tableau, stage structure and time step

### 7.1 Tableaux (main L2359–2409), `RKs = 2` by default (L43)

| RKs | Ae (explicit) | Ai (implicit) | b |
|---|---|---|---|
| 1 | [0] | [1] | [1] |
| 2 | [[0,0],[1/(2γ),0]] | [[γ,0],[1−γ,γ]] | (1−γ, γ) |
| 3 | [[0,0,0],[½,0,0],[½,½,0]] | [[¼,0,0],[0,¼,0],[⅓,⅓,⅓]] | (⅓,⅓,⅓) |

For RKs=2, γ = 1 − 1/√2. All tableaux are **stiffly accurate**: b equals the last row of Ai. The loop never uses `bRK`, and the solution after the last stage *is* u^{n+1}.

### 7.2 Stage structure (L2425–2494)
The structure is the "H(Q_E, Q_I)" form (Boscheri–Pareschi-style semi-implicit IMEX). A single stage derivative K_s = Hs(s) is stored and advanced with **both** tableaux:

```
Qe = uⁿ + dt Σ_{r<s} Ae(s,r) Hs(r)        // explicit-argument state (reconstruction, frozen coeffs, χ)
Qs = uⁿ + dt Σ_{r<s} Ai(s,r) Hs(r)        // implicit base
adt = Ai(s,s) dt
clampQ(Qe), clampQ(Qs)
TVDReconstruction(Qe) → Path_RusanovC_Update:  u_con = Qs + adt·L_E(Qe)
clamp α, m_k; reset_vanishing_phase
PressureSubsystem(u_con, Qs, Qe, …, h=adt)           // two linear solves (§3.5)
clamp_pressure_floor; reset_vanishing_phase
PressureRelaxation_Boscheri(u_con, mu)               // uses full dt, not adt
clamp_pressure_floor
Hs(s) = (u_con − Qs)/adt                              // includes explicit + implicit + relaxation + clamps
```

`vx_max` and `vy_max`, which set both the Rusanov S and the κ of the χ-Laplacian, are computed once per step in `caldt` from uⁿ. They are not updated per stage.

### 7.3 Time step, `caldt` (L1106)
- `Fullydt`: `a = |u|max + c_max`. Here c is the **Wood** speed, despite the "frozen" comment.
- `IMEXdt` (the default, L31): the acoustic CFL is used for `timestep < 5`, and the **material speed only** afterwards (`a = |u|max`).
- In both cases `dt = CFL·min(dx,dy)/max(a_x,a_y)` with **CFL = 0.1** (L44), then a growth cap `dt ≤ 1.2·dt_prev`.
- If |u|max = 0 after step 4, dt → ∞ and only the 1.2× cap bounds it.
- Final step: `dt = min(dt, time_stop − t)`.

---

## 8. Test cases, parameters and switches

### 8.1 Global configuration (L31–62)
- `test = testDodecaneShockTube`, `nx = 1000`, `ny = 1`, `ghost = 3`
- `RKs = 2`, `CFL = 0.1`, `Second_Order`, transmissive BCs
- `ALPHA_MIN = 1e-8`, `mu = 1e60`, `IMEXdt`

### 8.2 Test cases
In all cases the initial condition is written in primitive form (α1, ρ1, ρ2, u, v, α1p1, α2p2) (L843–969) and converted with `upri2ucon`.

| test | domain, t_end | EOS (phase1 / phase2) | IC |
|---|---|---|---|
| `testDodecaneShockTube` (L926, L2329) | [0,1], x_disc = 0.75, t = 473 µs | SG(1.025, 0) "vapour" / SG(2.35, 4e8) "liquid" | left: α1=1e-8, ρ1=2, ρ2=500, p=1e8; right: α1=1−1e-8, same ρ's, p=1e5; u=0 |
| `testshockbubble` (L851, L2306) | [0,1], t = 1.4 ms, 1D | SG(1.4,0) / SG(1.67,0) | shocked gas x≤0.05: ρ1=1.3333, u=0.3535√1e5 ≈ 111.8, p=1.5e5; light bubble 0.4<x≤0.6 (α1=1e-6, ρ2=0.1379); elsewhere ρ1=1, p=1e5 |
| `test2dahshockbubble` (L873, L2316) | [0,325]×[0,44.5] mm, t_end = 280 (non-dimensional) | SG(1.4,0) / SG(1.67,0) | air–helium; shock at x=50 (ρ=1.3764, u=0.394, p=1.5698); He bubble centre (175, 0), R=25, ρ2=0.138; pre-shock ρ=1, p=1 |

Notes on the test cases:
- The 2D bubble is centred on y=0 with y∈[0,44.5], so only a half-bubble is in the domain, and a reflective y BC would be needed. The defaults are transmissive and `nyCells = 1`, so this case **cannot actually run as 2D** without editing the constants.
- The local `ALPHA_MIN = 1e-6` shadows the global one inside that case.
- The dodecane EOS here is a simple SG pair. The variant uses the Le Métayer–Saurel water SG pair instead.

### 8.3 Enums and switches, live or dead

| symbol | line | status |
|---|---|---|
| `dtCalculation caltimestep` | 27/31 | live |
| `x_bc`, `y_bc` | 28–34 | live |
| `Test test` | 21/36 | live |
| `Reconstruction recon` | 22/45 | live |
| `Formula formula` | 20/32 | **unused** |
| `SolvePressure PressureSolver` | 26/35 | **unused** (two solves are hard-wired) |
| `shock_detector`, `shock_threshold_*`, `shock_kappa1`, `shock_threshold`, `add_Evisc` | 23–25, 47, 57–61 | **unused** |
| `GMREStol`, `GMRESiter`, `GMRESrest`, `add_restart` | 51–54 | used only by the dead `GMRES_Solver` |
| `Bx_p0`, `By_0p`, `B_max_*`, `ch_*`, `cp_*`, `ca_*`, `Bz_star` | various | leftovers from an MHD code; unused |

---

## 9. Bugs and hazards found

1. **Out-of-bounds write in reflective BCs.** `bound_cond_u` writes component `k=7` (L1007–1008, L1059–1060) while `NUM_VARS = 7`, which corrupts α of the neighbouring cell. It is latent because the tests use transmissive BCs. The variant has the same bug (O:1646).
2. **Garbage output.** In `main`, `shock_eta` (L2356) is never filled, but it is written to `shock_de.dat` from uninitialised memory, because `Array` is not zero-initialised.
3. **Asymmetric χ spreading.** The in-place sweep chains rightward but steps only one cell leftward per sweep (§6.1).
4. **Σ velocity correction double-counted in solve 2** (§3.5).
5. **Compact Laplacian vs wide-stencil momentum gradient** (§3.5): an odd–even decoupling risk.
6. **Periodic BC with direct assembly**: wrap-around couplings are not extracted (§3.4).
7. **SG-only routines inside a "generic EOS" code.** `reset_vanishing_phase` and `PressureRelaxation_Exact` use γ and π directly, and WB is exact only for SG (§4.5).
8. **First step has zero explicit dissipation** for quiescent initial conditions (S = vx_max = 0). After step 4 of `IMEXdt`, dt is set by |u|max alone and capped only by the 1.2× growth limiter.
9. **Minor code smells:**
   - `calkinetic` reads `con[5]` as vz but does not use it.
   - `calEnthalpy_phasic` results `h1_*` and `h2_*` are computed in `AX` and `PressureSubsystem` but never used (L1500–1513, L2012–2017).
   - The prototypes at L451 and L457 disagree with the definitions.
   - The comment at L463 describes an implicit α solve that does not exist.
10. **Hard-coded output paths** (`/Users/chenzefeng/Desktop/`).

---

## 10. Worth reusing / worth avoiding

**Reuse:**
- **The WB construction (§4):** Σ's α-gradient taken as the explicit PC increment, plus an enthalpy flux built from the *same* α Rusanov face flux, plus primitive (α_kp_k) reconstruction with a shared θ. It gives exact contact preservation for SG while keeping the mixture energy conservative. It is the cleanest idea in the file.
- **Kinetic-energy-only explicit energy rows** with dissipation on the m_kK jump only.
- **`e_linear_coeff` / `A_k, B_k`:** an EOS-agnostic linearisation that makes the pressure system linear in α_kp_k.
- **Putting A_k inside the divergence** of the artificial-viscosity term, so the phase sum is a flux difference and stays conservative for any Γ(ρ).
- **`PressureSolve_Direct`:** coloured probing of a matrix-free operator, a built-in self-check, symmetric Jacobi equilibration, SparseLU, then a preconditioned-GMRES fallback. It is robust for stiff liquid/gas systems in 1D and small 2D.
- **The `posTheta` positivity limiter** (linear constraints, closed-form θ).

**Avoid or fix:**
- Relaxation inside the RK stages, with its increment fed back through `Hs` (see the variant's evidence of pumping). Prefer Strang or step-end relaxation.
- Dead switches that suggest options which do not exist (`PressureSolver`, `add_Evisc`, the GMRES settings).
- Arithmetic face averages of C_k at material interfaces. The variant argues for pressure-upwinding (O:2222–2230).
- The q_k = α_kp_k unknown when α_k → 0, which gives 1/α amplification on the diagonal. The variant switched to p_k unknowns and replaced the row with p_k = p_j for vanishing phases.
- Non-zero-initialised `Array`, the reflective BC `k=7` bug, and the global-only Rusanov speed when trace phases are present (see the comments at O:1113–1129).

---

## 11. Variant `de47bc68-6-equation-diffuse-IMEX_WB.cpp` (3501 lines): structural diff

The core pipeline and solver are the same: PC-Rusanov explicit step, `PressureSubsystem` with two linear solves, `PressureSolve_Direct`, the WB construction and the same IMEX tableaux. The variant contains many fixes and documented failed experiments that the 2597-line file lacks. Several of its comments describe the 2597-line file's choices as already tried and reverted. That suggests the larger file may in fact be the more developed branch, so check which one is current.

| area | 2597-line file | variant |
|---|---|---|
| Relaxation module | `PressureRelaxation_Boscheri` (μ=1e60, linearised, per stage) | **`relax6.h` inlined** as `namespace r6` (O:21–371): SG thermo struct (γ, π, q, q′, c_v, covolume); `mech` (Pelanti IJMF 2022 §7.2.1, exponential semi-exact `α1 −= Δp/S·(1−e^{−μSdt})`, mixture-energy-consistent p from (59)); `therm` (§7.2.2, `1/ξ^T_k`, eq. (62)); `pT_coeffs` (23a–d); `chem` (§7.2.3, Gibbs relaxation with a **sign correction to printed eq. (72)**, complete-evaporation cap (83)). Driver `Relaxation_Pelanti` (O:2962) runs mech → therm → chem with generic-EOS `relax_mech_cell` (O:2916) and `mixture_pressure` (O:2899). Rates: `relax_mu = −1` (instantaneous), `relax_theta = 0`, `relax_nu = 0` (off); <0 means instantaneous, 0 off, >0 finite. There are guards: chemistry requires instantaneous thermal relaxation, and thermal/chemical need c_v>0. |
| Where relaxation is called | inside each RK stage, increment fed into `Hs` | **Strang half-steps** before and after the RK loop (O:~3270, O:3405). The comments record that per-stage relaxation caused pumping. |
| EOS | SG, Ideal, JWL, CochranChan | + `NASGEOS` (O:684), SG with (q, c_v, q′); base class gains `getRhoMax`, `getCv`, `getQ`, `getQprime` |
| Tests | 3 | + `testDetonationWater` (Pelanti 8.1.1, JWL + NASG, [−10,10] m, 0.9 ms; **default**), + `testWaterCavitation` (8.2.1). Dodecane uses the Le Métayer water SG pair. Defaults: `nx=500`, `RKs=1`. |
| Implicit unknown | q_k = α_kp_k | **p_k** (`QX = α·x`), initial guess p_kⁿ |
| Vanishing phase in matrix | none | row replaced by `A_j(p_k − p_j) = 0` with RHS 0 when α_k < 1e-6 (O:2048, O:2367, O:2717) |
| Elliptic face weight W_k | arithmetic mean of C_k | **pressure-upwinded** `Wup` (O:2271); `Hk` and `C_k` clamped ≥ 0 for ellipticity |
| Artificial viscosity | χ-gated, `κ = χ·vmax·dx` | **always on**, `κ = 0.5·vmax·dx` ("B-B (39)/(40)"), no χ (O:2305) |
| Use of χ | in the pressure Laplacian | **MUSCL flattener** `slope *= (1−χ)` for every variable except α (`use_shock_flattener`, O:548, O:1352–1367); detector called before reconstruction each stage |
| Robustness | — | `Array` zero-initialised (O:518); mass floor 1e-30·ρ instead of 1e-8·ρ; `enforce_alpha_floor` (NASG covolume α_k ≥ m_k b_k); `check_admissible` (ρ_kh_k > 0); generic-EOS `reset_vanishing_phase` |
| Removed or vestigial | — | `PressureRelaxation_Boscheri` and `PressureRelaxation_Exact` are gone; only the orphan header comment of the latter remains (O:3097–3102). `rescale_alpha` (O:1198) is left over from a reverted α^{n+1}-in-mass-flux attempt and is not called. `r6::relax_mech` is a stub that returns zeros (O:160). |
| Documented reverted experiments | — | α^{n+1} rescaling of mass fluxes (broke WB: 1.4e-7 → 9628); implicit mass correction (diverges at step 2); χ_k added into C_k; harmonic face averaging; B-B (36)–(38) global dissipation (broke WB); eliminating η_k from the elliptic coefficient |
| Shared bugs | reflective BC `k=7` | same (O:1646); same hard-coded output paths |
