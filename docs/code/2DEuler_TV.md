# `2DEuler_TV.cpp`: code notes

Source: `2DEuler_TV.cpp` (3275 lines), the 2D ideal-gas Euler solver with Toro–Vázquez-Cendón (TV) flux splitting and a semi-implicit IMEX-RK integrator.
Line numbers below refer to that file. The purpose of this note is to make the TV-split IMEX structure reusable for the 6-equation two-phase model.

---

## 1. Architecture

### 1.1 Global configuration (lines 19–75)

Everything is set by globals and enums at the top of the file. There is no input file.

| Switch | Line | Default | Meaning |
|---|---|---|---|
| `test` | 36 | `test2d1` | initial/boundary data (see §7) |
| `nxCells, nyCells` | 42–43 | 256 × 256 | `const int` |
| `ghost_cells` | 45 | 5 | enough for WENO5 |
| `RKs` | 47 | 2 | IMEX stages (1, 2, 3) |
| `CFL` | 49 | 0.3 | |
| `recon` | 51 | `Minmod` | `NoRecon, Minmod, VanLeer, ENO3, WENO5` |
| `full_step` | 52 | 10 | number of start-up steps that use the acoustic CFL |
| `add_pvisc` | 53 | `MaterialSpeedp` | adds an implicit pressure Laplacian to the p-system |
| `add_shock_detector` | 54 | `Yes_add` | weights that Laplacian by the shock sensor χ |
| `shock_kappa1` | 55 | 1e-3 | sensor threshold κ₁ |
| `GMREStol, GMRESiter` | 57–58 | 1e-10, 1e16 | Krylov tolerance (absolute) and max iterations |
| `U_update` | 61 | `DefinitionofU` | how E is rebuilt after the p-solve |
| `PressureSolver` | 41 | `all_correction` | `no_correction`/`once_correction`/`all_correction` (Picard pass) |
| `caltimestep` | 35 | `IMEXdt` | `Fullydt` (always acoustic) or `IMEXdt` (material after `full_step` steps) |
| `flux_diss` | 38 | `Global_Diss` | ν for the pressure Laplacian: global `vx_max·dx` or local `|u_i|·dx` |
| `vortex_Mach` | 33 | 0.1 | Gresho / BC pressure level |
| `formula` | 37 | unused | |

`gamma_plasma = 1.4` (line 73) is the ideal-gas γ. The "plasma" name and `Bx_p0`, `By_0p` (line 2965) are left over from an MHD code and are unused.

### 1.2 Data layout (lines 79–162)

- `Grid(nx,ny)`: a flat `double*` of size `(nx+2g)(ny+2g)·NUM_VARS` with `NUM_VARS=4`. Index: `NUM_VARS*(i + j*(nx+2g)) + k`, so the variables are interleaved per cell (AoS) and `i` is the fast index.
- `Array(nx,ny)`: the same layout for one scalar.
- `Tensor(nRK,nx,ny)`: a `Grid` per RK stage, used for the stage derivatives `Hs`.
- Conservative variables: `(ρ, ρu, ρv, E)`. Primitive variables: `(ρ, u, v, p)`.
- The interior is `i ∈ [g, nx+g)`, `j ∈ [g, ny+g)`. The cell centre is `x_0 + (i-g+0.5)dx`.
- None of these structs defines a copy constructor or assignment operator (the rule of three is violated), so copying one would cause a double free. They are never copied. `Grid`s are heap-allocated inside hot functions: `TVDReconstruction` (l. 514) and `ShockDetector` (l. 1669, 1674). The second one is called inside every mat-vec (see §3.5).
- Eigen is included (`Dense`, `SparseLU`, `IterativeLinearSolvers`, `unsupported/IterativeSolvers`). The `MatrixReplacement` wrapper (l. 165–228), `SpMat M_p` and `triplets_p` (l. 2091–2092) are dead code. The solve is matrix-free with a hand-written BiCGStab (§3.4).

### 1.3 Main loop (lines 2773–3113)

```
main: set domain/BC/time_stop per test (2778–2956), dx,dy (2959)
      u_initialization(u_pri) -> upri2ucon -> bound_cond_u            (2966–2970)
      set Butcher tableaux Ae, Ai, bRK                                 (2994–3044)
      do {
        dt = min(caldt(u_pri), time_stop - t); t += dt                 (3049–3051)
        u_con_n = u_con
        for rks in 0..RKs-1:                                            (3058)
          Qe = Qs = u_con_n
          Qe += Σ_{r<rks} Ae(rks,r) dt Hs(r);  Qs += Σ_{r<rks} Ai(rks,r) dt Hs(r)
          adt = Ai(rks,rks) dt
          TVDReconstruction(Qe -> w_xL,w_xR,w_yL,w_yR)                  (3076)
          Update(u_con = Qs - adt div F_conv(Qe))                       (3079)  explicit convective
          PressureSubsystem(u_con, Qs, Qe, ..., adt, rks)               (3083)  implicit pressure
          Hs(rks) = (u_con - Qs)/adt                                     (3089–3095)
        output (3101), ShockDetector diagnostic (3103)
      } while (t < time_stop)
      convergence-test L2 errors (3115–3158); dump .dat files (3160–3273)
```

`t` is advanced before the stages, and the stage times are never used because the system is autonomous.

### 1.4 Time step / CFL (`caldt`, lines 1538–1590)

```cpp
dt = CFL * 1 / (a_max_x / dx + lambda_y / dy);   // l.1588
```

- `Fullydt`, or `IMEXdt` while `timestep < full_step+1` (the first 11 steps): an acoustic estimate `a_max_x = |u|+c` is intended.
- `IMEXdt` after that: the material estimate `a_max_x = vx_max = max|u|`.
- The same loop sets `vx_max = max|u|` and `vy_max = max|v|`, which are global variables. They are the dissipation speeds of the explicit Rusanov flux (§2.2) and the ν of the pressure Laplacian. They are frozen for the whole step, from `u_pri` at time `tⁿ`.

**Bugs/caveats:**
1. `a_max_x = fmax(fabs(u)+a, fabs(lambda_x))` with `lambda_x = 0` (l. 1556, 1568) is **not a running max**. It holds the value of the *last* cell visited, the top-right interior corner. During the start-up steps the acoustic dt is therefore wrong.
2. `lambda_y` is always 0, so the y-direction speed never enters dt. `a_max_y` is computed but not used. The 2D CFL is effectively `max|u|·dt/dx ≤ CFL`.
3. If `max|u| = 0` (for example a fluid at rest after the start-up steps), dt = ∞ and is clipped to `time_stop - t`.
4. `fullydt` and `imexdt` are never computed (they stay 0), so `dt_fullydt.dat` and `dt_imexdt.dat` contain `inf`.

### 1.5 Boundary conditions

- `bound_cond_u` (l. 1366–1536) fills all `g` ghost layers.
  - Transmissive: zero-order copy.
  - Periodic.
  - Reflective: mirror. **Both** `ρu` and `ρv` are negated in each direction (l. 1389–1393, 1472–1476). This is a no-slip-like wall, not a slip wall: the tangential momentum should not change sign.
  - Dirichlet: hard-coded states for Gresho and RPContact only.
- `bound_cond_array` (l. 1592–1646) is the scalar version. Reflective negates the scalar, which is right for the normal momentum and wrong for p. That is why the p-solution of `x_Reflective` runs uses the Transmissive version (l. 2192–2196).
- `M_P(i,j)` (l. 1749–1800) maps a (possibly ghost) cell to an unknown index `(I-g) + nx(J-g)`. Periodic wraps around. Transmissive, Dirichlet and Reflective clamp to the boundary cell, which gives a homogeneous Neumann condition on p in the operator.
- Dirichlet in the p-operator (only used by the commented-out Gresho Dirichlet option, l. 2896): the ghost value is `-p_mirror` (unknown part, `get_x_unknown`, l. 1832–1857) plus `p_bc` (constant part, `get_x_const`, l. 2135–2141). The face value is then `p_bc/2`, not `p_bc`, which looks like a missing factor 2. The pressure Laplacian term also omits the constant part. After the solve, the ghosts are overwritten with `p_bc` (l. 2198–2224), which is inconsistent with the operator. Do not reuse this.

### 1.6 Output

- `printAllVariablesToVTK` (l. 2651–2771) writes ASCII VTK `RECTILINEAR_GRID` cell data: `rho, mom_X, mom_Y, U(=E), velX, velY, p`. The path is hard-coded under `/Users/chenzefeng/...` and encodes the switches (l. 2733).
  - Minor bugs: the X coordinates are shifted by `dx/2` and the Y coordinates lack the `y_0` offset (l. 2599–2609).
  - The `p/p_ref*p_epsilon` scaling is never applied in the VTK file because the condition `var-NUM_VARS+1 == 4` is never true (l. 2635).
- `IntermediatePlotIntermediatePlot` (l. 2476–2577) writes a file every 100 steps, at the listed fractions of `time_stop`, and at step 0.
- At the end, `.dat` dumps are written to `~/Desktop`: `RK_rho`, `RK_vx`, `RK_vy`, `RK_E`, `RK_p` (scaled), `RK_v_mag`, mid-line cuts and `shock_de`.

---

## 2. The TV split and the explicit convective subsystem

### 2.1 Split (lines 276–360)

For a normal x-flux, with `k = ½ρ|v|²` (`calkinetic`, l. 1656) and specific enthalpy `h = e + p/ρ` (`calenthalpy`, l. 1648):

| eq. | convective `f_convective` (l. 276) | pressure `f_pressure` (l. 293) |
|---|---|---|
| ρ | ρu | 0 |
| ρu | ρu² | p |
| ρv | ρuv | 0 |
| E | **k·u** | **ρu·h** = (ρe+p)u |

The y-fluxes (`g_convective`, `g_pressure`, l. 319–350) are analogous. This is exactly the Toro–Vázquez split: the convective system carries only the kinetic energy, and the pressure system carries the enthalpy flux. The pressure subsystem does not change ρ.

### 2.2 Explicit numerical flux (lines 417–447)

It is a **Rusanov (local Lax–Friedrichs) flux on the convective flux only, with a global material-speed dissipation**. It is not path-conservative.

```cpp
flux[k] = 0.5*(fL_c[k] + fR_c[k]) - 0.5*vx_max*(QR[k] - QL[k]);   // RusanovC_x, l.429
```

- The dissipation speed is `vx_max = max_cells |u|` (or `vy_max` in y). It is global, frozen per step, and contains no sound speed. This is consistent with the TV convective system, whose eigenvalues are 0, u, u, 2u, so the true bound is 2|u|; see caveats.
- The jump is applied to all four conserved variables, including E, so the internal-energy jumps are also diffused with `vx_max`.
- `f_pressure` is evaluated in `RusanovC_*`, but the result is unused. `Q_cellL/R` in `Update` are also unused.
- Caveat: the eigenvalues of the TV convective system in x are `(0, u, u, 2u)`, so `½·max|u|` does not bound `|2u|`. Rusanov stability therefore relies on CFL = 0.3. Using `2·max|u|` would be safer.
- With `max|u| = 0` (for example acoustic waves at rest) the explicit flux has no dissipation at all.

### 2.3 Update (lines 450–506)

```cpp
u_con(i,j,k) = Qs(i,j,k) - adt/dx*(fp[k]-fm[k]) - adt/dy*(gp[k]-gm[k]);   // l.502
```

The face states are reconstructed from **Qe**. The explicit increment is added to **Qs** with the factor `adt = a_ii·dt` (see §5).

### 2.4 Reconstruction (`TVDReconstruction`, lines 510–781)

- Reconstruction is done on **primitive** variables of `Qe` (l. 514–515). The results are converted back to conservative variables in place (l. 759–762), then `bound_cond_u` is applied to the face arrays (l. 777–780).
- Minmod (l. 532–562): when `test == convergencetest` the limiter is replaced by the unlimited central slope.
- Van Leer (l. 563–590): the form `c·slope_L` and `c·slope_R` is non-standard. Two different slopes are used for the left and right values, so this is not a single limited slope.
- ENO3 (l. 591–670) and WENO5-JS with ε = 1e-6 and power 2 (l. 671–753).
- x is reconstructed for `i ∈ [g-1, nx+g]`, but y only for interior `j`. The boundary-face states in y (and in x, which `bound_cond_u` overwrites) are ghost copies of interior face values, not true reconstructions. This is fine for periodic BCs and approximate for transmissive and reflective ones. For reflective BCs the mirrored `w_xR(g-1)` should come from `w_xL(g)`, not `w_xR(g)`.

---

## 3. The implicit pressure subsystem (`PressureSubsystem`, lines 2031–2468)

### 3.1 Inputs

- `u_con` holds the explicit result `Q* = (ρ*, m*, E*)`, where `ρ* = ρ^{n+1}` for this stage.
- `Qe` is the explicit stage state, used for linearisation.
- `h ≡ adt = a_ii·dt`.
- `m* = (rhovx_star, rhovy_star)` is copied with BCs applied (l. 2045–2088).

### 3.2 Linearised enthalpy at cells and faces

```cpp
double h00 = calenthalpy(con00_n) * Qe(i,j,0) / u_con(i,j,0);   // l.2115 (same in AX l.1860)
```

That is, `H_c = (ρh)(Qe)_c / ρ*_c = (ρe + p)^e / ρ^{n+1}`. The enthalpy *density* is frozen at the explicit state and divided by the new density, so that `H·m^{n+1} ≈ (ρe+p)^e u^{n+1}`. Face values are arithmetic means, `H_{i+1/2} = ½(H_i + H_{i+1})`. H is **not** updated during the Picard pass.

### 3.3 Discrete equations as implemented

Momentum update, done after the solve (l. 2227–2257). This is collocated and uses central face averages of p:

```
m_x,i^{n+1} = m*_x,i - h (p_{i+1} - p_{i-1}) / (2dx)
m_y,j^{n+1} = m*_y,j - h (p_{j+1} - p_{j-1}) / (2dy)
```

The energy update is conservative, with a face mass flux equal to the average of the cell momenta:

```
E^{n+1} = E* - h [ (H_{i+½} m̄_{i+½} - H_{i-½} m̄_{i-½})/dx + (y) ],   m̄_{i+½} = ½(m_i + m_{i+1})^{n+1}
```

The EOS is `E = p/(γ-1) + k`. Substituting gives the **pressure equation** (operator `Pressure_Subsystem_AX`, l. 1802–1936; RHS l. 2122):

```
p_i/(γ-1) - h² (Lap_x + Lap_y) - [η_i h/(γ-1)] (ν_x δ²_x p / dx² + ν_y δ²_y p / dy²)
   = E*_i - k_i^{(l)} - h [ (H_{i+½} m̄*_{i+½} - H_{i-½} m̄*_{i-½})/dx + (y) ]
```

with

```
G_i   = (p_{i+1} - p_{i-1})/(2dx)                     // gradP_i, cell-centred
Lap_x = [ H_{i+½} ½(G_i + G_{i+1}) - H_{i-½} ½(G_{i-1} + G_i) ] / dx
      = [ H_{i+½}(p_{i+2} + p_{i+1} - p_i - p_{i-1}) - H_{i-½}(p_{i+1} + p_i - p_{i-1} - p_{i-2}) ] / (4dx²)
```

Code:

```cpp
Ax(M_P(i,j)) = 1.0/(gamma_plasma - 1) * p_i - h*h * (Lap_x + Lap_y);          // l.1904
...
Ax(M_P(i,j)) -= shock_eta(i,j) * 1.0/(gamma_plasma - 1) * h * lap_p;            // l.1930
// lap_p = ν_x/dx²·(p_{i+1}-2p_i+p_{i-1}) + ν_y/dy²·(...),  ν_x = vx_max·dx (Global) or |u_i|dx (Local)
```

Stencil coefficients in x, with `c = 1/(4dx²)` and before the factor `-h²`:

| p_{i-2} | p_{i-1} | p_i | p_{i+1} | p_{i+2} |
|---|---|---|---|---|
| c·H₋ | c·(H₋ - H₊) | -c·(H₊ + H₋) | c·(H₊ - H₋) | c·H₊ |

So the operator is a **wide 9-point cross** (i±2 and j±2), not a compact 5-point stencil. The code comment "composed 5-point" means 5 points per direction. When H is constant, the p_{i±1} coefficients vanish and the operator becomes `H(p_{i+2} - 2p_i + p_{i-2})/(4dx²)`. This operator **decouples odd and even cells and has checkerboard null modes**. That is the reason the compact Laplacian `η·ν·δ²p` has to be added (§6).

The scheme is the exact Schur complement of the collocated momentum and energy updates. As a result, the discrete energy is conservative with flux `H m̄^{n+1}` even though E is rebuilt from p (§4).

The kinetic-energy linearisation (Picard on k) works as follows:
- First solve: `k^{(0)} = |m*|²/(2ρ*)` (l. 2122).
- If `PressureSolver == all_correction` (every stage) or `once_correction` (last stage only), a second solve uses `k^{(1)} = |m^{(1)}|²/(2ρ*)`, where `m^{(1)}` is the momentum after the first p-solve (l. 2295). The divergence term still uses `m*`, which is correct.
- `no_correction` gives a single solve.
- There is only one Picard correction. No convergence check is made.

### 3.4 Linear solver (`GMRES_Solver`, lines 1939–2028)

- Despite its name, it is a **hand-written, unpreconditioned BiCGStab**. It is matrix-free and calls `Pressure_Subsystem_AX` for every mat-vec (two per iteration).
- Stopping test: **absolute** `‖r‖₂ < GMREStol = 1e-10`, with `GMRESiter = 1e16`, which is effectively unlimited. There is no relative tolerance and no warning or iteration log. The `converged` and `breakdown` flags are computed but ignored. For large p (low-Mach Gresho at M = 1e-3, p ≈ 7e5, N = 65536) the absolute 1e-10 tolerance is close to round-off, so the solver can stall or loop for a very long time.
- Initial guess: `p_vect = RHS` (l. 2178, 2350). Since RHS ≈ p/(γ-1), the guess is wrong by a factor (γ-1) = 0.4. The second Picard solve also restarts from RHS instead of from p^{(1)}. Using `(γ-1)·RHS`, `pⁿ` or `p^{(1)}` would be better guesses.
- Dirichlet constants are moved to the RHS through `L_const` (l. 2143–2176). See the caveat in §1.5.
- The operator is non-symmetric (the H weights and the local ν), so BiCGStab is a reasonable choice. A Jacobi or ILU preconditioner, or an assembled SpMat with Eigen BiCGSTAB+IncompleteLUT, would be the obvious improvement.

### 3.5 Performance issue

`Pressure_Subsystem_AX` calls `ShockDetector(shock_eta, Qe)` on **every mat-vec** (l. 1804–1805). That means two `Grid` allocations, a full con→pri conversion and a 3×3 sweep each time, even though η depends only on `Qe`, which is fixed during the solve. Compute η once per stage and hoist it out of the operator. The enthalpies H are also recomputed in every mat-vec.

---

## 4. Update order after the pressure solve

1. First solve, giving p^{(1)}. Fill the ghosts (Transmissive if `x_Reflective`; Dirichlet overrides) (l. 2185–2224).
2. **Momentum**: `m^{(1)} = m* - h·∇_c p^{(1)}` with central face averages `½(p_i + p_{i±1})` (l. 2227–2257). Then `bound_cond_u`.
3. (Optional) second solve, giving p^{(2)}, with k from m^{(1)} (l. 2264–2354).
4. **Momentum is NOT updated again with p^{(2)}.** The stored momentum stays m^{(1)}.
5. **Energy** (l. 2399–2465):
   - `DefinitionofU` (default): `E = p^{final}/(γ-1) + |m^{(1)}|²/(2ρ)` (l. 2403). Algebraically this equals `E* - h·div(H·(m* - h∇p^{final}))`, so total energy is conservative in flux form, apart from the η-Laplacian (see §6). However, the momentum that defines the energy flux (m^{(2)}) is not the stored momentum (m^{(1)}). With `all_correction` this is an internal inconsistency. **Recommendation:** after the final p-solve, update momentum with the final p as well.
   - `ConservationofU` (l. 2407–2464): **broken**. `flux_xp`, `flux_xm`, `flux_yp` and `flux_ym` are never assigned: `fL` and `fR` are computed and then discarded, so E is updated with uninitialised values (undefined behaviour). It also rebuilds h from the new p instead of reusing H. Do not use it.
6. `bound_cond_u(u_con)`.

---

## 5. IMEX-RK scheme (lines 2994–3044, 3058–3098)

This is the **semi-implicit (partitioned / linearly-implicit) IMEX-RK** of Boscarino–Filbet–Russo. There is a single stage derivative per stage:

```
Qe^(i) = Qⁿ + dt Σ_{j<i} ã_ij K_j          (explicit/linearisation state)
Qs^(i) = Qⁿ + dt Σ_{j<i} a_ij K_j          (implicit predictor)
Q^(i)  = Qs^(i) + a_ii dt K_i,   K_i = -∇·F_conv(Qe^(i)) - ∇·F_p(p^(i); H(Qe^(i)), m^(i))
K_i    = Hs(i) = (Q^(i) - Qs^(i)) / (a_ii dt)        // l.3092
```

- Explicit terms are evaluated at `Qe` and scaled by `a_ii·dt`. They become "explicit" because they are evaluated on `Qe`, not because they use a separate K.
- `Hs` contains both the convective and the pressure contributions.
- There is no final `Qⁿ + dt Σ b_i K_i` assembly: `u_con` after the last stage is taken as `Q^{n+1}`. This is only valid because both implicit tableaux are **stiffly accurate** (the last row of `Ai` equals `b`). `bRK` is set but never used.

Tableaux:

- **RKs = 1**: `Ae = [0]`, `Ai = [1]`, b = 1. This is linearly-implicit backward Euler (explicit forward Euler for the convective part).
- **RKs = 2 (default)**, with γ = 1 - 1/√2 and c = 1/(2γ) ≈ 1.707:

  ```
  Ae = | 0   0 |     Ai = | γ    0 |     b = (1-γ, γ)
       | c   0 |          | 1-γ  γ |
  ```

  The second-order conditions hold: b·c_e = γ·c = ½, b·c_i = (1-γ)γ + γ = ½ and Σb = 1. The implicit part is L-stable, SDIRK and stiffly accurate. The scheme is type A (`a_11 = γ ≠ 0`), so stage 1 is already implicit with `Qe = Qs = Qⁿ`. Note that the explicit abscissa `c_2 = 1.707 > 1`.
- **RKs = 3**:

  ```
  Ae = | 0   0   0 |   Ai = | 1/4  0    0  |   b = (1/3,1/3,1/3)
       | ½   0   0 |        | 0    1/4  0  |
       | ½   ½   0 |        | 1/3  1/3  1/3|
  ```

  Here b·c_e = b·c_i = ½, but b·c_e² = 5/12 ≠ 1/3, so this scheme is also **second order** only.

`SpaceOrder` is printed but has no effect on time accuracy.

---

## 6. Shock capturing, all-Mach features

- **Shock sensor** (`ShockDetector`, l. 1667–1747; formulas labelled (4.40)–(4.41) in the comments):

  ```
  χ̃_i = min(1, max(0, -(h_min ∇·u + κ₁ c_min)/(κ₁ c_min))),  c_min over 3×3,  ∇·u central,  κ₁ = 1e-3
  χ_i  = max_{3×3} χ̃
  ```

  χ = 1 only in compressive cells (`h∇·u < -2κ₁c_min`). `No_add` gives χ ≡ 1.
- **Implicit pressure artificial viscosity** (`add_pvisc = MaterialSpeedp`, l. 1906–1932): `-(χ_i h/(γ-1))·ν·δ²p` with `ν = vx_max·dx` (Rusanov-like at the material speed). It is added to the **operator only**, which amounts to an implicit backward-Euler diffusion of ρe inside the p-equation. Because `χ_i` multiplies the Laplacian outside the divergence, the term is **not in flux form** where χ varies, and E conservation is lost slightly there. The term is also needed to damp the checkerboard modes of the wide operator (§3.3). With the sensor on, it vanishes in smooth or expanding regions, so smooth flows rely only on the odd–even-coupled wide stencil.
- **All-Mach**:
  - The CFL is material-based after `full_step` start-up steps, the explicit dissipation scales with |u|, and pressure waves are handled implicitly.
  - There is no low-Mach preconditioning of the explicit flux. None is needed, since the convective dissipation is O(|u|).
  - Low-Mach tests: Gresho (`vortex_Mach`), DoubleShearLayer (p = 1e5/γ), KH (p = 1e4/γ), CollidingPulses (Klein, M = 1/11, `p_ref`), RiemannLike (p = 1/ε, ε = 0.99).

---

## 7. Test cases (init l. 880–1363; domain/BC/time l. 2778–2956)

Primitive values are `(ρ, u, v, p)`.

| Test | Domain | BC | t_end | Data |
|---|---|---|---|---|
| testToroModified | [0,1]² | T | 0.2 | x<0.3: (1, 0.75, 0, 1); else (0.125, 0, 0, 0.1) |
| testToro1 (Sod) | [0,1]² | T | 0.25 | (1,0,0,1) / (0.125,0,0,0.1) at x=0.5 |
| testToro2 | [0,1]² | T | 0.15 | (1,-2,0,0.4) / (1,2,0,0.4) |
| testToro3 | [0,1]² | T | 0.012 | p 1000 / 0.01 |
| testToro4 | [0,1]² | T | 0.035 | p 0.01 / 100 |
| testToro5 | [0,1]² | T | 0.035 | (5.99924,19.5975,0,460.894)/(5.99242,-6.19633,0,46.095) |
| RPContact | [0,1]² | T (Dirichlet data hard-coded) | 0.5 | ρ 1000/0.01 at x=0.25, u=1, p=1e5 |
| RPLax / RPFVVEM3 | [0,1]² | T | 0.14 | (0.445,0.698,0,3.528)/(0.5,0,0,0.571) |
| RPTwoShocks | [0,1]² | T | 0.8 | (1,2,0,0.1)/(1,-2,0,0.1) |
| RPTwoRare / RPFVVEM2 | [0,1]² | T | 0.15 | (1,-1,0,0.4)/(1,1,0,0.4) |
| RPFVVEM1 | [0,1]² | T | 0.2 | Sod |
| RPFVVEM4 | [0,1]² | T | 0.012 | u=-19.59745, p erf-smoothed 1000→0.01 at x=0.8 |
| testExplosion | [0,2]² | T | 0.25 | circle r≤0.4 about (1,1): (1,0,0,1), else (0.125,0,0,0.1) |
| testSmoothAcousticWave | [0,4]² | P | 1 | ρ=1, u=0, p = 1+exp(-40 r²) centred at (2,2) |
| testDoubleShearLayer | [0,1]² | P | 1.8 | u = tanh(30(y-0.25)) or tanh(30(0.75-y)), v = 0.05 sin 2πx, p = 1e5/γ |
| testWoodwardColellaBlastWave | [0,1]² | R | 0.038 | p = 1000 / 0.01 / 100 |
| testGreshoVortex | [-0.5,0.5]² | P (Dirichlet commented) | 0.4π | standard Gresho, p0 = ρ/(γM²), M = `vortex_Mach` |
| testCollidingPulses | [-L,L]×[0,1], L=22 | P | 0.841 | Klein: M0=1/11, ρ=0.955+M0/2·2(1-cos), u=∓√γ(1-cos), p=(1+M0γ(1-cos))p_ref |
| testRiemannLike | [0,1]×[-0.5,0.5] | P | 0.05 | ε=0.99: u = 1∓ε/2 in bands, p = 1/ε |
| convergencetest | [0,10]² | P | 1 | Hu–Shu isentropic vortex β=5, (1,1,1,1) background, 3×3 Gauss cell averages; L2 errors printed |
| test2d1 (default) | [-0.5,0.5]² | T | 0.8 | four quadrants split at (0.3,0.3): (1.5,0,0,1.5), (0.5323,1.206,0,0.3), (0.138,1.206,1.206,0.029), (0.5323,0,1.206,0.3) |
| test2d2 | [-0.5,0.5]² | T | 0.25 | 1.1/0.5065, 0.8939 quadrants |
| test2d3, test2d4 | [-0.5,0.5]² | T | 0.23, 0.3 | vortex-sheet configurations (ρ 1/2/1/3, u ±0.75, v ±0.5, p 1) |
| KH_instability | [-1,1]² | P | 5.0 | ρ = 1-0.25 tanh(25(|y|-0.5)), u = -0.5 tanh(…), v = 0.01 sin2πx cos2πy, p = 1e4/γ |

T = transmissive, P = periodic, R = reflective. In `testRiemannLike` the variable is named `p_epsilon` and is used to scale the output.

---

## 8. Reuse and avoid list (for the 6-equation two-phase model)

**Reuse:**
- The **TV split** with the convective flux carrying `k·u` and the pressure flux carrying `(ρe+p)u`. For the two-phase model, put the phase masses, the volume fraction advection and the kinetic parts in the explicit system, and the pressure–work terms in the implicit system.
- The **semi-implicit IMEX-RK with one K per stage** (`Qe`/`Qs`/`Hs`, l. 3058–3098). It is simple and stage-generic. The RK2 tableau is stiffly accurate, so the last stage is the solution.
- The **frozen enthalpy density divided by the new density**, `H = (ρe+p)^e/ρ^{n+1}`. It keeps the p-equation linear while making `H·m^{n+1}` approximate `(ρe+p)u^{n+1}`.
- **Picard on the kinetic energy** with a fixed linear operator: only the RHS changes between passes.
- The **matrix-free Schur complement**, which yields a conservative energy flux by construction.

**Avoid / fix:**
- The **ideal-gas factor `1/(γ-1)` is hard-coded** in the operator (l. 1904, 1930), in `L_const`, and in the E rebuild (l. 2403). Generalise it to `∂(ρe)/∂p` and an EOS offset. For a stiffened gas, `ρe = (p+γπ)/(γ-1)`, so the `γπ/(γ-1)` term moves to the RHS. For the 6-equation model, use a mixture `Σ α_k/(γ_k-1)` after pressure relaxation.
- The **collocated central p-gradient with averaged face momenta** gives the wide, odd–even-decoupled operator. Prefer a compact face-based gradient (`m_{i+½} = m̄* - h(p_{i+1}-p_i)/dx` for the energy flux), or keep the χ-weighted Laplacian in flux form, `∇·(χν∇p)`.
- The second Picard p is used for E but not for momentum (§4). Update momentum with the final p.
- `ConservationofU` is broken (uninitialised fluxes).
- dt: the "max" is not a max, the y-direction is ignored, and the explicit dissipation `½·max|u|` does not bound the convective eigenvalue 2u.
- Reflective BCs negate the tangential momentum.
- Dirichlet p-BC constants are inconsistent.
- BiCGStab: absolute tolerance, infinite iteration cap, poor initial guess, no preconditioner, and the shock detector is recomputed in every mat-vec.
- Hard-coded output paths. Many dead declarations (`Roe_Flux`, `MatrixReplacement`, the Eigen sparse objects, `Formula`, `Bx/By` arrays, `bRK`).
