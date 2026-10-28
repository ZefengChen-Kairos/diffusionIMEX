# CLAUDE.md — latte/ (Robotic Latte Art: cup-liquid FVM solver)

Handoff notes for Claude Code working on this directory. The rest of the repository (`../src`, `../paper`, …)
is an unrelated compressible multiphase IMEX paper project: do not touch it from latte work.

## Project context
- Goal (user): latte-art physics simulation for a CVPR-level submission. Components: FVM cup-liquid model (this
  directory), MPM pour model (pitcher → milk stream), MPM–FVM coupling, parameter identification from video,
  pitcher 6-DoF action extraction from three-view video, forward reproduction, interactive web page.
- Specification: Notion page "V3 · 拉花物理模型：守恒接触、分层流动与逐级排错" §1–§17
  (spec id `v3_proposal_2026-10-28`). **Never read or use its §18 (development log)** — user instruction.
  Other project pages (V0.5, V0.5 diagnostics, V1, 三维倒奶 v0.1, video analysis, baselines) were reviewed;
  the review is in `docs/review_v3_notion.md`.
- CVPR 2027 deadlines (UK GMT): registration 2026-11-11, paper 2026-11-17, supplementary 2026-11-24.

## User decisions so far (2026-10-28)
- Do NOT reproduce the V0.5 baseline (G0, T01, T04 are `not_run`).
- The pour model is **MPM** (replaces the Notion choice of 3D VOF / OpenFOAM in 三维倒奶 v0.1).
- Work in stages; each stage needs a runnable minimal test before the next. Stage 1 = FVM cup-liquid solver (done).
- Not yet answered by the user (ask before relying on them): which 2 of the 5 three-view videos train;
  whether μ depends on composition; contact-line model; one-way vs two-way MPM–FVM coupling; the CVPR scope.

## Working rules from the user
- Precise reasoning; say "不确定" when unsure. No physics assumption that is not in Notion or confirmed by the
  user — if unavoidable, flag it (see `docs/stage1_design.md` §4, items A1'–A7).
- Parameters live in config files (`configs/…`); fixed seeds; every run records its resolved config.
- At the end of each stage report: what was done, verification numbers/figures, known issues. User writes Chinese.

## Layout and commands
```
latte/latte/{mesh,operators,transport,inlet,solver,remap,linsolve,runio,failures}.py
latte/latte/stage1/{t_inlet,t_transport,t_dynamics,suite,common}.py
latte/configs/stage1/suite.yaml       # all stage-1 physics/numerics/thresholds/test parameters
latte/docs/stage1_design.md           # discretisation, flagged assumptions, known issues
latte/reports/stage1/{README.md,RESULTS.md,summary.json,figures/}
latte/runs/stage1/<test>_<hash>/      # manifest.json, config_resolved.yaml, metrics.json, diagnostics/*.csv
```
Run from `latte/`: `python -m pytest -q tests/` (seconds); `python -m latte.stage1.suite` (~45 min, CPU);
`python -m latte.stage1.suite --only T15,T17` (partial rerun is merged into reports/stage1/summary.json).
Dependencies: numpy, scipy, pyamg, pyyaml, matplotlib, pytest. Lint with `python -m pyflakes latte tests`.

## Numerical design in one paragraph (details: docs/stage1_design.md)
Cartesian columns (staircase cup mask) × N terrain-following σ-layers over a flat bottom; free surface η = D(x,y,t).
Cell-centred momentum + exactly projected face volume fluxes (Basilisk centered.h structure); one common pressure
solve with Dirichlet p = σκ on the free surface; the pressure/viscosity operator is a sum of squares over faces and
edges including slope cross terms, hence SPD by construction. One set of relative face fluxes moves volume, the
components B/M (G derived), mass F_ρ = Σρ_j F_j and momentum; high-order flux = flux-limited LW (Superbee, one φ per
face = min over components) + selective joint Zalesak FCT. Explicit kinematic free-surface update with the stored
projected fluxes (exact GCL). Inlet: V3 §5.2 literal (prescribed normal velocity on the patch, Neumann pressure),
mesh-aligned sub-grid mapping with totals from the exact section integral.

## Open issues / next steps (priority order)
1. Inlet ↔ free-surface connection (V3 G4): the patch follows the input section, not the surrounding surface →
   steps at the patch edge, OUT_OF_SCOPE for fast jets. Needs a sealed rule before G4 (C01–C04).
2. Replace V3's synthetic reference inlet (Fr = 5, out of model range) with in-range synthetic parameters (ask user).
3. Crosswind diffusion of the unsplit transport (T06 disc / T07 rate ≈ 0.5) → corner transport (CTU) or a
   multi-dimensional limiter before pattern-level tests (G5).
4. Before variable viscosity / rheology: replace the vector-Laplacian viscous form by the full symmetric-strain
   operator (A3) and add the normal viscous stress to the free-surface condition (A2).
5. Contact line at the cup wall (T17b blocked); μ(composition) decision (A4).
6. Next stage per the user's plan: MPM pour model; the FVM consumes `inlet.SurfaceSource` (per-column
   volume/component/mass/momentum rates) built from samples {x, dA, j, v, C}, so an MPM adapter only has to
   produce those samples at the handoff section.
7. Re-run the full stage-1 suite on the current commit (results currently mix two code versions, see
   reports/stage1/README.md).
