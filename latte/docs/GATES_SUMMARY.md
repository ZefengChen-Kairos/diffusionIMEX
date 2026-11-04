# Gate summary

| gate | check | value | threshold | result |
|---|---|---|---|---|
| G0 | G0.1 static 1000 steps: max|P| (exact zero) | 0 | == 0 | PASS |
| G0 | G0.2a foam disc: max|u| bounded [m/s] | 0.006863 | < 0.5 | PASS |
| G0 | G0.2b foam disc: checkerboard index of u (kappa on) | 0.0009005 | < 0.05 | PASS |
| G0 | G0.3a pulse: max per-step rel. error [mass, milk, tag, mom] | [3.057e-16, 5.996e-16, 3.769e-16, 4.568e-16] | < 1e-6 | PASS |
| G0 | G0.3b pulse: cumulative mass / milk closure | [4.164e-16, 5.692e-16] | < 1e-9 | PASS |
| G0 | G0.4 pitcher moves with Q=0: fields bitwise unchanged | True | bitwise | PASS |
| G0 | G0.5a mirror: coarse M, Phi, P_x mirror to rounding (rel. max diff) | [7.463e-15, 2.861e-14] | < 1e-12 | PASS |
| G0 | G0.5b mirror: white masks mirror IoU; fine max|dY| (rounding grown by the BFECC clamp) | [1, 7.955e-12] | IoU >= 0.999, dY < 1e-4 | PASS |
| G0 | G0.7 positivity M >= 0.3 rho_c h0 (leaf, 2.5 s): violations | 0 | 0 | PASS |
| G0 | G0.9 fine Y*M_interp vs coarse Phi (leaf 2.5 s): total rel. diff, restricted L1 rel. diff | [0.02089, 0.05585] | total < 5 % (diagnostic) | PASS |
| G0 | G0.8 numba == numpy reference (80 steps): rel max|dU|, max|dY| | [3.344e-16, 5.829e-16] | < 1e-12 | PASS |
| G0 | G0.6 refinement 128->192 big heart: white IoU, class (notch/lobes) unchanged | [0.9944, True] | IoU >= 0.85 & same class | PASS |
| G1 | G1.1 A(t) non-decreasing during the low dwell [0.8, 5.2] s (min dA, mm^2) | 59.25 | >= -0.12 | PASS |
| G1 | G1.2 mean dA/dt [mm^2/s] for Q = 6, 9, 12 ml/s strictly increasing | [272.5, 348.2, 431] | increasing | PASS |
| G1 | G1.3 impact hole: alpha(x_p) < 0.5 and lateral brown width 3-10 mm (alpha, mm) | [0.8389, 0] | alpha<0.5, 3..10 mm | **FAIL** |
| G1 | G1.4 horseshoe: white centroid ahead of x_p (dy mm) | 6.584 | > 0 | PASS |
| G1 | G1.5 post-stop velocity e-folding time of the white region [s] | 0.1789 | 0.1..0.5 | PASS |
| G1 | G1.6 creep: boundary displacement 1 s -> 6 s after stop, max / mean [mm] | [0, 0] | max < 2 mm | PASS |
| G1 | G1.7 after the cut: rear notch, two lobes, tip >= 5 mm, connected | {lobes: 2, notch: True, tip_len: 0.0115, ncomp: 2, width_max: 0.0545, cut_line: 0.35, cut_rows: 117} | notch & lobes>=2 & tip>=5mm & ncomp<=2 | PASS |
| G1 | G1.8 control (same pass at h=1.2 cm): brown cut line share drops below half of G1.7 | {lobes: 2, notch: True, tip_len: 0.0128, ncomp: 1, width_max: 0.056, cut_line: 0, cut_rows: 152} | cut_line < 0.5 x heart | PASS |
| G2 | G2.1 wobble amplitude 0: number of stripes on the axis | 0 | == 0 | PASS |
| G2 | G2.2 a=9 mm, f=1.5 Hz: number of stripes on the axis | 11 | >= 4 | PASS |
| G2 | G2.3 median stripe spacing increases with v_rel/f (f = 2, 1.5, 1 Hz) [mm] | [2.25, 3, 4.25] | increasing | PASS |
| G2 | G2.5 share of tracked stripes that bend backwards at x = +-8 mm (forward-convex arcs) | [0.7273, 11] | >= 0.6 of >= 3 tracked | PASS |
| G2 | G2.4 stripe contrast loss from stop+0.5 s to stop+10.5 s | [0.3765, 0.3738, 0.007251] | < 30 % | PASS |
| G3 | G3.1 retained milk mass ratio low (h=1.5 cm) / high (h=10 cm), same path | 444.2 | >= 10 | PASS |
| G3 | G3.2a high pass drags the old white: |d centroid(A)| [mm] | 0.05375 | > 0.5 mm | **FAIL** |
| G3 | G3.2b high pass adds ~no white: retained mass relative to the low pass | 0.002251 | < 0.1 | PASS |
| G3 | G3.3 milk removed by E is exactly the Phi ledger change: |dPhi - (D - E + R - X_wall)| / E | 3.771e-14 | < 1e-9 | PASS |
| G3 | G3.4a C_t = 0: along-path drag flux < 50 % of the reference | [2.39e-08, 1.504e-05] | < 50 % | PASS |
| G3 | G3.4b C_E = 0: inward (collapse) flux < 50 % of the reference (signed: < 0 means outward) | [-1.598e-07, 4.833e-07] | < 50 % | PASS |
| G4 | G4.1 dot B pushes A forward: front edge of A, t=3.55 -> 7.75 s [mm] | 8.25 | > 2 mm | PASS |
| G4 | G4.2 A is a crescent around B: angular coverage of A seen from B's centroid [deg] | 360 | > 180 | PASS |
| G4 | G4.3 brown band between B (centroid) and A's front edge on the axis: min alpha, width [mm] | [0.5057, 1] | alpha < 0.6, width > 0 | PASS |
| G4 | G4.4 control (flow runs out after part A): front edge of A over the same interval [mm], ratio | [0.5, 0.06061] | < 20 % of G4.1 | PASS |
| G4 | G4.5 tag A ledger: surface + bulk (incl. E-subducted, minus R) = injected | 2.246e-15 | < 1e-9 | PASS |
| G5 | G5.1 rib count along a line beside the axis | 6 | >= 5 | PASS |
| G5 | G5.2 left/right symmetry: mirror IoU of the white mask | 0.9887 | >= 0.75 | PASS |
| G5 | G5.4 midline vein: share of transverse profiles with a dip at x=0 | 1 | >= 0.5 | PASS |
| G5 | G5.3 rib spacing in the retreat zone before the pull-through: reference vs retreat x1.5 [mm] | [4.625, 6.625] | increases | PASS |
| G6 | G6.1a batch A is pressed open, not erased and redrawn: surface A mass kept 5.2 -> 10.6 s; components of A material (Y_A > 0.2) | [0.9475, 1] | kept >= 80 %, 1..2 pieces | PASS |
| G6 | G6.1b tag A ledger closure | 2.561e-15 | < 1e-9 | PASS |
| G6 | G6.2 white bands on the axis before the pull-through (one per part A, B, C) | 2 | == 3 | **FAIL** |
| G7 | G7.1 wing / neck / head white areas [mm^2] | [2030, 207.4, 295.6] | each > 20 mm^2 | PASS |
| G7 | G7.2 neck width < wing rib spacing [mm] | [4.075, 8] | neck < rib | PASS |
| G7 | G7.3 same physics parameters as all other patterns (digest) | 1a6b8d7baeea | == shared | PASS |
| G8 | G8.1-heart-0 heart +-10 % perturbation: class kept, white IoU vs nominal | [True, 0.9355] | class & IoU >= 0.7 | PASS |
| G8 | G8.1-heart-1 heart +-10 % perturbation: class kept, white IoU vs nominal | [True, 0.9401] | class & IoU >= 0.7 | PASS |
| G8 | G8.1-leaf-0 leaf +-10 % perturbation: class kept, white IoU vs nominal | [True, 0.8947] | class & IoU >= 0.7 | PASS |
| G8 | G8.1-leaf-1 leaf +-10 % perturbation: class kept, white IoU vs nominal | [True, 0.9511] | class & IoU >= 0.7 | PASS |
| G8 | G8.2 dt halved (CFL 0.35 -> 0.175): heart white IoU | 0.9992 | >= 0.9 | PASS |
| G8 | G8.3 cost reported (s per physical s, peak RSS MB) | 435.7 | report | PASS |
| CF | CF1a colour-only paint (tag field) mid-run: M, P, Phi and Y bitwise unchanged | True | bitwise | PASS |
| CF | CF1b injection with no volume / momentum input (D, E, tau off): max |P| after 3 s | 0 | == 0 | PASS |
| CF | CF2 height label x3 with the same (Q, v_j): output bitwise identical | True | bitwise | PASS |
| CF | CF3 R = 0 (bare sink): white area far from the cut line (|x|>10 mm) shrinks [fraction] | 0.2505 | > 10 % | PASS |
| G9 | G9-heart heart: white IoU vs the rectified video frame, pour axes aligned (axis tolerance -6 deg); IoU with the simulation upside down; white fraction sim / video | [0.6552, 0.5906, 0.4387, 0.4129] | report | PASS |
| G9 | G9-push_heart push_heart: white IoU vs the rectified video frame, pour axes aligned (axis tolerance -6 deg); IoU with the simulation upside down; white fraction sim / video | [0.398, 0.4014, 0.5109, 0.2985] | report | PASS |
| G9 | G9-layered_heart layered_heart: white IoU vs the rectified video frame, pour axes aligned (axis tolerance -6 deg); IoU with the simulation upside down; white fraction sim / video | [0.657, 0.6292, 0.4819, 0.4999] | report | PASS |
| G9 | G9-tulip tulip: white IoU vs the rectified video frame, pour axes aligned (axis tolerance -2 deg); IoU with the simulation upside down; white fraction sim / video | [0.5119, 0.4691, 0.5568, 0.3748] | report | PASS |
| G9 | G9-leaf leaf: white IoU vs the rectified video frame, pour axes aligned (axis tolerance +6 deg); IoU with the simulation upside down; white fraction sim / video | [0.6733, 0.6313, 0.5207, 0.5875] | report | PASS |
| G9 | G9-swan swan: white IoU vs the rectified video frame, pour axes aligned (axis tolerance -4 deg); IoU with the simulation upside down; white fraction sim / video | [0.5263, 0.53, 0.5577, 0.3767] | report | PASS |
