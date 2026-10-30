# Gate summary

| gate | check | value | threshold | result |
|---|---|---|---|---|
| G0 | G0.1 static 1000 steps: max|P| (exact zero) | 0 | == 0 | PASS |
| G0 | G0.2a foam disc: max|u| bounded [m/s] | 0.006863 | < 0.5 | PASS |
| G0 | G0.2b foam disc: checkerboard index of u (kappa on) | 0.0009005 | < 0.05 | PASS |
| G0 | G0.3a pulse: max per-step rel. error [mass, milk, tag, mom] | [2.163e-16, 2.655e-16, 2.397e-16, 8.865e-16] | < 1e-6 | PASS |
| G0 | G0.3b pulse: cumulative mass / milk closure | [2.208e-15, 5.206e-16] | < 1e-9 | PASS |
| G0 | G0.4 pitcher moves with Q=0: fields bitwise unchanged | True | bitwise | PASS |
| G0 | G0.5a mirror: coarse M, Phi, P_x mirror to rounding (rel. max diff) | [3.479e-14, 6.335e-14] | < 1e-12 | PASS |
| G0 | G0.5b mirror: white masks mirror IoU; fine max|dY| (rounding grown by the BFECC clamp) | [1, 3.227e-06] | IoU >= 0.999, dY < 1e-4 | PASS |
| G0 | G0.7 positivity M >= 0.3 rho_c h0 (leaf, 2.5 s): violations | 0 | 0 | PASS |
| G0 | G0.9 fine Y*M_interp vs coarse Phi (leaf 2.5 s): total rel. diff, restricted L1 rel. diff | [0.01366, 0.05612] | total < 5 % (diagnostic) | PASS |
| G0 | G0.8 numba == numpy reference (80 steps): rel max|dU|, max|dY| | [6.386e-16, 1.776e-15] | < 1e-12 | PASS |
| G0 | G0.6 refinement 128->192 big heart: white IoU, class (notch/lobes) unchanged | [0.9918, True] | IoU >= 0.85 & same class | PASS |
| G1 | G1.1 A(t) non-decreasing during the low dwell [0.8, 5.2] s (min dA, mm^2) | 60.62 | >= -0.12 | PASS |
| G1 | G1.2 mean dA/dt [mm^2/s] for Q = 6, 9, 12 ml/s strictly increasing | [262.9, 346.9, 390.1] | increasing | PASS |
| G1 | G1.3 impact hole: alpha(x_p) < 0.5 and lateral brown width 3-10 mm (alpha, mm) | [0.9358, 0] | alpha<0.5, 3..10 mm | **FAIL** |
| G1 | G1.4 horseshoe: white centroid ahead of x_p (dy mm) | 7.862 | > 0 | PASS |
| G1 | G1.5 post-stop velocity e-folding time of the white region [s] | 0.1782 | 0.1..0.5 | PASS |
| G1 | G1.6 creep: boundary displacement 1 s -> 6 s after stop, max / mean [mm] | [0, 0] | max < 2 mm | PASS |
| G1 | G1.7 after the cut: rear notch, two lobes, tip >= 5 mm, connected | {lobes: 2, notch: True, tip_len: 0.0198, ncomp: 1, width_max: 0.0405, cut_line: 0.414, cut_rows: 133} | notch & lobes>=2 & tip>=5mm & ncomp<=2 | PASS |
| G1 | G1.8 control (same pass at h=1.2 cm): brown cut line share drops below half of G1.7 | {lobes: 2, notch: True, tip_len: 0.024, ncomp: 1, width_max: 0.048, cut_line: 0, cut_rows: 190} | cut_line < 0.5 x heart | PASS |
| G2 | G2.1 wobble amplitude 0: number of stripes on the axis | 0 | == 0 | PASS |
| G2 | G2.2 a=9 mm, f=1.5 Hz: number of stripes on the axis | 12 | >= 4 | PASS |
| G2 | G2.3 median stripe spacing increases with v_rel/f (f = 2, 1.5, 1 Hz) [mm] | [2.5, 3, 4.125] | increasing | PASS |
| G2 | G2.5 share of tracked stripes that bend backwards at x = +-8 mm (forward-convex arcs) | [0.5455, 11] | >= 0.6 of >= 3 tracked | **FAIL** |
| G2 | G2.4 stripe contrast loss from stop+0.5 s to stop+10.5 s | [0.3415, 0.3405, 0.003046] | < 30 % | PASS |
| G3 | G3.1 retained milk mass ratio low (h=1.5 cm) / high (h=10 cm), same path | 449 | >= 10 | PASS |
| G3 | G3.2a high pass drags the old white: |d centroid(A)| [mm] | 1.111 | > 0.5 mm | PASS |
| G3 | G3.2b high pass adds ~no white: retained mass relative to the low pass | 0.002227 | < 0.1 | PASS |
| G3 | G3.3 milk removed by E is exactly the Phi ledger change: |dPhi - (D - E + R)| / E | 5.91e-16 | < 1e-9 | PASS |
| G3 | G3.4a C_t = 0: along-path drag flux < 50 % of the reference | [6.845e-08, 5.515e-05] | < 50 % | PASS |
| G3 | G3.4b C_E = 0: inward (collapse) flux < 50 % of the reference | [6.115e-08, 7.016e-06] | < 50 % | PASS |
| G4 | G4.1 dot B pushes A forward: front edge of A, t=1.9 -> 3.6 s [mm] | 3.25 | > 2 mm | PASS |
| G4 | G4.2 A is a crescent around B: angular coverage of A seen from B's centroid [deg] | 140 | > 180 | **FAIL** |
| G4 | G4.3 brown band between A and B on the axis: min alpha, width [mm] | [0.04092, 1.25] | alpha < 0.6, width > 0 | PASS |
| G4 | G4.4 control (no dot B): front edge of A over the same interval [mm], ratio | [0, 0] | < 20 % of G4.1 | PASS |
| G4 | G4.5 tag A ledger: surface + bulk (incl. E-subducted, minus R) = injected | 1.211e-14 | < 1e-9 | PASS |
| G5 | G5.1 rib count along a line beside the axis | 8 | >= 5 | PASS |
| G5 | G5.2 left/right symmetry: mirror IoU of the white mask | 0.8659 | >= 0.75 | PASS |
| G5 | G5.4 midline vein: share of transverse profiles with a dip at x=0 | 1 | >= 0.5 | PASS |
| G5 | G5.3 rib spacing before the cut: reference vs retreat x1.5 [mm] | [5.607, 6.792] | increases | PASS |
| G6 | G6.1a batch A is pressed open, not erased and redrawn: surface A mass kept 2.5 -> 4.9 s; components of A material (Y_A > 0.2) | [1.402, 1] | kept >= 80 %, 1..2 pieces | PASS |
| G6 | G6.1b tag A ledger closure | 3.132e-15 | < 1e-9 | PASS |
| G6 | G6.2 white bands on the axis before the cut (one per injection segment) | 3 | == 3 | PASS |
| G7 | G7.1 wing / neck / head white areas [mm^2] | [1530, 241.9, 215.1] | each > 20 mm^2 | PASS |
| G7 | G7.2 neck width < wing rib spacing [mm] | [3.889, 5.303] | neck < rib | PASS |
| G7 | G7.3 same physics parameters as all other patterns (digest) | 4512db61f9c2 | == shared | PASS |
| G8 | G8.1-heart-0 heart +-10 % perturbation: class kept, white IoU vs nominal | [True, 0.9319] | class & IoU >= 0.7 | PASS |
| G8 | G8.1-heart-1 heart +-10 % perturbation: class kept, white IoU vs nominal | [True, 0.9509] | class & IoU >= 0.7 | PASS |
| G8 | G8.1-leaf-0 leaf +-10 % perturbation: class kept, white IoU vs nominal | [True, 0.8504] | class & IoU >= 0.7 | PASS |
| G8 | G8.1-leaf-1 leaf +-10 % perturbation: class kept, white IoU vs nominal | [True, 0.8327] | class & IoU >= 0.7 | PASS |
| G8 | G8.2 dt halved (CFL 0.35 -> 0.175): heart white IoU | 0.9976 | >= 0.9 | PASS |
| G8 | G8.3 cost reported (s per physical s, peak RSS MB) | 463 | report | PASS |
| CF | CF1a colour-only paint (tag field) mid-run: M, P, Phi and Y bitwise unchanged | True | bitwise | PASS |
| CF | CF1b injection with no volume / momentum input (D, E, tau off): max |P| after 3 s | 0 | == 0 | PASS |
| CF | CF2 height label x3 with the same (Q, v_j): output bitwise identical | True | bitwise | PASS |
| CF | CF3 R = 0 (bare sink): white area far from the cut line (|x|>10 mm) shrinks [fraction] | 0.9881 | > 10 % | PASS |
