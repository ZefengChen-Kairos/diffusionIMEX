# Finite-rate relaxation: backward-Euler and exponential closures, rate sweeps

Scripts behind `docs/figures/expo_relax.png` (Pelanti 2022 Sec. 8.2.1, 5000 cells; centre pressure and L1 error of p
against the exact fans as functions of theta and nu, for `relaxform=0` (backward Euler, the scheme) and `relaxform=1`
(frozen-coefficient exponential), the infinite rate through `unify=1`).

- `jobs_expo.txt`, `runexpo.sh`: the runs (outputs in `run/expo/`, prefixes `E1_` to `E4_`), built from `src/SPB6_PC_IMEX.cpp`.
- `exact_vals.py` -> `exact_vals.json`: exact centre values (5-eq., 4-eq. pT, pT->pTg fan) from `tests/exact_relaxed.py`.
- `ana.py`: tables of the runs; `l1.py`: L1 errors; `fig.py`: the figure.

The 0D comparison of the two closure forms (order of LSDIRK2) is `tests/plot_relax_forms.py`.
