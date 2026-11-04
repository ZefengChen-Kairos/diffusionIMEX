# Benchmark (phase one: CPU reference)

## Environment

```
{
 "python": "3.13.17",
 "numpy": "2.5.3",
 "machine": "x86_64",
 "cpus": 4,
 "nvidia_smi": "not available (FileNotFoundError)",
 "warp": "not installed (ModuleNotFoundError); GPU backend is phase two"
}
```

## Speed

| backend | N | ms/step | steps/s | s per physical s | speedup vs numpy |
|---|---|---|---|---|---|
| numpy | 128 | 65.77 | 15.2 | 65.4 | 1.00x |
| numba | 128 | 17.33 | 57.7 | 17.6 | 3.80x |
| warp (GPU) | - | phase two | | | |

## Hot spots, numpy backend (cumulative share of step time)

| operator | what | calls | cum. s | share |
|---|---|---|---|---|
| fine_advect | fine material BFECC (numpy) | 25 | 4.166 | 63.4 % |
| rhs | coarse fluxes (MUSCL + pressure + kappa), x2 per step (RK2) | 200 | 1.811 | 27.5 % |
| fine_sources | fine source mixing in windows | 100 | 0.243 | 3.7 % |
| prims | primitive variables / pressure | 300 | 0.152 | 2.3 % |
| sources | coarse D/E/R/tau_j with windowed kernels | 100 | 0.072 | 1.1 % |
| kernel | kernel evaluation (exp on windows) | 500 | 0.071 | 1.1 % |
| dt_cfl | CFL time step | 100 | 0.063 | 1.0 % |
| check | conservation assertion | 100 | 0.027 | 0.4 % |

## Hot spots, numba backend (cumulative share of step time)

| operator | what | calls | cum. s | share |
|---|---|---|---|---|
| _fine_advect_nb | fine material BFECC (numba kernels) | 75 | 2.060 | 39.7 % |
| _rhs_nb | coarse fluxes (numba kernel), x2 per step | 600 | 1.442 | 27.8 % |
| fine_sources | fine source mixing in windows | 300 | 0.729 | 14.0 % |
| prims | primitive variables / pressure | 900 | 0.459 | 8.8 % |
| kernel | kernel evaluation (exp on windows) | 1500 | 0.255 | 4.9 % |
| sources | coarse D/E/R/tau_j with windowed kernels | 300 | 0.249 | 4.8 % |
| dt_cfl | CFL time step | 300 | 0.196 | 3.8 % |
| check | conservation assertion | 300 | 0.085 | 1.6 % |

peak RSS: 227 MB

