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
| numpy | 128 | 90.14 | 11.1 | 114.0 | 1.00x |
| numba | 128 | 29.49 | 33.9 | 37.3 | 3.06x |
| warp (GPU) | - | phase two | | | |

## Hot spots, numpy backend (cumulative share of step time)

| operator | what | calls | cum. s | share |
|---|---|---|---|---|
| fine_advect | fine material BFECC (numpy) | 25 | 5.254 | 58.3 % |
| rhs | coarse fluxes (MUSCL + pressure + kappa), x2 per step (RK2) | 200 | 2.252 | 25.0 % |
| fine_sources | fine source mixing in windows | 100 | 1.108 | 12.3 % |
| kernel | kernel evaluation (exp on windows) | 500 | 0.192 | 2.1 % |
| prims | primitive variables / pressure | 300 | 0.180 | 2.0 % |
| sources | coarse D/E/R/tau_j with windowed kernels | 100 | 0.117 | 1.3 % |
| dt_cfl | CFL time step | 100 | 0.074 | 0.8 % |
| check | conservation assertion | 100 | 0.029 | 0.3 % |

## Hot spots, numba backend (cumulative share of step time)

| operator | what | calls | cum. s | share |
|---|---|---|---|---|
| _fine_advect_nb | fine material BFECC (numba kernels) | 75 | 3.340 | 37.8 % |
| fine_sources | fine source mixing in windows | 300 | 2.387 | 27.0 % |
| _rhs_nb | coarse fluxes (numba kernel), x2 per step | 600 | 2.009 | 22.7 % |
| kernel | kernel evaluation (exp on windows) | 1500 | 0.605 | 6.8 % |
| prims | primitive variables / pressure | 900 | 0.532 | 6.0 % |
| sources | coarse D/E/R/tau_j with windowed kernels | 300 | 0.341 | 3.9 % |
| dt_cfl | CFL time step | 300 | 0.209 | 2.4 % |
| check | conservation assertion | 300 | 0.086 | 1.0 % |

peak RSS: 279 MB
