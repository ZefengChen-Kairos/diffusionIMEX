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
| numpy | 128 | 86.48 | 11.6 | 108.0 | 1.00x |
| numba | 128 | 26.68 | 37.5 | 33.3 | 3.24x |
| warp (GPU) | - | phase two | | | |

## Hot spots, numpy backend (cumulative share of step time)

| operator | what | calls | cum. s | share |
|---|---|---|---|---|
| fine_advect | fine material BFECC (numpy) | 25 | 4.996 | 57.8 % |
| rhs | coarse fluxes (MUSCL + pressure + kappa), x2 per step (RK2) | 200 | 2.164 | 25.0 % |
| fine_sources | fine source mixing in windows | 100 | 1.103 | 12.8 % |
| kernel | kernel evaluation (exp on windows) | 500 | 0.185 | 2.1 % |
| prims | primitive variables / pressure | 300 | 0.172 | 2.0 % |
| sources | coarse D/E/R/tau_j with windowed kernels | 100 | 0.104 | 1.2 % |
| dt_cfl | CFL time step | 100 | 0.068 | 0.8 % |
| check | conservation assertion | 100 | 0.029 | 0.3 % |

## Hot spots, numba backend (cumulative share of step time)

| operator | what | calls | cum. s | share |
|---|---|---|---|---|
| _fine_advect_nb | fine material BFECC (numba kernels) | 75 | 2.793 | 34.9 % |
| fine_sources | fine source mixing in windows | 300 | 2.524 | 31.5 % |
| _rhs_nb | coarse fluxes (numba kernel), x2 per step | 600 | 1.634 | 20.4 % |
| kernel | kernel evaluation (exp on windows) | 1500 | 0.632 | 7.9 % |
| prims | primitive variables / pressure | 900 | 0.497 | 6.2 % |
| sources | coarse D/E/R/tau_j with windowed kernels | 300 | 0.321 | 4.0 % |
| dt_cfl | CFL time step | 300 | 0.201 | 2.5 % |
| check | conservation assertion | 300 | 0.084 | 1.0 % |

peak RSS: 258 MB
