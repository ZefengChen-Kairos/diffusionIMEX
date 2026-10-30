"""Phase-two GPU backend (NVIDIA Warp) -- interface placeholder.

Phase one delivers the fp64 numpy reference (``core/solver_np.py``, backend="numpy") and a
numba backend whose kernels (``core/kernels_nb.py``) are written as plain loop nests over
cells and faces.  Each of them maps one-to-one onto a ``@wp.kernel``:

    kernels_nb.rhs_kernel           -> flux_x_kernel + flux_y_kernel + wall_kernel (dim=(B, N, N))
    kernels_nb.trace_kernel         -> trace_kernel (shared back-trace for all material fields)
    kernels_nb.sample_kernel        -> sample_kernel
    kernels_nb.bfecc_finish_kernel  -> bfecc_finish_kernel
    Sim.sources / fine_sources      -> source_window_kernel launched on the footprint window only

Planned (section 6): fp32 fields + fp64 ledger accumulators, batch dimension B for patterns /
parameter populations, CUDA-graph capture of a fixed-dt segment, wp.Tape recording for dI/dchi.
The equivalence target is the numpy reference (IoU >= 0.99, relative L2 <= 1e-3 at t = 2 s and
at the end).  Nothing here is used by the phase-one gates.
"""


def available():
    try:
        import warp  # noqa: F401
        return True
    except Exception:
        return False


class WarpSim:
    def __init__(self, *a, **k):
        raise NotImplementedError("The Warp backend is phase two; use core.solver_np.Sim "
                                  "(backend='numpy' reference or 'numba').")
