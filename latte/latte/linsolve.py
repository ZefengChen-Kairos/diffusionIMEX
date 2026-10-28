"""SPD linear solves with true-residual reporting (V3 Sec. 7.3: keep real residuals,
iteration counts and failure status)."""
import numpy as np
from scipy import sparse
from scipy.sparse import linalg as spla

from .failures import SimFailure


class SPDSolver:
    """method: 'direct' (sparse LU), 'amg' (pyamg smoothed aggregation as CG preconditioner,
    hierarchy reused while it stays effective), 'jacobi' (diagonally preconditioned CG),
    'auto' (direct below `direct_max` unknowns, else amg)."""

    def __init__(self, method="auto", rtol=1e-12, maxiter=2000, direct_max=4000,
                 reuse_iters=40, fail_code="PRESSURE_FAILED", accept_rtol=None):
        self.method = method
        self.rtol = rtol
        self.maxiter = maxiter
        self.direct_max = direct_max
        self.reuse_iters = reuse_iters
        self.fail_code = fail_code
        self.accept_rtol = accept_rtol if accept_rtol is not None else 100 * rtol
        self._ml = None
        self._ml_shape = None

    def solve(self, A, b, x0=None):
        x, info = self.solve_many(A, [b], [x0])
        return x[0], info

    def solve_many(self, A, bs, x0s=None):
        """Solve A x = b for several right-hand sides with one factorisation / hierarchy."""
        A = sparse.csr_matrix(A)
        n = A.shape[0]
        method = self.method
        if method == "auto":
            method = "direct" if n <= self.direct_max else "amg"
        x0s = x0s if x0s is not None else [None] * len(bs)
        info = {"method": method, "iters": 0, "rel_res": 0.0, "abs_res": 0.0}
        xs = []
        lu = None
        M = None
        if method == "direct":
            lu = spla.splu(A.tocsc())
        elif method == "amg":
            import pyamg
            if self._ml is None or self._ml_shape != A.shape:
                self._ml = pyamg.smoothed_aggregation_solver(A, symmetry="symmetric", max_coarse=500)
                self._ml_shape = A.shape
            M = self._ml.aspreconditioner(cycle="V")
        elif method == "jacobi":
            d = A.diagonal()
            M = spla.LinearOperator(A.shape, matvec=lambda v: v / d)
        else:
            raise ValueError(method)
        for b, x0 in zip(bs, x0s):
            bnorm = float(np.linalg.norm(b))
            if bnorm == 0.0:
                xs.append(np.zeros(n))
                continue
            if lu is not None:
                x = lu.solve(b)
            else:
                it = [0]

                def cb(_):
                    it[0] += 1
                x, flag = spla.cg(A, b, x0=x0, rtol=self.rtol, atol=0.0, maxiter=self.maxiter, M=M, callback=cb)
                info["iters"] += it[0]
                if method == "amg" and it[0] > self.reuse_iters:
                    self._ml = None          # rebuild the hierarchy next time
            r = b - A @ x
            rel = float(np.linalg.norm(r) / bnorm)
            info["rel_res"] = max(info["rel_res"], rel)
            info["abs_res"] = max(info["abs_res"], float(np.abs(r).max()))
            if not np.all(np.isfinite(x)) or rel > self.accept_rtol:
                raise SimFailure(self.fail_code, f"linear solve not converged (rel residual {rel:.3e})",
                                 method=method, iters=info["iters"], rel_res=rel)
            xs.append(x)
        return xs, info
