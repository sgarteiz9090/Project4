"""Discrete adjoint for an implicit-Euler, Newton-solved transient problem.

State equation, one residual per time step n = 1..N (T_0 is the fixed initial field):

    R_n(T_n, T_{n-1}, d) = 0

Functional of the whole history (for example a KS aggregate over space and time):

    J(d) = Phi(T_0, ..., T_N, d)

Lagrangian  L = Phi + sum_n lam_n^T R_n.  Setting dL/dT_n = 0 gives the backward sweep

    A_n^T lam_n = -(dPhi/dT_n)^T - B_{n+1}^T lam_{n+1},      lam_{N+1} = 0,
    A_n = dR_n/dT_n   (the Newton tangent at the converged step),
    B_n = dR_n/dT_{n-1},

and the total derivative is

    dJ/dd = dPhi/dd + sum_n lam_n^T dR_n/dd.

The sweep runs backward through the load history (the "path / time dependence"
part of the rubric). Several functionals share one factorization of A_n per step.

The problem object `prob` must provide (see ainote.md, "Solver API"):
    prob.residual(d, T_new, T_old, n)   -> (ndof,) array, complex-step safe in d
    prob.jac_T(d, T_new, T_old, n)      -> sparse (ndof, ndof)
    prob.jac_Told(d, T_new, T_old, n)   -> sparse (ndof, ndof)
    optional prob.vjp_d(d, T_new, T_old, n, lam) -> (nvar,) or (nvar, nfun) for lam (ndof, nfun)
    optional prob.vjp_d_all(d, T_hist, lam_hist) -> (nvar,) = sum_n lam_n^T dR_n/dd (preferred)
"""
from __future__ import annotations

import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

CS_STEP = 1e-30  # complex-step size; exact to machine precision, no subtraction error


def dR_dd_complex_step(prob, d, T_new, T_old, n):
    """Dense dR_n/dd (ndof, nvar) by complex step, one residual call per variable."""
    d = np.asarray(d, dtype=float)
    cols = []
    for j in range(d.size):
        dc = d.astype(complex)
        dc[j] += 1j * CS_STEP
        cols.append(np.imag(prob.residual(dc, T_new, T_old, n)) / CS_STEP)
    return np.column_stack(cols)


def _vjp_d(prob, d, T_new, T_old, n, Lam):
    """Lam^T dR_n/dd for Lam of shape (ndof, nfun); returns (nvar, nfun)."""
    if hasattr(prob, "vjp_d"):
        out = np.asarray(prob.vjp_d(d, T_new, T_old, n, Lam))
        return out.reshape(len(d), Lam.shape[1])
    return dR_dd_complex_step(prob, d, T_new, T_old, n).T @ Lam


def adjoint_gradients(prob, d, T_hist, dPhi_dT, dPhi_dd=None, return_lambdas=False):
    """Gradients of several history functionals by one backward adjoint sweep.

    T_hist  : (N+1, ndof) converged temperatures, T_hist[0] = initial field.
    dPhi_dT : (nfun, N+1, ndof) partial derivatives of each functional w.r.t. T_n
              (row 0 is ignored: T_0 does not depend on d).
    dPhi_dd : (nfun, nvar) explicit partials, or None.
    Returns grads (nfun, nvar) and timing info.
    """
    d = np.asarray(d, dtype=float)
    dPhi_dT = np.asarray(dPhi_dT)
    nfun, Np1, ndof = dPhi_dT.shape
    N = Np1 - 1
    grads = np.zeros((nfun, d.size)) if dPhi_dd is None else np.array(dPhi_dd, dtype=float)
    lam_next = np.zeros((ndof, nfun))
    lams = [None] * (N + 1)
    use_all = hasattr(prob, "vjp_d_all")
    lam_hist = np.zeros((nfun, N + 1, ndof)) if use_all else None   # row 0 unused
    t0 = time.perf_counter()
    t_vjp = 0.0
    for n in range(N, 0, -1):
        rhs = -dPhi_dT[:, n, :].T
        if n < N:
            B_next = prob.jac_Told(d, T_hist[n + 1], T_hist[n], n + 1)
            rhs = rhs - B_next.T @ lam_next
        A = prob.jac_T(d, T_hist[n], T_hist[n - 1], n)
        lu = spla.splu(sp.csc_matrix(A.T))
        lam = lu.solve(np.ascontiguousarray(rhs))
        lam = lam.reshape(ndof, nfun)
        if use_all:
            lam_hist[:, n, :] = lam.T
        else:
            tv = time.perf_counter()
            grads += _vjp_d(prob, d, T_hist[n], T_hist[n - 1], n, lam).T
            t_vjp += time.perf_counter() - tv
        lam_next = lam
        if return_lambdas:
            lams[n] = lam
    if use_all:
        # one vectorized pass over all steps: sum_n lam_n^T dR_n/dd, per functional
        tv = time.perf_counter()
        for k in range(nfun):
            grads[k] += np.asarray(prob.vjp_d_all(d, T_hist, lam_hist[k])).ravel()
        t_vjp += time.perf_counter() - tv
    info = {"time_total": time.perf_counter() - t0, "time_vjp": t_vjp, "n_steps": N}
    if return_lambdas:
        info["lambdas"] = lams
    return grads, info
