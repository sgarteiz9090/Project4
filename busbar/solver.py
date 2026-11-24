"""Implicit Euler in time, full Newton on the coupled (psi, T, V_in) system at every step."""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np
import scipy.sparse.linalg as spla

from .model import BusbarModel


class NewtonError(RuntimeError):
    pass


@dataclass
class Solution:
    t: np.ndarray
    U: np.ndarray            # (N+1, ndof)
    newton: list = field(default_factory=list)


def initial_state(model: BusbarModel):
    u = np.zeros(model.ndof)
    u[model.nn:2 * model.nn] = model.p.T_amb
    return u


def newton_step(model, u_old, dt, n, u_guess=None, tol=None, maxit=None):
    """Solve R_n(u, u_old) = 0. Stops when the update is at round-off level (relative 1e-13 per
    field) or the residual falls below tol * initial; one extra iteration is then taken."""
    p = model.p
    tol = p.newton_tol if tol is None else tol
    maxit = p.newton_maxit if maxit is None else maxit
    u = np.array(u_old if u_guess is None else u_guess, dtype=float)
    nn = model.nn
    hist, ref, extra = [], None, False
    for _ in range(maxit):
        R = model.residual(u, u_old, dt, n)
        rn = float(np.linalg.norm(R))
        hist.append(rn)
        ref = max(rn, 1e-300) if ref is None else ref
        if extra:
            return u, hist
        du = spla.spsolve(model.jacobian(u, u_old, dt, n).tocsc(), -R)
        u = u + du
        sp_ = np.max(np.abs(du[:nn])) <= 1e-13 * max(np.max(np.abs(u[:nn])), 1e-300)
        sT = np.max(np.abs(du[nn:2 * nn])) <= 1e-13 * np.max(np.abs(u[nn:2 * nn]))
        if (sp_ and sT) or rn <= tol * ref:
            extra = True
    raise NewtonError(f"Newton did not converge at step {n}: {hist}")


def solve_transient(model: BusbarModel, u0=None) -> Solution:
    s = model.sched
    N = s.n_steps
    U = np.empty((N + 1, model.ndof))
    U[0] = initial_state(model) if u0 is None else u0
    hist = [[]]
    for n in range(1, N + 1):
        U[n], h = newton_step(model, U[n - 1], s.dt[n], n)
        hist.append(h)
    return Solution(t=s.t.copy(), U=U, newton=hist)
