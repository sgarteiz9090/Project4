"""Implicit-Euler time stepping with a full Newton solve at every step."""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np
import scipy.sparse.linalg as spla

from .loads import Schedule, brake_schedule  # noqa: F401  (re-exported: solver.brake_schedule)
from .model import ThermalModel


class NewtonError(RuntimeError):
    pass


@dataclass
class Solution:
    t: np.ndarray                 # (N+1,)
    T: np.ndarray                 # (N+1, n_nodes) K
    newton: list = field(default_factory=list)   # per step: list of residual 2-norms


def newton_step(model: ThermalModel, T_old, dt, n, T_guess=None, tol=None, maxit=None):
    """Solve R_n(T, T_old) = 0. Converged when the update is below 1e-12 |T| (machine level) or the
    residual falls below tol * reference; one extra iteration is always taken after the residual
    test is met so the state is converged to round-off (needed for adjoint gradient checks)."""
    p = model.p
    tol = p.newton_tol if tol is None else tol
    maxit = p.newton_maxit if maxit is None else maxit
    T = np.array(T_old if T_guess is None else T_guess, dtype=float)
    hist = []
    ref = None
    extra = False
    for _ in range(maxit):
        R = model.residual(T, T_old, dt, n)
        rn = float(np.linalg.norm(R))
        hist.append(rn)
        if ref is None:
            ref = max(rn, 1e-300)
        if extra:
            return T, hist
        J = model.jacobian(T, dt, n)
        dT = spla.spsolve(J.tocsc(), -R)
        T = T + dT
        small_update = np.max(np.abs(dT)) <= 1e-12 * np.max(np.abs(T))
        if small_update or rn <= tol * ref:
            extra = True
            if small_update:
                hist.append(float(np.linalg.norm(model.residual(T, T_old, dt, n))))
                return T, hist
    raise NewtonError(f"Newton did not converge at step {n}: residual history {hist}")


def solve_transient(model: ThermalModel, sched: Schedule | None = None, T0=None) -> Solution:
    if sched is not None and sched is not model.sched:
        model.sched = sched
    s = model.sched
    N = s.n_steps
    T = np.empty((N + 1, model.nn))
    T[0] = model.p.T_init if T0 is None else T0
    hist = [[]]
    for n in range(1, N + 1):
        T[n], h = newton_step(model, T[n - 1], s.dt[n], n)
        hist.append(h)
    return Solution(t=s.t.copy(), T=T, newton=hist)
