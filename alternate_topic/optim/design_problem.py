"""Scaled optimization problem: objective, aggregated constraints, adjoint gradients.

Scaling (handout section 5): optimize x = d / d0, objective f = m(d) / m(d0),
constraints g_k = KS_k - 1 <= 0 with each KS built on T / T_limit, so every
quantity the optimizer sees is O(1).
"""
from __future__ import annotations

import time

import numpy as np

from .adjoint import adjoint_gradients
from .functionals import ks_history, true_max


class DesignProblem:
    """Wraps a solver object with the ainote.md API.

    constraints: list of dicts {name, nodes_attr, limit_C, rho or rho_per_K}. KS is built on
    a = T_C / limit_C, so rho (dimensionless) = rho_per_K * limit_C. nodes_attr names the
    attribute of the forward result (or the solver) that holds the node indices.
    """

    def __init__(self, solver, constraints, offset=273.15, steps=None):
        self.solver = solver
        self.cons = constraints
        self.offset = offset          # kelvin -> deg C if the solver works in K
        self.steps = steps
        self.d0 = np.asarray(solver.x0, dtype=float).copy()
        self.m0 = solver.mass(self.d0)
        self._cache_x = None
        self.n_forward = 0
        self.n_adjoint = 0
        self.history = []             # one entry per objective evaluation

    # --- helpers ---------------------------------------------------------
    def d_of(self, x):
        return np.asarray(x, dtype=float) * self.d0

    def _nodes(self, sol, attr):
        if isinstance(sol, dict) and attr in sol:
            return np.asarray(sol[attr])
        if hasattr(sol, attr):
            return np.asarray(getattr(sol, attr))
        return np.asarray(getattr(self.solver, attr))

    def _hist(self, sol):
        return np.asarray(sol["T"] if isinstance(sol, dict) else sol.T)

    def evaluate(self, x, need_grad=True):
        x = np.asarray(x, dtype=float)
        if self._cache_x is not None and np.array_equal(x, self._cache_x) and \
                (not need_grad or "grad_g" in self._cache):
            return self._cache
        d = self.d_of(x)
        t0 = time.perf_counter()
        sol = self.solver.forward(d)
        t_fwd = time.perf_counter() - t0
        self.n_forward += 1
        T = self._hist(sol)
        g, dT, tmax = [], [], []
        for c in self.cons:
            nodes = self._nodes(sol, c["nodes_attr"])
            rho = c["rho"] if "rho" in c else c.get("rho_per_K", 0.1) * c["limit_C"]
            ks, dks = ks_history(T, nodes, c["limit_C"], rho=rho,
                                 offset=self.offset, steps=self.steps)
            g.append(ks - 1.0)
            dT.append(dks)
            tmax.append(true_max(T, nodes, offset=self.offset, steps=self.steps)[0])
        out = {"x": x.copy(), "d": d, "f": self.solver.mass(d) / self.m0,
               "mass": self.solver.mass(d), "g": np.array(g), "Tmax": np.array(tmax),
               "T": T, "sol": sol, "t_forward": t_fwd, "_dT": dT}
        if need_grad:
            grads, info = adjoint_gradients(self.solver, d, T, np.array(dT))
            self.n_adjoint += 1
            out["grad_g"] = grads * self.d0          # chain rule d/dx = d0 * d/dd
            out["grad_f"] = self.solver.dmass_dd(d) * self.d0 / self.m0
            out["t_adjoint"] = info["time_total"]
        self._cache_x, self._cache = x.copy(), out
        return out

    # --- scipy callbacks --------------------------------------------------
    def f(self, x):
        return self.evaluate(x, need_grad=False)["f"]

    def df(self, x):
        d = self.d_of(x)
        return self.solver.dmass_dd(d) * self.d0 / self.m0

    def c_ineq(self, x):      # SLSQP convention: c(x) >= 0
        return -self.evaluate(x)["g"]

    def dc_ineq(self, x):
        return -self.evaluate(x)["grad_g"]

    def bounds(self):
        return list(zip(self.solver.lb / self.d0, self.solver.ub / self.d0))
