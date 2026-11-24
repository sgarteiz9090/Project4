"""Scaled design problem for the busbar: mass objective, three constraints, adjoint gradients.

    x = d / d0,  f = m / m0,
    g_bar   = KS over all nodes and all time levels of (T_C / 105 C)        - 1 <= 0
    g_joint = KS over the joint-face nodes and all time levels of (T_C / 70 C) - 1 <= 0
    g_R     = R(steady continuous rating) / R_max                            - 1 <= 0
The history functionals are differentiated by the generic backward sweep `busbar.adjoint`
(written first for the brake-disc study); the state is the full coupled vector u_n.
"""
from __future__ import annotations

import time
import numpy as np

from .adjoint import adjoint_gradients
from .ks import ks_history, true_max
from .params import Params, PHI_SCALE

C0 = 273.15


def constraint_spec(p: Params | None = None):
    p = Params() if p is None else p
    return [
        {"name": "peak_bar_T", "kind": "ks", "dofs": "temp_dofs", "limit": p.T_bar_max - C0,
         "limit_C": p.T_bar_max - C0, "unit": "C", "rho_per_K": p.ks_rho, "label": "Peak busbar temperature (insulation)"},
        {"name": "peak_joint_T", "kind": "ks", "dofs": "joint_dofs", "limit": p.T_joint_max - C0,
         "limit_C": p.T_joint_max - C0, "unit": "C", "rho_per_K": p.ks_rho, "label": "Peak joint / cell-terminal temperature"},
        {"name": "resistance", "kind": "resistance", "limit": p.R_max * 1e6, "limit_C": p.R_max * 1e6, "unit": "micro-ohm",
         "label": "Terminal-to-terminal resistance at the continuous rating"},
    ]


class BusbarSolver:
    """Binds busbar.api with a fixed physics flag and parameter set ."""
    nvar = 24

    def __init__(self, nonlinear=True, params: Params | None = None):
        from . import api
        self.api, self.p, self.nonlinear = api, params, nonlinear
        self.kw = {"nonlinear": nonlinear, "p": params}
        self.x0 = np.asarray(api.x0(params), float)
        lb, ub = api.bounds(params)
        self.lb, self.ub = np.asarray(lb, float), np.asarray(ub, float)

    def forward(self, d):
        return self.api.forward(d, **self.kw)

    def residual(self, d, u_new, u_old, n):
        return self.api.residual(d, u_new, u_old, n, **self.kw)

    def jac_T(self, d, u_new, u_old, n):
        return self.api.jac_T(d, u_new, u_old, n, **self.kw)

    def jac_Told(self, d, u_new, u_old, n):
        return self.api.jac_Told(d, u_new, u_old, n, **self.kw)

    def vjp_d_all(self, d, U, lam):
        return self.api.vjp_d_all(d, U, lam, **self.kw)

    def mass(self, d):
        return self.api.mass(d, self.p)

    def dmass_dd(self, d):
        return self.api.dmass_dd(d, self.p)


class BusbarProblem:
    def __init__(self, solver: BusbarSolver, constraints=None):
        self.solver = solver
        self.cons = constraints or constraint_spec(solver.p)
        self.d0 = solver.x0.copy()
        self.m0 = solver.mass(self.d0)
        self._cache_x, self._cache = None, None
        self.n_forward = self.n_adjoint = 0

    def d_of(self, x):
        return np.asarray(x, float) * self.d0

    def evaluate(self, x, need_grad=True):
        x = np.asarray(x, float)
        if self._cache_x is not None and np.array_equal(x, self._cache_x) and \
                (not need_grad or "grad_g" in self._cache):
            return self._cache
        d = self.d_of(x)
        t0 = time.perf_counter()
        fw = self.solver.forward(d)
        t_fwd = time.perf_counter() - t0
        self.n_forward += 1
        U = fw.U
        g, dU, vmax = [], [], []
        for c in self.cons:
            if c["kind"] == "ks":
                dofs = getattr(fw, c["dofs"])
                ks, dks = ks_history(U, dofs, c["limit"], rho=c["rho_per_K"] * c["limit"], offset=C0)
                g.append(ks - 1.0); dU.append(dks)
                vmax.append(true_max(U, dofs, offset=C0)[0])
            else:
                I1 = fw.sched.I[1]
                R = PHI_SCALE * U[1, fw.iV] / I1 * 1e6           # micro-ohm, steady rating step
                g.append(R / c["limit"] - 1.0)
                dr = np.zeros_like(U); dr[1, fw.iV] = PHI_SCALE / I1 * 1e6 / c["limit"]
                dU.append(dr); vmax.append(R)
        m = self.solver.mass(d)
        out = {"x": x.copy(), "d": d, "f": m / self.m0, "mass": m, "g": np.array(g),
               "Tmax": np.array(vmax), "T": U, "sol": fw, "t_forward": t_fwd}
        if need_grad:
            grads, info = adjoint_gradients(self.solver, d, U, np.array(dU))
            self.n_adjoint += 1
            out["grad_g"] = grads * self.d0
            out["grad_f"] = self.solver.dmass_dd(d) * self.d0 / self.m0
            out["t_adjoint"] = info["time_total"]
        self._cache_x, self._cache = x.copy(), out
        return out

    def f(self, x):
        return self.evaluate(x, need_grad=False)["f"]

    def df(self, x):
        return self.solver.dmass_dd(self.d_of(x)) * self.d0 / self.m0

    def c_ineq(self, x):
        return -self.evaluate(x)["g"]

    def dc_ineq(self, x):
        return -self.evaluate(x)["grad_g"]

    def bounds(self):
        return list(zip(self.solver.lb / self.d0, self.solver.ub / self.d0))
