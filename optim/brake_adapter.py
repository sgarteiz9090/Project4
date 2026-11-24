"""Binds brakeopt.api (WS1) to the interface optim.adjoint / optim.design_problem expect.

The radiation / kT switches are fixed per adapter instance, so every call inside one
gradient evaluation uses the same physics (ainote.md requirement).
"""
from __future__ import annotations

import numpy as np


class BrakeSolver:
    nvar = 24

    def __init__(self, radiation=True, kT=True, params=None):
        from brakeopt import api   # imported lazily so the toy tests run without it
        self.api = api
        self.kw = {"radiation": radiation, "kT": kT}
        if params is not None:
            self.kw["p"] = params
        self.radiation, self.kT = radiation, kT
        self.x0 = np.asarray(api.x0(), dtype=float)
        lb, ub = api.bounds()
        self.lb, self.ub = np.asarray(lb, float), np.asarray(ub, float)

    def forward(self, d):
        return self.api.forward(d, **self.kw)

    def residual(self, d, T_new, T_old, n):
        return self.api.residual(d, T_new, T_old, n, **self.kw)

    def jac_T(self, d, T_new, T_old, n):
        return self.api.jac_T(d, T_new, T_old, n, **self.kw)

    def jac_Told(self, d, T_new, T_old, n):
        return self.api.jac_Told(d, T_new, T_old, n, **self.kw)

    def vjp_d(self, d, T_new, T_old, n, lam):
        lam = np.asarray(lam)
        if lam.ndim == 1:
            return self.api.vjp_d(d, T_new, T_old, n, lam, **self.kw)
        return np.column_stack([self.api.vjp_d(d, T_new, T_old, n, lam[:, k], **self.kw)
                                for k in range(lam.shape[1])])

    def vjp_d_all(self, d, T_hist, lam_hist):
        return self.api.vjp_d_all(d, T_hist, lam_hist, **self.kw)

    def mass(self, d):
        return float(self.api.mass(d))

    def dmass_dd(self, d):
        return np.asarray(self.api.dmass_dd(d), dtype=float)
