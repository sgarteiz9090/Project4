"""Small stand-in problem with the same API as the brake-disc solver.

A 1D rod of 24 segments (design variable = thickness of each segment) heated by
pulsed frictional flux at x = 0, with temperature-dependent conductivity k(T),
convection and radiation on its faces, implicit Euler in time, Newton per step.
It exists so the adjoint, gradient checks and optimizer can be developed and
unit-tested before (and independently of) the 3D solver. Nothing in the report
comes from this model.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

SIGMA = 5.670374419e-8


class ToyRod:
    nvar = 24

    def __init__(self, radiation=True, kT=True, n_steps=60, dt=0.5):
        self.radiation, self.kT = radiation, kT
        self.L, self.w = 0.12, 0.02          # length, width [m]
        self.rho, self.c = 7200.0, 500.0     # cast iron
        self.k0, self.beta = 50.0, -4e-4     # k(T) = k0 (1 + beta (T - 293))
        self.h, self.eps, self.Tinf = 60.0, 0.6, 293.15
        self.n_steps, self.dt = n_steps, dt
        self.ne = self.nvar
        self.nn = self.ne + 1
        self.le = self.L / self.ne
        self.x0 = np.full(self.nvar, 0.008)
        self.lb = np.full(self.nvar, 0.003)
        self.ub = np.full(self.nvar, 0.02)
        self.friction_nodes = np.arange(0, 3)
        self.hub_nodes = np.arange(self.nn - 3, self.nn)

    # loads: on/off braking pulses
    def q_in(self, n):
        t = n * self.dt
        return 4.0e5 if (t % 10.0) < 5.0 else 0.0   # W/m^2 on the end face

    def _k(self, T):
        return self.k0 * (1 + self.beta * (T - 293.15)) if self.kT else self.k0 + 0 * T

    def _dk(self, T):
        return self.k0 * self.beta + 0 * T if self.kT else 0 * T

    def _node_area_perim(self, d):
        A = d * self.w                        # element cross-section
        P = 2 * (self.w + d)                  # element perimeter
        Mn = np.zeros(self.nn, dtype=np.result_type(d, float))
        Sn = np.zeros_like(Mn)
        Mn[:-1] += 0.5 * A * self.le; Mn[1:] += 0.5 * A * self.le
        Sn[:-1] += 0.5 * P * self.le; Sn[1:] += 0.5 * P * self.le
        Aend = A[0]
        return A, Mn, Sn, Aend

    def residual(self, d, T_new, T_old, n):
        A, Mn, Sn, Aend = self._node_area_perim(d)
        T = T_new
        R = self.rho * self.c * Mn * (T - T_old) / self.dt
        Tm = 0.5 * (T[:-1] + T[1:])
        ke = self._k(Tm) * A / self.le
        flux = ke * (T[1:] - T[:-1])
        R[:-1] = R[:-1] - flux
        R[1:] = R[1:] + flux
        R = R + self.h * Sn * (T - self.Tinf)
        if self.radiation:
            R = R + self.eps * SIGMA * Sn * (T**4 - self.Tinf**4)
        R[0] = R[0] - self.q_in(n) * Aend
        return R

    def jac_T(self, d, T_new, T_old, n):
        A, Mn, Sn, _ = self._node_area_perim(d)
        T = T_new
        Tm = 0.5 * (T[:-1] + T[1:])
        ke = self._k(Tm) * A / self.le
        dke = 0.5 * self._dk(Tm) * A / self.le    # d ke / d T_i (same for both ends)
        dT = T[1:] - T[:-1]
        diag = self.rho * self.c * Mn / self.dt + self.h * Sn
        if self.radiation:
            diag = diag + 4 * self.eps * SIGMA * Sn * T**3
        rows, cols, vals = [], [], []
        i = np.arange(self.ne); j = i + 1
        # R_i -= flux, R_j += flux ; flux = ke (T_j - T_i)
        dflux_di = -ke + dke * dT
        dflux_dj = ke + dke * dT
        for r, sgn in ((i, -1.0), (j, 1.0)):
            rows += [r, r]; cols += [i, j]; vals += [sgn * dflux_di, sgn * dflux_dj]
        rows.append(np.arange(self.nn)); cols.append(np.arange(self.nn)); vals.append(diag)
        return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(self.nn, self.nn))

    def jac_Told(self, d, T_new, T_old, n):
        _, Mn, _, _ = self._node_area_perim(d)
        return sp.diags(-self.rho * self.c * Mn / self.dt)

    def forward(self, d, tol=1e-13, maxit=30):
        T = np.full(self.nn, self.Tinf)
        hist = [T.copy()]
        newton_iters = []
        for n in range(1, self.n_steps + 1):
            T_old = hist[-1]
            T = T_old.copy()
            for it in range(1, maxit + 1):
                R = self.residual(d, T, T_old, n)
                dT = spla.spsolve(self.jac_T(d, T, T_old, n).tocsc(), R)
                T = T - dT
                if np.linalg.norm(dT) < tol * np.linalg.norm(T):
                    break
            newton_iters.append(it)
            hist.append(T.copy())
        return {"T": np.array(hist), "newton_iters": newton_iters}

    def mass(self, d):
        return float(self.rho * np.sum(d * self.w * self.le))

    def dmass_dd(self, d):
        return np.full(len(d), self.rho * self.w * self.le)
