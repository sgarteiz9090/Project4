"""Functional facade used by the adjoint / optimizer workstream (WS2+3). Contract: ainote.md.

All functions take the design vector d (24,) in metres. `radiation` and `kT` switch the two
nonlinear terms (see Params). Step index n runs 1..N (step from t[n-1] to t[n]).
Partial derivatives with respect to d are computed by complex step through the mesh generator and
the residual (exact to round-off; the residual code is complex safe).
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import numpy as np

from .functionals import mass as _mass
from .loads import Schedule
from .mesh import Mesh, build_mesh
from .model import ThermalModel
from .params import Params, initial_design, design_bounds, N_DESIGN
from .solver import solve_transient

CS_STEP = 1e-30
BASE = Params()


@dataclass
class Forward:
    T: np.ndarray             # (N+1, n_nodes) K
    t: np.ndarray             # (N+1,)
    dt: np.ndarray            # (N+1,), dt[0] = nan
    friction_nodes: np.ndarray
    hub_nodes: np.ndarray
    mesh: Mesh
    model: ThermalModel
    sched: Schedule
    newton: list


def _params(radiation, kT, p):
    p = BASE if p is None else p
    return p.replace(radiation=bool(radiation), kT=bool(kT))


@lru_cache(maxsize=16)
def _model_cached(dkey, radiation, kT, pkey):
    p = _params(radiation, kT, _PCACHE[pkey])
    d = np.frombuffer(dkey, dtype=float)
    return ThermalModel(build_mesh(d, p), p)


_PCACHE: dict = {}


def _pkey(p):
    p = BASE if p is None else p
    key = repr(sorted(p.to_dict().items()))
    _PCACHE[key] = p
    return key


def model_for(d, radiation=True, kT=True, p: Params | None = None) -> ThermalModel:
    d = np.asarray(d)
    if np.iscomplexobj(d):
        q = _params(radiation, kT, p)
        return ThermalModel(build_mesh(d, q), q)
    return _model_cached(np.ascontiguousarray(d, dtype=float).tobytes(), bool(radiation), bool(kT), _pkey(p))


def forward(d, radiation=True, kT=True, p: Params | None = None) -> Forward:
    model = model_for(d, radiation, kT, p)
    sol = solve_transient(model)
    m = model.mesh
    return Forward(T=sol.T, t=sol.t, dt=model.sched.dt, friction_nodes=m.face_nodes("friction"),
                   hub_nodes=m.face_nodes("hub"), mesh=m, model=model, sched=model.sched, newton=sol.newton)


def residual(d, T_new, T_old, n, radiation=True, kT=True, p=None):
    model = model_for(d, radiation, kT, p)
    return model.residual(T_new, T_old, model.sched.dt[n], n)


def jac_T(d, T_new, T_old, n, radiation=True, kT=True, p=None):
    model = model_for(d, radiation, kT, p)
    return model.jacobian(np.asarray(T_new, dtype=float), model.sched.dt[n], n)


def jac_Told(d, T_new, T_old, n, radiation=True, kT=True, p=None):
    model = model_for(d, radiation, kT, p)
    return -model.mass_matrix() / model.sched.dt[n]


@lru_cache(maxsize=4)
def _cs_models(dkey, radiation, kT, pkey):
    d = np.frombuffer(dkey, dtype=float)
    q = _params(radiation, kT, _PCACHE[pkey])
    out = []
    for j in range(N_DESIGN):
        dc = d.astype(complex)
        dc[j] += 1j * CS_STEP
        out.append(ThermalModel(build_mesh(dc, q), q))
    return out


def _cs(d, radiation, kT, p):
    return _cs_models(np.ascontiguousarray(d, dtype=float).tobytes(), bool(radiation), bool(kT), _pkey(p))


def vjp_d(d, T_new, T_old, n, lam, radiation=True, kT=True, p=None):
    """lam^T dR_n/dd, (24,)."""
    out = np.empty(N_DESIGN)
    for j, mc in enumerate(_cs(d, radiation, kT, p)):
        out[j] = np.imag(lam @ mc.residual(T_new, T_old, mc.sched.dt[n], n)) / CS_STEP
    return out


def vjp_d_all(d, T_hist, lam_hist, radiation=True, kT=True, p=None):
    """sum_{n=1..N} lam_n^T dR_n/dd, (24,). T_hist, lam_hist: (N+1, n_nodes); lam_hist[0] unused."""
    T_hist = np.asarray(T_hist, dtype=float)
    lam = np.asarray(lam_hist, dtype=float)[1:]
    N = T_hist.shape[0] - 1
    steps = np.arange(1, N + 1)
    out = np.empty(N_DESIGN)
    for j, mc in enumerate(_cs(d, radiation, kT, p)):
        R = mc.residual(T_hist[1:], T_hist[:-1], mc.sched.dt[1:], steps)
        out[j] = np.imag(np.sum(lam * R)) / CS_STEP
    return out


def mass(d, p=None):
    q = BASE if p is None else p
    return float(_mass(build_mesh(np.asarray(d, dtype=float), q), q))


def dmass_dd(d, p=None):
    q = BASE if p is None else p
    d = np.asarray(d, dtype=float)
    g = np.empty(N_DESIGN)
    for j in range(N_DESIGN):
        dc = d.astype(complex)
        dc[j] += 1j * CS_STEP
        g[j] = np.imag(_mass(build_mesh(dc, q), q)) / CS_STEP
    return g


def x0(p=None):
    return initial_design(BASE if p is None else p)


def bounds(p=None):
    return design_bounds(BASE if p is None else p)
