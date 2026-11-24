"""Functional facade for the adjoint and the optimizer (same contract as brakeopt.api).

State u_n = [psi, T, V_in] per time level (see busbar.model). `nonlinear` switches every
temperature dependence (sigma(T), k(T), h(dT), radiation) at once; use the same flag for every
call inside one gradient evaluation. Partial derivatives w.r.t. d use complex step through the
mesh generator and the residual (exact to round-off).
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import numpy as np

from .mesh import Mesh, build_mesh
from .model import BusbarModel, Schedule
from .params import Params, initial_design, design_bounds, N_DESIGN
from .solver import solve_transient

CS_STEP = 1e-30
BASE = Params()
_PCACHE: dict = {}


@dataclass
class Forward:
    U: np.ndarray            # (N+1, ndof)
    T: np.ndarray            # alias of U (the optimizer's history array is the full state)
    t: np.ndarray
    dt: np.ndarray
    temp_dofs: np.ndarray    # state indices of all nodal temperatures
    joint_dofs: np.ndarray   # state indices of the temperatures on the joint faces
    iV: int
    mesh: Mesh
    model: BusbarModel
    sched: Schedule
    newton: list


def _pkey(p):
    p = BASE if p is None else p
    key = repr(sorted(p.to_dict().items()))
    _PCACHE[key] = p
    return key


def _params(nonlinear, p):
    return (BASE if p is None else p).replace(nonlinear=bool(nonlinear))


@lru_cache(maxsize=16)
def _model_cached(dkey, nonlinear, pkey):
    q = _params(nonlinear, _PCACHE[pkey])
    return BusbarModel(build_mesh(np.frombuffer(dkey, dtype=float), q), q)


def model_for(d, nonlinear=True, p=None):
    d = np.asarray(d)
    if np.iscomplexobj(d):
        q = _params(nonlinear, p)
        return BusbarModel(build_mesh(d, q), q)
    return _model_cached(np.ascontiguousarray(d, dtype=float).tobytes(), bool(nonlinear), _pkey(p))


def forward(d, nonlinear=True, p=None) -> Forward:
    model = model_for(d, nonlinear, p)
    sol = solve_transient(model)
    m, nn = model.mesh, model.nn
    joints = np.unique(np.r_[m.face_nodes("contact_in"), m.face_nodes("contact_out")])
    return Forward(U=sol.U, T=sol.U, t=sol.t, dt=model.sched.dt, temp_dofs=np.arange(nn, 2 * nn),
                   joint_dofs=joints + nn, iV=model.iV, mesh=m, model=model, sched=model.sched,
                   newton=sol.newton)


def residual(d, u_new, u_old, n, nonlinear=True, p=None):
    model = model_for(d, nonlinear, p)
    return model.residual(u_new, u_old, model.sched.dt[n], n)


def jac_T(d, u_new, u_old, n, nonlinear=True, p=None):
    model = model_for(d, nonlinear, p)
    return model.jacobian(np.asarray(u_new, float), np.asarray(u_old, float), model.sched.dt[n], n)


def jac_Told(d, u_new, u_old, n, nonlinear=True, p=None):
    model = model_for(d, nonlinear, p)
    return -model.mass_matrix() / model.sched.dt[n]


@lru_cache(maxsize=4)
def _cs_models(dkey, nonlinear, pkey):
    d = np.frombuffer(dkey, dtype=float)
    q = _params(nonlinear, _PCACHE[pkey])
    out = []
    for j in range(N_DESIGN):
        dc = d.astype(complex)
        dc[j] += 1j * CS_STEP
        out.append(BusbarModel(build_mesh(dc, q), q))
    return out


def _cs(d, nonlinear, p):
    return _cs_models(np.ascontiguousarray(d, dtype=float).tobytes(), bool(nonlinear), _pkey(p))


def vjp_d(d, u_new, u_old, n, lam, nonlinear=True, p=None):
    out = np.empty(N_DESIGN)
    for j, mc in enumerate(_cs(d, nonlinear, p)):
        out[j] = np.imag(lam @ mc.residual(u_new, u_old, mc.sched.dt[n], n)) / CS_STEP
    return out


def vjp_d_all(d, U_hist, lam_hist, nonlinear=True, p=None):
    """sum_{n=1..N} lam_n^T dR_n/dd (24,)."""
    U = np.asarray(U_hist, dtype=float)
    lam = np.asarray(lam_hist, dtype=float)[1:]
    N = U.shape[0] - 1
    steps = np.arange(1, N + 1)
    out = np.empty(N_DESIGN)
    for j, mc in enumerate(_cs(d, nonlinear, p)):
        R = mc.residual(U[1:], U[:-1], mc.sched.dt[1:], steps)
        out[j] = np.imag(np.sum(lam * R)) / CS_STEP
    return out


def _mass_of(d, q):
    from brakeopt.fem import HexGeometry
    m = build_mesh(d, q)
    return q.rho * np.sum(HexGeometry(m.X, m.hexes).wdetJ) * m.multiplicity


def mass(d, p=None):
    return float(np.real(_mass_of(np.asarray(d, float), BASE if p is None else p)))


def dmass_dd(d, p=None):
    q = BASE if p is None else p
    d = np.asarray(d, dtype=float)
    g = np.empty(N_DESIGN)
    for j in range(N_DESIGN):
        dc = d.astype(complex); dc[j] += 1j * CS_STEP
        g[j] = np.imag(_mass_of(dc, q)) / CS_STEP
    return g


def x0(p=None):
    return initial_design(BASE if p is None else p)


def bounds(p=None):
    return design_bounds(BASE if p is None else p)
