"""SLSQP driver with per-iteration history and KKT read-out (active set, multipliers)."""
from __future__ import annotations

import time

import numpy as np
from scipy.optimize import minimize, nnls


def run_slsqp(prob, x_start, maxiter=200, ftol=1e-7, label="start"):
    hist = []

    def cb(xk):
        e = prob.evaluate(xk)
        hist.append({"iter": len(hist) + 1, "f": float(e["f"]), "mass_kg": float(e["mass"]),
                     "g": e["g"].tolist(), "Tmax_C": e["Tmax"].tolist(),
                     "max_violation": float(max(0.0, e["g"].max()))})

    e0 = prob.evaluate(np.asarray(x_start, float))
    hist.append({"iter": 0, "f": float(e0["f"]), "mass_kg": float(e0["mass"]),
                 "g": e0["g"].tolist(), "Tmax_C": e0["Tmax"].tolist(),
                 "max_violation": float(max(0.0, e0["g"].max()))})
    nf0, na0 = prob.n_forward, prob.n_adjoint
    t0 = time.perf_counter()
    res = minimize(prob.f, x_start, jac=prob.df, method="SLSQP", bounds=prob.bounds(),
                   constraints=[{"type": "ineq", "fun": prob.c_ineq, "jac": prob.dc_ineq}],
                   options={"maxiter": maxiter, "ftol": ftol}, callback=cb)
    wall = time.perf_counter() - t0
    kkt = kkt_readout(prob, res.x)
    return {"label": label, "x_start": list(map(float, x_start)), "x_opt": res.x.tolist(),
            "d_opt_m": prob.d_of(res.x).tolist(), "success": bool(res.success),
            "message": str(res.message), "nit": int(res.nit), "wall_s": wall,
            "n_forward": prob.n_forward - nf0, "n_adjoint": prob.n_adjoint - na0,
            "slsqp_multipliers": _jsonable(getattr(res, "multipliers", None)),
            "history": hist, "kkt": kkt}


def _jsonable(m):
    if m is None:
        return None
    if isinstance(m, dict):
        return {k: np.asarray(v).tolist() for k, v in m.items()}
    return np.asarray(m).tolist()


def kkt_readout(prob, x, act_tol=1e-3, bnd_tol=1e-6):
    """Active set and multipliers from grad f + sum mu_k grad g_k + bound terms = 0.

    Multipliers are fitted by nonnegative least squares on the active set, so they
    can be read even if the optimizer's own estimates are unavailable. In scaled
    units mu_k = d(f)/d(relaxation of g_k): relaxing limit k by 1 % lowers the mass
    ratio by about 0.01 mu_k, i.e. 0.01 mu_k m0 kg.
    """
    e = prob.evaluate(x)
    gf, G, g = e["grad_f"], e["grad_g"], e["g"]
    lb, ub = np.array(prob.bounds()).T
    active_c = [k for k in range(len(g)) if g[k] > -act_tol]
    at_lb = np.where(x <= lb + bnd_tol * np.maximum(1, abs(lb)))[0]
    at_ub = np.where(x >= ub - bnd_tol * np.maximum(1, abs(ub)))[0]
    n = len(x)
    cols, labels = [], []
    for k in active_c:
        cols.append(G[k]); labels.append(("con", k))
    for i in at_lb:                       # lower bound lb - x <= 0, gradient -e_i
        v = np.zeros(n); v[i] = -1; cols.append(v); labels.append(("lb", int(i)))
    for i in at_ub:
        v = np.zeros(n); v[i] = 1; cols.append(v); labels.append(("ub", int(i)))
    if cols:
        A = np.column_stack(cols)
        mu, resid = nnls(A, -gf)
    else:
        mu, resid = np.zeros(0), float(np.linalg.norm(gf))
    con_mu = {prob.cons[k]["name"]: float(mu[i]) for i, (t, k) in enumerate(labels) if t == "con"}
    stat = float(resid / max(np.linalg.norm(gf), 1e-300))
    return {"g": g.tolist(), "Tmax_C": e["Tmax"].tolist(),
            "active_constraints": [prob.cons[k]["name"] for k in active_c],
            "constraint_multipliers": con_mu,
            "kg_per_percent_relaxation": {k: v * 0.01 * prob.m0 for k, v in con_mu.items()},
            "lower_bounds_active": at_lb.tolist(), "upper_bounds_active": at_ub.tolist(),
            "bound_multipliers": {f"{t}{k}": float(mu[i]) for i, (t, k) in enumerate(labels)
                                  if t != "con"},
            "stationarity_residual_rel": stat,
            "max_violation": float(max(0.0, g.max()))}
