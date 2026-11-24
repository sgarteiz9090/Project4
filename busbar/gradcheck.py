"""Gradient verification: adjoint vs finite differences, Taylor test, cost comparison."""
from __future__ import annotations

import time

import numpy as np


def _values(prob, d):
    x = d / prob.d0
    e = prob.evaluate(x, need_grad=False)
    return np.concatenate([[e["f"]], e["g"]])


def fd_table(prob, x, var_indices, h_rel=1e-5):
    """Central differences of [f, g_1..g_m] w.r.t. selected design variables.

    Returns list of rows {var, fun, adjoint, fd, rel_err}. Derivatives are w.r.t.
    the physical variable d_j (per metre), which is what the report quotes.
    """
    e = prob.evaluate(x, need_grad=True)
    d = prob.d_of(x)
    adj = np.vstack([e["grad_f"], e["grad_g"]]) / prob.d0     # back to d/dd
    names = ["mass/m0"] + [c["name"] for c in prob.cons]
    rows = []
    for j in var_indices:
        h = h_rel * abs(d[j])
        dp, dm = d.copy(), d.copy()
        dp[j] += h; dm[j] -= h
        fd = (_values(prob, dp) - _values(prob, dm)) / (2 * h)
        for k, nm in enumerate(names):
            a = adj[k, j]
            scale = np.abs(adj[k]).max()            # size of this gradient
            rel = abs(a - fd[k]) / abs(fd[k]) if fd[k] != 0 else float("inf")
            rows.append({"var": int(j), "fun": nm, "adjoint": float(a), "fd": float(fd[k]),
                         "rel_err": float(rel),
                         "err_vs_grad_norm": float(abs(a - fd[k]) / scale) if scale else 0.0})
    prob.evaluate(x, need_grad=True)   # restore cache
    return rows


def taylor_test(prob, x, v=None, eps_list=None, seed=0):
    """Central-difference directional derivative error vs step size.

    err(eps) = | (g(d + eps v) - g(d - eps v)) / (2 eps) - grad g . v |
    should fall as eps^2 until round-off. Also returns the first-order remainder
    | g(d + eps v) - g(d) - eps grad g . v | which falls as eps^2 too.
    v is a random unit direction in scaled space x (so every variable moves).
    """
    rng = np.random.default_rng(seed)
    if v is None:
        v = rng.standard_normal(len(x)); v /= np.linalg.norm(v)
    if eps_list is None:
        eps_list = np.logspace(-1, -8, 15)
    e0 = prob.evaluate(x, need_grad=True)
    G = np.vstack([e0["grad_f"], e0["grad_g"]])
    base = np.concatenate([[e0["f"]], e0["g"]])
    dirder = G @ v
    out = {"eps": list(map(float, eps_list)), "dirder_adjoint": dirder.tolist(),
           "central_err": [], "remainder": []}
    for eps in eps_list:
        p = _values(prob, prob.d_of(x + eps * v))
        m = _values(prob, prob.d_of(x - eps * v))
        out["central_err"].append((np.abs((p - m) / (2 * eps) - dirder)).tolist())
        out["remainder"].append((np.abs(p - base - eps * dirder)).tolist())
    prob.evaluate(x, need_grad=True)
    return out


def cost_comparison(prob, x):
    """Wall time: one forward + one adjoint sweep vs central FD over all variables."""
    prob._cache_x = None
    t0 = time.perf_counter(); e = prob.evaluate(x, need_grad=False); t_fwd = time.perf_counter() - t0
    prob._cache_x = None
    t0 = time.perf_counter(); e = prob.evaluate(x, need_grad=True); t_both = time.perf_counter() - t0
    n = len(x)
    return {"n_var": n, "n_constraints": len(prob.cons),
            "t_forward_s": t_fwd, "t_forward_plus_adjoint_s": t_both,
            "t_adjoint_only_s": e["t_adjoint"],
            "t_fd_central_est_s": 2 * n * t_fwd,
            "t_fd_forward_est_s": (n + 1) * t_fwd,
            "speedup_vs_central": 2 * n * t_fwd / t_both}
