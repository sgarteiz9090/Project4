"""WS2+3 stage of `python scripts/run_all.py`: sensitivities and optimization.

    from optim.run_ws23 import run
    blocks = run(quick=False)        # dict with keys "sensitivities", "optimization",
                                     # "designs_optimized_d", "designs_optimized_linear_d"

Also runnable alone (writes out/ws23.json and the figures):
    python -m optim.run_ws23 [--quick]
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

from .design_problem import DesignProblem
from .gradcheck import cost_comparison, fd_table, taylor_test
from .optimize import run_slsqp

N_ST = 12


def constraint_spec(p=None):
    if p is None:
        from brakeopt.params import Params
        p = Params()
    return [
        {"name": "peak_face_T", "nodes_attr": "friction_nodes",
         "limit_C": p.T_face_max - 273.15, "rho_per_K": p.ks_rho},
        {"name": "peak_hub_T", "nodes_attr": "hub_nodes",
         "limit_C": p.T_hub_max - 273.15, "rho_per_K": p.ks_rho},
    ]


def make_problem(radiation=True, kT=True, cons=None):
    from .brake_adapter import BrakeSolver
    return DesignProblem(BrakeSolver(radiation=radiation, kT=kT), cons or constraint_spec())


def var_label(j):
    return f"vane[{j}]" if j < N_ST else f"cheek[{j - N_ST}]"


def sensitivities(prob, x, check_vars, label):
    rows = fd_table(prob, x, check_vars)
    for r in rows:
        r["var_label"] = var_label(r["var"])
        r["var_type"] = "vane thickness" if r["var"] < N_ST else "cheek thickness"
    tay = taylor_test(prob, x)
    e = prob.evaluate(x)
    return {"label": label, "radiation": prob.solver.radiation, "kT": prob.solver.kT,
            "functionals": ["mass/m0"] + [c["name"] for c in prob.cons],
            "values": [e["f"]] + e["g"].tolist(), "Tmax_C": e["Tmax"].tolist(),
            "fd_table": rows, "taylor": tay,
            "max_rel_err_vs_grad_norm": max(r["err_vs_grad_norm"] for r in rows),
            "gradient_x": {"mass/m0": e["grad_f"].tolist(),
                           **{c["name"]: e["grad_g"][k].tolist()
                              for k, c in enumerate(prob.cons)}}}


def _design_summary(prob, x):
    e = prob.evaluate(x)
    d = prob.d_of(x)
    return {"d_m": d.tolist(),
            "vane_thickness_mm": (1e3 * d[:N_ST]).tolist(),
            "cheek_thickness_mm": (1e3 * d[N_ST:]).tolist(),
            "mass_kg": float(e["mass"]),
            "constraints": [{"name": c["name"], "ks_value_C": float((1 + e["g"][k]) * c["limit_C"]),
                             "true_max_C": float(e["Tmax"][k]), "limit_C": c["limit_C"],
                             "pass": bool(e["Tmax"][k] <= c["limit_C"] + 1e-9)}
                            for k, c in enumerate(prob.cons)]}


def optimization(quick=False):
    out = {"algorithm": "SciPy SLSQP", "scaling": "x = d/d0; f = m/m0; g_k = KS_k(T/T_lim) - 1 <= 0",
           "starts": []}
    prob = make_problem(True, True)
    lb, ub = np.array(prob.bounds()).T
    starts = {"baseline": np.ones(len(lb)), "all_max": ub.copy(),
              "midpoint": 0.5 * (lb + ub)}
    if quick:
        starts = {"baseline": starts["baseline"]}
    best = None
    for name, xs in starts.items():
        r = run_slsqp(prob, xs, label=name, maxiter=40 if quick else 200)
        r["summary"] = _design_summary(prob, np.array(r["x_opt"]))
        out["starts"].append(r)
        feas = r["kkt"]["max_violation"] < 1e-4
        if feas and (best is None or r["history"][-1]["f"] < best["history"][-1]["f"]):
            best = r
    best = best or min(out["starts"], key=lambda r: r["kkt"]["max_violation"])
    out["best_start"] = best["label"]
    out["spread_mass_kg"] = [s["summary"]["mass_kg"] for s in out["starts"]]
    out["initial"] = _design_summary(prob, np.ones(len(lb)))
    x_nl = np.array(best["x_opt"])

    # N2 comparison: optimize the linearized model, then evaluate it with full physics
    lin = make_problem(False, False)
    rl = run_slsqp(lin, best["x_start"], label="linear_model", maxiter=40 if quick else 200)
    x_lin = np.array(rl["x_opt"])
    rl["summary_linear_model"] = _design_summary(lin, x_lin)
    rl["summary_full_physics"] = _design_summary(prob, x_lin)
    out["linear"] = rl
    out["n2_comparison"] = {
        "mass_nonlinear_opt_kg": best["summary"]["mass_kg"],
        "mass_linear_opt_kg": rl["summary_linear_model"]["mass_kg"],
        "linear_opt_checked_with_full_physics": rl["summary_full_physics"]["constraints"],
        "d_nonlinear_mm": (1e3 * prob.d_of(x_nl)).tolist(),
        "d_linear_mm": (1e3 * lin.d_of(x_lin)).tolist()}
    if not quick:
        out["robust"] = robust_block(x_nl, prob)
    return out, prob.d_of(x_nl), lin.d_of(x_lin)


def robust_block(x_nominal_opt, prob_nominal):
    """Bonus: worst case over sampled scenarios; nominal optimum re-checked, then re-optimized."""
    from .robust import RobustProblem
    rp = RobustProblem()
    check_nominal = _design_summary(rp, x_nominal_opt)
    r = run_slsqp(rp, x_nominal_opt, label="robust", maxiter=200)
    r["summary"] = _design_summary(rp, np.array(r["x_opt"]))
    r["nominal_optimum_under_scenarios"] = check_nominal
    r["scenarios"] = rp.scenario_params
    r["mass_penalty_kg"] = r["summary"]["mass_kg"] - _design_summary(prob_nominal, x_nominal_opt)["mass_kg"]
    return r


def run(quick=False, fig_dir=None):
    t0 = time.perf_counter()
    check_vars = [0, 4, 8, 11, 12, 16, 20, 23] if not quick else [0, 11, 12, 23]
    sens = {}
    p_nl = make_problem(True, True)
    x0 = np.ones(24)
    sens["nonlinear"] = sensitivities(p_nl, x0, check_vars, "radiation + k(T)")
    sens["linear"] = sensitivities(make_problem(False, False), x0, check_vars[:4], "linearized")
    sens["cost"] = cost_comparison(p_nl, x0)
    opt, d_nl, d_lin = optimization(quick)
    blocks = {"sensitivities": sens, "optimization": opt,
              "designs_optimized_d": d_nl.tolist(), "designs_optimized_linear_d": d_lin.tolist(),
              "wall_s": time.perf_counter() - t0}
    if fig_dir:
        make_figures(blocks, fig_dir)
    return blocks


def make_figures(blocks, fig_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(fig_dir, exist_ok=True)
    tay = blocks["sensitivities"]["nonlinear"]["taylor"]
    names = blocks["sensitivities"]["nonlinear"]["functionals"]
    eps = np.array(tay["eps"]); err = np.array(tay["central_err"])
    fig, ax = plt.subplots(figsize=(5.5, 4))
    for k, nm in enumerate(names[1:], start=1):
        ax.loglog(eps, np.maximum(err[:, k], 1e-18), "o-", label=nm)
    ref = err[0, 1] * (eps / eps[0]) ** 2
    ax.loglog(eps, ref, "k--", lw=1, label=r"slope 2 ($\epsilon^2$)")
    ax.set_xlabel(r"step $\epsilon$ (scaled design units)")
    ax.set_ylabel(r"|central FD $-$ adjoint $\nabla g\cdot v$|")
    ax.set_title("Taylor test (random direction v)")
    ax.legend(); ax.grid(True, which="both", alpha=.3); fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "taylor_test.png"), dpi=150); plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.5, 4))
    for s in blocks["optimization"]["starts"]:
        h = s["history"]
        ax.plot([q["iter"] for q in h], [q["mass_kg"] for q in h], "o-", ms=3, label=s["label"])
    ax.set_xlabel("SLSQP iteration"); ax.set_ylabel("rotor mass (kg)")
    ax.set_title("Convergence history"); ax.legend(); ax.grid(alpha=.3); fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "opt_history.png"), dpi=150); plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    b = run(quick=a.quick, fig_dir="out/fig")
    os.makedirs("out", exist_ok=True)
    with open("out/ws23.json", "w") as f:
        json.dump(b, f, indent=1, default=float)
    print(json.dumps({k: b["optimization"]["n2_comparison"][k] for k in
                      ("mass_nonlinear_opt_kg", "mass_linear_opt_kg")}, indent=1))
