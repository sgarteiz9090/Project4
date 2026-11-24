"""One command regenerates busbar/site/assets/results.json, the figures and every number on the
busbar site:

    python -m busbar.run_all            # full run (~1 h on a laptop)
    python -m busbar.run_all --quick    # smoke run

Stages: verification ladder, initial design, gradient checks (adjoint vs FD, Taylor test, cost),
SLSQP from three starts, N2 comparison (constant-property model) with attribution, design under
uncertainty (three current/ambient scenarios), and the design blocks the 3D viewer reads.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import subprocess
import sys
import time

import numpy as np

from optim.gradcheck import cost_comparison, fd_table, taylor_test
from optim.optimize import run_slsqp
from . import api
from .ladder import run_ladder
from .params import Params, initial_design, station_x, split_design, N_STATIONS, PHI_SCALE
from .problem import BusbarProblem, BusbarSolver, constraint_spec

HERE = os.path.dirname(os.path.abspath(__file__))
ASSETS = os.path.join(HERE, "site", "assets")
C0 = 273.15


def _r(x, nd):
    return np.round(np.asarray(x, dtype=float), nd).tolist()


def var_label(j):
    return f"width[{j}]" if j < N_STATIONS else f"thickness[{j - N_STATIONS}]"


def make_problem(nonlinear=True, p=None):
    return BusbarProblem(BusbarSolver(nonlinear, p), constraint_spec(p))


# ------------------------------------------------------------------ design block for the viewer
def design_block(d, p: Params | None = None, label="", max_frames=80):
    p = Params() if p is None else p
    d = np.asarray(d, float)
    fw = api.forward(d, True, p)
    M, m, U = fw.model, fw.mesh, fw.U
    ps, T, V = M.split(U)
    TC = T - C0
    jn = fw.joint_dofs - M.nn
    s = fw.sched
    R = PHI_SCALE * V / s.I * 1e6
    w, t = split_design(d)
    cons = [
        {"name": "peak_bar_T", "label": "Peak busbar temperature (insulation)", "value_C": round(float(TC[1:].max()), 2),
         "true_max_C": round(float(TC[1:].max()), 2), "limit_C": round(p.T_bar_max - C0, 2),
         "pass": bool(TC[1:].max() <= p.T_bar_max - C0 + 1e-9)},
        {"name": "peak_joint_T", "label": "Peak joint / cell-terminal temperature",
         "value_C": round(float(TC[1:, jn].max()), 2), "true_max_C": round(float(TC[1:, jn].max()), 2),
         "limit_C": round(p.T_joint_max - C0, 2), "pass": bool(TC[1:, jn].max() <= p.T_joint_max - C0 + 1e-9)},
        {"name": "resistance", "label": "Resistance at the continuous rating", "value": round(float(R[1]), 3),
         "limit": round(p.R_max * 1e6, 3), "unit": "µΩ", "sense": "<=", "pass": bool(R[1] <= p.R_max * 1e6 + 1e-9)},
    ]
    # current density magnitude per element at the pulse peak (tooltip / viewer)
    n_pk = int(np.argmax(TC[1:].max(1))) + 1
    g = M.geo
    gp = np.einsum("egai,ea->egi", g.dNdx, ps[n_pk][m.hexes]) * PHI_SCALE
    Tg = T[n_pk][m.hexes] @ g.N.T
    jmag = (M.sigma(Tg) * np.sqrt(np.sum(gp * gp, -1))).mean(1) / 1e6       # A/mm^2
    shown = np.arange(1, len(s.t))
    idx = shown[np.unique(np.round(np.linspace(0, len(shown) - 1, min(max_frames, len(shown)))).astype(int))]
    return {
        "label": label,
        "d": {"width_mm": _r(1e3 * w, 4), "thickness_mm": _r(1e3 * t, 4), "station_x_mm": _r(1e3 * station_x(p), 3)},
        "mass_kg": round(float(api.mass(d, p)), 5),
        "constraints": cons,
        "mesh": {"nodes": _r(np.real(m.X), 6), "hexes": m.hexes.tolist(),
                 "faces": {k: v.tolist() for k, v in m.faces.items()},
                 "symmetry": "half model: mirror y = 0"},
        "temperature": {"unit": "C", "times": _r(s.t[idx], 2), "frames": [_r(TC[i], 2) for i in idx]},
        "current_density": {"unit": "A/mm^2", "step": n_pk, "time": float(s.t[n_pk]), "elem": _r(jmag, 3)},
        "history": {"times": _r(s.t[1:], 2), "peak_bar_C": _r(TC[1:].max(1), 3),
                    "peak_joint_C": _r(TC[1:, jn].max(1), 3), "current_A": _r(s.I[1:], 1),
                    "resistance_uohm": _r(R[1:], 4), "mean_C": _r(TC[1:].mean(1), 3)},
        "power_W": {k: round(2 * v, 4) for k, v in M.heat_flows(U[1]).items()},
        "model": {"nonlinear": True, "n_nodes": m.n_nodes, "n_elem": m.n_elem, "n_steps": s.n_steps,
                  "ndof": M.ndof},
    }


# ------------------------------------------------------------------ sensitivities
def sensitivities(prob, x, check_vars, label):
    rows = fd_table(prob, x, check_vars)
    for r in rows:
        r["var_label"] = var_label(r["var"])
        r["var_type"] = "width" if r["var"] < N_STATIONS else "thickness"
    tay = taylor_test(prob, x)
    e = prob.evaluate(x)
    return {"label": label, "nonlinear": prob.solver.nonlinear,
            "functionals": ["mass/m0"] + [c["name"] for c in prob.cons],
            "values": [e["f"]] + e["g"].tolist(), "Tmax_C": e["Tmax"].tolist(), "fd_table": rows, "taylor": tay,
            "max_rel_err_vs_grad_norm": max(r["err_vs_grad_norm"] for r in rows)}


def summary(prob, x):
    e = prob.evaluate(x)
    d = prob.d_of(x)
    return {"d_m": d.tolist(), "width_mm": (1e3 * d[:N_STATIONS]).tolist(),
            "thickness_mm": (1e3 * d[N_STATIONS:]).tolist(), "mass_kg": float(e["mass"]),
            "constraints": [{"name": c["name"], "ks_or_value": float((1 + e["g"][k]) * c["limit"]),
                             "true_value": float(e["Tmax"][k]), "limit": c["limit"], "unit": c["unit"],
                             "pass": bool(e["Tmax"][k] <= c["limit"] + 1e-9)} for k, c in enumerate(prob.cons)]}


def n2_attribution(d):
    """The design evaluated with full physics and with each temperature dependence frozen in turn."""
    rows = []
    variants = [("full physics", {}), ("constant-property model (all frozen at 70 C)", {"nonlinear": False}),
                ("sigma frozen at 70 C", {"alpha": 0.0, "sigma0": Params().sigma0 / (1 + Params().alpha * 50.0)}),
                ("natural convection h frozen at h(30 K)", {"nc_delta": 30.0, "nc_coef": Params().nc_coef * 0.0 + Params().nc_coef}),
                ]
    for label, kw in variants[:2]:
        p = Params().replace(**{k: v for k, v in kw.items() if k != "nonlinear"})
        pr = make_problem(kw.get("nonlinear", True), p)
        e = pr.evaluate(np.asarray(d) / pr.d0, need_grad=False)
        rows.append({"model": label, "peak_bar_C": float(e["Tmax"][0]), "peak_joint_C": float(e["Tmax"][1]),
                     "resistance_uohm": float(e["Tmax"][2])})
    # sigma(T) alone frozen: alpha = 0 with sigma0 set to sigma(70 C)
    p0 = Params()
    p = p0.replace(alpha=0.0, sigma0=p0.sigma0 / (1 + p0.alpha * (p0.T_lin - p0.T_ref)), T_ref=p0.T_ref)
    pr = make_problem(True, p)
    e = pr.evaluate(np.asarray(d) / pr.d0, need_grad=False)
    rows.append({"model": "only sigma frozen at sigma(70 C)", "peak_bar_C": float(e["Tmax"][0]),
                 "peak_joint_C": float(e["Tmax"][1]), "resistance_uohm": float(e["Tmax"][2])})
    return rows


SCENARIOS = [("nominal", {}), ("+10 % peak current (1320 A)", {"I_peak": 1.1 * Params().I_peak}),
             ("+5 K pack air and terminals", {"T_amb": Params().T_amb + 5, "T_term": Params().T_term + 5})]


class RobustProblem:
    """Worst case over sampled scenarios: every constraint in every scenario (bonus)."""

    def __init__(self):
        self.subs = []
        for label, kw in SCENARIOS:
            p = Params().replace(**kw)
            cons = [dict(c, name=f"{c['name']} | {label}") for c in constraint_spec(p)]
            self.subs.append(BusbarProblem(BusbarSolver(True, p), cons))
        b = self.subs[0]
        self.d0, self.m0, self.solver = b.d0, b.m0, b.solver
        self.cons = [c for s in self.subs for c in s.cons]

    n_forward = property(lambda self: sum(s.n_forward for s in self.subs))
    n_adjoint = property(lambda self: sum(s.n_adjoint for s in self.subs))
    d_of = lambda self, x: self.subs[0].d_of(x)

    def evaluate(self, x, need_grad=True):
        es = [s.evaluate(x, need_grad) for s in self.subs]
        out = dict(es[0])
        out["g"] = np.concatenate([e["g"] for e in es])
        out["Tmax"] = np.concatenate([e["Tmax"] for e in es])
        if need_grad:
            out["grad_g"] = np.vstack([e["grad_g"] for e in es])
        return out

    f = BusbarProblem.f
    df = BusbarProblem.df
    c_ineq = BusbarProblem.c_ineq
    dc_ineq = BusbarProblem.dc_ineq
    bounds = BusbarProblem.bounds


MEANING = {
    "peak_bar_T": "kg saved per K the insulation limit is relaxed",
    "peak_joint_T": "kg saved per K the joint limit is relaxed: the joint heat (contact resistance) and the heat "
                    "conducted back from the span must go into the cell terminal",
    "resistance": "kg saved per micro-ohm of extra resistance allowed: copper cross-section bought for efficiency",
}
UNIT = {"peak_bar_T": "kg/K", "peak_joint_T": "kg/K", "resistance": "kg/µΩ"}


def optimization(quick=False):
    out = {"algorithm": "SciPy SLSQP (analytic adjoint gradients)",
           "scaling": "x = d/d0; f = m/m0; g_k = KS_k(T/T_lim) - 1 <= 0; g_R = R/R_max - 1 <= 0", "starts": []}
    prob = make_problem(True)
    lb, ub = np.array(prob.bounds()).T
    starts = {"baseline": np.ones(24), "all_max": ub.copy(), "midpoint": 0.5 * (lb + ub)}
    if quick:
        starts = {"baseline": starts["baseline"]}
    best = None
    for name, xs in starts.items():
        print(f"  SLSQP start {name}", flush=True)
        r = run_slsqp(prob, xs, label=name, maxiter=30 if quick else 150)
        r["summary"] = summary(prob, np.array(r["x_opt"]))
        out["starts"].append(r)
        if r["kkt"]["max_violation"] < 1e-4 and (best is None or r["history"][-1]["f"] < best["history"][-1]["f"]):
            best = r
    best = best or min(out["starts"], key=lambda r: r["kkt"]["max_violation"])
    out["best_start"] = best["label"]
    out["spread_mass_kg"] = [s["summary"]["mass_kg"] for s in out["starts"]]
    out["initial"] = summary(prob, np.ones(24))
    x_nl = np.array(best["x_opt"])
    print("  N2 comparison", flush=True)
    lin = make_problem(False)
    rl = run_slsqp(lin, np.ones(24), label="constant_property_model", maxiter=30 if quick else 150)
    x_lin = np.array(rl["x_opt"])
    rl["summary_linear_model"] = summary(lin, x_lin)
    rl["summary_full_physics"] = summary(prob, x_lin)
    out["linear"] = rl
    out["n2_comparison"] = {"mass_nonlinear_opt_kg": best["summary"]["mass_kg"],
                            "mass_linear_opt_kg": rl["summary_linear_model"]["mass_kg"],
                            "linear_opt_checked_with_full_physics": rl["summary_full_physics"]["constraints"],
                            "d_nonlinear_mm": (1e3 * prob.d_of(x_nl)).tolist(),
                            "d_linear_mm": (1e3 * lin.d_of(x_lin)).tolist(),
                            "attribution_linear_opt": n2_attribution(lin.d_of(x_lin)),
                            "attribution_nonlinear_opt": n2_attribution(prob.d_of(x_nl))}
    if not quick:
        print("  robust design", flush=True)
        rp = RobustProblem()
        chk = summary(rp, x_nl)
        rr = run_slsqp(rp, x_nl, label="robust", maxiter=150)
        rr["summary"] = summary(rp, np.array(rr["x_opt"]))
        rr["nominal_optimum_under_scenarios"] = chk
        rr["scenarios"] = [{"label": l, **{k: float(v) for k, v in kw.items()}} for l, kw in SCENARIOS]
        rr["mass_penalty_kg"] = rr["summary"]["mass_kg"] - best["summary"]["mass_kg"]
        out["robust"] = rr
    return out, prob.d_of(x_nl), lin.d_of(x_lin)


def website_views(sens, opt):
    gc = []
    for key in ("nonlinear", "linear"):
        for r in sens[key]["fd_table"]:
            gc.append({"function": r["fun"], "variable": r["var_label"], "model": sens[key]["label"],
                       "adjoint": r["adjoint"], "fd": r["fd"], "rel_error": r["rel_err"],
                       "err_vs_grad_norm": r["err_vs_grad_norm"]})
    sens["gradient_check"] = gc
    tay = sens["nonlinear"]["taylor"]
    err = np.array(tay["central_err"])
    sens["taylor_plot"] = {"eps": tay["eps"], **{nm: err[:, k].tolist() for k, nm in
                                                 enumerate(sens["nonlinear"]["functionals"])}}
    c = sens["cost"]
    sens["cost_view"] = {"n_vars": c["n_var"], "forward_s": c["t_forward_s"], "adjoint_s": c["t_adjoint_only_s"],
                         "fd_s": c["t_fd_central_est_s"], "speedup": c["speedup_vs_central"]}
    for s in opt["starts"] + [opt["linear"]] + ([opt["robust"]] if "robust" in opt else []):
        s["name"], s["iterations"], s["converged"] = s["label"], s["nit"], s["success"]
        s["final_mass_kg"] = s["history"][-1]["mass_kg"]
        for h in s["history"]:
            h["peak_bar_C"], h["peak_joint_C"], h["resistance_uohm"] = h["Tmax_C"][:3]
            for k in ("g_face_C", "g_hub_C"):
                h.pop(k, None)
    best = next(s for s in opt["starts"] if s["label"] == opt["best_start"])
    k = best["kkt"]
    opt["active"] = [{"constraint": nm, "multiplier": k["kg_per_K_relaxation"][nm], "unit": UNIT.get(nm, ""),
                      "multiplier_scaled": mu, "meaning": MEANING.get(nm, "")}
                     for nm, mu in k["constraint_multipliers"].items()]
    opt["active_bounds"] = {"lower": [var_label(i) for i in k["lower_bounds_active"]],
                            "upper": [var_label(i) for i in k["upper_bounds_active"]]}


# ------------------------------------------------------------------ figures
def figures(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figs = {}
    conv = out["convergence"]
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    for key, lab in (("mms_space", "MMS coupled 3D, T (V7a)"), ("mms_space_psi", "MMS coupled 3D, potential (V7b)"),
                     ("coupled_1d", "sigma(T) bar (V4)")):
        if key in conv:
            c = conv[key]
            ax[0].loglog(c["h_rel"], c["err_L2"], "o-", label=f"{lab}: slope {c['slope']:.2f}")
    ax[0].set_xlabel("relative element size h/h0"); ax[0].set_ylabel("RMS nodal error")
    ax[0].legend(fontsize=7); ax[0].set_title("Spatial convergence")
    for key, lab in (("natconv_time", "natural convection lump (V5)"), ("radiation_time", "radiation lump (V6)")):
        c = conv[key]
        ax[1].loglog(c["dt"], c["err"], "o-", label=f"{lab}: slope {c['slope']:.2f}")
    ax[1].set_xlabel("time step dt (s)"); ax[1].set_ylabel("temperature error (K)")
    ax[1].legend(fontsize=7); ax[1].set_title("Temporal convergence (implicit Euler)")
    fig.tight_layout(); fig.savefig(os.path.join(ASSETS, "fig_convergence.png"), dpi=150); plt.close(fig)
    figs["convergence"] = "assets/fig_convergence.png"
    r = conv["newton"]["residuals"]
    fig, ax = plt.subplots(figsize=(4.5, 3.4))
    ax.semilogy(range(len(r)), r, "o-"); ax.set_xlabel("Newton iteration"); ax.set_ylabel("residual norm")
    ax.set_title("Newton, first pulse step (V9)"); fig.tight_layout()
    fig.savefig(os.path.join(ASSETS, "fig_newton.png"), dpi=150); plt.close(fig)
    figs["newton"] = "assets/fig_newton.png"
    fig, ax = plt.subplots(figsize=(7, 3.6))
    for key, blk in out["designs"].items():
        h = blk["history"]
        ax.plot(h["times"], h["peak_bar_C"], label=f"{key}: busbar peak")
        ax.plot(h["times"], h["peak_joint_C"], "--", label=f"{key}: joint")
    p = Params()
    ax.axhline(p.T_bar_max - C0, color="k", lw=0.8, ls=":"); ax.axhline(p.T_joint_max - C0, color="k", lw=0.8, ls=":")
    ax.set_xlabel("time after the steady continuous rating (s)"); ax.set_ylabel("temperature (C)")
    ax.legend(fontsize=7); fig.tight_layout()
    fig.savefig(os.path.join(ASSETS, "fig_history.png"), dpi=150); plt.close(fig)
    figs["history"] = "assets/fig_history.png"
    if "sensitivities" in out:
        tay = out["sensitivities"]["nonlinear"]["taylor"]
        names = out["sensitivities"]["nonlinear"]["functionals"]
        eps, err = np.array(tay["eps"]), np.array(tay["central_err"])
        fig, ax = plt.subplots(figsize=(5.5, 4))
        for k, nm in enumerate(names):
            ax.loglog(eps, np.maximum(err[:, k], 1e-18), "o-", label=nm)
        ax.loglog(eps, err[0, 1] * (eps / eps[0]) ** 2, "k--", lw=1, label="slope 2")
        ax.set_xlabel("step eps (scaled design units)"); ax.set_ylabel("|central FD - adjoint grad . v|")
        ax.set_title("Taylor test (random direction)"); ax.legend(fontsize=7); ax.grid(True, which="both", alpha=.3)
        fig.tight_layout(); fig.savefig(os.path.join(ASSETS, "taylor_test.png"), dpi=150); plt.close(fig)
        figs["taylor"] = "assets/taylor_test.png"
        fig, ax = plt.subplots(figsize=(5.5, 4))
        for s in out["optimization"]["starts"]:
            h = s["history"]
            ax.plot([q["iter"] for q in h], [q["mass_kg"] for q in h], "o-", ms=3, label=s["label"])
        ax.set_xlabel("SLSQP iteration"); ax.set_ylabel("busbar mass (kg)"); ax.legend(); ax.grid(alpha=.3)
        fig.tight_layout(); fig.savefig(os.path.join(ASSETS, "opt_history.png"), dpi=150); plt.close(fig)
        figs["opt_history"] = "assets/opt_history.png"
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
        xs = 1e3 * station_x(p)
        for key, blk in out["designs"].items():
            ax[0].plot(xs, blk["d"]["width_mm"], "o-", ms=3, label=key)
            ax[1].plot(xs, blk["d"]["thickness_mm"], "o-", ms=3, label=key)
        ax[0].set_ylabel("width (mm)"); ax[1].set_ylabel("thickness (mm)")
        for a in ax:
            a.set_xlabel("x along the bar (mm)"); a.axvspan(0, 1e3 * p.l_c, color="0.9"); a.axvspan(1e3 * (p.L - p.l_c), 1e3 * p.L, color="0.9")
        ax[1].legend(fontsize=7); fig.tight_layout()
        fig.savefig(os.path.join(ASSETS, "fig_designs.png"), dpi=150); plt.close(fig)
        figs["designs"] = "assets/fig_designs.png"
    return figs


def _git():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE, text=True).strip()
    except Exception:
        return "unknown"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--skip-opt", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(ASSETS, exist_ok=True)
    t0 = time.perf_counter()
    p = Params()
    out = {"meta": {"generated": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
                    "git_commit": _git(), "command": "python -m busbar.run_all" + (" --quick" if a.quick else ""),
                    "topic": "busbar"},
           "params": p.to_dict()}
    print("[1/4] verification ladder", flush=True)
    rows, conv = run_ladder(include_slow=not a.quick)
    out["verification"], out["convergence"] = rows, conv
    for r in rows:
        print(f"  {r['id']:5s} {'PASS' if r['pass'] else 'FAIL'} rel={r['rel_error']:.2e} {r['name']}", flush=True)
    print("[2/4] initial design", flush=True)
    designs = {"initial": design_block(initial_design(p), p, "initial (uniform 30 x 5 mm bar)")}
    if not a.skip_opt:
        print("[3/4] sensitivities + optimization", flush=True)
        check = [0, 4, 8, 11, 12, 16, 20, 23] if not a.quick else [0, 12]
        sens = {"nonlinear": sensitivities(make_problem(True), np.ones(24), check, "full physics"),
                "linear": sensitivities(make_problem(False), np.ones(24), check[:4], "constant-property model")}
        sens["cost"] = cost_comparison(make_problem(True), np.ones(24))
        opt, d_nl, d_lin = optimization(a.quick)
        website_views(sens, opt)
        out["sensitivities"], out["optimization"] = sens, opt
        print("[4/4] optimized design blocks", flush=True)
        designs["optimized"] = design_block(d_nl, p, "optimized (full physics)")
        designs["optimized_linear"] = design_block(d_lin, p, "optimized with the constant-property model, "
                                                             "evaluated with full physics")
        if "robust" in opt:
            designs["optimized_robust"] = design_block(opt["robust"]["d_opt_m"], p,
                                                       "robust optimum (3 scenarios), nominal scenario shown")
    out["designs"] = designs
    out["figures"] = figures(out)
    out["meta"]["wall_time_s"] = round(time.perf_counter() - t0, 1)
    path = os.path.join(ASSETS, "results.json")
    with open(path, "w") as f:
        json.dump(out, f, separators=(",", ":"), default=float)
    print(f"wrote {path} ({os.path.getsize(path) / 1e6:.2f} MB) in {out['meta']['wall_time_s']} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
