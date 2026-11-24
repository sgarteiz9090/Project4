"""One command that regenerates docs/assets/results.json, the figures, and every number on the site.

    python scripts/run_all.py            # full run
    python scripts/run_all.py --quick    # fast smoke run (skips the slow mesh study, one SLSQP start)

Stages: verification ladder + convergence figures (WS1), initial design (WS1),
sensitivities + optimization (WS2+3, optim.run_ws23), optimized designs (evaluated by WS1 code).
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from brakeopt.params import Params, initial_design  # noqa: E402
from brakeopt.results import design_block  # noqa: E402
from brakeopt.loads import energy_into_one_face, face_power_factor  # noqa: E402
from verification.ladder import run_ladder  # noqa: E402

ASSETS = os.path.join(ROOT, "docs", "assets")


def _git_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        return "unknown"


def convergence_figures(conv):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figs = {}
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    for key, lab in (("mms_space", "MMS, sector mesh (V5)"), ("kirchhoff", "k(T) slab (V3)")):
        c = conv[key]
        h = np.array(c.get("h_rel", c.get("h")))
        h = h / h[0]
        ax[0].loglog(h, c["err_L2"], "o-", label=f"{lab}: slope {c['slope']:.2f}")
    ax[0].loglog([1, 0.25], [1, 1 / 16], "k--", lw=0.8, label="slope 2")
    ax[0].set_xlabel("relative element size h/h0")
    ax[0].set_ylabel("L2 temperature error (K)")
    ax[0].legend(fontsize=7)
    ax[0].set_title("Spatial convergence")
    for key, lab in (("mms_time", "MMS transient (V6)"), ("radiation_time", "radiation plate (V4)")):
        c = conv[key]
        ax[1].loglog(c["dt"], c["err"], "o-", label=f"{lab}: slope {c['slope']:.2f}")
    ax[1].set_xlabel("time step dt (s)")
    ax[1].set_ylabel("temperature error (K)")
    ax[1].legend(fontsize=7)
    ax[1].set_title("Temporal convergence (implicit Euler)")
    fig.tight_layout()
    figs["convergence"] = "assets/fig_convergence.png"
    fig.savefig(os.path.join(ASSETS, "fig_convergence.png"), dpi=150)
    plt.close(fig)
    r = conv["newton"]["residuals"]
    fig, ax = plt.subplots(figsize=(4.5, 3.4))
    ax.semilogy(range(len(r)), r, "o-")
    ax.set_xlabel("Newton iteration")
    ax.set_ylabel("residual norm (W)")
    ax.set_title(f"Newton, braking step {conv['newton']['step']} (V8)")
    fig.tight_layout()
    fig.savefig(os.path.join(ASSETS, "fig_newton.png"), dpi=150)
    plt.close(fig)
    figs["newton"] = "assets/fig_newton.png"
    return figs


def history_figure(designs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 3.6))
    for key, blk in designs.items():
        h = blk["history"]
        ax.plot(h["times"], h["peak_face_C"], label=f"{key}: friction face")
        ax.plot(h["times"], h["peak_hub_C"], "--", label=f"{key}: hub")
    any_blk = next(iter(designs.values()))
    lim = {c["name"]: c.get("limit_C") for c in any_blk["constraints"]}
    ax.axhline(lim["peak_face_T"], color="k", lw=0.8, ls=":")
    ax.axhline(lim["peak_hub_T"], color="k", lw=0.8, ls=":")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("peak temperature (C)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(ASSETS, "fig_history.png"), dpi=150)
    plt.close(fig)
    return "assets/fig_history.png"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--skip-opt", action="store_true", help="skip the WS2+3 stage")
    args = ap.parse_args(argv)
    os.makedirs(ASSETS, exist_ok=True)
    t0 = time.perf_counter()
    p = Params()
    out = {"meta": {"generated": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
                    "git_commit": _git_commit(),
                    "command": "python scripts/run_all.py" + (" --quick" if args.quick else "")},
           "params": p.to_dict()}
    out["params"]["derived"] = {"face_power_fraction": face_power_factor(p),
                                "energy_per_face_J": energy_into_one_face(p)}
    print("[1/4] verification ladder ...", flush=True)
    rows, conv = run_ladder(include_slow=not args.quick)
    out["verification"] = rows
    out["convergence"] = conv
    out["figures"] = convergence_figures(conv)
    for r in rows:
        print(f"  {r['id']:5s} {'PASS' if r['pass'] else 'FAIL'} rel={r['rel_error']:.2e} {r['name']}")
    print("[2/4] initial design ...", flush=True)
    designs = {"initial": design_block(initial_design(p), p, "initial (baseline)")}
    if not args.skip_opt:
        print("[3/4] sensitivities + optimization (WS2+3) ...", flush=True)
        try:
            from optim.run_ws23 import run as run_ws23
        except ImportError as e:  # WS2+3 not merged yet
            print("  WS2+3 stage unavailable:", e)
        else:
            blocks = run_ws23(quick=args.quick, fig_dir=ASSETS)
            out["sensitivities"] = blocks["sensitivities"]
            out["optimization"] = blocks["optimization"]
            print("[4/4] evaluating optimized designs ...", flush=True)
            designs["optimized"] = design_block(blocks["designs_optimized_d"], p,
                                                "optimized (radiation + k(T))")
            designs["optimized_linear"] = design_block(
                blocks["designs_optimized_linear_d"], p,
                "optimized with the linearized model, evaluated with full physics")
    out["designs"] = designs
    out["figures"]["history"] = history_figure(designs)
    out["meta"]["wall_time_s"] = round(time.perf_counter() - t0, 1)
    path = os.path.join(ASSETS, "results.json")
    with open(path, "w") as f:
        json.dump(out, f, separators=(",", ":"), default=float)
    print(f"wrote {path} ({os.path.getsize(path) / 1e6:.2f} MB) in {out['meta']['wall_time_s']} s")
    n_fail = sum(not r["pass"] for r in rows)
    if n_fail:
        print(f"WARNING: {n_fail} verification test(s) failed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
