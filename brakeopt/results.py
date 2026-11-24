"""Builders for docs/assets/results.json blocks (schema in ainote.md)."""
from __future__ import annotations

import numpy as np

from .functionals import mass, peak_temperatures
from .mesh import build_mesh
from .model import ThermalModel
from .params import Params, station_radii, split_design
from .solver import solve_transient

C0 = 273.15


def _r(x, nd):
    return np.round(np.asarray(x, dtype=float), nd).tolist()


def design_block(d, p: Params | None = None, label: str = "", max_frames: int = 80, sol=None, model=None):
    p = Params() if p is None else p
    d = np.asarray(d, dtype=float)
    if model is None:
        model = ThermalModel(build_mesh(d, p), p)
    m = model.mesh
    if sol is None:
        sol = solve_transient(model)
    s = model.sched
    tv, tc = split_design(d)
    tf, th = peak_temperatures(m, sol.T)
    fr, hb = m.face_nodes("friction"), m.face_nodes("hub")
    idx = np.unique(np.round(np.linspace(0, len(sol.t) - 1, min(max_frames, len(sol.t)))).astype(int))
    cons = [
        {"name": "peak_face_T", "label": "Peak friction-surface temperature", "value_C": round(tf - C0, 2),
         "limit_C": round(p.T_face_max - C0, 2), "pass": bool(tf <= p.T_face_max + 1e-9)},
        {"name": "peak_hub_T", "label": "Peak hub / bell-mounting temperature", "value_C": round(th - C0, 2),
         "limit_C": round(p.T_hub_max - C0, 2), "pass": bool(th <= p.T_hub_max + 1e-9)},
    ]
    lo_ok = bool(np.all(d >= np.r_[[p.tv_min] * 12, [p.tc_min] * 12] - 1e-12))
    cons.append({"name": "min_thickness", "label": "Minimum vane / cheek thickness (bounds)",
                 "value_mm": round(1e3 * min(tv.min(), tc.min()), 3), "limit_mm": round(1e3 * min(p.tv_min, p.tc_min), 3),
                 "pass": lo_ok})
    return {
        "label": label,
        "d": {"vane_thickness_mm": _r(1e3 * tv, 4), "cheek_thickness_mm": _r(1e3 * tc, 4),
              "station_radius_mm": _r(1e3 * station_radii(p), 3)},
        "mass_kg": round(float(mass(m, p, model)), 4),
        "constraints": cons,
        "mesh": {"nodes": _r(m.X, 6), "hexes": m.hexes.tolist(), "n_vanes": p.n_vanes,
                 "faces": {k: v.tolist() for k, v in m.faces.items()},
                 "symmetry": "quarter sector: mirror z=0 and theta=0 (y=0), then rotate by 2*pi/n_vanes, n_vanes times"},
        "temperature": {"unit": "C", "times": _r(sol.t[idx], 3), "frames": [_r(sol.T[i] - C0, 1) for i in idx]},
        "history": {"times": _r(sol.t, 3), "peak_face_C": _r(sol.T[:, fr].max(1) - C0, 2),
                    "peak_hub_C": _r(sol.T[:, hb].max(1) - C0, 2), "mean_C": _r(sol.T.mean(1) - C0, 2),
                    "speed_kmh": _r(3.6 * s.v, 2), "braking": s.braking.tolist()},
        "model": {"radiation": p.radiation, "kT": p.kT, "n_nodes": m.n_nodes, "n_elem": m.n_elem,
                  "n_steps": s.n_steps},
    }
