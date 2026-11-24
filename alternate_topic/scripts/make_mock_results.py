"""Write docs/assets/mock_results.json: an invented results file that follows the
results.json schema in ainote.md, so the website can be built and tested before
the real pipeline exists.

Every number in it is made up. The file sets meta.mock = true and the site shows
a red "mock data" banner whenever it is used. The site only falls back to this
file when docs/assets/results.json (written by scripts/run_all.py) is missing,
or when the page is opened with ?mock.

    python scripts/make_mock_results.py
"""
from __future__ import annotations

import datetime as dt
import json
import pathlib

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "assets" / "mock_results.json"

# Geometry roughly as in brakeopt/params.py (copied, not imported, so this script
# keeps working whatever state the package is in).
R_I, R_O, R_PAD, H_V, N_VANES = 0.110, 0.178, 0.128, 0.012, 48
RHO = 7100.0
T_FACE_MAX, T_HUB_MAX = 700.0, 300.0
STATIONS = np.linspace(R_I, R_O, 12)
HEX_FACES = np.array([[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [2, 3, 7, 6], [0, 4, 7, 3], [1, 2, 6, 5]])


def quarter_sector(tv, tc, nr=16, n_side=3, n_vane=1, nz_vane=2, nz_cheek=2):
    """Quarter sector: theta in [0, pi/N], z in [0, h_v/2 + t_c(r)], hex8 (VTK order)."""
    half = np.pi / N_VANES
    r = np.linspace(R_I, R_O, nr + 1)
    tvr, tcr = np.interp(r, STATIONS, tv), np.interp(r, STATIONS, tc)
    idx, X = {}, []

    def nid(i, j, k):
        if (i, j, k) not in idx:
            a = min(tvr[i] / (2 * r[i]), 0.8 * half)
            th = np.concatenate([np.linspace(0, a, n_vane + 1), np.linspace(a, half, n_side + 1)[1:]])[j]
            z = np.concatenate([np.linspace(0, H_V / 2, nz_vane + 1), H_V / 2 + np.linspace(0, tcr[i], nz_cheek + 1)[1:]])[k]
            X.append([r[i] * np.cos(th), r[i] * np.sin(th), z])
            idx[(i, j, k)] = len(X) - 1
        return idx[(i, j, k)]

    hexes = []
    for i in range(nr):
        for j in range(n_vane + n_side):
            for k in range(nz_vane + nz_cheek):
                if k < nz_vane and j >= n_vane:
                    continue  # air channel
                hexes.append([nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                              nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)])
    return np.array(X), np.array(hexes)


def face_sets(X, hexes, tc):
    faces = hexes[:, HEX_FACES].reshape(-1, 4)
    key = np.sort(faces, axis=1)
    _, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
    bnd = faces[cnt[inv.ravel()] == 1]
    c = X[bnd].mean(axis=1)
    r = np.hypot(c[:, 0], c[:, 1])
    th = np.arctan2(c[:, 1], c[:, 0])
    ztop = H_V / 2 + np.interp(r, STATIONS, tc)
    zr = X[bnd][:, :, 2]
    rr = np.hypot(X[bnd][:, :, 0], X[bnd][:, :, 1])
    sets = {"friction": [], "face_free": [], "channel": [], "rim": [], "hub": []}
    for q, cq, rq, thq, zt, zq, rqs in zip(bnd, c, r, th, ztop, zr, rr):
        if np.all(zq < 1e-9) or np.all(np.abs(X[q][:, 1]) < 1e-9) or abs(thq - np.pi / N_VANES) < 1e-3 and np.ptp(np.arctan2(X[q][:, 1], X[q][:, 0])) < 1e-6:
            continue  # symmetry planes
        if np.all(np.abs(rqs - R_O) < 1e-9):
            sets["rim"].append(q)
        elif np.all(np.abs(rqs - R_I) < 1e-9):
            sets["hub"].append(q)
        elif abs(cq[2] - zt) < 1e-6:
            sets["friction" if rq >= R_PAD else "face_free"].append(q)
        else:
            sets["channel"].append(q)
    return {k: np.array(v).tolist() for k, v in sets.items()}


def schedule():
    t, v = [0.0], [250.0]
    t0 = 0.0
    for _ in range(3):
        for s in np.linspace(0.3, 3.2, 10):
            t.append(t0 + s); v.append(250 - 170 * s / 3.2)
        for s in np.linspace(1.2, 12, 9):
            t.append(t0 + 3.2 + s); v.append(80 + 170 * s / 12)
        t0 += 15.2
    return np.array(t), np.array(v)


def temperatures(X, t, v, tv, tc, friction_nodes):
    r = np.hypot(X[:, 0], X[:, 1])
    z = X[:, 2]
    tcr, tvr = np.interp(r, STATIONS, tc), np.interp(r, STATIONS, tv)
    depth = H_V / 2 + tcr - z
    thin = (0.009 / tcr) ** 0.9 * (0.005 / tvr) ** 0.15
    band = np.clip((r - R_I) / (R_PAD - R_I), 0, 1) ** 1.5 * (0.8 + 0.35 * np.exp(-((r - 0.155) / 0.02) ** 2))
    braking = np.r_[False, np.diff(v) < 0]
    T = np.zeros((len(t), len(X)))
    bulk, flash, prev = 0.0, 0.0, 0.0
    for n, tn in enumerate(t):
        d = tn - prev; prev = tn
        if braking[n]:
            bulk += 22.0 * d; flash = 210 * min(1, (tn % 15.2) / 3.2) ** 0.5
        else:
            bulk *= np.exp(-d / 40.0); flash *= np.exp(-d / 1.5)
        T[n] = 200 + thin * band * (bulk + flash * np.exp(-depth / 0.003)) + 1.6 * thin * bulk * (0.35 + 0.65 * (r - R_I) / (R_O - R_I))
    return 200 + 0.575 * (T - 200)   # calibrated so the mock optimum sits on the 700 C limit


def design(tv_mm, tc_mm, t, v):
    tv, tc = np.asarray(tv_mm) * 1e-3, np.asarray(tc_mm) * 1e-3
    X, hexes = quarter_sector(tv, tc)
    fs = face_sets(X, hexes, tc)
    fr_nodes = np.unique(np.array(fs["friction"]))
    hub_nodes = np.unique(np.array(fs["hub"]))
    T = temperatures(X, t, v, tv, tc, fr_nodes)
    # volume by tetra split is overkill for a mock: use the analytic sector volume
    rr = np.linspace(R_I, R_O, 400)
    area = 2 * np.pi * rr * (H_V / 2 * 0 + 2 * np.interp(rr, STATIONS, tc)) + N_VANES * np.interp(rr, STATIONS, tv) * H_V
    mass = RHO * np.trapezoid(area, rr) if hasattr(np, "trapezoid") else RHO * np.trapz(area, rr)
    pf, ph = T[:, fr_nodes].max(axis=1), T[:, hub_nodes].max(axis=1)
    cons = [
        dict(name="peak_face_T", label="Peak friction-face temperature", value_C=float(pf.max() + 4.0), true_max_C=float(pf.max()),
             limit_C=T_FACE_MAX),
        dict(name="peak_hub_T", label="Peak hub / bell-mount temperature", value_C=float(ph.max() + 2.0), true_max_C=float(ph.max()),
             limit_C=T_HUB_MAX),
    ]
    for c in cons:
        c["pass"] = bool(c["value_C"] <= c["limit_C"])
    return {
        "d": {"vane_thickness_mm": list(map(float, tv_mm)), "cheek_thickness_mm": list(map(float, tc_mm)),
              "station_radius_mm": np.round(STATIONS * 1e3, 3).tolist()},
        "mass_kg": float(mass),
        "constraints": cons,
        "mesh": {"nodes": np.round(X, 6).tolist(), "hexes": hexes.tolist(), "n_vanes": N_VANES,
                 "symmetry": "quarter sector: mirror z=0 and theta=0, then rotate n_vanes times",
                 "faces": fs},
        "temperature": {"unit": "C", "times": np.round(t, 3).tolist(), "frames": np.round(T, 1).tolist()},
        "history": {"times": np.round(t, 3).tolist(), "peak_face_C": np.round(pf, 2).tolist(),
                    "peak_hub_C": np.round(ph, 2).tolist(), "speed_kmh": np.round(v, 2).tolist()},
    }


def main():
    rng = np.random.default_rng(7)
    t, v = schedule()
    s = np.linspace(0, 1, 12)
    init = design(np.full(12, 5.0), np.full(12, 9.0), t, v)
    opt = design(np.round(2 + 2.6 * np.exp(-((s - 0.6) / 0.3) ** 2), 3), np.round(6.3 + 2.4 * np.exp(-((s - 0.62) / 0.3) ** 2), 3), t, v)
    lin = design(np.round(2 + 2.2 * np.exp(-((s - 0.6) / 0.3) ** 2), 3), np.round(5.3 + 2.0 * np.exp(-((s - 0.62) / 0.3) ** 2), 3), t, v)
    m0, m1 = init["mass_kg"], opt["mass_kg"]

    def hist(ma, mb, n, g0):
        k = np.arange(n)
        return {"iter": k.tolist(), "mass_kg": (mb + (ma - mb) * np.exp(-k / 3)).tolist(),
                "g_face_C": (700 + g0 * np.exp(-k / 2.5)).tolist(), "g_hub_C": (230 + 20 * np.exp(-k / 3)).tolist(),
                "kkt": (10.0 ** (-0.5 * k)).tolist()}

    eps = np.logspace(-1, -8, 15)
    res = {
        "meta": {"mock": True, "generated": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                 "git_commit": "mock", "command": "python scripts/run_all.py",
                 "note": "MOCK DATA from scripts/make_mock_results.py; every number is invented."},
        "params": {"r_i": {"value": R_I, "unit": "m"}, "r_o": {"value": R_O, "unit": "m"}, "n_vanes": {"value": N_VANES, "unit": "-"},
                   "tv_min": {"value": 0.002, "unit": "m"}, "tc_min": {"value": 0.005, "unit": "m"},
                   "tv_max": {"value": 0.008, "unit": "m"}, "tc_max": {"value": 0.012, "unit": "m"},
                   "T_face_max": {"value": 973.15, "unit": "K"}, "T_hub_max": {"value": 573.15, "unit": "K"}},
        "designs": {"initial": init, "optimized": opt, "optimized_linear": lin},
        "verification": [
            dict(id="V1", name="Steady 1D slab, fixed T both ends", ours=412.0, exact=412.0, rel_error=3e-15, pass_=True, note="linear profile"),
            dict(id="V2", name="Steady 1D slab, convective end", ours=388.512, exact=388.515, rel_error=8e-6, pass_=True, note=""),
            dict(id="V3", name="MMS steady 3D, distorted mesh, k(T)", ours=2.01, exact=2.0, rel_error=5e-3, pass_=True, note="L2 slope"),
            dict(id="V4", name="MMS transient (implicit Euler)", ours=0.99, exact=1.0, rel_error=1e-2, pass_=True, note="slope in dt"),
            dict(id="V5", name="Energy balance over braking cycle", ours=1.0003, exact=1.0, rel_error=3e-4, pass_=True, note="in / (stored+out)"),
            dict(id="V6", name="Nonlinear 1D slab k(T) (Kirchhoff)", ours=611.2, exact=611.25, rel_error=8e-5, pass_=True, note="Newton quadratic"),
            dict(id="V7", name="Rigid rotation invariance of peak T", ours=943.1, exact=943.1, rel_error=1e-13, pass_=True, note=""),
        ],
        "convergence": {"mms_space": {"h": [0.02, 0.01, 0.005, 0.0025], "err_L2": [3.1e-2, 7.8e-3, 1.95e-3, 4.9e-4], "slope": 2.0},
                        "mms_time": {"dt": [0.4, 0.2, 0.1, 0.05], "err": [2.2e-2, 1.1e-2, 5.5e-3, 2.75e-3], "slope": 1.0},
                        "newton": {"residuals": [1.0, 4e-2, 2e-5, 5e-11, 3e-16]}},
        "sensitivities": {
            "gradient_check": [dict(function=f, variable=var, adjoint=a, fd=a * (1 + e), rel_error=abs(e))
                               for f in ("mass", "g_face", "g_hub")
                               for var, a, e in (("t_v,4 (vane)", float(rng.normal()), float(rng.uniform(1e-8, 1e-6))),
                                                 ("t_c,9 (cheek)", float(rng.normal()), float(rng.uniform(1e-8, 1e-6))))],
            "taylor": {"eps": eps.tolist(), "err_fd": (3 * eps).tolist(), "err_adjoint": (5 * eps ** 2 + 1e-12).tolist()},
            "cost": {"n_vars": 24, "forward_s": 4.1, "adjoint_s": 5.3, "fd_s": 25 * 4.1},
        },
        "optimization": {
            "algorithm": "SLSQP (scipy.optimize.minimize)",
            "scaling": "x = d / d_ref; objective m / m_0; constraints g / g_lim - 1",
            "starts": [dict(name="A: baseline (5 mm vanes, 9 mm cheeks)", history=hist(m0, m1, 15, -150), final_mass_kg=m1, converged=True),
                       dict(name="B: thick (8 mm, 12 mm)", history=hist(m0 * 1.3, m1 * 1.003, 19, -260), final_mass_kg=m1 * 1.003, converged=True)],
            "active": [dict(constraint="peak_face_T", multiplier=0.0042, unit="kg/°C",
                            meaning="1 °C more allowed on the rubbing face would save about 4 g"),
                       dict(constraint="t_v ≥ 2 mm (outer 3 stations)", multiplier=0.11, unit="kg/mm",
                            meaning="the casting limit, not heat, sets the outer vanes")],
        },
    }
    # json does not allow a key named "pass" via dict(pass=...), so rename here
    for row in res["verification"]:
        row["pass"] = row.pop("pass_")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(res, separators=(",", ":")))
    print(f"wrote {OUT.relative_to(ROOT)} ({OUT.stat().st_size / 1e6:.2f} MB, MOCK)",
          f"m0={m0:.2f} m*={m1:.2f}", [round(c["value_C"]) for c in opt["constraints"]], [round(c["value_C"]) for c in init["constraints"]])


if __name__ == "__main__":
    main()
