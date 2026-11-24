"""Verification ladder for the busbar solver.  python -m busbar.ladder

Every row compares the solver with an exact answer (or an exact discrete identity / an exact
rate of convergence). The joints are moved to the end faces (p.contacts_at_ends) for the 1D
problems, so that they have closed-form solutions.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from .fem import HexGeometry, QuadGeometry
from .mesh import build_mesh
from .model import BusbarModel, Schedule, EXPOSED
from .params import Params, initial_design, PHI_SCALE, SIGMA_SB
from .solver import solve_transient, newton_step, initial_state

C0 = 273.15


def _row(id_, name, ours, exact, rel, tol, note="", **kw):
    r = {"id": id_, "name": name, "ours": float(ours), "exact": float(exact), "rel_error": float(rel),
         "tol": tol, "pass": bool(rel <= tol), "note": note}
    r.update(kw)
    return r


def _steady_sched(I, n=1, dt=1e15):
    return Schedule(t=np.arange(n + 1) * dt, dt=np.r_[np.nan, np.full(n, dt)], I=np.full(n + 1, I),
                    pulse=np.zeros(n + 1, bool))


def _uniform(p, w=0.03, t=0.005):
    return np.r_[np.full(12, w), np.full(12, t)]


def _p1d(**kw):
    base = dict(contacts_at_ends=True, nc_coef=0.0, emissivity=0.0, nx_c=2, nx_m=10, ny=2, nz=1)
    base.update(kw)
    return Params().replace(**base)


# ---------------------------------------------------------------- V1 resistance
def v1_resistance():
    p = _p1d(nonlinear=False)
    w, t = 0.03, 0.005
    M = BusbarModel(build_mesh(_uniform(p, w, t), p), p, _steady_sched(500.0))
    U = solve_transient(M).U
    A = w * t
    sig = p.sigma0 / (1 + p.alpha * (p.T_lin - p.T_ref))
    R_ex = 2 * p.rho_c / A + p.L / (sig * A)
    R = M.resistance(U[1], 1)
    return _row("V1", "Resistance of a uniform bar with end joints: 2 rho_c/A + L/(sigma A)", R * 1e6, R_ex * 1e6,
                abs(R - R_ex) / R_ex, 1e-12, "micro-ohm; linear potential is in the element space, so exact")


# ---------------------------------------------------------------- V2 invariance
def v2_invariance():
    p = Params()
    d = initial_design(p) * (1 + 0.2 * np.sin(np.arange(24)))
    sch = _steady_sched(p.I_base)
    M = BusbarModel(build_mesh(d, p), p, sch)
    U1 = solve_transient(M).U[1]
    m2 = build_mesh(d, p)
    rng = np.random.default_rng(3)
    Q, _ = np.linalg.qr(rng.standard_normal((3, 3)))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    m2.X = m2.X @ Q.T + np.array([0.3, -1.2, 2.0])
    U2 = solve_transient(BusbarModel(m2, p, sch)).U[1]
    err = np.max(np.abs(U2 - U1)[M.nn:2 * M.nn]) / np.ptp(U1[M.nn:2 * M.nn])
    return _row("V2", "Rigid rotation + translation of the mesh (tapered bar, full physics)",
                U2[M.nn:2 * M.nn].max() - C0, U1[M.nn:2 * M.nn].max() - C0, err, 1e-9,
                "peak T in C; error relative to the temperature range")


# ---------------------------------------------------------------- V3 1D Joule + Robin
def v3_joule_1d():
    p = _p1d(nonlinear=False)
    w, t, I = 0.03, 0.005, 800.0
    M = BusbarModel(build_mesh(_uniform(p, w, t), p), p, _steady_sched(I))
    U = solve_transient(M).U[1]
    T = U[M.nn:2 * M.nn]
    sig = p.sigma0 / (1 + p.alpha * (p.T_lin - p.T_ref))
    k = p.k0 * (1 + p.k_beta * (p.T_lin - p.T_ref))
    j = I / (w * t)
    q = j * j / sig
    T_end = p.T_term + (q * p.L / 2 + p.joint_heat_share * p.rho_c * j * j) / p.h_term
    x = np.real(M.mesh.X[:, 0])
    Tex = T_end + q * x * (p.L - x) / (2 * k)
    err = np.max(np.abs(T - Tex)) / np.ptp(Tex)
    return _row("V3", "1D Joule heating with joint heat and Robin ends: parabolic T(x)", T.max() - C0, Tex.max() - C0,
                err, 1e-9, "peak T in C; nodally exact for linear elements in 1D (max nodal error / range)")


# ---------------------------------------------------------------- V4 coupled sigma(T) 1D (nonlinear benchmark)
def v4_coupled_1d(levels=(1, 2, 4)):
    errs, hs, peaks = [], [], []
    for f in levels:
        p = _p1d(nonlinear=True, k_beta=0.0, h_term=3000.0).refined(f)
        w, t, I = 0.03, 0.004, 500.0
        M = BusbarModel(build_mesh(_uniform(p, w, t), p), p, _steady_sched(I))
        U = solve_transient(M).U[1]
        T = U[M.nn:2 * M.nn]
        A = w * t
        j = I / A
        a = j * j / p.sigma0
        lam = np.sqrt(a * p.alpha / p.k0)
        L, h = p.L, p.h_term
        Cc = (h * (p.T_ref - 1 / p.alpha - p.T_term) - p.joint_heat_share * p.rho_c * j * j) / \
             (p.k0 * lam * np.sin(lam * L / 2) - h * np.cos(lam * L / 2))
        x = np.real(M.mesh.X[:, 0])
        Tex = p.T_ref - 1 / p.alpha + Cc * np.cos(lam * (x - L / 2))
        errs.append(np.sqrt(np.mean((T - Tex) ** 2)))
        hs.append(1.0 / f)
        peaks.append((T.max(), Tex.max()))
    slope = np.polyfit(np.log(hs), np.log(errs), 1)[0]
    Tlin = None
    return [_row("V4a", "Two-way coupled sigma(T) bar, closed-form cosine profile (Newton)",
                 peaks[-1][0] - C0, peaks[-1][1] - C0, abs(peaks[-1][0] - peaks[-1][1]) / (peaks[-1][1] - 313.15), 1e-3,
                 "peak T in C on the finest mesh; error relative to the temperature rise"),
            _row("V4b", "Same benchmark: convergence order of the RMS nodal error", slope, 2.0, abs(slope - 2) / 2, 0.1,
                 "expected order 2", h_rel=hs, err_L2=errs)]


# ---------------------------------------------------------------- lumped cooling (V5 natural convection, V6 radiation)
def _lump(p, dT0, t_end, dt):
    p = p.replace(contacts_at_ends=True, h_term=0.0, k0=1e7, k_beta=0.0, nx_c=1, nx_m=2, ny=1, nz=1)
    n = int(round(t_end / dt))
    sch = Schedule(t=np.arange(n + 1) * dt, dt=np.r_[np.nan, np.full(n, dt)], I=np.zeros(n + 1),
                   pulse=np.zeros(n + 1, bool))
    M = BusbarModel(build_mesh(_uniform(p), p), p, sch)
    u0 = initial_state(M)
    u0[M.nn:2 * M.nn] = p.T_amb + dT0
    U = solve_transient(M, u0).U
    A = sum(np.sum(M.fgeo[nm].wdA) for nm in EXPOSED if nm in M.fgeo)
    Cth = p.rho * p.c * M.volume()
    return U[-1, M.nn:2 * M.nn].mean(), A, Cth, M


def v5_natural_convection(dts=(8.0, 4.0, 2.0, 1.0)):
    p = Params().replace(emissivity=0.0, nc_delta=1e-6)
    dT0, t_end = 60.0, 400.0
    res = [_lump(p, dT0, t_end, dt) for dt in dts]
    _, A, Cth, M = res[0]
    a = M.C_nc * A / Cth
    ex = p.T_amb + (dT0 ** -0.25 + 0.25 * a * t_end) ** -4
    errs = [abs(r[0] - ex) for r in res]
    slope = np.polyfit(np.log(dts), np.log(errs), 1)[0]
    rich = 2 * res[-1][0] - res[-2][0]
    return [_row("V5a", "Lumped bar cooled by natural convection h ~ dT^1/4: exact dT(t) = (dT0^-1/4 + a t/4)^-4",
                 rich - C0, ex - C0, abs(rich - ex) / dT0, 1e-4, "T at 400 s in C (Richardson from dt = 2, 1 s)"),
            _row("V5b", "Same: order of the implicit-Euler time error", slope, 1.0, abs(slope - 1), 0.05,
                 "expected order 1", dt=list(dts), err=errs)]


def v6_radiation(dts=(8.0, 4.0, 2.0, 1.0)):
    p = Params().replace(nc_coef=0.0)
    dT0, t_end = 120.0, 400.0
    res = [_lump(p, dT0, t_end, dt) for dt in dts]
    _, A, Cth, _ = res[0]
    a = p.T_amb
    b = p.emissivity * SIGMA_SB * A / Cth
    G = lambda T: (1 / (2 * a * a)) * (np.log((T - a) / (T + a)) / (2 * a) - np.arctan(T / a) / a)
    T0 = a + dT0
    ex = brentq(lambda T: G(T) - G(T0) + b * t_end, a + 1e-6, T0)
    errs = [abs(r[0] - ex) for r in res]
    slope = np.polyfit(np.log(dts), np.log(errs), 1)[0]
    rich = 2 * res[-1][0] - res[-2][0]
    return [_row("V6a", "Lumped bar cooled by radiation eps sigma (T^4 - T_amb^4): exact implicit T(t)",
                 rich - C0, ex - C0, abs(rich - ex) / dT0, 1e-4, "T at 400 s in C (Richardson from dt = 2, 1 s)"),
            _row("V6b", "Same: order of the implicit-Euler time error", slope, 1.0, abs(slope - 1), 0.05,
                 "expected order 1", dt=list(dts), err=errs)]


# ---------------------------------------------------------------- V7 MMS, full coupled 3D
class _MMS:
    def __init__(self, p):
        self.p = p
        self.L, self.W, self.H = p.L, 0.04, 0.008

    def fields(self, X):
        """Returns T, grad T, lap T, psi, grad psi, lap psi at points X (..., 3)."""
        x, y, z = X[..., 0], X[..., 1], X[..., 2]
        L, W, H = self.L, self.W, self.H
        a, b = np.pi / L, np.pi / W
        T = 340.0 + 25 * np.sin(a * x) + 6 * np.cos(b * y) + 4 * (z / H) ** 2
        gT = np.stack([25 * a * np.cos(a * x), -6 * b * np.sin(b * y), 8 * z / H ** 2], -1)
        lT = -25 * a * a * np.sin(a * x) - 6 * b * b * np.cos(b * y) + 8 / H ** 2
        ps = 30 * (1 - x / L) + 3 * np.sin(a * x) * y * z / (W * H)
        gp = np.stack([-30 / L + 3 * a * np.cos(a * x) * y * z / (W * H), 3 * np.sin(a * x) * z / (W * H),
                       3 * np.sin(a * x) * y / (W * H)], -1)
        lp = -3 * a * a * np.sin(a * x) * y * z / (W * H)
        return T, gT, lT, ps, gp, lp


def v7_mms(levels=(1, 2, 4)):
    errs_T, errs_p, hs = [], [], []
    for f in levels:
        p = Params().refined(f)
        d = initial_design(p) * (1 + 0.15 * np.cos(np.arange(24)))
        mesh = build_mesh(d, p)
        M = BusbarModel(mesh, p, _steady_sched(0.0))
        mms = _MMS(p)
        s = PHI_SCALE
        g = M.geo
        T, gT, lT, ps, gp, lp = mms.fields(g.xg)
        sig, dsig, k, dk = M.sigma(T), M.dsigma(T), M.k(T), M.dk(T)
        f_e = -s * (dsig * np.sum(gT * gp, -1) + sig * lp)
        f_T = -(dk * np.sum(gT * gT, -1) + k * lT) - s * s * sig * np.sum(gp * gp, -1)
        F = np.zeros(M.ndof)
        h = mesh.hexes
        np.add.at(F, h, np.einsum("eg,ga->ea", f_e * g.wdetJ, g.N))
        np.add.at(F, h + M.nn, np.einsum("eg,ga->ea", f_T * g.wdetJ, g.N))
        V_star = 35.0
        I_half = 0.0
        for name, q in mesh.faces.items():
            if len(q) == 0:
                continue
            fg = QuadGeometry(mesh.X, q)
            Tq, gTq, _, pq, gpq, _ = mms.fields(fg.xg)
            nrm = fg.normal
            fe = s * M.sigma(Tq) * np.sum(gpq * nrm, -1)
            fT = M.k(Tq) * np.sum(gTq * nrm, -1)
            if name in EXPOSED:
                fT = fT + M.exposed_flux(Tq)
            elif name.startswith("contact"):
                Vt = V_star if name == "contact_in" else 0.0
                j = s * (pq - Vt) / p.rho_c
                fe = fe + j
                fT = fT + p.h_term * (Tq - p.T_term) - p.joint_heat_share * p.rho_c * j * j
                if name == "contact_in":
                    I_half += -np.sum(j * fg.wdA)
            np.add.at(F, q, np.einsum("fg,fg,ga->fa", fe, fg.wdA, fg.N))
            np.add.at(F, q + M.nn, np.einsum("fg,fg,ga->fa", fT, fg.wdA, fg.N))
        M.sched = _steady_sched(2 * I_half)
        M.extra_load = lambda n: F
        Xn = np.real(mesh.X)
        Tn, _, _, pn, _, _ = mms.fields(Xn)
        u0 = np.r_[pn, Tn, V_star] * 0.9 + np.r_[0 * pn, 0.1 * Tn, 0]
        U = solve_transient(M, u0).U[1]
        errs_T.append(np.sqrt(np.mean((U[M.nn:2 * M.nn] - Tn) ** 2)))
        errs_p.append(np.sqrt(np.mean((U[:M.nn] - pn) ** 2)))
        hs.append(1.0 / f)
    sT = np.polyfit(np.log(hs), np.log(errs_T), 1)[0]
    sp_ = np.polyfit(np.log(hs), np.log(errs_p), 1)[0]
    return [_row("V7a", "MMS, coupled 3D on a tapered mesh, all nonlinear terms: order of the T error", sT, 2.0,
                 abs(sT - 2) / 2, 0.1, "RMS nodal error, meshes x1 x2 x4", h_rel=hs, err_L2=errs_T),
            _row("V7b", "Same MMS: order of the potential error", sp_, 2.0, abs(sp_ - 2) / 2, 0.1,
                 "RMS nodal error of psi (mV)", h_rel=hs, err_L2=errs_p)]


# ---------------------------------------------------------------- V8 energy balance
def v8_energy():
    p = Params()
    M = BusbarModel(build_mesh(initial_design(p) * (1 + 0.1 * np.sin(np.arange(24))), p), p)
    U = solve_transient(M).U
    hf = M.heat_flows(U[1])
    VI = PHI_SCALE * U[1, M.iV] * 0.5 * M.sched.I[1]
    elec = hf["joule_bulk"] + hf["joule_joints_total"]
    rows = [_row("V8a", "Electrical power V_in I = Joule heat in copper + joints (steady rating)", elec, VI,
                 abs(elec - VI) / VI, 1e-10, "W, half model")]
    Mm = M.mass_matrix()
    one = np.zeros(M.ndof); one[M.nn:2 * M.nn] = 1.0
    bal, scale = 0.0, 0.0
    for n in range(2, M.sched.n_steps + 1):
        h = M.heat_flows(U[n])
        src = h["joule_bulk"] + h["joule_joints_into_bar"]
        out = sum(v for k_, v in h.items() if k_.startswith("out_"))
        stored = one @ (Mm @ (U[n] - U[n - 1]))
        bal += stored - M.sched.dt[n] * (src - out)
        scale += M.sched.dt[n] * src
    rows.append(_row("V8b", "Transient energy: stored + lost = generated over pulse and recovery", 0.0, 0.0,
                     abs(bal) / scale, 1e-10, "J imbalance / J generated"))
    return rows


# ---------------------------------------------------------------- V9 Newton, V10 mesh, V11 tangent, V12 faces
def v9_newton():
    p = Params()
    M = BusbarModel(build_mesh(initial_design(p), p), p)
    U = solve_transient(M).U
    _, h = newton_step(M, U[1], M.sched.dt[2], 2)
    r = np.array(h)
    orders = [np.log(r[i + 1] / r[i]) / np.log(r[i] / r[i - 1]) for i in range(1, len(r) - 1)
              if r[i + 1] < 1e-2 * r[i]]
    order = max(orders) if orders else 2.0
    return _row("V9", "Newton on the coupled pulse step: quadratic convergence", order, 2.0, max(0.0, 2 - order) / 2,
                0.15, "observed order; residuals " + ", ".join(f"{x:.1e}" for x in r), residuals=r.tolist(), step=2)


def v10_mesh(levels=(1, 2, 3)):
    peaks, Rs = [], []
    for f in levels:
        p = Params().refined(f)
        M = BusbarModel(build_mesh(initial_design(p), p), p)
        U = solve_transient(M).U
        peaks.append(U[1:, M.nn:2 * M.nn].max() - C0)
        Rs.append(M.resistance(U[1], 1) * 1e6)
    rich = peaks[-1] + (peaks[-1] - peaks[-2]) / ((levels[-1] / levels[-2]) ** 2 - 1)
    return _row("V10", "Mesh convergence of the peak temperature (x1, x2, x3 meshes)", peaks[0], rich,
                abs(peaks[0] - rich) / (rich - 40.0), 0.02,
                "peaks " + ", ".join(f"{x:.2f}" for x in peaks) + " C; R " + ", ".join(f"{x:.3f}" for x in Rs) +
                " micro-ohm; error of the optimizer mesh vs Richardson, relative to the rise", peaks=peaks, R=Rs)


def v11_tangent():
    p = Params()
    M = BusbarModel(build_mesh(initial_design(p) * (1 + 0.1 * np.sin(np.arange(24))), p), p)
    U = solve_transient(M).U
    n = 20
    J1 = M.jacobian(U[n], U[n - 1], M.sched.dt[n])
    J2 = M.jacobian(U[n], U[n - 1], M.sched.dt[n], cs_volume=True)
    err = abs(J1 - J2).max() / abs(J2).max()
    rng = np.random.default_rng(0)
    v = rng.standard_normal(M.ndof) * np.r_[np.full(M.nn, 1.0), np.full(M.nn, 1.0), 1.0]
    eps = 1e-6
    fd = (M.residual(U[n] + eps * v, U[n - 1], M.sched.dt[n], n) -
          M.residual(U[n] - eps * v, U[n - 1], M.sched.dt[n], n)) / (2 * eps)
    e2 = np.max(np.abs(J1 @ v - fd)) / np.max(np.abs(fd))
    return [_row("V11a", "Analytic Newton tangent vs complex step (volume terms)", err, 0.0, err, 1e-12, "max |dJ| / max |J|"),
            _row("V11b", "Newton tangent vs central differences of the residual (random direction)", e2, 0.0, e2, 1e-6,
                 "max error / max |J v|")]


def v12_faces():
    p = Params()
    mesh = build_mesh(initial_design(p), p)
    from .fem import HEX_FACES
    cnt = {}
    for e in mesh.hexes:
        for loc in HEX_FACES.values():
            key = tuple(sorted(e[loc]))
            cnt[key] = cnt.get(key, 0) + 1
    ext = {k for k, v in cnt.items() if v == 1}
    sets = [tuple(sorted(f)) for q in mesh.faces.values() for f in q]
    ok = len(sets) == len(set(sets)) == len(ext) and set(sets) == ext
    return _row("V12", "Every exterior face belongs to exactly one boundary set", len(sets), len(ext),
                0.0 if ok else 1.0, 0.0, "face counts")


def run_ladder(include_slow=True):
    rows = [v1_resistance(), v2_invariance(), v3_joule_1d()]
    rows += v4_coupled_1d()
    rows += v5_natural_convection()
    rows += v6_radiation()
    rows += v7_mms((1, 2, 4) if include_slow else (1, 2))
    rows += v8_energy()
    rows.append(v9_newton())
    if include_slow:
        rows.append(v10_mesh())
    rows += v11_tangent()
    rows.append(v12_faces())
    conv = {}
    for r in rows:
        if r["id"] == "V7a":
            conv["mms_space"] = {"h_rel": r["h_rel"], "err_L2": r["err_L2"], "slope": r["ours"]}
        if r["id"] == "V7b":
            conv["mms_space_psi"] = {"h_rel": r["h_rel"], "err_L2": r["err_L2"], "slope": r["ours"]}
        if r["id"] == "V4b":
            conv["coupled_1d"] = {"h_rel": r["h_rel"], "err_L2": r["err_L2"], "slope": r["ours"]}
        if r["id"] == "V5b":
            conv["natconv_time"] = {"dt": r["dt"], "err": r["err"], "slope": r["ours"]}
        if r["id"] == "V6b":
            conv["radiation_time"] = {"dt": r["dt"], "err": r["err"], "slope": r["ours"]}
        if r["id"] == "V9":
            conv["newton"] = {"residuals": r["residuals"], "step": r["step"]}
        if r["id"] == "V10":
            conv["mesh"] = {"peaks_C": r["peaks"], "R_micro_ohm": r["R"]}
    slim = [{k: v for k, v in r.items() if k not in ("h_rel", "err_L2", "dt", "err", "residuals", "peaks", "R")}
            for r in rows]
    return slim, conv


if __name__ == "__main__":
    import sys
    rows, _ = run_ladder(include_slow="--quick" not in sys.argv)
    for r in rows:
        print(f"{r['id']:5s} {'PASS' if r['pass'] else 'FAIL'} ours={r['ours']:.6g} exact={r['exact']:.6g} "
              f"rel={r['rel_error']:.2e}  {r['name']}")
