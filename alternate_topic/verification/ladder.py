"""Verification ladder for the brake-disc thermal solver (WS1).

Every test runs the *production* ThermalModel / Newton solver (brakeopt.model, brakeopt.solver),
either on the real sector mesh or on box meshes built with the same Mesh class. Each test returns
a row {"id", "name", "ours", "exact", "rel_error", "pass", "note"} for the report table; the
convergence studies also return their error-vs-h / error-vs-dt data for plots.

  V1  steady 1D slab, convection both ends, constant k         exact linear profile (round-off)
  V2  V1 under an arbitrary 3D rotation + translation          identical temperatures (invariance)
  V3  nonlinear k(T) slab (Kirchhoff transform)                 exact flux and profile; O(h^2)
  V4  radiation-cooled plate (lumped, T_amb = 0)                exact T(t); implicit Euler O(dt)
  V5  MMS steady 3D on the distorted sector mesh, k(T)+rad      L2 error slope 2 in h
  V6  MMS transient (all nonlinear terms)                       error slope 1 in dt
  V7  energy balance over the braking sequence                  heat in = stored + lost; KE budget
  V8  Newton convergence on the real problem                    quadratic (order ~2)
  V9  mesh convergence of the peak friction-face temperature    Richardson order ~2 (slow)
  V10 boundary closure: every exterior face in exactly one set  counts match
"""
from __future__ import annotations

import numpy as np

from brakeopt.fem import HEX_FACES
from brakeopt.loads import Schedule, energy_into_one_face
from brakeopt.mesh import Mesh, build_mesh
from brakeopt.model import ThermalModel, RAD_SETS
from brakeopt.params import Params, SIGMA_SB, initial_design
from brakeopt.solver import solve_transient, newton_step

PHYS_SETS = ("friction", "face_free", "channel", "rim", "hub")
STEADY_DT = 1e30


def row(id_, name, ours, exact, tol, note="", rel=None):
    if rel is None:
        rel = abs(ours - exact) / max(abs(exact), 1e-300)
    return {"id": id_, "name": name, "ours": float(ours), "exact": float(exact), "rel_error": float(rel),
            "tol": tol, "pass": bool(rel <= tol), "note": note}


# ----------------------------------------------------------------------------- helpers
def box_mesh(L, n, mapping, Rot=None, shift=(0, 0, 0)) -> Mesh:
    """Box [0,Lx]x[0,Ly]x[0,Lz] of n hexes; mapping side -> physics set ('x-': 'hub', ...)."""
    nx, ny, nz = n
    I, J, K = np.meshgrid(np.arange(nx + 1), np.arange(ny + 1), np.arange(nz + 1), indexing="ij")
    ijk = np.stack([I.ravel(), J.ravel(), K.ravel()], 1)
    X = ijk * (np.array(L) / np.array(n))
    nid = lambda i, j, k: (i * (ny + 1) + j) * (nz + 1) + k
    off = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]])
    E = np.array([(i, j, k) for i in range(nx) for j in range(ny) for k in range(nz)])
    c = E[:, None, :] + off[None]
    hexes = nid(c[..., 0], c[..., 1], c[..., 2])
    sides = {"x-": (E[:, 0] == 0, "xi-"), "x+": (E[:, 0] == nx - 1, "xi+"),
             "y-": (E[:, 1] == 0, "eta-"), "y+": (E[:, 1] == ny - 1, "eta+"),
             "z-": (E[:, 2] == 0, "zeta-"), "z+": (E[:, 2] == nz - 1, "zeta+")}
    faces = {s: np.zeros((0, 4), int) for s in PHYS_SETS}
    for side, (mask, key) in sides.items():
        q = hexes[mask][:, HEX_FACES[key]]
        tgt = mapping.get(side, "sym_" + side)
        faces[tgt] = np.vstack([faces[tgt], q]) if tgt in faces else q
    if Rot is not None:
        X = X @ np.asarray(Rot).T
    X = X + np.asarray(shift)
    return Mesh(X=X, hexes=hexes, faces=faces, node_rtz=X.copy(), node_ijk=ijk, elem_ijk=E,
                dX_dd=np.zeros((len(X), 3, 24)), sector_multiplicity=1, n_vanes=1)


def const_schedule(n_steps, dt, h_face=0.0, h_chan=0.0):
    t = np.arange(n_steps + 1) * dt
    z = np.zeros(n_steps + 1)
    return Schedule(t=t, dt=np.r_[np.nan, np.full(n_steps, dt)], v=z, braking=z.astype(bool),
                    P_face=z.copy(), h_face=np.full(n_steps + 1, h_face), h_chan=np.full(n_steps + 1, h_chan))


def steady(model, T_guess):
    T, hist = newton_step(model, np.full(model.nn, T_guess, dtype=float), STEADY_DT, 1)
    return T, hist


def random_rotation(seed=7):
    q, r = np.linalg.qr(np.random.default_rng(seed).standard_normal((3, 3)))
    q = q * np.sign(np.diag(r))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


def mms_load(model, fsrc, gflux, n):
    """Nodal load: int N f dV + sum_sets int N g dA, f(x)->(e,g), g(name, xg, normal)->(f,g)."""
    geo = model.geo
    F = model.S_hex @ np.einsum("eg,ga->ea", fsrc(geo.xg) * geo.wdetJ, geo.N).ravel()
    for name, fg in model.fgeo.items():
        if len(model.mesh.faces[name]) == 0:
            continue
        g = gflux(name, fg.xg, fg.normal, n)
        F = F + model.S_face[name] @ np.einsum("fg,ga->fa", g * fg.wdA, fg.N).ravel()
    return F


# ----------------------------------------------------------------------------- V1, V2
def v1_v2_slab():
    L = 0.05
    p = Params().replace(k_beta=0.0, h_hub=500.0, T_bell=300.0, T_amb=800.0, T_init=500.0)
    h2 = 200.0
    k = p.k0
    q = (p.T_amb - p.T_bell) / (1 / p.h_hub + L / k + 1 / h2)
    T0 = p.T_bell + q / p.h_hub
    exact = lambda x: T0 + q * x / k
    mapping = {"x-": "hub", "x+": "channel"}
    sched = const_schedule(1, STEADY_DT, h_chan=h2)
    m = box_mesh((L, 0.01, 0.01), (6, 2, 2), mapping)
    model = ThermalModel(m, p, sched)
    T, _ = steady(model, 500.0)
    err = np.max(np.abs(T - exact(m.X[:, 0]))) / (p.T_amb - p.T_bell)
    out = [row("V1", "Steady slab, convection both ends (exact linear profile)",
               T[np.argmax(m.X[:, 0])], exact(L), 1e-10,
               f"max nodal error / dT = {err:.1e}; q = {q:.1f} W/m^2", rel=err)]
    Rot = random_rotation()
    m2 = box_mesh((L, 0.01, 0.01), (6, 2, 2), mapping, Rot=Rot, shift=(0.3, -0.2, 0.1))
    T2, _ = steady(ThermalModel(m2, p, sched), 500.0)
    inv = np.max(np.abs(T2 - T)) / (p.T_amb - p.T_bell)
    out.append(row("V2", "Invariance: same slab under arbitrary 3D rotation + translation",
                   T2.max(), T.max(), 1e-10, f"max nodal difference / dT = {inv:.1e}", rel=inv))
    return out


# ----------------------------------------------------------------------------- V3
def v3_kirchhoff(ns=(4, 8, 16, 32)):
    """Steady slab with k(T) = k0 (1 + beta (T - Tr)), near-Dirichlet ends (h = 1e10)."""
    L = 0.02
    p = Params().replace(h_hub=1e10, T_bell=400.0, T_amb=1000.0, k_beta=-5e-4)
    k0, b, Tr = p.k0, p.k_beta, p.T_kref
    U = lambda T: k0 * ((T - Tr) + 0.5 * b * (T - Tr) ** 2)
    Uinv = lambda u: Tr + (-1 + np.sqrt(1 + 2 * b * u / k0)) / b
    u1, u2 = U(p.T_bell), U(p.T_amb)
    q_exact = (u2 - u1) / L
    Tex = lambda x: Uinv(u1 + (u2 - u1) * x / L)
    hs, errs, qs = [], [], []
    for n in ns:
        m = box_mesh((L, 0.004, 0.004), (n, 1, 1), {"x-": "hub", "x+": "channel"})
        model = ThermalModel(m, p, const_schedule(1, STEADY_DT, h_chan=1e10))
        T, _ = steady(model, 700.0)
        geo = model.geo
        Tg = np.einsum("ga,ea->eg", geo.N, T[m.hexes])
        e2 = np.sum(geo.wdetJ * (Tg - Tex(geo.xg[..., 0])) ** 2) / np.sum(geo.wdetJ)
        hs.append(L / n)
        errs.append(float(np.sqrt(e2)))
        fl = model.heat_flows(T, 1)
        qs.append(fl["conv_hub"] / (0.004 * 0.004))    # heat leaving through the cold end
    slope = np.polyfit(np.log(hs), np.log(errs), 1)[0]
    rows = [row("V3a", "Nonlinear k(T) slab: heat flux vs Kirchhoff exact solution", qs[-1], q_exact, 1e-4,
                f"finest mesh n = {ns[-1]}"),
            row("V3b", "Nonlinear k(T) slab: L2 error convergence order", slope, 2.0, 0.1,
                f"L2 error {errs[0]:.2e} -> {errs[-1]:.2e} K")]
    return rows, {"h": hs, "err_L2": errs, "slope": float(slope)}


# ----------------------------------------------------------------------------- V4
def v4_radiation(Ns=(25, 50, 100, 200, 400)):
    Lz, T0, t_end = 0.002, 1000.0, 250.0
    p = Params().replace(k0=1e5, k_beta=0.0, T_amb=0.0, T_init=T0)
    a = p.emissivity * SIGMA_SB / (p.rho * p.c * Lz)
    T_exact = (T0 ** -3 + 3 * a * t_end) ** (-1 / 3)
    dts, errs, Tend = [], [], []
    for N in Ns:
        m = box_mesh((0.01, 0.01, Lz), (1, 1, 1), {"z+": "rim"})
        model = ThermalModel(m, p, const_schedule(N, t_end / N, h_face=0.0))
        sol = solve_transient(model)
        Tn = sol.T[-1].mean()
        dts.append(t_end / N)
        errs.append(abs(Tn - T_exact))
        Tend.append(Tn)
    slope = np.polyfit(np.log(dts), np.log(errs), 1)[0]
    # Richardson-extrapolated end temperature (first order)
    T_rich = 2 * Tend[-1] - Tend[-2]
    rows = [row("V4a", "Radiation-cooled plate: T(t_end) vs exact (Richardson, dt -> 0)", T_rich, T_exact, 1e-5,
                f"T_exact = {T_exact:.3f} K after {t_end:.0f} s; raw finest dt error {errs[-1]:.2e} K"),
            row("V4b", "Radiation-cooled plate: implicit-Euler order in dt", slope, 1.0, 0.05,
                f"error {errs[0]:.2e} -> {errs[-1]:.2e} K")]
    return rows, {"dt": dts, "err": [float(e) for e in errs], "slope": float(slope)}


# ----------------------------------------------------------------------------- V5 MMS steady on sector
def _mms_fields(p):
    A, B, C, T0 = 100.0, 80.0, 60.0, 650.0
    a, b, c = 40.0, 200.0, 150.0

    def T(x):
        return T0 + A * np.sin(a * x[..., 0]) * np.cos(b * x[..., 1]) + B * np.cos(c * x[..., 2]) + C * x[..., 0] / 0.1

    def grad(x):
        gx = A * a * np.cos(a * x[..., 0]) * np.cos(b * x[..., 1]) + C / 0.1
        gy = -A * b * np.sin(a * x[..., 0]) * np.sin(b * x[..., 1])
        gz = -B * c * np.sin(c * x[..., 2])
        return np.stack([gx, gy, gz], -1)

    def lap(x):
        return (-A * (a ** 2 + b ** 2) * np.sin(a * x[..., 0]) * np.cos(b * x[..., 1])
                - B * c ** 2 * np.cos(c * x[..., 2]))
    return T, grad, lap


def v5_mms_space(factors=(1, 2, 3, 4)):
    p0 = Params()
    Tf, gradf, lapf = _mms_fields(p0)
    hs, errs = [], []
    for fac in factors:
        p = p0.refined(fac)
        m = build_mesh(initial_design(p), p)
        model = ThermalModel(m, p)
        n_load = 10   # take h_face, h_chan of a braking step; P_face is switched off below
        model.sched.P_face[:] = 0.0

        def fsrc(x):
            k, dk = model.k_of_T(Tf(x))
            g = gradf(x)
            return -(k * lapf(x) + dk * np.sum(g * g, -1))

        F = mms_load(model, fsrc, lambda nm, x, nr, n: _g(model, nm, x, nr, n, Tf, gradf), n_load)
        model.extra_load = lambda n: F
        T, _ = newton_step(model, np.full(model.nn, 650.0), STEADY_DT, n_load)
        geo = model.geo
        Tg = np.einsum("ga,ea->eg", geo.N, T[m.hexes])
        e2 = np.sum(geo.wdetJ * (Tg - Tf(geo.xg)) ** 2) / np.sum(geo.wdetJ)
        hs.append(1.0 / fac)
        errs.append(float(np.sqrt(e2)))
    slope = np.polyfit(np.log(hs[1:]), np.log(errs[1:]), 1)[0]
    rows = [row("V5", "MMS steady 3D on sector mesh (k(T), radiation, convection): L2 order in h", slope, 2.0, 0.15,
                f"L2 error {errs[0]:.2e} -> {errs[-1]:.2e} K over refinement x{factors[0]}..x{factors[-1]}")]
    return rows, {"h_rel": hs, "err_L2": errs, "slope": float(slope)}


def _g(model, name, x, nrm, n, Tf, gradf):
    """MMS boundary flux g = BC(T_ex) + k(T_ex) grad T_ex . n."""
    Tq = Tf(x)
    k, _ = model.k_of_T(Tq)
    kgn = k * np.sum(gradf(x) * nrm, -1)
    if name.startswith("sym"):
        return kgn
    return model.bc_flux(name, Tq[None], np.array([n]))[0] + kgn


# ----------------------------------------------------------------------------- V6 MMS transient
def v6_mms_time(Ns=(10, 20, 40, 80, 160)):
    Lx, tau, t_end = 0.02, 2.0, 2.0
    p = Params().replace(T_init=500.0)

    def Tf(x, t):
        return 500.0 + 200.0 * (x[..., 0] / Lx) * (1 - np.exp(-t / tau)) + 50.0 * np.sin(t / tau)

    def dTdt(x, t):
        return 200.0 * (x[..., 0] / Lx) * np.exp(-t / tau) / tau + 50.0 * np.cos(t / tau) / tau

    def gradf(x, t):
        g = np.zeros(x.shape)
        g[..., 0] = 200.0 / Lx * (1 - np.exp(-t / tau))
        return g

    mapping = {"x-": "hub", "x+": "channel", "z+": "rim"}
    dts, errs = [], []
    for N in Ns:
        dt = t_end / N
        m = box_mesh((Lx, 0.004, 0.004), (4, 1, 1), mapping)
        model = ThermalModel(m, p, const_schedule(N, dt, h_face=80.0, h_chan=150.0))
        loads = {}
        for n in range(1, N + 1):
            t = n * dt

            def fsrc(x, t=t):
                k, dk = model.k_of_T(Tf(x, t))
                g = gradf(x, t)
                return p.rho * p.c * dTdt(x, t) - dk * np.sum(g * g, -1)

            loads[n] = mms_load(model, fsrc,
                                lambda nm, x, nr, n, t=t: _g(model, nm, x, nr, n, lambda y: Tf(y, t), lambda y: gradf(y, t)),
                                n)
        model.extra_load = lambda n: loads[n]
        sol = solve_transient(model, T0=Tf(m.X, 0.0))
        errs.append(float(np.max(np.abs(sol.T[-1] - Tf(m.X, t_end)))))
        dts.append(dt)
    slope = np.polyfit(np.log(dts), np.log(errs), 1)[0]
    rows = [row("V6", "MMS transient (k(T), radiation, convection): implicit-Euler order in dt", slope, 1.0, 0.05,
                f"max nodal error at t_end {errs[0]:.2e} -> {errs[-1]:.2e} K")]
    return rows, {"dt": dts, "err": errs, "slope": float(slope)}


# ----------------------------------------------------------------------------- V7, V8 on the real problem
def v7_v8_real(p: Params | None = None):
    p = Params() if p is None else p
    m = build_mesh(initial_design(p), p)
    model = ThermalModel(m, p)
    sol = solve_transient(model)
    s = model.sched
    M = model.mass_matrix()
    E_in = sum(s.dt[n] * model.heat_flows(sol.T[n], n)["in"] for n in range(1, s.n_steps + 1))
    E_out = sum(s.dt[n] * sum(v for k, v in model.heat_flows(sol.T[n], n).items() if k != "in")
                for n in range(1, s.n_steps + 1))
    E_store = float(np.ones(model.nn) @ (M @ (sol.T[-1] - sol.T[0])))
    E_budget = energy_into_one_face(p) * model.sector_share
    imbalance = abs(E_in - E_store - E_out) / E_in
    rows = [row("V7a", "Heat input = braking kinetic-energy budget (J, quarter sector)", E_in, E_budget, 1e-12,
                "sum of dt * frictional power vs partition x front bias x KE / (faces x sectors)"),
            row("V7b", "Energy balance over the braking sequence: stored + lost vs input (J)", E_store + E_out, E_in, 1e-9,
                f"stored {E_store:.1f} J, lost {E_out:.1f} J", rel=imbalance)]
    # Newton order from the step with the longest history (during braking)
    n_star = max(range(1, s.n_steps + 1), key=lambda n: len(sol.newton[n]))
    r = np.array(sol.newton[n_star])
    r = r[r > 1e-11 * r[0]]
    orders = np.log(r[2:] / r[1:-1]) / np.log(r[1:-1] / r[:-2]) if len(r) >= 3 else np.array([np.nan])
    order = float(orders[-1])
    rows.append(row("V8", "Newton convergence order on the braking problem", order, 2.0, 0.25,
                    "residual norms " + ", ".join(f"{x:.1e}" for x in sol.newton[n_star])))
    return rows, {"newton": {"step": int(n_star), "residuals": [float(x) for x in sol.newton[n_star]]},
                  "energy": {"in": E_in, "stored": E_store, "lost": E_out}}, (model, sol)


def v9_mesh_convergence(factors=(1, 2, 3)):
    from brakeopt.functionals import peak_temperatures
    p0 = Params()
    peaks = []
    for f in factors:
        p = p0.refined(f)
        m = build_mesh(initial_design(p), p)
        sol = solve_transient(ThermalModel(m, p))
        peaks.append(peak_temperatures(m, sol.T)[0])
    e1, e2 = peaks[0] - peaks[1], peaks[1] - peaks[2]
    r = factors[1] / factors[0]
    order = float(np.log(abs(e1 / e2)) / np.log(factors[2] / factors[1])) if len(factors) == 3 else np.nan
    # Richardson with the observed order on the two finest meshes
    q = factors[2] / factors[1]
    T_inf = peaks[2] + (peaks[2] - peaks[1]) / (q ** order - 1)
    rows = [row("V9", "Mesh convergence: peak friction-face T, coarse mesh vs Richardson limit (K)", peaks[0], T_inf,
                0.02, f"peaks {', '.join(f'{x - 273.15:.1f}' for x in peaks)} C for refinement {factors}; "
                      f"observed order {order:.2f}")]
    return rows, {"factors": list(factors), "peak_face_K": peaks, "order": order, "richardson_K": float(T_inf)}


def v10_closure(p: Params | None = None):
    p = Params() if p is None else p
    m = build_mesh(initial_design(p), p)
    allf = np.vstack([m.hexes[:, f] for f in HEX_FACES.values()])
    key = np.sort(allf, axis=1)
    _, cnt = np.unique(key, axis=0, return_counts=True)
    n_ext = int(np.sum(cnt == 1))
    n_sets = int(sum(len(v) for v in m.faces.values()))
    return [row("V10", "Boundary closure: exterior faces vs faces in boundary sets", n_sets, n_ext, 0.0,
                f"{n_ext} exterior faces; sets: " + ", ".join(f"{k} {len(v)}" for k, v in m.faces.items()))]


def run_ladder(include_slow=True):
    rows, conv = [], {}
    rows += v1_v2_slab()
    r, conv["kirchhoff"] = v3_kirchhoff(); rows += r
    r, conv["radiation_time"] = v4_radiation(); rows += r
    r, conv["mms_space"] = v5_mms_space(); rows += r
    r, conv["mms_time"] = v6_mms_time(); rows += r
    r, extra, _ = v7_v8_real(); rows += r
    conv.update(extra)
    if include_slow:
        r, conv["mesh_convergence"] = v9_mesh_convergence(); rows += r
    rows += v10_closure()
    return rows, conv


if __name__ == "__main__":
    rows, conv = run_ladder(include_slow=False)
    for r_ in rows:
        print(f"{r_['id']:5s} {'PASS' if r_['pass'] else 'FAIL'}  ours={r_['ours']:.6g} exact={r_['exact']:.6g} "
              f"rel={r_['rel_error']:.2e}  {r_['name']} | {r_['note']}")
