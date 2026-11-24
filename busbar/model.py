"""Two-way coupled electro-thermal busbar model: residual, exact tangent, current schedule.

State of one time level (ndof = 2 n_nodes + 1):
    u = [psi (n_nodes, potential in mV), T (n_nodes, K), V_in (1, inlet terminal potential in mV)]
Strong form (phi = PHI_SCALE * psi):
    div(sigma(T) grad phi) = 0
    rho c dT/dt - div(k(T) grad T) = sigma(T) |grad phi|^2
Boundary conditions
    joint faces (lap joints, z = 0):  j_n = (phi - V_t) / rho_c   (outward current density)
                                       -k dT/dn = h_term (T - T_term) - share * rho_c j_n^2
        inlet V_t = V_in (unknown, fixed by the total current), outlet V_t = 0
    exposed faces: insulated electrically; -k dT/dn = h(T) (T - T_amb) + eps sigma_SB (T^4 - T_amb^4),
        h(T) = C (|T - T_amb|^2 + delta^2)^(1/8),  C = nc_coef / nc_length^(1/4)
    y = 0: mirror symmetry (natural).
Global row: (half model)  -int_inlet j_n dA - I(t)/2 = 0.
Implicit Euler for step n: R_n(u_n, u_{n-1}) = 0 (only T carries a time derivative).

The tangent dR/du is exact: the volume terms are differentiated analytically, the boundary terms
by complex step over each face's 8 (9 with V_in) local unknowns (exact to round-off); the
verification ladder checks the analytic part against complex step. The residual is
complex-step safe in the node coordinates (for dR/dd) and in the state (for the tangent).
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import scipy.sparse as sp

from brakeopt.fem import HexGeometry, QuadGeometry
from .mesh import Mesh
from .params import Params, SIGMA_SB, PHI_SCALE

EXPOSED = ("bottom", "top", "side", "ends")
CS = 1e-30


@dataclass
class Schedule:
    t: np.ndarray        # (N+1,) t[1] = 0 is the end of the steady step
    dt: np.ndarray       # (N+1,), dt[0] = nan
    I: np.ndarray        # (N+1,) total current during step n (A)
    pulse: np.ndarray    # (N+1,) bool, step in the boost pulse

    @property
    def n_steps(self):
        return len(self.t) - 1


def current_schedule(p: Params) -> Schedule:
    n1 = int(round(p.t_pulse / p.dt_pulse))
    n2 = int(round(p.t_after / p.dt_after))
    dt = np.r_[np.nan, p.dt_steady, np.full(n1, p.dt_pulse), np.full(n2, p.dt_after)]
    t = np.r_[-p.dt_steady, 0.0, np.cumsum(dt[2:])]
    I = np.r_[p.I_base, p.I_base, np.full(n1, p.I_peak), np.full(n2, p.I_base)]
    pulse = np.r_[False, False, np.ones(n1, bool), np.zeros(n2, bool)]
    return Schedule(t=t, dt=dt, I=I, pulse=pulse)


class BusbarModel:
    def __init__(self, mesh: Mesh, p: Params, sched: Schedule | None = None):
        self.mesh, self.p = mesh, p
        self.sched = sched if sched is not None else current_schedule(p)
        self.nn = mesh.n_nodes
        self.ndof = 2 * self.nn + 1
        self.iV = 2 * self.nn
        self.geo = HexGeometry(mesh.X, mesh.hexes)
        g = self.geo
        self.Me = p.rho * p.c * np.einsum("eg,ga,gb->eab", g.wdetJ, g.N, g.N)
        h = mesh.hexes
        self.vdofs = np.hstack([h, h + self.nn])                       # (e, 16)
        self.fgeo, self.fdofs = {}, {}
        for name, q in mesh.faces.items():
            if name.startswith("sym") or len(q) == 0:
                continue
            self.fgeo[name] = QuadGeometry(mesh.X, q)
            dofs = np.hstack([q, q + self.nn])
            if name == "contact_in":
                dofs = np.hstack([dofs, np.full((len(q), 1), self.iV)])
            self.fdofs[name] = dofs
        self.extra_load = None      # verification: callable n -> (ndof,) load added to the RHS
        self.C_nc = p.nc_coef / p.nc_length ** 0.25
        self._M = None

    # ------------------------------------------------------------------ material laws
    def sigma(self, T):
        p = self.p
        Te = T if p.nonlinear else p.T_lin + 0.0 * T
        return p.sigma0 / (1.0 + p.alpha * (Te - p.T_ref))

    def dsigma(self, T):
        p = self.p
        if not p.nonlinear:
            return 0.0 * T
        return -p.sigma0 * p.alpha / (1.0 + p.alpha * (T - p.T_ref)) ** 2

    def dk(self, T):
        p = self.p
        return (p.k0 * p.k_beta if p.nonlinear else 0.0) + 0.0 * T

    def k(self, T):
        p = self.p
        Te = T if p.nonlinear else p.T_lin + 0.0 * T
        return p.k0 * (1.0 + p.k_beta * (Te - p.T_ref))

    def exposed_flux(self, T):
        """Outward heat flux on exposed faces, W/m^2 (natural convection + radiation)."""
        p = self.p
        es = p.emissivity * SIGMA_SB
        if p.nonlinear:
            dT = T - p.T_amb
            h = self.C_nc * (dT * dT + p.nc_delta ** 2) ** 0.125
            return h * dT + es * (T ** 4 - p.T_amb ** 4)
        dl = p.dT_lin_h
        h = self.C_nc * (dl * dl + p.nc_delta ** 2) ** 0.125
        Tl = p.T_lin_rad
        return h * (T - p.T_amb) + es * (Tl ** 4 - p.T_amb ** 4 + 4 * Tl ** 3 * (T - Tl))

    # ------------------------------------------------------------------ local residuals
    def _vol_local(self, ue, Toe, dt, g=None, Me=None):
        """ue (Z, e, 16) = [psi_a, T_a]; Toe (Z, e, 8); dt (Z,). Returns (Z, e, 16)."""
        g = self.geo if g is None else g
        Me = self.Me if Me is None else Me
        s = PHI_SCALE
        ps, Te = ue[..., :8], ue[..., 8:]
        Tg = np.einsum("ga,Zea->Zeg", g.N, Te)
        gT = np.einsum("egai,Zea->Zegi", g.dNdx, Te)
        gp = np.einsum("egai,Zea->Zegi", g.dNdx, ps)
        sg = self.sigma(Tg) * g.wdetJ[None]
        kg = self.k(Tg) * g.wdetJ[None]
        r_p = s * np.einsum("egai,Zegi->Zea", g.dNdx, sg[..., None] * gp)
        joule = s * s * sg * np.einsum("Zegi,Zegi->Zeg", gp, gp)
        r_T = np.einsum("egai,Zegi->Zea", g.dNdx, kg[..., None] * gT)
        r_T = r_T + np.einsum("eab,Zeb->Zea", Me, Te - Toe) / dt[:, None, None]
        r_T = r_T - np.einsum("ga,Zeg->Zea", g.N, joule)
        return np.concatenate([r_p, r_T], axis=-1)

    def _face_local(self, name, uf):
        """uf (Z, f, 8 or 9) = [psi_a, T_a, (V_in)]. Returns same shape."""
        p, fg = self.p, self.fgeo[name]
        s = PHI_SCALE
        ps, Tf = uf[..., :4], uf[..., 4:8]
        Tq = np.einsum("ga,Zfa->Zfg", fg.N, Tf)
        if name in EXPOSED:
            qT = self.exposed_flux(Tq)
            r_T = np.einsum("Zfg,fg,ga->Zfa", qT, fg.wdA, fg.N)
            return np.concatenate([0.0 * ps, r_T], axis=-1)
        pq = np.einsum("ga,Zfa->Zfg", fg.N, ps)
        Vt = uf[..., 8:9] if name == "contact_in" else 0.0
        j = s * (pq - Vt) / p.rho_c                       # outward current density (A/m^2)
        qT = p.h_term * (Tq - p.T_term) - p.joint_heat_share * p.rho_c * j * j
        r_p = np.einsum("Zfg,fg,ga->Zfa", j, fg.wdA, fg.N)
        r_T = np.einsum("Zfg,fg,ga->Zfa", qT, fg.wdA, fg.N)
        out = [r_p, r_T]
        if name == "contact_in":
            out.append(-np.einsum("Zfg,fg->Zf", j, fg.wdA)[..., None])
        return np.concatenate(out, axis=-1)

    # ------------------------------------------------------------------ global residual
    def residual(self, u, u_old, dt, n):
        """R_n(u_n, u_{n-1}); u, u_old (ndof,) or (B, ndof); dt, n scalars or (B,)."""
        u = np.asarray(u)
        single = u.ndim == 1
        U = np.atleast_2d(u)
        B = U.shape[0]
        Uo = np.broadcast_to(np.atleast_2d(u_old), U.shape)
        n = np.broadcast_to(np.atleast_1d(n), (B,))
        dt = np.broadcast_to(np.atleast_1d(np.asarray(dt, dtype=float)), (B,))
        dtype = np.result_type(U, self.mesh.X)
        R = np.zeros((B, self.ndof), dtype=dtype)
        To = Uo[:, self.nn:2 * self.nn][:, self.mesh.hexes]
        re = self._vol_local(U[:, self.vdofs], To, dt)
        self._scatter_add(R, self.vdofs, re)
        for name, dofs in self.fdofs.items():
            self._scatter_add(R, dofs, self._face_local(name, U[:, dofs]))
        R[:, self.iV] -= 0.5 * self.sched.I[n]
        if self.extra_load is not None:
            R = R - np.stack([self.extra_load(int(m)) for m in n])
        return R[0] if single else R

    @staticmethod
    def _scatter_add(R, dofs, vals):
        flat = dofs.ravel()
        for b in range(R.shape[0]):
            R[b] += np.bincount(flat, weights=np.real(vals[b]).ravel(), minlength=R.shape[1])
            if np.iscomplexobj(R):
                R[b] += 1j * np.bincount(flat, weights=np.imag(vals[b]).ravel(), minlength=R.shape[1])

    # ------------------------------------------------------------------ tangent
    def _local_jac(self, f, ul, *args):
        """Complex-step Jacobian of a local residual f(ul (Z, m, a)) -> (m, a, a)."""
        m, a = ul.shape
        Z = np.repeat(ul[None].astype(complex), a, axis=0)            # (a, m, a)
        Z[np.arange(a), :, np.arange(a)] += 1j * CS
        r = f(Z, *args)                                                # (a, m, a_out)
        return np.transpose(np.imag(r) / CS, (1, 2, 0))               # (m, out, in)

    def vol_tangent(self, u, dt):
        """Analytic element tangent (e, 16, 16) of the volume terms."""
        g, s = self.geo, PHI_SCALE
        ue = u[self.vdofs]
        ps, Te = ue[:, :8], ue[:, 8:]
        Tg = Te @ g.N.T
        G = g.dNdx                                               # (e, g, 8, 3)
        gT = np.einsum("egai,ea->egi", G, Te)
        gp = np.einsum("egai,ea->egi", G, ps)
        w = g.wdetJ
        sg, dsg = self.sigma(Tg) * w, self.dsigma(Tg) * w
        kg, dkg = self.k(Tg) * w, self.dk(Tg) * w
        GG = np.einsum("egai,egbi->egab", G, G)                 # grad N_a . grad N_b
        Gp = np.einsum("egai,egi->ega", G, gp)                  # grad N_a . grad psi
        GT = np.einsum("egai,egi->ega", G, gT)
        N = g.N
        K = np.zeros((len(Te), 16, 16))
        K[:, :8, :8] = s * np.einsum("eg,egab->eab", sg, GG)
        K[:, :8, 8:] = s * np.einsum("eg,ega,gb->eab", dsg, Gp, N)
        K[:, 8:, 8:] = (np.einsum("eg,egab->eab", kg, GG) + np.einsum("eg,ega,gb->eab", dkg, GT, N)
                        + self.Me / dt
                        - s * s * np.einsum("ga,eg,eg,gb->eab", N, dsg, np.sum(gp * gp, -1), N))
        K[:, 8:, :8] = -2 * s * s * np.einsum("ga,eg,egb->eab", N, sg, Gp)
        return K

    def jacobian(self, u, u_old, dt, n=None, cs_volume=False):
        """Exact dR_n/du_n (CSR). Real u only. Volume terms analytic (cs_volume=True: complex step,
        used by the tangent test); face terms by complex step over the local unknowns."""
        u = np.asarray(u, dtype=float)
        if cs_volume:
            To = np.asarray(u_old, dtype=float)[self.nn:2 * self.nn][self.mesh.hexes]
            Ke = self._local_jac(lambda Z: self._vol_local(Z, np.broadcast_to(To, Z.shape[:2] + (8,)),
                                                           np.full(Z.shape[0], dt)), u[self.vdofs])
        else:
            Ke = self.vol_tangent(u, dt)
        rows = [np.repeat(self.vdofs, 16, axis=1).ravel()]
        cols = [np.tile(self.vdofs, (1, 16)).ravel()]
        vals = [Ke.ravel()]
        for name, dofs in self.fdofs.items():
            Kf = self._local_jac(lambda Z, nm=name: self._face_local(nm, Z), u[dofs])
            a = dofs.shape[1]
            rows.append(np.repeat(dofs, a, axis=1).ravel()); cols.append(np.tile(dofs, (1, a)).ravel())
            vals.append(Kf.ravel())
        return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(self.ndof, self.ndof))

    def mass_matrix(self):
        """Heat-capacity block in the full state layout; dR_n/du_{n-1} = -mass_matrix()/dt_n."""
        if self._M is None:
            h = self.mesh.hexes + self.nn
            self._M = sp.csr_matrix((np.real(self.Me).ravel(), (np.repeat(h, 8, 1).ravel(),
                                                                np.tile(h, (1, 8)).ravel())),
                                    shape=(self.ndof, self.ndof))
        return self._M

    # ------------------------------------------------------------------ diagnostics
    def volume(self):
        return np.sum(self.geo.wdetJ)

    def split(self, u):
        u = np.asarray(u)
        return u[..., :self.nn], u[..., self.nn:2 * self.nn], u[..., self.iV]

    def resistance(self, u, n):
        """Terminal-to-terminal resistance V_in / I (Ohm)."""
        return PHI_SCALE * self.split(u)[2] / self.sched.I[n]

    def heat_flows(self, u):
        """Half-model power balance (W): Joule in the copper, joint heat, heat leaving each face set."""
        p, g, s = self.p, self.geo, PHI_SCALE
        ps, T, V = self.split(np.asarray(u, float))
        Te, pe = T[self.mesh.hexes], ps[self.mesh.hexes]
        Tg = Te @ g.N.T
        gp = np.einsum("egai,ea->egi", g.dNdx, pe)
        out = {"joule_bulk": float(np.sum(s * s * self.sigma(Tg) * np.sum(gp * gp, -1) * g.wdetJ))}
        jh = 0.0
        for name, fg in self.fgeo.items():
            q = self.mesh.faces[name]
            Tq = T[q] @ fg.N.T
            if name in EXPOSED:
                out["out_" + name] = float(np.sum(self.exposed_flux(Tq) * fg.wdA))
            else:
                Vt = V if name == "contact_in" else 0.0
                j = s * (ps[q] @ fg.N.T - Vt) / p.rho_c
                jh += float(np.sum(p.rho_c * j * j * fg.wdA))
                out["out_" + name] = float(np.sum(p.h_term * (Tq - p.T_term) * fg.wdA))
        out["joule_joints_total"] = jh
        out["joule_joints_into_bar"] = p.joint_heat_share * jh
        return out
