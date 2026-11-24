"""Transient nonlinear heat conduction in the disc sector: residual and exact Newton tangent.

Weak form of  rho c dT/dt = div(k(T) grad T)  with boundary conditions
  friction faces : -k dT/dn = h_face (T - T_amb) + eps sigma (T^4 - T_amb^4) - q(r, t)
  free face, rim : -k dT/dn = h_face (T - T_amb) + eps sigma (T^4 - T_amb^4)
  vane channels  : -k dT/dn = h_chan (T - T_amb)
  hub (r = r_i)  : -k dT/dn = h_hub (T - T_bell)
  symmetry planes theta = 0, theta = pi/n_vanes, z = 0 : -k dT/dn = 0  (natural)
Implicit Euler for step n (t[n-1] -> t[n]):
  R_n = M (T_n - T_{n-1}) / dt_n + K(T_n) T_n + H_n (T_n - T_amb) + Rad(T_n) + Hub (T_n - T_bell) - Q_n
with the frictional flux q(r, t) = P_sector(t) r / int(r dA) (uniform pad pressure: q proportional
to sliding speed, i.e. to r), normalised on the discrete face so the input energy is exact.

All residual code is complex-step safe in the node coordinates. The residual accepts a batch of
states (B, n_nodes) with step indices n (B,) so that WS2 can evaluate all steps at once.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from .fem import HexGeometry, QuadGeometry
from .loads import Schedule, brake_schedule
from .mesh import Mesh
from .params import Params, SIGMA_SB

RAD_SETS = ("friction", "face_free", "rim")
FACE_CONV_SETS = ("friction", "face_free", "rim")


class ThermalModel:
    def __init__(self, mesh: Mesh, p: Params, sched: Schedule | None = None):
        self.mesh, self.p = mesh, p
        self.sched = sched if sched is not None else brake_schedule(p)
        X, hexes = mesh.X, mesh.hexes
        self.nn, self.ne = mesh.n_nodes, mesh.n_elem
        self.geo = HexGeometry(X, hexes)
        g = self.geo
        self.Me = p.rho * p.c * np.einsum("eg,ga,gb->eab", g.wdetJ, g.N, g.N)
        self.S_hex = self._scatter(hexes)
        self._hex_rows = np.repeat(hexes, 8, axis=1).ravel()
        self._hex_cols = np.tile(hexes, (1, 8)).ravel()
        self.fgeo, self.S_face, self._frows, self._fcols = {}, {}, {}, {}
        for name, quads in mesh.faces.items():
            self.fgeo[name] = QuadGeometry(X, quads)
            self.S_face[name] = self._scatter(quads)
            self._frows[name] = np.repeat(quads, 4, axis=1).ravel()
            self._fcols[name] = np.tile(quads, (1, 4)).ravel()
        # frictional flux shape: q = P_sector * r / C on projected area of the friction faces
        fg = self.fgeo["friction"]
        rq = np.sqrt(fg.xg[..., 0] ** 2 + fg.xg[..., 1] ** 2)
        C = np.sum(fg.wdA_z * rq)
        self.Fq = self.S_face["friction"] @ np.einsum("fg,ga->fa", fg.wdA_z * rq / C, fg.N).ravel()
        self.sector_share = 1.0 / (2 * p.n_vanes)     # quarter sector covers 1/(2 n_vanes) of a face
        self._M = None

    # ------------------------------------------------------------------ helpers
    def _scatter(self, conn):
        m, a = conn.shape
        return sp.csr_matrix((np.ones(m * a), (conn.ravel(), np.arange(m * a))), shape=(self.nn, m * a))

    def k_of_T(self, T):
        p = self.p
        if p.kT:
            return p.k0 * (1.0 + p.k_beta * (T - p.T_kref)), p.k0 * p.k_beta + 0.0 * T
        kc = p.k0 * (1.0 + p.k_beta * (p.T_init - p.T_kref))
        return kc + 0.0 * T, 0.0 * T

    def rad_flux(self, T):
        """Radiative flux (W/m^2) and its derivative w.r.t. T."""
        p = self.p
        es = p.emissivity * SIGMA_SB
        if p.radiation:
            return es * (T ** 4 - p.T_amb ** 4), 4 * es * T ** 3
        Ti = p.T_init
        return es * (Ti ** 4 - p.T_amb ** 4 + 4 * Ti ** 3 * (T - Ti)), 4 * es * Ti ** 3 + 0.0 * T

    def face_coeffs(self, n):
        """Convection coefficient per face set for step index/indices n."""
        s = self.sched
        n = np.atleast_1d(n)
        h = {name: s.h_face[n] for name in FACE_CONV_SETS}
        h["channel"] = s.h_chan[n]
        h["hub"] = np.full(n.shape, self.p.h_hub)
        return h

    def sink(self, name):
        return self.p.T_bell if name == "hub" else self.p.T_amb

    # ------------------------------------------------------------------ residual
    def residual(self, T, T_old, dt, n):
        """R_n(T_n, T_{n-1}); T, T_old (n_nodes,) or (B, n_nodes); dt, n scalars or (B,)."""
        T = np.asarray(T)
        single = T.ndim == 1
        T2 = np.atleast_2d(T)
        To2 = np.broadcast_to(np.atleast_2d(T_old), T2.shape)
        B = T2.shape[0]
        n = np.broadcast_to(np.atleast_1d(n), (B,))
        dt = np.broadcast_to(np.atleast_1d(np.asarray(dt, dtype=float)), (B,))
        g, hexes = self.geo, self.mesh.hexes
        Te = T2[:, hexes]
        Tg = np.einsum("ga,Zea->Zeg", g.N, Te)
        grad = np.einsum("egai,Zea->Zegi", g.dNdx, Te)
        k, _ = self.k_of_T(Tg)
        qg = (k * g.wdetJ[None])[..., None] * grad
        re = np.einsum("egai,Zegi->Zea", g.dNdx, qg)
        re = re + np.einsum("eab,Zeb->Zea", self.Me, Te - To2[:, hexes]) / dt[:, None, None]
        R = (self.S_hex @ re.reshape(B, -1).T).T
        h = self.face_coeffs(n)
        for name, fg in self.fgeo.items():
            quads = self.mesh.faces[name]
            Tq = np.einsum("ga,Zfa->Zfg", fg.N, T2[:, quads])
            val = h[name][:, None, None] * (Tq - self.sink(name))
            if name in RAD_SETS:
                val = val + self.rad_flux(Tq)[0]
            rf = np.einsum("Zfg,fg,ga->Zfa", val, fg.wdA, fg.N)
            R = R + (self.S_face[name] @ rf.reshape(B, -1).T).T
        P = self.sched.P_face[n] * self.sector_share
        R = R - P[:, None] * self.Fq[None, :]
        return R[0] if single else R

    # ------------------------------------------------------------------ tangent
    def jacobian(self, T, dt, n):
        """Exact dR_n/dT_n (scipy CSR). Real T only."""
        g, hexes = self.geo, self.mesh.hexes
        Te = T[hexes]
        Tg = Te @ g.N.T
        grad = np.einsum("egai,ea->egi", g.dNdx, Te)
        k, dk = self.k_of_T(Tg)
        Ke = np.einsum("egai,egbi->eab", (k * g.wdetJ)[..., None, None] * g.dNdx, g.dNdx)
        Ke += np.einsum("egai,egi,gb->eab", (dk * g.wdetJ)[..., None, None] * g.dNdx, grad, g.N)
        Ke += self.Me / dt
        rows, cols, vals = [self._hex_rows], [self._hex_cols], [np.real(Ke).ravel()]
        h = self.face_coeffs(n)
        for name, fg in self.fgeo.items():
            quads = self.mesh.faces[name]
            Tq = T[quads] @ fg.N.T
            c = h[name][0] + 0.0 * Tq
            if name in RAD_SETS:
                c = c + self.rad_flux(Tq)[1]
            Kf = np.einsum("fg,ga,gb->fab", c * fg.wdA, fg.N, fg.N)
            rows.append(self._frows[name]); cols.append(self._fcols[name]); vals.append(np.real(Kf).ravel())
        return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(self.nn, self.nn))

    def mass_matrix(self):
        """Consistent heat-capacity matrix M (CSR); dR_n/dT_{n-1} = -M / dt_n."""
        if self._M is None:
            self._M = sp.csr_matrix((np.real(self.Me).ravel(), (self._hex_rows, self._hex_cols)),
                                    shape=(self.nn, self.nn))
        return self._M

    # ------------------------------------------------------------------ diagnostics
    def heat_flows(self, T, n):
        """Sector heat flows at state T for step n, W: dict with 'in' and each outflow."""
        out = {"in": float(np.real(self.sched.P_face[n] * self.sector_share * self.Fq.sum()))}
        h = self.face_coeffs(n)
        for name, fg in self.fgeo.items():
            Tq = T[self.mesh.faces[name]] @ fg.N.T
            out["conv_" + name] = float(np.sum(h[name][0] * (Tq - self.sink(name)) * fg.wdA))
            if name in RAD_SETS:
                out["rad_" + name] = float(np.sum(self.rad_flux(Tq)[0] * fg.wdA))
        return out

    def volume(self):
        return np.sum(self.geo.wdetJ)
