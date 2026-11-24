"""Structured hex mesh of one quarter vane-pitch sector of a vented brake disc.

Coordinates: cylindrical (r, theta, z) mapped to Cartesian (x, y, z) = (r cos theta, r sin theta, z).
  theta in [0, pi/n_vanes]   theta = 0 is the mid-plane of a vane (mirror symmetry)
  z     in [0, h_v/2 + t_c]  z = 0 is the mid-plane of the disc between the two cheeks (mirror symmetry)
Solid regions:
  vane : 0 <= theta <= theta_v(r) = t_v(r) / (2 r),  0 <= z <= h_v/2
  cheek: all theta,                                   h_v/2 <= z <= h_v/2 + t_c(r)
The structured index grid (i, j, k) has nodes along r, theta, z. Elements with j in the channel
columns and k in the vane layers are air and are not created. The topology never changes with the
design; only node coordinates move (shape parameterization).
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .fem import HEX_FACES
from .params import Params, N_STATIONS, station_radii, split_design


def hat_matrix(r_nodes: np.ndarray, r_st: np.ndarray) -> np.ndarray:
    """Piecewise-linear interpolation matrix Phi (len(r_nodes), len(r_st)): f(r) = Phi @ f_st."""
    Phi = np.zeros((len(r_nodes), len(r_st)))
    for s in range(len(r_st)):
        e = np.zeros(len(r_st))
        e[s] = 1.0
        Phi[:, s] = np.interp(r_nodes, r_st, e)
    return Phi


@dataclass
class Mesh:
    X: np.ndarray            # (n_nodes, 3) Cartesian coordinates, m (complex if d is complex)
    hexes: np.ndarray        # (n_elem, 8)
    faces: dict              # name -> (n_face, 4) node ids, outward cyclic order
    node_rtz: np.ndarray     # (n_nodes, 3) (r, theta, z)
    node_ijk: np.ndarray     # (n_nodes, 3) structured indices
    elem_ijk: np.ndarray     # (n_elem, 3)
    dX_dd: np.ndarray        # (n_nodes, 3, 24) exact derivative of X w.r.t. d
    sector_multiplicity: int
    n_vanes: int

    @property
    def n_nodes(self):
        return self.X.shape[0]

    @property
    def n_elem(self):
        return self.hexes.shape[0]

    def face_nodes(self, name):
        return np.unique(self.faces[name])


def radial_nodes(p: Params) -> np.ndarray:
    a = np.linspace(p.r_i, p.r_pad_in, p.nr_inner + 1)
    b = np.linspace(p.r_pad_in, p.r_o, p.nr_outer + 1)
    return np.r_[a, b[1:]]


def _topology(p: Params):
    nr = p.nr_inner + p.nr_outer
    ntv, ntc, nzv, nzc = p.nth_vane, p.nth_chan, p.nz_vane, p.nz_cheek
    nt, nz = ntv + ntc, nzv + nzc

    def solid(j, k):
        return (j < ntv) or (k >= nzv)

    elems = [(i, j, k) for i in range(nr) for j in range(nt) for k in range(nz) if solid(j, k)]
    elem_ijk = np.array(elems, dtype=int)
    # local node offsets in VTK order (xi -> i, eta -> j, zeta -> k)
    off = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                    [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]])
    corner = elem_ijk[:, None, :] + off[None, :, :]                   # (e, 8, 3)
    flat = (corner[..., 0] * (nt + 1) + corner[..., 1]) * (nz + 1) + corner[..., 2]
    used, inv = np.unique(flat.ravel(), return_inverse=True)
    hexes = inv.reshape(flat.shape)
    node_ijk = np.stack([used // ((nt + 1) * (nz + 1)), (used // (nz + 1)) % (nt + 1), used % (nz + 1)], 1)

    def local_faces(mask, key):
        f = np.array(HEX_FACES[key])
        return hexes[mask][:, f]

    ei, ej, ek = elem_ijk.T
    faces = {
        "friction": local_faces((ek == nz - 1) & (ei >= p.nr_inner), "zeta+"),
        "face_free": local_faces((ek == nz - 1) & (ei < p.nr_inner), "zeta+"),
        "channel": np.vstack([local_faces((ej == ntv - 1) & (ek < nzv), "eta+"),
                              local_faces((ek == nzv) & (ej >= ntv), "zeta-")]),
        "rim": local_faces(ei == nr - 1, "xi+"),
        "hub": local_faces(ei == 0, "xi-"),
        # symmetry planes (adiabatic; used for display and verification only)
        "sym_theta0": local_faces(ej == 0, "eta-"),
        "sym_thetamax": local_faces((ej == nt - 1) & (ek >= nzv), "eta+"),
        "sym_z0": local_faces((ek == 0) & (ej < ntv), "zeta-"),
    }
    return hexes, faces, node_ijk, elem_ijk


def build_mesh(d, p: Params) -> Mesh:
    """Mesh for design d (24,), metres. d may be complex (complex-step derivatives)."""
    d = np.asarray(d)
    if d.shape != (2 * N_STATIONS,):
        raise ValueError(f"design vector must have shape (24,), got {d.shape}")
    hexes, faces, node_ijk, elem_ijk = _topology(p)
    ntv, ntc, nzv, nzc = p.nth_vane, p.nth_chan, p.nz_vane, p.nz_cheek
    r_nodes = radial_nodes(p)
    Phi = hat_matrix(r_nodes, station_radii(p))
    dv, dc = split_design(d)
    tv = Phi @ dv                      # vane thickness at each radial node line
    tc = Phi @ dc                      # cheek thickness
    i, j, k = node_ijk.T
    r = r_nodes[i]
    th_max = np.pi / p.n_vanes
    th_v = tv[i] / (2.0 * r)
    if np.any(np.real(th_v) >= th_max):
        raise ValueError("vane thicker than the vane pitch")
    # theta of each node and its derivative w.r.t. theta_v
    in_vane = j <= ntv
    s_v = j / ntv
    s_c = (j - ntv) / ntc
    theta = np.where(in_vane, s_v * th_v, th_v + s_c * (th_max - th_v))
    dth_dthv = np.where(in_vane, s_v, 1.0 - s_c)
    # z of each node and its derivative w.r.t. t_c
    in_vz = k <= nzv
    z = np.where(in_vz, (k / nzv) * (p.h_v / 2), p.h_v / 2 + ((k - nzv) / nzc) * tc[i])
    dz_dtc = np.where(in_vz, 0.0, (k - nzv) / nzc)

    X = np.stack([r * np.cos(theta), r * np.sin(theta), z], axis=1)
    # exact derivative dX/dd (evaluated at the real part)
    thr = np.real(theta)
    dX = np.zeros((len(r), 3, 2 * N_STATIONS))
    dth_ddv = (dth_dthv / (2.0 * r))[:, None] * Phi[i]                  # (n, 12)
    dX[:, 0, :N_STATIONS] = (-r * np.sin(thr))[:, None] * dth_ddv
    dX[:, 1, :N_STATIONS] = (r * np.cos(thr))[:, None] * dth_ddv
    dX[:, 2, N_STATIONS:] = dz_dtc[:, None] * Phi[i]
    rtz = np.stack([r + 0 * theta, theta, z], axis=1)
    return Mesh(X=X, hexes=hexes, faces=faces, node_rtz=rtz, node_ijk=node_ijk, elem_ijk=elem_ijk,
                dX_dd=dX, sector_multiplicity=4 * p.n_vanes, n_vanes=p.n_vanes)


def full_disc_copies(X: np.ndarray, n_vanes: int):
    """Generator of the 4*n_vanes rigid copies of the quarter sector that make the full disc (display)."""
    for m in range(n_vanes):
        for sz in (1, -1):
            for sth in (1, -1):
                Y = X.copy()
                Y[:, 1] *= sth
                Y[:, 2] *= sz
                a = 2 * np.pi * m / n_vanes
                c, s = np.cos(a), np.sin(a)
                yield np.stack([c * Y[:, 0] - s * Y[:, 1], s * Y[:, 0] + c * Y[:, 1], Y[:, 2]], 1)
