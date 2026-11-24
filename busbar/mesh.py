"""Structured hex mesh of half of a straight busbar with variable width and thickness.

Coordinates: x along the bar (0..L), y across the width (mirror plane y = 0, the bar is
y in [0, w(x)/2]), z through the thickness (0..t(x); the joint faces are on z = 0).
w(x) and t(x) are piecewise-linear between 12 equally spaced stations. Only node coordinates
depend on the design; the topology is fixed, so shape gradients are smooth.

Face sets: 'contact_in' (z=0, x<=l_c), 'contact_out' (z=0, x>=L-l_c), 'bottom' (rest of z=0),
'top' (z=t), 'side' (y=w/2), 'ends' (x=0 and x=L), 'sym_y' (y=0).
With p.contacts_at_ends (verification) the joints are the end faces instead.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from brakeopt.fem import HEX_FACES
from .params import Params, N_STATIONS, station_x, split_design


def hat_matrix(x_nodes, x_st):
    Phi = np.zeros((len(x_nodes), len(x_st)))
    for s in range(len(x_st)):
        e = np.zeros(len(x_st)); e[s] = 1.0
        Phi[:, s] = np.interp(x_nodes, x_st, e)
    return Phi


@dataclass
class Mesh:
    X: np.ndarray
    hexes: np.ndarray
    faces: dict
    node_ijk: np.ndarray
    x_nodes: np.ndarray
    multiplicity: int = 2          # half model

    @property
    def n_nodes(self):
        return self.X.shape[0]

    @property
    def n_elem(self):
        return self.hexes.shape[0]

    def face_nodes(self, name):
        return np.unique(self.faces[name])


def x_grid(p: Params):
    a = np.linspace(0.0, p.l_c, p.nx_c + 1)
    b = np.linspace(p.l_c, p.L - p.l_c, p.nx_m + 1)
    c = np.linspace(p.L - p.l_c, p.L, p.nx_c + 1)
    return np.r_[a, b[1:], c[1:]]


def build_mesh(d, p: Params) -> Mesh:
    d = np.asarray(d)
    w_st, t_st = split_design(d)
    xs = x_grid(p)
    nx, ny, nz = len(xs) - 1, p.ny, p.nz
    Phi = hat_matrix(xs, station_x(p))
    w = Phi @ w_st
    t = Phi @ t_st
    I, J, K = np.meshgrid(np.arange(nx + 1), np.arange(ny + 1), np.arange(nz + 1), indexing="ij")
    I, J, K = I.ravel(), J.ravel(), K.ravel()
    nid = lambda i, j, k: (i * (ny + 1) + j) * (nz + 1) + k
    dtype = np.result_type(d, float)
    X = np.empty((I.size, 3), dtype=dtype)
    X[:, 0] = xs[I]
    X[:, 1] = (J / ny) * 0.5 * w[I]
    X[:, 2] = (K / nz) * t[I]
    ei, ej, ek = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    ei, ej, ek = ei.ravel(), ej.ravel(), ek.ravel()
    hexes = np.stack([nid(ei, ej, ek), nid(ei + 1, ej, ek), nid(ei + 1, ej + 1, ek), nid(ei, ej + 1, ek),
                      nid(ei, ej, ek + 1), nid(ei + 1, ej, ek + 1), nid(ei + 1, ej + 1, ek + 1),
                      nid(ei, ej + 1, ek + 1)], axis=1)

    def fset(mask, local):
        return hexes[mask][:, HEX_FACES[local]]

    xc = 0.5 * (xs[ei] + xs[ei + 1])
    bot = ek == 0
    if p.contacts_at_ends:
        faces = {"contact_in": fset(ei == 0, "xi-"), "contact_out": fset(ei == nx - 1, "xi+"),
                 "bottom": fset(bot, "zeta-"), "ends": np.zeros((0, 4), int)}
    else:
        cin = bot & (xc < p.l_c)
        cout = bot & (xc > p.L - p.l_c)
        faces = {"contact_in": fset(cin, "zeta-"), "contact_out": fset(cout, "zeta-"),
                 "bottom": fset(bot & ~cin & ~cout, "zeta-"),
                 "ends": np.vstack([fset(ei == 0, "xi-"), fset(ei == nx - 1, "xi+")])}
    faces["top"] = fset(ek == nz - 1, "zeta+")
    faces["side"] = fset(ej == ny - 1, "eta+")
    faces["sym_y"] = fset(ej == 0, "eta-")
    return Mesh(X=X, hexes=hexes, faces=faces, node_ijk=np.stack([I, J, K], 1), x_nodes=xs)
