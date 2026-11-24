"""Finite element building blocks: 8-node trilinear hexahedra and 4-node bilinear boundary quads.

Everything here is written with plain numpy arithmetic so that node coordinates may be complex
(complex-step derivatives with respect to the design pass straight through).
"""
from __future__ import annotations

import numpy as np

# Local node coordinates of the 8-node hex (VTK / Abaqus C3D8 ordering)
HEX_NODES = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                      [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], dtype=float)

# Local faces of the hex, each listed counter-clockwise when seen from outside (outward normal)
HEX_FACES = {
    "xi-": [0, 4, 7, 3], "xi+": [1, 2, 6, 5],
    "eta-": [0, 1, 5, 4], "eta+": [3, 7, 6, 2],
    "zeta-": [0, 3, 2, 1], "zeta+": [4, 5, 6, 7],
}

QUAD_NODES = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], dtype=float)


def gauss_1d(n: int):
    x, w = np.polynomial.legendre.leggauss(n)
    return x, w


def hex_shape(xi: np.ndarray):
    """Shape functions and local derivatives at points xi (g, 3) -> N (g, 8), dN (g, 8, 3)."""
    xi = np.atleast_2d(xi)
    a = HEX_NODES
    f = 1.0 + xi[:, None, :] * a[None, :, :]          # (g, 8, 3)
    N = 0.125 * f[..., 0] * f[..., 1] * f[..., 2]
    dN = np.empty(xi.shape[:1] + (8, 3))
    dN[..., 0] = 0.125 * a[:, 0] * f[..., 1] * f[..., 2]
    dN[..., 1] = 0.125 * a[:, 1] * f[..., 0] * f[..., 2]
    dN[..., 2] = 0.125 * a[:, 2] * f[..., 0] * f[..., 1]
    return N, dN


def quad_shape(s: np.ndarray):
    """Bilinear quad: N (g, 4), dN (g, 4, 2)."""
    s = np.atleast_2d(s)
    a = QUAD_NODES
    f = 1.0 + s[:, None, :] * a[None, :, :]
    N = 0.25 * f[..., 0] * f[..., 1]
    dN = np.empty(s.shape[:1] + (4, 2))
    dN[..., 0] = 0.25 * a[:, 0] * f[..., 1]
    dN[..., 1] = 0.25 * a[:, 1] * f[..., 0]
    return N, dN


def hex_rule(n: int = 2):
    x, w = gauss_1d(n)
    pts = np.array([[a, b, c] for c in x for b in x for a in x])
    wts = np.array([wa * wb * wc for wc in w for wb in w for wa in w])
    return pts, wts


def quad_rule(n: int = 2):
    x, w = gauss_1d(n)
    pts = np.array([[a, b] for b in x for a in x])
    wts = np.array([wa * wb for wb in w for wa in w])
    return pts, wts


def det3(J):
    """Determinant of (..., 3, 3), complex safe."""
    return (J[..., 0, 0] * (J[..., 1, 1] * J[..., 2, 2] - J[..., 1, 2] * J[..., 2, 1])
            - J[..., 0, 1] * (J[..., 1, 0] * J[..., 2, 2] - J[..., 1, 2] * J[..., 2, 0])
            + J[..., 0, 2] * (J[..., 1, 0] * J[..., 2, 1] - J[..., 1, 1] * J[..., 2, 0]))


def inv3(J, detJ=None):
    """Inverse of (..., 3, 3) via the adjugate, complex safe."""
    if detJ is None:
        detJ = det3(J)
    C = np.empty_like(J)
    C[..., 0, 0] = J[..., 1, 1] * J[..., 2, 2] - J[..., 1, 2] * J[..., 2, 1]
    C[..., 0, 1] = J[..., 0, 2] * J[..., 2, 1] - J[..., 0, 1] * J[..., 2, 2]
    C[..., 0, 2] = J[..., 0, 1] * J[..., 1, 2] - J[..., 0, 2] * J[..., 1, 1]
    C[..., 1, 0] = J[..., 1, 2] * J[..., 2, 0] - J[..., 1, 0] * J[..., 2, 2]
    C[..., 1, 1] = J[..., 0, 0] * J[..., 2, 2] - J[..., 0, 2] * J[..., 2, 0]
    C[..., 1, 2] = J[..., 0, 2] * J[..., 1, 0] - J[..., 0, 0] * J[..., 1, 2]
    C[..., 2, 0] = J[..., 1, 0] * J[..., 2, 1] - J[..., 1, 1] * J[..., 2, 0]
    C[..., 2, 1] = J[..., 0, 1] * J[..., 2, 0] - J[..., 0, 0] * J[..., 2, 1]
    C[..., 2, 2] = J[..., 0, 0] * J[..., 1, 1] - J[..., 0, 1] * J[..., 1, 0]
    return C / detJ[..., None, None]


class HexGeometry:
    """Gauss-point geometry of a hex mesh: shape functions, physical gradients, weights."""

    def __init__(self, X, hexes, order: int = 2):
        pts, wts = hex_rule(order)
        self.N, dNdxi = hex_shape(pts)                 # (g, 8), (g, 8, 3)
        Xe = X[hexes]                                   # (e, 8, 3)
        J = np.einsum("gai,eaj->egij", dNdxi, Xe)       # J_ij = dx_j / dxi_i
        detJ = det3(J)
        if np.any(np.real(detJ) <= 0):
            raise ValueError("inverted or degenerate hex element (detJ <= 0)")
        invJ = inv3(J, detJ)
        self.dNdx = np.einsum("egij,gaj->egai", invJ, dNdxi)   # (e, g, 8, 3)
        self.wdetJ = detJ * wts[None, :]                         # (e, g)
        self.xg = np.einsum("ga,eai->egi", self.N, Xe)           # Gauss points (e, g, 3)


class QuadGeometry:
    """Gauss-point geometry of boundary quads (nodes in cyclic order)."""

    def __init__(self, X, quads, order: int = 2):
        pts, wts = quad_rule(order)
        self.N, dN = quad_shape(pts)                    # (g, 4), (g, 4, 2)
        Xf = X[quads]                                   # (f, 4, 3)
        ts = np.einsum("ga,fai->fgi", dN[..., 0], Xf)
        tt = np.einsum("ga,fai->fgi", dN[..., 1], Xf)
        n = np.cross(ts, tt)                            # area-weighted normal (f, g, 3)
        dA = np.sqrt(np.sum(n * n, axis=-1))            # complex-safe norm
        self.wdA = dA * wts[None, :]                    # (f, g)
        self.wdA_z = n[..., 2] * wts[None, :]           # projected area on the z = const plane (signed)
        self.normal = n / dA[..., None]                 # unit outward normal (f, g, 3)
        self.xg = np.einsum("ga,fai->fgi", self.N, Xf)  # (f, g, 3)
