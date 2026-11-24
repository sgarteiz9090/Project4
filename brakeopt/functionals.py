"""Objective and constraint functionals.

Mass is the objective. The two temperature constraints are the peak friction-surface temperature
and the peak hub/mounting-face temperature over the whole braking sequence. WS2+3 owns the smooth
(KS) aggregation used in the optimizer; `ks_max` here is a reference implementation and the exact
peaks are what the report shows next to the KS values.
"""
from __future__ import annotations

import numpy as np

from .mesh import Mesh
from .params import Params


def mass(mesh: Mesh, p: Params, model=None):
    """Full-disc rotor mass, kg (quarter-sector volume x 4 n_vanes x rho). Complex safe."""
    if model is not None:
        vol = model.volume()
    else:
        from .fem import HexGeometry
        vol = np.sum(HexGeometry(mesh.X, mesh.hexes).wdetJ)
    return p.rho * vol * mesh.sector_multiplicity


def ks_max(x, rho):
    """Kreisselmeier-Steinhauser smooth maximum (upper bound on max(x), within log(len(x))/rho)."""
    x = np.ravel(x)
    m = np.max(np.real(x))
    return m + np.log(np.sum(np.exp(rho * (x - m)))) / rho


def friction_nodes(mesh: Mesh):
    return mesh.face_nodes("friction")


def hub_nodes(mesh: Mesh):
    return mesh.face_nodes("hub")


def peak_temperatures(mesh: Mesh, T_hist):
    """Exact peaks (K) over space and time: (peak friction-face T, peak hub T)."""
    return float(np.max(T_hist[:, friction_nodes(mesh)])), float(np.max(T_hist[:, hub_nodes(mesh)]))


def constraints(mesh: Mesh, sol, p: Params):
    """name -> (value K, limit K), constraint form value <= limit (exact peaks)."""
    tf, th = peak_temperatures(mesh, sol.T)
    return {"peak_face_T": (tf, p.T_face_max), "peak_hub_T": (th, p.T_hub_max)}
