"""Adjoint machinery on the 1D stand-in model (fast; independent of the 3D solver)."""
import numpy as np

from optim.design_problem import DesignProblem
from optim.gradcheck import fd_table, taylor_test
from optim.toy_problem import ToyRod

CONS = [{"name": "friction", "nodes_attr": "friction_nodes", "limit_C": 400, "rho": 100},
        {"name": "hub", "nodes_attr": "hub_nodes", "limit_C": 150, "rho": 100}]


def _prob(**kw):
    rod = ToyRod(n_steps=30, **kw)
    rod.hub_nodes = [3, 4, 5]      # close enough to the heated end to warm up within 15 s
    return DesignProblem(rod, CONS)


def test_adjoint_matches_fd_nonlinear():
    p = _prob()
    x = np.linspace(0.8, 1.3, 24)
    rows = fd_table(p, x, [0, 7, 15, 23])
    assert max(r["err_vs_grad_norm"] for r in rows) < 1e-5


def test_adjoint_matches_fd_linear():
    p = _prob(radiation=False, kT=False)
    rows = fd_table(p, np.ones(24), [1, 12])
    assert max(r["err_vs_grad_norm"] for r in rows) < 1e-5


def test_taylor_second_order():
    p = _prob()
    t = taylor_test(p, np.ones(24), eps_list=[1e-1, 1e-2])
    e1, e2 = np.array(t["central_err"])[:, 1]     # friction constraint
    assert e2 < e1 / 50                            # ~ eps^2: factor 100 per decade
