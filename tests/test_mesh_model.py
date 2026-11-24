import numpy as np
import pytest

from brakeopt import api
from brakeopt.mesh import build_mesh
from brakeopt.model import ThermalModel
from brakeopt.params import Params, initial_design, design_bounds


@pytest.fixture(scope="module")
def setup():
    p = Params()
    d = initial_design(p)
    m = build_mesh(d, p)
    return p, d, m, ThermalModel(m, p)


def test_volume_matches_analytic(setup):
    p, d, m, model = setup
    th = np.pi / p.n_vanes
    v_cheek = 0.5 * th * (p.r_o ** 2 - p.r_i ** 2) * p.tc0
    v_vane = 0.5 * p.tv0 * (p.r_o - p.r_i) * 0.5 * p.h_v
    assert model.volume() == pytest.approx(v_cheek + v_vane, rel=2e-3)


def test_dX_dd_matches_finite_difference(setup):
    p, d, m, _ = setup
    rng = np.random.default_rng(1)
    v = rng.standard_normal(24) * 1e-3
    h = 1e-4
    fd = (build_mesh(d + h * v, p).X - build_mesh(d - h * v, p).X) / (2 * h)
    assert np.allclose(fd, m.dX_dd @ v, rtol=1e-6, atol=1e-12)


def test_extreme_designs_build(setup):
    p = setup[0]
    for d in design_bounds(p):
        ThermalModel(build_mesh(d, p), p)


@pytest.mark.parametrize("radiation,kT", [(True, True), (False, False)])
def test_jacobian_matches_finite_difference(setup, radiation, kT):
    p, d, m, _ = setup
    q = p.replace(radiation=radiation, kT=kT)
    model = ThermalModel(m, q)
    rng = np.random.default_rng(0)
    T = 600 + 200 * rng.random(model.nn)
    To = 500 + 100 * rng.random(model.nn)
    n, dt = 5, model.sched.dt[5]
    J = model.jacobian(T, dt, n)
    v = rng.standard_normal(model.nn)
    h = 1e-4
    fd = (model.residual(T + h * v, To, dt, n) - model.residual(T - h * v, To, dt, n)) / (2 * h)
    assert np.linalg.norm(J @ v - fd) <= 1e-7 * np.linalg.norm(fd)


def test_conduction_conserves_energy(setup):
    """Sum of the residual = net heat flow (conduction terms cancel exactly)."""
    p, d, m, model = setup
    rng = np.random.default_rng(2)
    T = 500 + 300 * rng.random(model.nn)
    n = 3
    R = model.residual(T, T, 1.0, n)
    f = model.heat_flows(T, n)
    out = sum(v for k, v in f.items() if k != "in")
    assert R.sum() == pytest.approx(out - f["in"], rel=1e-10)


def test_api_complex_step_partials():
    d = api.x0()
    fw = api.forward(d)
    n = 12
    lam = np.random.default_rng(3).standard_normal(fw.T.shape[1])
    g = api.vjp_d(d, fw.T[n], fw.T[n - 1], n, lam)
    for j in (0, 7, 13, 23):
        e = np.zeros(24)
        e[j] = 1e-7
        fd = lam @ (api.residual(d + e, fw.T[n], fw.T[n - 1], n) - api.residual(d - e, fw.T[n], fw.T[n - 1], n)) / 2e-7
        assert g[j] == pytest.approx(fd, rel=1e-5, abs=1e-6 * np.abs(g).max())
    gm = api.dmass_dd(d)
    for j in (2, 14):
        e = np.zeros(24)
        e[j] = 1e-7
        assert gm[j] == pytest.approx((api.mass(d + e) - api.mass(d - e)) / 2e-7, rel=1e-6)
    assert np.allclose(api.jac_Told(d, fw.T[n], fw.T[n - 1], n).toarray(),
                       -fw.model.mass_matrix().toarray() / fw.dt[n])
