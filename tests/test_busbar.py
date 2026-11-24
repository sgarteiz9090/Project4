"""Fallback topic (busbar): fast checks of the coupled electro-thermal model."""
import numpy as np
import pytest

from busbar import ladder
from busbar.mesh import build_mesh
from busbar.model import BusbarModel
from busbar.params import Params, initial_design


def test_mesh_volume_and_faces():
    p = Params()
    m = build_mesh(initial_design(p), p)
    from busbar.fem import HexGeometry
    vol = np.sum(HexGeometry(m.X, m.hexes).wdetJ) * 2
    assert vol == pytest.approx(p.L * p.w0 * p.t0, rel=1e-12)
    assert ladder.v12_faces()["pass"]


def test_exact_resistance_and_joule_profile():
    assert ladder.v1_resistance()["pass"]
    assert ladder.v3_joule_1d()["pass"]


def test_tangent_and_energy():
    for r in ladder.v11_tangent() + ladder.v8_energy():
        assert r["pass"], r


@pytest.mark.slow
def test_full_ladder():
    rows, _ = ladder.run_ladder(include_slow=False)
    assert all(r["pass"] for r in rows), [r for r in rows if not r["pass"]]


@pytest.mark.slow
def test_adjoint_vs_fd():
    from busbar.gradcheck import fd_table
    from busbar.problem import BusbarProblem, BusbarSolver
    pr = BusbarProblem(BusbarSolver(True))
    rows = fd_table(pr, np.ones(24), [0, 13])
    assert max(r["err_vs_grad_norm"] for r in rows) < 1e-6


def test_busbar_site_keys_and_results():
    import json, pathlib, re
    site = pathlib.Path(__file__).resolve().parents[1] / "docs"
    md = (site / "index.md").read_text()
    js = (site / "assets" / "js" / "report.js").read_text()
    for kind in ("table", "plot", "calc"):
        for key in set(re.findall(rf'data-{kind}="([a-z_A-Z0-9]+)"', md)):
            assert re.search(rf"\n\s+{key}[(:]", js), f"data-{kind}={key} used in index.md but not defined in report.js"
    p = site / "assets" / "results.json"
    if not p.exists():
        pytest.skip("busbar results.json not generated yet")
    r = json.loads(p.read_text())
    for key, d in r["designs"].items():
        m = d["mesh"]
        n = len(m["nodes"])
        assert all(len(h) == 8 and max(h) < n for h in m["hexes"]), key
        fr = d["temperature"]["frames"]
        assert len(fr) == len(d["temperature"]["times"]) and all(len(f) == n for f in fr), key
        for c in d["constraints"]:
            assert "pass" in c and "name" in c, c
