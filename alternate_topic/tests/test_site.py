"""Checks that the report site and results.json fit together (no browser needed)."""
import json
import pathlib
import re

DOCS = pathlib.Path(__file__).resolve().parents[1] / "docs"


def _results():
    for name in ("results.json", "mock_results.json"):
        p = DOCS / "assets" / name
        if p.exists():
            return json.loads(p.read_text())
    raise AssertionError("no results file in docs/assets")


def test_every_table_and_plot_in_index_is_defined():
    md = (DOCS / "index.md").read_text()
    js = (DOCS / "assets" / "js" / "report.js").read_text()
    for kind in ("table", "plot", "calc"):
        for key in set(re.findall(rf'data-{kind}="([a-z_]+)"', md)):
            assert re.search(rf"\n\s+{key}[(:]", js), f"data-{kind}={key} used in index.md but not defined in report.js"


def test_results_have_what_the_viewer_needs():
    r = _results()
    for key, d in r["designs"].items():
        m = d["mesh"]
        n = len(m["nodes"])
        assert all(len(h) == 8 and max(h) < n for h in m["hexes"]), key
        assert m["n_vanes"] > 0
        fr = d["temperature"]["frames"]
        assert len(fr) == len(d["temperature"]["times"]) and all(len(f) == n for f in fr), key
        assert isinstance(d["mass_kg"], (int, float))
        for c in d["constraints"]:
            assert "pass" in c and "name" in c, c


def test_mock_results_are_flagged():
    p = DOCS / "assets" / "mock_results.json"
    if p.exists():
        assert json.loads(p.read_text())["meta"]["mock"] is True
