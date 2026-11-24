"""The verification ladder must pass (fast rungs here; V9 mesh study is slow)."""
import pytest

from verification import ladder


def _check(rows):
    for r in rows:
        assert r["pass"], r


def test_v1_v2_slab_and_invariance():
    _check(ladder.v1_v2_slab())


def test_v3_kirchhoff_nonlinear_k():
    _check(ladder.v3_kirchhoff()[0])


def test_v4_radiation_plate():
    _check(ladder.v4_radiation()[0])


def test_v5_mms_space():
    _check(ladder.v5_mms_space(factors=(1, 2, 3))[0])


def test_v6_mms_time():
    _check(ladder.v6_mms_time()[0])


def test_v7_v8_energy_and_newton():
    _check(ladder.v7_v8_real()[0])


def test_v10_closure():
    _check(ladder.v10_closure())


@pytest.mark.slow
def test_v9_mesh_convergence():
    _check(ladder.v9_mesh_convergence()[0])
