import numpy as np
from brakeopt.params import Params, initial_design, design_bounds, N_DESIGN


def test_design_vector_shape_and_bounds():
    p = Params()
    d0 = initial_design(p)
    lo, hi = design_bounds(p)
    assert d0.shape == (N_DESIGN,) == lo.shape == hi.shape
    assert np.all(lo <= d0) and np.all(d0 <= hi)
