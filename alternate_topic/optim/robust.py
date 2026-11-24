"""Design under uncertainty (bonus): several sampled braking scenarios.

Each scenario s (initial speed, front brake bias) has its own transient solve and its
own adjoint sweep; the robust design must satisfy every constraint in every scenario,
g_{k,s}(d) <= 0, i.e. the worst case over the sampled set. Gradients are exact because
each g_{k,s} is an ordinary adjoint-differentiated history functional.
"""
from __future__ import annotations

import numpy as np

from .design_problem import DesignProblem

SCENARIOS = [  # (label, change in initial speed km/h, change in front bias) vs Params()
    ("nominal", 0.0, 0.0),
    ("+20 km/h initial speed", 20.0, 0.0),
    ("+0.05 front brake bias", 0.0, 0.05),
]


class RobustProblem:
    def __init__(self, scenarios=SCENARIOS, radiation=True, kT=True):
        from brakeopt.params import Params
        from .brake_adapter import BrakeSolver
        from .run_ws23 import constraint_spec
        self.scenarios = scenarios
        self.subs = []
        self.scenario_params = []
        for label, dv, dbias in scenarios:
            p0 = Params()
            p = p0.replace(v_high=p0.v_high + dv / 3.6, front_bias=p0.front_bias + dbias)
            self.scenario_params.append({"label": label, "v_high_kmh": p.v_high * 3.6,
                                         "front_bias": p.front_bias})
            cons = [dict(c, name=f"{c['name']} | {label}") for c in constraint_spec(p)]
            self.subs.append(DesignProblem(BrakeSolver(radiation, kT, params=p), cons))
        self.base = self.subs[0]
        self.d0, self.m0, self.solver = self.base.d0, self.base.m0, self.base.solver
        self.cons = [c for s in self.subs for c in s.cons]

    @property
    def n_forward(self):
        return sum(s.n_forward for s in self.subs)

    @property
    def n_adjoint(self):
        return sum(s.n_adjoint for s in self.subs)

    def d_of(self, x):
        return self.base.d_of(x)

    def evaluate(self, x, need_grad=True):
        es = [s.evaluate(x, need_grad) for s in self.subs]
        out = dict(es[0])
        out["g"] = np.concatenate([e["g"] for e in es])
        out["Tmax"] = np.concatenate([e["Tmax"] for e in es])
        if need_grad:
            out["grad_g"] = np.vstack([e["grad_g"] for e in es])
        return out

    f = DesignProblem.f
    df = DesignProblem.df
    c_ineq = DesignProblem.c_ineq
    dc_ineq = DesignProblem.dc_ineq
    bounds = DesignProblem.bounds
