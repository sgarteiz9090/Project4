"""Aggregated temperature constraints over space and time (KS function).

For a set of nodes S and all time steps n = 1..N, with a_{i,n} = T_{i,n} / T_limit
(temperatures in the same units as the limit, here degrees C),

    KS(a) = a_max + (1/rho) ln sum_{i in S, n} exp(rho (a_{i,n} - a_max)),
    a_max <= KS <= a_max + ln(|S| N) / rho,

and the scaled constraint is g = KS - 1 <= 0. KS is a smooth upper bound on the true
maximum over the braking history; the report shows the true maximum next to it.
Its partial derivative w.r.t. each a_{i,n} is the softmax weight
w_{i,n} = exp(rho (a_{i,n} - KS)), which sums to 1.
"""
from __future__ import annotations

import numpy as np


def ks_history(T_hist, nodes, T_limit, rho=100.0, offset=0.0, steps=None):
    """KS aggregate of (T - offset)/T_limit over `nodes` and time steps.

    T_hist : (N+1, ndof) temperatures; offset converts to the limit's units
             (offset=273.15 if T_hist is in kelvin and the limit is in deg C).
    steps  : indices of time levels to include (default 1..N).
    Returns (value, dKS_dT) with dKS_dT of shape (N+1, ndof).
    """
    T_hist = np.asarray(T_hist)
    N = T_hist.shape[0] - 1
    steps = np.arange(1, N + 1) if steps is None else np.asarray(steps)
    nodes = np.asarray(nodes)
    a = (T_hist[np.ix_(steps, nodes)] - offset) / T_limit
    amax = a.max()
    s = np.exp(rho * (a - amax)).sum()
    ks = amax + np.log(s) / rho
    w = np.exp(rho * (a - ks))  # softmax weights, sum to 1
    dT = np.zeros_like(T_hist, dtype=float)
    dT[np.ix_(steps, nodes)] = w / T_limit
    return ks, dT


def true_max(T_hist, nodes, offset=0.0, steps=None):
    T_hist = np.asarray(T_hist)
    N = T_hist.shape[0] - 1
    steps = np.arange(1, N + 1) if steps is None else np.asarray(steps)
    sub = T_hist[np.ix_(steps, np.asarray(nodes))] - offset
    k = np.unravel_index(np.argmax(sub), sub.shape)
    return float(sub[k]), int(steps[k[0]]), int(np.asarray(nodes)[k[1]])
