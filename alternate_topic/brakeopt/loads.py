"""Braking sequence: time grid, vehicle speed, frictional heat input and convection coefficients.

The car does `n_stops` stops from v_high to v_low at constant deceleration, re-accelerating back to
v_high (linearly in time) between stops. Brake power of the whole car while braking is M a v(t).
A fraction front_bias goes to the front axle, split over 2 discs, each with 2 rubbing faces, and a
fraction `partition` of it enters the disc (the rest goes into the pads).

Per step n (from t[n-1] to t[n]) the heat input uses the *step-averaged* power, so the discrete
energy input equals the kinetic-energy budget exactly for any time step. Convection coefficients
use the speed at the end of the step (implicit Euler).
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .params import Params


@dataclass
class Schedule:
    t: np.ndarray          # (N+1,) s
    dt: np.ndarray         # (N+1,) s, dt[0] = nan
    v: np.ndarray          # (N+1,) m/s, vehicle speed at t[n]
    braking: np.ndarray    # (N+1,) bool: step n (ending at t[n]) is a braking step
    P_face: np.ndarray     # (N+1,) W, step-averaged heat power into ONE rubbing face of ONE disc
    h_face: np.ndarray     # (N+1,) W/(m^2 K), friction faces, free face, rim
    h_chan: np.ndarray     # (N+1,) W/(m^2 K), vane channels

    @property
    def n_steps(self):
        return len(self.t) - 1


def segments(p: Params):
    """List of (duration, v_start, v_end, is_braking, dt)."""
    t_brake = (p.v_high - p.v_low) / p.decel
    seg = []
    for s in range(p.n_stops):
        seg.append((t_brake, p.v_high, p.v_low, True, p.dt_brake))
        if s < p.n_stops - 1:
            seg.append((p.t_accel, p.v_low, p.v_high, False, p.dt_accel))
    if p.t_cool_end > 0:
        seg.append((p.t_cool_end, p.v_low, p.v_low, False, p.dt_accel))
    return seg


def h_face_of_speed(v, p: Params):
    """Average rotating-disc heat transfer coefficient (turbulent), Nu = 0.015 Re_w^0.8."""
    omega = np.asarray(v) / p.tyre_radius
    Re = omega * p.r_o ** 2 / p.nu_air
    h = 0.015 * Re ** 0.8 * p.k_air / p.r_o
    return np.maximum(h, p.h_min)


def h_chan_of_speed(v, p: Params):
    """Vane-channel coefficient, Nu = 0.023 (1 + (d_h/l)^0.67) Re^0.8 Pr^(1/3), baseline channel."""
    omega = np.asarray(v) / p.tyre_radius
    r_mean = 0.5 * (p.r_i + p.r_o)
    u = p.c_pump * omega * r_mean
    w = 2 * np.pi * r_mean / p.n_vanes - p.tv0      # baseline channel width (h independent of d)
    dh = 2 * w * p.h_v / (w + p.h_v)
    L = p.r_o - p.r_i
    Re = u * dh / p.nu_air
    Nu = 0.023 * (1 + (dh / L) ** 0.67) * Re ** 0.8 * p.Pr_air ** (1 / 3)
    return np.maximum(Nu * p.k_air / dh, p.h_min)


def face_power_factor(p: Params) -> float:
    """Fraction of the car's braking power that enters one rubbing face of one front disc."""
    return p.partition * p.front_bias / 2.0 / 2.0


def brake_schedule(p: Params) -> Schedule:
    ts, vs, br, Pf = [0.0], [p.v_high], [False], [0.0]
    t0 = 0.0
    fac = face_power_factor(p)
    for dur, va, vb, is_brake, dt in segments(p):
        n = max(1, int(round(dur / dt)))
        tt = np.linspace(0.0, dur, n + 1)[1:]
        vv = va + (vb - va) * tt / dur
        vprev = np.r_[va, vv[:-1]]
        ts += list(t0 + tt)
        vs += list(vv)
        br += [is_brake] * n
        # step-averaged power M a (v_{n-1}+v_n)/2 (v is linear in t within the segment)
        Pf += list(fac * p.vehicle_mass * p.decel * 0.5 * (vprev + vv)) if is_brake else [0.0] * n
        t0 += dur
    t = np.array(ts)
    v = np.array(vs)
    dt = np.r_[np.nan, np.diff(t)]
    return Schedule(t=t, dt=dt, v=v, braking=np.array(br), P_face=np.array(Pf),
                    h_face=h_face_of_speed(v, p), h_chan=h_chan_of_speed(v, p))


def energy_into_one_face(p: Params) -> float:
    """Exact frictional energy entering one rubbing face over the whole sequence, J."""
    ke = 0.5 * p.vehicle_mass * (p.v_high ** 2 - p.v_low ** 2) * p.n_stops
    return face_power_factor(p) * ke
