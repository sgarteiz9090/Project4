"""Physical, geometric and numerical parameters of the brake-disc model.

All values are SI (m, s, kg, K, W). Every number that ends up in the report should come from here
so it is written once and copied into results.json["params"] by scripts/run_all.py.

Sources (cited in the report, section 2/3):
- Grey cast iron properties: rho ~ 7100-7200 kg/m^3, c ~ 460-550 J/kgK, k ~ 50 W/mK at room
  temperature falling to ~35 W/mK near 700 C (ASM Handbook Vol. 1, gray iron; Limpert, *Brake Design
  and Safety*, 3rd ed., SAE 2011, ch. 3).
- Peak rubbing-surface temperature limit ~700 C for grey-iron race rotors (thermal cracking, pad fade);
  Limpert ch. 3, Day, *Braking of Road Vehicles* (2014) ch. 3.
- Rotating-disc convection: Dorfman (1963) / Limpert ch. 3, average Nu = 0.015 Re_w^0.8 (turbulent).
- Vane-channel convection: Limpert ch. 3, Nu = 0.023 (1 + (d_h/l)^0.67) Re^0.8 Pr^(1/3).
- Emissivity of oxidised cast iron ~0.55-0.65 (Incropera, *Fundamentals of Heat and Mass Transfer*, App. A).
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
import numpy as np

N_STATIONS = 12          # radial stations per variable type
N_DESIGN = 2 * N_STATIONS


@dataclass
class Params:
    # ---------------- geometry (full disc) ----------------
    r_i: float = 0.110          # inner (mounting / bell) radius, m
    r_o: float = 0.178          # outer radius, m (356 mm race rotor)
    r_pad_in: float = 0.128     # inner radius of the pad rubbing band, m
    h_v: float = 0.012          # vane (channel) height, m (axial gap between the two cheeks)
    n_vanes: int = 48           # number of straight radial vanes

    # ---------------- design variable bounds ----------------
    tv_min: float = 0.002       # min vane thickness, m (casting limit)
    tv_max: float = 0.008       # max vane thickness, m (channel must stay open)
    tc_min: float = 0.005       # min cheek thickness, m (pad wear allowance + strength)
    tc_max: float = 0.012       # max cheek thickness, m
    tv0: float = 0.005          # baseline vane thickness, m
    tc0: float = 0.009          # baseline cheek thickness, m

    # ---------------- material: grey cast iron ----------------
    rho: float = 7100.0         # kg/m^3
    c: float = 500.0            # J/(kg K)
    k0: float = 54.0            # W/(m K) at T_kref
    k_beta: float = -5.0e-4     # 1/K, k(T) = k0 (1 + beta (T - T_kref))
    T_kref: float = 293.15      # K
    emissivity: float = 0.60

    # ---------------- vehicle / braking ----------------
    vehicle_mass: float = 1300.0    # kg (car + driver)
    front_bias: float = 0.65        # fraction of braking energy on the front axle
    partition: float = 0.90         # fraction of friction heat entering the disc (rest to pad)
    tyre_radius: float = 0.33       # m
    v_high: float = 250.0 / 3.6     # m/s, speed at start of each stop
    v_low: float = 80.0 / 3.6       # m/s, speed at end of each stop
    decel: float = 1.5 * 9.81       # m/s^2
    t_accel: float = 12.0           # s, re-acceleration (cool-down) between stops
    n_stops: int = 3
    t_cool_end: float = 0.0         # s, extra cool-down at constant v_low after the last stop
    dt_brake: float = 0.1           # s, time step while braking
    dt_accel: float = 0.5           # s, time step while accelerating / cooling

    # ---------------- environment / convection ----------------
    T_amb: float = 308.15           # K, air and radiation surroundings (35 C)
    T_init: float = 473.15          # K, disc temperature at the start of the sequence (200 C)
    T_bell: float = 308.15          # K, sink temperature behind the hub/bell joint
    h_hub: float = 1000.0           # W/(m^2 K), conductance of the bell/hub joint at r = r_i
    h_min: float = 10.0             # W/(m^2 K), natural-convection floor at zero speed
    k_air: float = 0.040            # W/(m K) at ~400-500 K film temperature
    nu_air: float = 3.2e-5          # m^2/s
    Pr_air: float = 0.69
    c_pump: float = 0.30            # channel air speed = c_pump * omega * r_mean (vane pumping)

    # ---------------- constraints ----------------
    T_face_max: float = 973.15      # K (700 C)
    T_hub_max: float = 573.15       # K (300 C) bell / bearing protection
    ks_rho: float = 0.1             # 1/K, KS aggregation parameter (owned by WS2+3)

    # ---------------- numerics ----------------
    nr_inner: int = 2       # radial elements between r_i and r_pad_in
    nr_outer: int = 8       # radial elements between r_pad_in and r_o
    nth_vane: int = 2       # angular elements across the half vane
    nth_chan: int = 2       # angular elements across the half channel
    nz_vane: int = 2        # axial elements over the half vane height
    nz_cheek: int = 2       # axial elements through the cheek
    newton_tol: float = 1e-10   # relative residual tolerance
    newton_maxit: int = 25

    # ---------------- model switches (N2 nonlinearity) ----------------
    radiation: bool = True   # False: radiation linearised about T_init (h_rad = 4 eps sigma T_init^3)
    kT: bool = True          # False: k frozen at k(T_init)

    def to_dict(self) -> dict:
        return asdict(self)

    def refined(self, factor: int) -> "Params":
        """Copy with every element count multiplied by `factor` (mesh-convergence studies)."""
        q = Params(**asdict(self))
        for name in ("nr_inner", "nr_outer", "nth_vane", "nth_chan", "nz_vane", "nz_cheek"):
            setattr(q, name, getattr(self, name) * factor)
        return q

    def replace(self, **kw) -> "Params":
        q = Params(**asdict(self))
        for key, val in kw.items():
            if not hasattr(q, key):
                raise AttributeError(key)
            setattr(q, key, val)
        return q


SIGMA_SB = 5.670374419e-8   # W/(m^2 K^4)


def station_radii(p: Params) -> np.ndarray:
    """Radii of the 12 design stations, equally spaced from r_i to r_o."""
    return np.linspace(p.r_i, p.r_o, N_STATIONS)


def design_bounds(p: Params) -> tuple[np.ndarray, np.ndarray]:
    lo = np.r_[np.full(N_STATIONS, p.tv_min), np.full(N_STATIONS, p.tc_min)]
    hi = np.r_[np.full(N_STATIONS, p.tv_max), np.full(N_STATIONS, p.tc_max)]
    return lo, hi


def initial_design(p: Params) -> np.ndarray:
    return np.r_[np.full(N_STATIONS, p.tv0), np.full(N_STATIONS, p.tc0)]


def split_design(d):
    """(vane thickness at stations, cheek thickness at stations)."""
    d = np.asarray(d)
    return d[:N_STATIONS], d[N_STATIONS:]
