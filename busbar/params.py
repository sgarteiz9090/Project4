"""Physical, geometric and numerical parameters of the busbar model (SI units, K inside).

Sources (cited in the report):
- Annealed copper (C11000 / Cu-ETP): resistivity 1.7241e-8 Ohm m at 20 C (100 % IACS) and
  temperature coefficient 0.00393 /K (IEC 60028, International Standard of Resistance for Copper);
  k ~ 398 W/mK at 300 K falling ~2 % per 100 K; rho = 8960 kg/m^3; c = 385 J/kgK
  (Incropera, Fundamentals of Heat and Mass Transfer, App. A, Table A.1).
- Bolted copper joints: contact resistance of the order of 5-20 micro-ohm for a 25 x 30 mm lap
  joint, i.e. specific contact resistance ~1e-9 to 1e-8 Ohm m^2 (Braunovic, Konchits, Myshkin,
  Electrical Contacts: Fundamentals, Applications and Technology, CRC 2007, ch. 2-3).
- Simplified natural convection for air, laminar: h = 1.42 (dT/L)^(1/4) W/m^2K
  (Holman, Heat Transfer, 10th ed., Table 7-2).
- Powder-coated / epoxy-insulated busbar emissivity ~0.85-0.9 (Incropera Table A.11, paints).
- Temperature limits: 105 C insulation class of the busbar sleeve / heat-shrink (UL 224, 105 C
  rating); 70 C at the joint face, which is bolted to the 40 C module terminal: with the ~10 K drop
  of the terminal path this keeps the cell tabs within the 55-60 C maximum operating temperature of
  typical Li-ion cell data sheets. Resistance budget 52 micro-ohm (about 8 W per busbar at 400 A).
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import numpy as np

N_STATIONS = 12
N_DESIGN = 2 * N_STATIONS
SIGMA_SB = 5.670374419e-8
PHI_SCALE = 1e-3          # state psi is the potential in millivolts: phi = PHI_SCALE * psi


@dataclass
class Params:
    # ---------------- geometry ----------------
    L: float = 0.300            # busbar length, m (module-to-module link)
    l_c: float = 0.025          # lap-joint overlap length at each end, m
    w_min: float = 0.012        # design bounds: width, m
    w_max: float = 0.050
    t_min: float = 0.0015       # thickness, m
    t_max: float = 0.008
    w0: float = 0.030           # baseline (uniform) 30 x 5 mm bar
    t0: float = 0.005

    # ---------------- copper ----------------
    rho: float = 8960.0
    c: float = 385.0
    sigma0: float = 5.80e7      # S/m at T_ref (100 % IACS)
    alpha: float = 3.93e-3      # 1/K, sigma(T) = sigma0 / (1 + alpha (T - T_ref))
    k0: float = 398.0           # W/mK at T_ref
    k_beta: float = -2.0e-4     # 1/K, k(T) = k0 (1 + k_beta (T - T_ref))
    T_ref: float = 293.15
    emissivity: float = 0.85

    # ---------------- joints / terminals ----------------
    rho_c: float = 7.5e-9       # specific contact resistance, Ohm m^2 (10 micro-ohm on 25 x 30 mm)
    joint_heat_share: float = 0.5   # share of the contact Joule heat entering the busbar
    h_term: float = 250.0       # W/m^2K, conductance from the joint face into the cell terminal
    T_term: float = 313.15      # K, battery terminal / cell temperature (40 C)

    # ---------------- environment ----------------
    T_amb: float = 313.15       # K, pack air (40 C)
    nc_coef: float = 1.42       # h = nc_coef (|dT| / nc_length)^(1/4)
    nc_length: float = 0.05     # m, characteristic length (bar mounted on edge)
    nc_delta: float = 0.05      # K, smoothing of |dT| at dT = 0 (keeps h differentiable)

    # ---------------- load (current) profile ----------------
    I_base: float = 400.0       # A, continuous current
    I_peak: float = 1200.0       # A, boost pulse (full-power acceleration / fast-charge peak)
    t_pulse: float = 60.0       # s
    t_after: float = 240.0      # s, recovery at I_base after the pulse
    dt_pulse: float = 2.0
    dt_after: float = 6.0
    dt_steady: float = 1.0e9    # step 1: (pseudo) steady state at I_base

    # ---------------- constraints ----------------
    T_bar_max: float = 378.15   # K (105 C) insulation rating, peak over the whole history
    T_joint_max: float = 343.15 # K (70 C) at the joint faces (cell terminal), whole history
    R_max: float = 52.0e-6      # Ohm, terminal-to-terminal resistance at the continuous rating
    ks_rho: float = 1.0         # 1/K, KS aggregation parameter

    # ---------------- numerics ----------------
    nx_c: int = 3               # elements along each joint overlap
    nx_m: int = 16              # elements along the free span
    ny: int = 3                 # elements across the half width
    nz: int = 2                 # elements through the thickness
    newton_tol: float = 1e-11
    newton_maxit: int = 30

    # ---------------- N2 switches ----------------
    nonlinear: bool = True      # False: constant-property ("hand calculation") model:
    T_lin: float = 293.15       # K, sigma and k frozen at data-sheet values (20 C)
    dT_lin_h: float = 30.0      # K, natural-convection h frozen at h(30 K)
    T_lin_rad: float = 343.15   # K, radiation linearized about 70 C
    contacts_at_ends: bool = False   # verification only: joints on the end faces x=0, x=L

    def to_dict(self):
        return asdict(self)

    def replace(self, **kw):
        q = Params(**asdict(self))
        for k, v in kw.items():
            if not hasattr(q, k):
                raise AttributeError(k)
            setattr(q, k, v)
        return q

    def refined(self, f: int):
        return self.replace(nx_c=self.nx_c * f, nx_m=self.nx_m * f, ny=self.ny * f, nz=self.nz * f)


def station_x(p: Params):
    return np.linspace(0.0, p.L, N_STATIONS)


def design_bounds(p: Params):
    lo = np.r_[np.full(N_STATIONS, p.w_min), np.full(N_STATIONS, p.t_min)]
    hi = np.r_[np.full(N_STATIONS, p.w_max), np.full(N_STATIONS, p.t_max)]
    return lo, hi


def initial_design(p: Params):
    return np.r_[np.full(N_STATIONS, p.w0), np.full(N_STATIONS, p.t0)]


def split_design(d):
    d = np.asarray(d)
    return d[:N_STATIONS], d[N_STATIONS:]
