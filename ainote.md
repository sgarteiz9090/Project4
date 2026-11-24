# ainote.md: shared status board

Every workstream (WS) keeps its own section below up to date: **Done / In progress / Next / Blocked on**.
Read the whole file before you start a work session, and update your section before you push.
The integrator (WS1) merges each workstream branch into `main` once it merges cleanly and
`python -m pytest -q` passes on the merged result. When you finish a piece, add a 2-3 sentence **plain-language note** to your section (the guide uses it) and record any bug the verification caught (AI-use statement). Do not edit another workstream's section,
except to answer a "Blocked on" item (sign it, e.g. "WS1: answered, see Interfaces").

Deadline: GitHub Pages site final on **Nov 25, 11:59 PM America/Phoenix** (assumed).

## Workstreams and branches

| WS | Scope | Branch (write yours here) | Owner thread |
|----|-------|---------------------------|--------------|
| WS1 | Repo scaffold, mesh, transient nonlinear solver, verification ladder, integration/merging | `claude/project-thread-3ve605` | Build brake disc heat solver |
| WS2+3 | Discrete adjoint (backward in time), KS aggregation, FD gradient checks, Taylor test, cost comparison; SLSQP optimization, >=2 starts, history, multipliers, with vs without nonlinearity, optional uncertainty bonus | `gradients-optimizer` | adjoint + optimizer thread |
| WS4 | Website: GitHub Pages report (10 sections), three.js viewer (8 requirements) | _fill in_ | |
| WS5 | GUIDE.md beginner guide (each WS also adds notes to its own GUIDE section) | `claude/project-thread-c068hs` | guide thread |

Rules for every branch:
- Branch off the latest `main`; merge `main` into your branch (do not rebase shared history) before asking for integration.
- `python -m pytest -q` must pass. Slow tests get `@pytest.mark.slow` (skipped by default, run with `-m slow`).
- Only `scripts/run_all.py` writes `docs/assets/results.json`. Never hand-edit numbers on the site.
- Units: SI inside the code (m, s, K, W). Temperatures are **Kelvin internally**; convert to Celsius only for display/JSON fields that say `_C`.

---

## Interfaces (contracts between workstreams)

### Design vector `d` (24 entries, metres)
```
d[0:12]  = full vane thickness t_v at 12 radial stations, inner -> outer (m), bounds [2 mm, 8 mm]
d[12:24] = friction-cheek thickness t_c at 12 radial stations, inner -> outer (m), bounds [5 mm, 12 mm]
```
Stations are equally spaced from inner radius `r_i` to outer radius `r_o` (see `brakeopt/params.py`,
`station_radii(p)`). Between stations the thickness is piecewise-linear in r.
`brakeopt.params.design_bounds(p)` returns `(lower, upper)` arrays; `initial_design(p)` returns a baseline.

### Geometry / mesh (`brakeopt/mesh.py`)
One vane pitch of a vented disc with `n_vanes` vanes, reduced by two mirror symmetries to a
**quarter sector**: theta in [0, pi/n_vanes] (theta=0 is the vane mid-plane), z in [0, h_v/2 + t_c(r)]
(z=0 is the disc mid-plane). Vane solid: theta*r <= t_v(r)/2 and z <= h_v/2. Cheek solid: h_v/2 <= z <= h_v/2 + t_c(r).
Air channel elements are simply not created. 8-node trilinear hexes, true cylindrical geometry.

```python
mesh = build_mesh(d, p)          # Mesh object
mesh.X          # (n_nodes, 3) node coordinates, m (Cartesian x, y, z)
mesh.hexes      # (n_elem, 8) int, node ids, standard VTK/Abaqus C3D8 ordering
mesh.faces[name]  # (n_face, 4) int, boundary quads, outward-oriented, for each set:
                  # 'friction' (top face z=h_v/2+t_c, r >= r_pad_in), 'face_free' (top face, r < r_pad_in),
                  # 'channel' (all faces wetted by vane-channel air), 'rim' (r=r_o), 'hub' (r=r_i)
mesh.node_rtz   # (n_nodes, 3) (r, theta, z) of each node
mesh.dX_dd      # (n_nodes, 3, 24) exact derivative of node coordinates w.r.t. design variables
mesh.sector_multiplicity  # = 4 * n_vanes  (copies of the quarter sector in the full disc)
```
The mesh topology (connectivity, face sets) is identical for every `d`; only `X` changes. That is what
makes shape gradients smooth: dR/dd = dR/dX . dX/dd.

### Physics model (`brakeopt/model.py`)
Implicit Euler in time, Newton in each step. Residual of step n (state `T` = nodal temperatures in K):
```
R_n(T_n, T_{n-1}, X) = M(X) (T_n - T_{n-1})/dt_n + K(T_n, X) T_n
                     + Hconv_n(X) (T_n - T_air)  + Rad(T_n, X)  - Q_n(X) = 0
```
- `K(T)` uses k(T) of grey cast iron evaluated at Gauss points (nonlinear).
- `Rad` = eps*sigma*(T^4 - T_amb^4) on friction, face_free and rim faces (nonlinear).
- `Hconv_n`: convection coefficients depend on time (vehicle speed) but not on T or d.
- `Q_n`: frictional heat flux on the friction faces during braking (energy from deceleration).
```python
model = ThermalModel(mesh, p)            # precomputes geometry (B matrices, detJ, face Jacobians)
R = model.residual(T, T_old, dt, n)      # (n_nodes,)   n = step index into the load schedule
J = model.jacobian(T, dt, n)             # scipy.sparse CSR, dR/dT_n  (exact Newton tangent)
Mdt = model.mass_matrix() / dt           # dR_n/dT_{n-1} = -M/dt  (needed by the adjoint)
```
`model.residual` is written with plain numpy ops so it also accepts **complex** `X`/`T` (complex step
works for partial derivatives dR/dX, which WS2 can use for dR/dd).

### Transient solve (`brakeopt/solver.py`)
```python
sched = brake_schedule(p)                  # Schedule: times t[0..N], dt[1..N], speed v(t), heat flux, h(t)
sol = solve_transient(model, sched)        # Solution
sol.T        # (N+1, n_nodes) temperature history in K, sol.T[0] = initial temperature
sol.t        # (N+1,) times in s
sol.newton   # list of per-step lists of residual norms (for the quadratic-convergence plot)
```

### Functional facade for WS2/WS3 (`brakeopt/api.py`), agreed with the adjoint+optimizer thread
All functions take the design vector `d` (24,) in metres. `radiation` and `kT` switch the two nonlinear
terms off (`kT=False` freezes k at k(T_ref)); use the same flags in every call of one gradient evaluation.
```python
from brakeopt import api
fw = api.forward(d, radiation=True, kT=True)   # Forward object
fw.T              # (N+1, n_nodes) K, fw.T[0] = initial state
fw.t, fw.dt       # (N+1,) times; dt[n] = t[n]-t[n-1] for n>=1 (dt[0] = nan)
fw.friction_nodes # node ids on the rubbing faces (constraint 1)
fw.hub_nodes      # node ids on the hub/bell mounting face r = r_i (constraint 2)
fw.mesh, fw.model, fw.sched, fw.newton   # full objects if needed
api.residual(d, T_new, T_old, n, radiation=True, kT=True)   # R_n (n_nodes,), complex-step safe in d
api.jac_T(d, T_new, T_old, n, ...)      # sparse dR_n/dT_new (exact Newton tangent)
api.jac_Told(d, T_new, T_old, n, ...)   # sparse dR_n/dT_old = -M/dt_n
api.vjp_d(d, T_new, T_old, n, lam, ...) # lam^T dR_n/dd  (24,)
api.vjp_d_all(d, T_hist, lam_hist, ...) # sum over n=1..N of lam_n^T dR_n/dd (24,), vectorized
api.mass(d), api.dmass_dd(d)            # kg, kg/m (24,)
api.x0(), api.bounds()                  # baseline design, (lower, upper)
```
Step index `n` runs 1..N and refers to the step from `t[n-1]` to `t[n]`.
WS2+3 owns KS aggregation of the two temperature constraints and results.json keys
`"sensitivities"` and `"optimization"` (and `designs.optimized`, `designs.optimized_linear`).

### Functionals (`brakeopt/functionals.py`)
```python
mass(mesh, p)                  # full-disc rotor mass, kg (objective)
g_face(model, sol, p)          # KS-aggregated peak friction-surface temperature over space AND time, K
g_hub(model, sol, p)           # KS-aggregated peak hub/mounting temperature (r = r_i face), K
constraints(model, sol, p)     # dict name -> (value, limit), constraint form value <= limit
```
KS: `KS(x) = x_max + (1/rho) * log( sum_w w_i exp(rho (x_i - x_max)) )` over nodal values of the face
set at every time step (w = 1, rho = `p.ks_rho`, in 1/K). Exact peaks are also reported for comparison.
The functionals are written so WS2 can get `dg/dT_n` (partials w.r.t. each step's state) and `dg/dX`.

### `docs/assets/results.json` schema (written only by `scripts/run_all.py`)
```jsonc
{
  "meta": {"generated": "ISO time", "git_commit": "...", "command": "python scripts/run_all.py"},
  "params": { ... all physical parameters with units ... },
  "designs": {
    "initial":   DESIGN,
    "optimized": DESIGN,             // WS3
    "optimized_linear": DESIGN       // WS3: optimum without the nonlinear terms (N2 comparison)
  },
  "verification": [ {"id": "V1", "name": "...", "ours": 1.0, "exact": 1.0, "rel_error": 1e-9,
                     "pass": true, "note": "..."} ],
  "convergence": {"mms_space": {"h": [...], "err_L2": [...], "slope": 2.0},
                  "mms_time": {"dt": [...], "err": [...], "slope": 1.0},
                  "newton": {"residuals": [...]}},
  "sensitivities": { ... WS2+3: check table, taylor test, cost comparison ... },
  "optimization": { ... WS2+3: per-start history, multipliers, active set ... }
}
DESIGN = {
  "d": {"vane_thickness_mm": [12], "cheek_thickness_mm": [12], "station_radius_mm": [12]},
  "mass_kg": 0.0,
  "constraints": [ {"name": "peak_face_T", "value_C": 0.0, "limit_C": 700.0, "pass": true} , ...],
  "mesh": {"nodes": [[x,y,z], ...] /* m, quarter sector */, "hexes": [[8 ints], ...],
           "n_vanes": 36, "symmetry": "quarter sector: mirror z=0 and theta=0, then rotate n_vanes times"},
  "temperature": {"unit": "C", "times": [...], "frames": [[per-node T], ...]   /* <= 80 frames */},
  "history": {"times": [...], "peak_face_C": [...], "peak_hub_C": [...], "speed_kmh": [...]}
}
```
Add new keys freely; do not rename or remove existing ones without a note here.

---

## WS1: Solver, verification, integration (branch `claude/project-thread-3ve605`)
**Done**
- Repo scaffold: package layout, requirements, pytest setup, this file, `brakeopt/params.py` (all physical parameters + design vector helpers).

**In progress**
- Mesh generator, ThermalModel (k(T) + radiation + convection + braking flux), Newton/implicit Euler.

**Next**
- Verification ladder V1-V7 (see plan below), `scripts/run_all.py` writing initial design to results.json.

**Blocked on**
- nothing

Verification ladder plan (WS1):
| ID | Test | Exact answer |
|----|------|--------------|
| V1 | Steady 1D slab, fixed T both ends, constant k | linear profile, exact to round-off |
| V2 | Steady 1D slab with convection (Robin) end | Analytical T(x) and heat flux |
| V3 | MMS steady 3D on distorted sector mesh with k(T) | O(h^2) L2 error slope |
| V4 | MMS transient (implicit Euler) | O(dt) error slope |
| V5 | Global energy balance over braking cycle | heat in = stored + convected + radiated |
| V6 | Nonlinear benchmark: 1D slab with k(T)=k0(1+bT) (Kirchhoff transform exact solution) and radiation-cooled lump | exact T profile; Newton quadratic rate |
| V7 | Invariance: rotate/translate mesh, mesh-independence of peak T | identical results |
| V8 | (WS2+3) adjoint vs central FD, every constraint, both variable types; Taylor test | |

## WS2+3: Adjoint gradients + optimization (branch `gradients-optimizer`)
_section owned by the adjoint + optimizer thread_

## WS4: Website
_not started_

## WS5: Guide (GUIDE.md)

- **Done:** First draft of GUIDE.md from the plan: what is being built, each piece in plain language (heat PDE, finite elements, time stepping, Newton/nonlinearity, KS aggregation, verification ladder, adjoint, SLSQP, KKT multipliers, results.json + website), glossary, rubric map, AI-use statement skeleton.
- **In progress:** Nothing.
- **Next:** Update each "Status" box and fill the *TBD* numbers as the other workstreams report here (solver/ladder results, gradient checks, optimum and multipliers, viewer). Fill the AI-use table and "bug the verification caught".
- **Blocked on:** Other workstreams' plain-language notes and results.
- **Needs from others:** (1) Final mesh sizes and braking scenario numbers. (2) Ladder table with errors. (3) One bug the ladder or gradient check caught. (4) Active constraints and multiplier values at the optimum. (5) With vs without nonlinearity comparison.
