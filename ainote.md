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
| WS2+3 | Discrete adjoint (backward in time), KS aggregation, FD gradient checks, Taylor test, cost comparison; SLSQP optimization, >=2 starts, history, multipliers, with vs without nonlinearity, optional uncertainty bonus | `claude/project-thread-1qs3ox` | adjoint + optimizer thread |
| WS4 | Website: GitHub Pages report (10 sections), three.js viewer (8 requirements) | `website` | Build the project report website |
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
- Repo scaffold: package layout, requirements, pytest setup, this file, `brakeopt/params.py`.
- `brakeopt/mesh.py` (quarter-sector hex mesh, exact dX/dd), `fem.py`, `loads.py` (braking schedule,
  step-averaged heat input, convection correlations), `model.py` (residual + exact tangent, k(T),
  radiation, batch- and complex-safe), `solver.py` (implicit Euler + Newton), `functionals.py`, `api.py`.
- Baseline (default mesh, 231 nodes, 120 hexes, 144 steps): forward solve ~1.7 s, mass 9.25 kg,
  peak friction face ~683 C, peak hub ~199 C (limits 700 C / 250 C). Newton converges quadratically
  (residual 1e2 -> 1e-1 -> 1e-7 -> 1e-13). `vjp_d_all` ~1.2 s.
- Scenario: 900 kg car, 60 % front, 3 stops 250->80 km/h at 1.5 g, 12 s re-acceleration between,
  disc starts at 150 C. All in `Params`.

**Plain-language note (for GUIDE.md)**
The solver chops a 1/192 slice of the disc (half of one vane and the cheek above it; the rest of the
disc is copies of this slice by symmetry) into small brick elements and tracks the temperature at
every brick corner, tenth of a second by tenth of a second, through three hard stops. At each time
step the equations are nonlinear (iron conducts worse when hot, and a hot disc glows away heat as
T^4), so Newton's method guesses, measures how wrong the guess is, and corrects, usually in 3-4 tries.

**Verification ladder results (default mesh; all PASS; `python -m verification.ladder`)**
| ID | Test | Result |
|----|------|--------|
| V1 | Steady slab, convection both ends, exact linear profile | max error 1e-16 of dT |
| V2 | Same slab, arbitrary 3D rotation + translation | identical to 3e-16 |
| V3a/b | Nonlinear k(T) slab vs Kirchhoff exact solution | flux rel. error 4e-7; L2 order 2.00 |
| V4a/b | Radiation-cooled plate (T^4) vs exact T(t) | rel. error 2e-6 (Richardson); order in dt 0.99 |
| V5 | MMS steady 3D on the curved sector mesh, k(T)+radiation+convection | L2 order 1.98 |
| V6 | MMS transient, all nonlinear terms | order in dt 1.00 |
| V7a/b | Heat input = kinetic-energy budget; stored + lost = input | exact; imbalance 2e-16 |
| V8 | Newton on the braking problem | quadratic (order ~2.2): 2e2, 6e-1, 2e-6, 3e-13 |
| V9 | Mesh convergence of peak face T (x1, x2, x3 meshes) | 683.2, 686.9, 687.6 C; coarse mesh 4.7 K (0.5 %) below the Richardson limit 687.8 C (optimizer mesh is slightly unconservative) |
| V10 | Every exterior face belongs to exactly one boundary set | 184 = 184 |

Bug caught by the ladder (for the AI-use statement): in V3 the agent's first version measured the
heat flux leaving the cold end with the wrong sign, so V3a failed (-1.29e6 vs +1.29e6 W/m^2); the
check against the exact Kirchhoff flux caught it. Separately, the first batched residual reused the
same einsum index letter for the time-step batch and the element node, which crashed on first run.

**Plain-language note: verification (for GUIDE.md)**
Before trusting the solver we fed it problems whose answers are known exactly: a plain wall with
heat flowing through it, the same wall turned in a random direction, a wall whose conductivity
changes with temperature, a plate cooling only by glowing, and "manufactured" answers where we pick
the temperature field first and work out what heating would produce it. The solver reproduces them,
and its error shrinks at the textbook rate as the bricks get smaller (4x less error per halving) and
as the time step gets smaller (2x less per halving). Energy is conserved exactly: all the braking
energy put in is accounted for as heat stored in the iron plus heat lost to air and radiation.

**Outputs**: `scripts/run_all.py` writes `docs/assets/results.json` (`verification`, `convergence`,
`designs.initial`, `figures`) and `docs/assets/fig_convergence.png`, `fig_newton.png`, `fig_history.png`.
Mesh JSON also includes `faces` (boundary sets, incl. `sym_*` planes) to help the viewer draw surfaces.

**In progress**
- Integration of branches (checking remote branches periodically).

**Next**
- Re-run `scripts/run_all.py` full once WS2+3 optimization is final.

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

## WS2+3: Adjoint gradients + optimization (branch `claude/project-thread-1qs3ox`)
_section owned by the adjoint + optimizer thread_

Code lives in `optim/` (does not touch `brakeopt/`):
- `optim/adjoint.py`: generic discrete adjoint, backward sweep over implicit-Euler steps, all functionals share one LU of A_n^T per step; uses `vjp_d_all` if present, else `vjp_d`, else complex step on `residual`.
- `optim/functionals.py`: KS over (node set x all time steps) of a = T_C / limit_C; rho_dimensionless = rho_per_K * limit_C (rho_per_K defaults to `p.ks_rho`). True maxima reported alongside.
- `optim/design_problem.py`: scaled problem x = d/d0, f = m/m0, g_k = KS_k - 1 <= 0, with caching.
- `optim/gradcheck.py`: FD table, Taylor test, adjoint-vs-FD cost. `optim/optimize.py`: SLSQP + history + KKT read-out (NNLS multipliers on the active set, plus SLSQP's own).
- `optim/brake_adapter.py`: binds `brakeopt.api` (radiation/kT flags fixed per instance).
- `optim/toy_problem.py`: 1D nonlinear transient stand-in with the same API, used only by `tests/test_adjoint_toy.py`.

**Done**
- Adjoint derivation + implementation; verified on the 1D stand-in (nonlinear k(T) + radiation): adjoint vs central FD agree to ~1e-10 relative, Taylor test falls as eps^2 over 5 decades.
- Bug the check caught (for the AI-use statement): the first Newton stopping test (residual relative to ||rho c T||) stopped after 1 iteration, so FD and adjoint disagreed at 4e-7 with a *flat* Taylor curve. Fix: stop on the Newton update, tol 1e-13. **WS1: please make sure `newton_tol` is tight enough for gradient checks (converge to ~1e-12 relative in the update), otherwise FD checks will flatten out.**

**In progress**
- Waiting for `brakeopt/api.py` to be real; then: gradient-check table (both variable types, both constraints, rad on/off), Taylor plot, cost comparison.

**Next**
- SLSQP from >= 2 starts (baseline, all-max), history per iteration, active set + multipliers in kg per % limit relaxation; linear (radiation=False, kT=False) optimum for the N2 comparison; optional uncertainty bonus (sampled initial speed / front bias, KS over scenarios).
- `optim/run_ws23.py` will return the `"sensitivities"`, `"optimization"`, `designs.optimized`, `designs.optimized_linear` blocks for `scripts/run_all.py` to merge in.

**Blocked on**
- `brakeopt/api.py` implementation (WS1).

**Plain-language note (for GUIDE.md)**
- The *adjoint* answers "if I make vane 5 a little thicker, how much hotter or cooler does the hottest point on the disc get during the whole braking sequence?" for all 24 dimensions at once. It does that by running the heat simulation once forward in time, then running one extra, cheaper "backward" calculation from the last stop to the first, which tracks how much each moment's temperatures influence the final peak. Checking against brute force (nudging each dimension and re-running the simulation) proves the numbers are right.

## WS4: Website (branch `website`)
**Done**
- `docs/_config.yml`, `docs/_layouts/default.html` (Jekyll 3.10, MathJax 3.2.2 and three.js 0.160.0 pinned), `docs/assets/css/site.css` (light/dark, phone width).
- `docs/index.md`: the report with all 10 handout sections + Reproducibility. **No number is typed in it**: `docs/assets/js/report.js` fills every value, table and plot from `results.json` at page load (spans `class="rj" data-k=/data-calc=/data-param=`, `<div class="rj-table" data-table=...>`, `<div class="rj-plot" data-plot=...>`). Missing keys show "pending", nothing breaks. So `python scripts/run_all.py` alone regenerates every number on the site; no extra build step.
- `docs/assets/js/viewer.js`: 3D viewer meeting the 8 requirements: initial vs optimized side by side with one shared camera; quarter sector mirrored (y, z) and rotated `n_vanes` times; draws `mesh.faces` sets except `sym*`; single-hue "Oranges" scale with colour bar and limit marks; over-limit faces hatched (shader), so not colour-only; faces subdivided with the bilinear face shape functions; play/scrub heating animation through the stops with a peak-T/speed chart; hover/tap tooltip (element, surface, radius, T now, peak, margin); metrics panel with pass/fail; dropdown to put `optimized_linear` on the right (N2 comparison).
- `docs/assets/js/data.js`: loads `assets/results.json`; falls back to `assets/mock_results.json` (red MOCK banner) if it is missing, or with `?mock` in the URL.
- `scripts/make_mock_results.py` writes `docs/assets/mock_results.json` (never `results.json`).
- `tests/test_site.py`: every table/plot/calc key used in index.md exists in report.js; results.json has what the viewer needs.
- Tested with Jekyll 3.10 + headless Chromium at 1280 px and 390 px against both the mock and WS1's real initial-design results.json: loads in about 1 s, no horizontal scroll, hover works, hatching shows on an over-limit design.

**Plain-language note (for GUIDE.md)**: The website is the report. It is a normal Markdown page, but when it opens in a browser a small script reads `results.json` and writes every number, table and plot into the page, so nothing can go stale. The 3D viewer takes the small slice of the disc the solver computed, mirrors and copies it into the whole disc, colours it by temperature, and lets you play the three stops to see where the iron gets hot and how much hotter the lighter disc runs.

**In progress**
- nothing (waiting for real optimized results to check the page against them).

**Next**
- Re-check the live page once `designs.optimized`, `optimized_linear`, `sensitivities`, `optimization` land; adjust tables if their shapes differ from below.
- After the first merge to main: Settings → Pages → Deploy from a branch → `main` / `/docs` (Sivan has to click this once), then check math and the viewer on the live URL.

**Blocked on / requests**
- WS2+3: the report reads these shapes (add keys freely; anything else falls back to a generic table):
  - `sensitivities.gradient_check`: list of `{function, variable, adjoint, fd, rel_error}`
  - `sensitivities.taylor`: `{eps: [...], <any series name>: [...]}` (log-log; e.g. `err_fd`, `err_adjoint`)
  - `sensitivities.cost`: `{n_vars, forward_s, adjoint_s, fd_s}`
  - `optimization.algorithm`, `optimization.scaling` (strings)
  - `optimization.starts`: list of `{name, iterations, final_mass_kg, converged, history: {iter: [...], mass_kg: [...], g_face_C: [...]}}`
  - `optimization.active`: list of `{constraint, multiplier, unit, meaning}` (meaning = one-line engineering interpretation)
  - `designs.optimized` / `designs.optimized_linear`: same shape as `designs.initial` (incl. `history`, `mesh.faces`).
- WS1: optional `true_max_C` next to `value_C` in the two temperature constraints (KS vs true peak is shown when present).
- AI-use statement (section 10): each WS please list one bug the verification caught in your section; I will quote them.

## WS5: Guide (GUIDE.md)

- **Done:** First draft of GUIDE.md from the plan: what is being built, each piece in plain language (heat PDE, finite elements, time stepping, Newton/nonlinearity, KS aggregation, verification ladder, adjoint, SLSQP, KKT multipliers, results.json + website), glossary, rubric map, AI-use statement skeleton.
- **In progress:** Nothing.
- **Next:** Update each "Status" box and fill the *TBD* numbers as the other workstreams report here (solver/ladder results, gradient checks, optimum and multipliers, viewer). Fill the AI-use table and "bug the verification caught".
- **Blocked on:** Other workstreams' plain-language notes and results.
- **Needs from others:** (1) Final mesh sizes and braking scenario numbers. (2) Ladder table with errors. (3) One bug the ladder or gradient check caught. (4) Active constraints and multiplier values at the optimum. (5) With vs without nonlinearity comparison.
