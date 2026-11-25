# ainote.md: shared status board

Project 4 topic: **EV battery busbar** (handout menu "Busbar or resistive heater"). Switched from the off-menu
brake disc on Nov 24 because no instructor approval arrived; the brake disc is archived, complete and runnable,
in `alternate_topic/` (its own status board is `alternate_topic/ainote.md`).

Deadline: GitHub Pages site final on **Nov 25, 11:59 PM America/Phoenix**.

Rules: `python -m pytest -q` must pass before merging to `main`. Only `scripts/run_all.py` writes
`docs/assets/results.json`; never hand-type numbers on the site. SI inside the code, Kelvin internally.

## Layout
- `busbar/`: params, mesh, model (coupled residual + tangent), solver, api, problem (constraints), ladder,
  run_all; generic pieces copied from the brake-disc study: fem, adjoint, ks, gradcheck, optimize.
- `docs/`: report (`index.md`, numbers filled from results.json by `assets/js/report.js`), viewer (`assets/js/viewer.js`).
- `tests/test_busbar.py`.

## Handout checklist (busbar)
| Requirement | Where |
|---|---|
| Objective + >= 2 constraints as PDE functionals, >= 20 variables | mass; bar T, joint T (KS over space and time), resistance; 24 vars (report s.2) |
| Own FE code, verification ladder, all kinds of tests | V1-V12 (report s.4) |
| Adjoint for every constraint, FD check per constraint and variable type, Taylor test, cost vs FD | report s.5 |
| Gradient optimizer, scaling, several starts, KKT active set and multipliers | SLSQP, 3 starts (report s.6) |
| Nonlinear component with/without comparison | N2 (report s.7) |
| Bonus fully nonlinear (+5) | Newton + implicit differentiation, nonlinear benchmark V4-V6 |
| Bonus time dependence (+5) | current pulse, adjoint swept backward, FD-checked |
| Bonus uncertainty (+5) | robust design, 3 scenarios (report s.9) |
| Website + 8 viewer requirements, math renders | `docs/` |
| One command reproduces everything | `python scripts/run_all.py` |
| Beginner guide for Sivan | `GUIDE.md` (busbar); brake disc guide in `alternate_topic/GUIDE.md` |
| AI-use statement with a bug the verification caught | report s.10; "what I checked myself" filled in by Sivan |

## Busbar: solver, verification, adjoint, optimization, site (owner: busbar thread)

**Why this topic.** The brake disc is off-menu and needed approval by Nov 4. Of the menu items, the
"Busbar or resistive heater" and the "Heat sink" both earn the fully-nonlinear (+5) and time-dependent (+5)
bonuses with our transient nonlinear heat solver; the busbar was chosen as the more impressive one
(two-way coupled electro-thermal multiphysics: sigma(T), Joule heating, contact resistance, monolithic Newton).

**Done** (numbers from `docs/assets/results.json`, written by `python scripts/run_all.py`, ~25 min)
- Coupled model `busbar/model.py`: state [potential, T, V_in]; bolted lap joints with contact resistance and a
  conductance into the 40 C terminal; natural convection h ~ dT^(1/4) and radiation on exposed faces; k(T).
  400 A steady (as step 1), 60 s at 1200 A, 240 s at 400 A. 24 variables: width and thickness at 12 stations.
- Verification `python -m busbar.ladder`: V1-V12 all pass (exact resistance 2e-16; 1D Joule profile nodally exact;
  closed-form two-way coupled sigma(T) benchmark, order 2; exact natural-convection and radiation cooling laws, order 1;
  3D coupled MMS order ~2; energy identities 1e-16; Newton quadratic; tangent vs complex step 4e-16; mesh study).
- Adjoint (reuses `optim/adjoint.py` unchanged) vs central FD: max error 4e-10 of the gradient norm; 14x faster than FD.
- SLSQP, 3 starts all reach 351.2 g from 403.2 g (-12.9 %). Initial bar fails the joint limit and the resistance
  budget; optimum passes all. Active: joint temperature only (0.011 kg/K). Widths at the 50 mm bound; thickness
  largest next to the joints (bar acts as a fin for the joint).
- N2: hand-calculation model (20 C sigma, k; fixed h) gives 322 g but fails the resistance budget with full physics
  (54.3 vs 52 micro-ohm); attribution shows sigma(T) is the nonlinearity that matters.
- Uncertainty bonus: robust design over (+10 % current, +5 K ambient) is 420.8 g (+70 g); nominal optimum fails the
  joint limit in both off-nominal scenarios.
- Site in `docs/` (report with all 10 handout sections, three.js viewer meeting the 8 requirements), checked in headless Chromium at 1280 and 390 px.

## Guide (GUIDE.md) (owner: guide thread, branch `claude/project-thread-c068hs`)
- **Done:** Busbar beginner guide as `GUIDE.md`: every piece in plain language, ladder V1-V12, gradient checks and cost, optimum and multiplier, hand-calculation (N2) comparison, robust bonus, glossary, rubric map, AI-use draft (section 8, Sivan's "what I checked myself" line filled in), 30-second summary (section 9). Numbers checked against `docs/assets/results.json` from 49526e2.
- **Next:** Re-check if results change.
- **Blocked on:** nothing.
