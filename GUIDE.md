# The beginner's guide to this project

*For someone who only half-follows the lectures and wants to know what was built, why, and how each piece earns points.*

> **Status of this guide (2026-11-24).** Updated from what the workstreams actually built and reported in [`ainote.md`](ainote.md), and checked against the final full run of `scripts/run_all.py` (`docs/assets/results.json`). The website always reads its numbers from `docs/assets/results.json`, so if this guide and the website ever disagree, **the website is right** and this guide is stale.

---

## 0. The whole project in one paragraph

We designed the **brake disc of a race car**: the spinning iron disc that the brake pads squeeze. We want it **as light as possible**, because every kilogram that spins and bounces with the wheel hurts acceleration and grip. But if we take too much iron away, the disc **gets too hot** during hard braking and cracks or loses grip ("brake fade"), and the heat soaks into the hub and **cooks the wheel bearing**. So the question is: *where exactly can we remove iron, and where must we keep it?* With 24 dimensions to tune, guessing by hand is hopeless. We wrote a program that (1) simulates how heat flows through the disc during three hard stops, (2) works out, for all 24 dimensions at once, "if I make this a bit thinner, how much hotter does the disc get?", and (3) gives that information to an optimizer, which keeps removing weight until the temperature limit is just barely met.

**The headline result:** the disc went from **9.25 kg to 7.20 kg (22 % lighter)** and still keeps the rubbing surface under 700 °C and the hub under 220 °C. If you leave out the nonlinear physics (the fact that hot iron conducts heat worse), the optimizer produces a 6.85 kg disc that actually reaches **722.8 °C** and would crack. So the nonlinear part of the model is what keeps the design safe.

---

## 1. The object: what a vented brake disc looks like

A race disc is not a solid plate. It is two flat rings ("**cheeks**", the parts the pads rub on) joined by many thin walls ("**vanes**"). Air is flung outward through the gaps between the vanes, like a centrifugal fan, which cools the disc from the inside.

```
   Side cut through the disc (one slice)        Looking at it from the front

      outer edge (r = 178 mm)                        . - ~ ~ ~ - .
   ┌──────────────┐  ← cheek (pad rubs here)     .'   \  |  /    '.
   │              │                             /  \   \ | /   /   \
   │  air channel │  ← vane between cheeks     |  --  ( hub )  --   |   ← 48 vanes
   │   (12 mm)    │                             \  /   / | \   \   /      like
   └──────────────┘  ← cheek (pad rubs here)     '.   /  |  \    .'      spokes
      inner edge (r = 110 mm) → hub / "bell"        ' - ~ ~ ~ - '
```

Our disc: 356 mm diameter, 48 straight radial vanes, grey cast iron, a 12 mm air gap between the cheeks.

**Symmetry trick.** All 48 vanes are identical, and each vane is also mirror-symmetric (left/right, and top cheek/bottom cheek). So we simulate only **1/192 of the disc**: half of one vane and the one cheek above it. The rest of the disc is copies of that slice. This is called *symmetry reduction*. It makes the computation roughly 200 times cheaper and loses no accuracy. The website's 3D viewer mirrors and copies the slice back into a full disc for display.

---

## 2. The design problem, written like an engineer

The handout says every problem has the form *"minimize something, subject to limits"*. Ours:

| Piece | In plain words | In our project |
|---|---|---|
| **Design variables** `d` (the knobs we turn) | The 24 numbers the optimizer is allowed to change | Vane thickness at 12 radii from inner to outer edge (2 to 8 mm), plus cheek thickness at the same 12 radii (5 to 12 mm). Starting design: all vanes 5 mm, all cheeks 9 mm |
| **Objective** `f` (what we minimize) | The thing we want small | Mass of the whole disc |
| **Constraint 1** | Don't overheat the rubbing surface | Peak friction-face temperature ≤ 700 °C. Above this, grey cast iron cracks and pads fade (Limpert, *Brake Design and Safety*; Day, *Braking of Road Vehicles*) |
| **Constraint 2** | Don't cook the wheel bearing | Peak temperature where the disc bolts to the hub ≤ 220 °C. That is the minimum dropping point (the temperature where it turns liquid and runs out) of NLGI GC-LB wheel-bearing grease (ASTM D4950), and it is also below the 230 °C minimum boiling point of DOT 4 brake fluid (FMVSS 116) |
| **Constraint 3** | Don't make parts impossibly thin | The thickness bounds above: 2 mm vanes are about the casting limit; 5 mm cheeks leave pad-wear allowance and strength |
| **The "state"** `u` | What the physics tells us once the knobs are set | The temperature at every point of the slice at every moment of the braking run |

**Why 12 + 12?** The handout requires at least 20 design variables, so that the clever gradient method matters. Letting thickness vary from inner to outer edge lets the optimizer put iron where the heat goes in (the band the pad rubs on) and remove it where it does no work.

**The braking scenario.** A 900 kg race car (60 % of the braking on the front wheels) makes **three stops from 250 to 80 km/h at 1.5 g**, with 12 s of re-acceleration in between, like three braking zones of a lap. The disc starts at 150 °C, and the air is 35 °C. 90 % of the friction heat goes into the disc and the rest into the pad. One stop is easy. Repeated stops are what overheat discs, because the heat keeps building up.

---

## 3. The big picture: how the pieces connect

```
          ┌───────────────────────────────────────────────────────────────┐
          │                     THE OPTIMIZATION LOOP                      │
          │                                                                │
  start   │   24 thicknesses ──► [1] BUILD MESH ──► [2] HEAT SIMULATION    │
  design ─┼─►      d                (geometry)        (temperatures over   │
          │       ▲                                     the whole braking   │
          │       │                                     run, ~1.8 s)       │
          │       │                                         │              │
          │       │                                         ▼              │
          │  [5] OPTIMIZER  ◄── mass, max temps ◄── [3] SCORE THE DESIGN   │
          │   (SLSQP: "move   + their gradients           │                │
          │    the 24 knobs        ▲                      ▼                │
          │    this way")          └──────────── [4] ADJOINT: gradients    │
          │                                      of every score w.r.t.     │
          │                                      all 24 knobs (~2.4 s)     │
          └───────────────────────────────────────────────────────────────┘
                     │ when nothing improves any more
                     ▼
          [6] RESULTS FILE (results.json) ──► [7] WEBSITE with 3D viewer
```

Before any of this ran on the real disc, piece [2] had to pass the **verification ladder** (section 4.6). The handout's tip says it best: *"A wrong gradient of a wrong model optimizes beautifully."*

Where each piece lives in the repository:

| Folder / file | What's in it | Workstream | Rubric points it mainly earns |
|---|---|---|---|
| `brakeopt/` | Mesh, heat solver, Newton, functionals | WS1 | Formulation 10, Analysis & verification 20, Nonlinear 15, fully-nonlinear bonus +5 |
| `verification/ladder.py` | The 10 verification tests (13 checks) | WS1 | Analysis & verification 20 |
| `optim/` | Adjoint, KS, gradient checks, SLSQP, with/without nonlinearity, uncertainty | WS2+3 | Sensitivity 20, Optimization 15, Nonlinear 15, time-dependence bonus +5, uncertainty bonus +5 |
| `docs/` | The website: report (`index.md`), 3D viewer, `results.json` | WS4 | Website 10, Presentation 5 |
| `scripts/run_all.py` | The one command that regenerates every number | WS1 | Reproducibility 5 |

---

## 4. Piece by piece

### 4.1 The physics: the heat equation

**What it is.** One rule governs how temperature moves through solid iron: *heat flows from hot to cold, faster when the iron conducts well and the temperature difference is steep, and a chunk of iron warms up according to how much heat goes in compared with how much it can store.* Written as math:

$$\rho c \frac{\partial T}{\partial t} = \nabla\cdot\big(k(T)\,\nabla T\big)$$

Read it as: *(how fast a point warms up) × (how much heat it takes to warm it) = (net heat flowing into that point)*. `ρ` is density (7100 kg/m³), `c` is heat capacity (500 J/kg·K), and `k` is thermal conductivity (54 W/m·K when cold, dropping as the iron heats up).

**Where the heat comes from and where it goes.** This is what the "boundary conditions" describe:

```
   pads rubbing  ──►  heat INTO the cheek faces   (the car's kinetic energy turned into heat;
                                                   90 % goes to the disc, 10 % to the pad)
   air in vanes  ◄──  heat OUT by convection      (h × (T − T_air), h from textbook formulas
                                                   for air pumped through the vane channels)
   air on faces  ◄──  heat OUT by convection      (formula for a spinning disc)
   glowing hot   ◄──  heat OUT by radiation       (εσ(T⁴ − T_air⁴); matters a lot at 600+ °C)
   hub bolts     ◄──  heat OUT into the hub/bell  (a "contact conductance" at the inner edge)
```

**Analogy.** Think of the disc as a sponge. The pads pour water (heat) onto its two faces during each stop. The water spreads through the sponge (conduction) and drips out of every surface (convection and radiation). More sponge (more iron) holds more water before it overflows (overheats), but it weighs more.

**Why this matters for the goal.** Our constraints are temperatures, so we need the temperature everywhere, at every moment. This equation is the "PDE" in "PDE-constrained optimization".

### 4.2 Finite elements: turning the equation into numbers

**What it is.** A computer can't solve the equation over a smooth curved shape directly. So we chop the slice of disc into small bricks ("**hexahedral elements**"), store a temperature at each brick corner ("**node**"), and approximate the temperature inside each brick by blending its 8 corner values. The heat equation then turns into a big set of ordinary equations, one per node:

$$M\,\dot{T} + K(T)\,T = F(t)$$

- `M` (the "capacity matrix"): how much heat each node can store.
- `K` (the "conductivity matrix"): how strongly each node is connected to its neighbours.
- `F`: the heat coming in from the pads (and going out to the air).

**What we actually used:** 120 bricks and 231 nodes in the optimizer's mesh. The bricks follow the true curved shape of the disc.

**The mesh changes with the design.** When the optimizer changes a thickness, the bricks stretch or squash, but the number of bricks and how they connect stay the same. That is how the design variables reach the physics, and it keeps the gradients smooth: the code knows *exactly* how far each node moves per millimetre of thickness change.

**Sparse matrices.** Each node only touches its neighbours, so `K` and `M` are almost all zeros. We store only the non-zeros (`scipy.sparse`), as the handout requires. The solver is written by us, not ANSYS or FEniCS.

### 4.3 Stepping through time (implicit Euler)

The braking run lasts about a minute, so we march forward in small time steps: 0.1 s while braking and 0.5 s while accelerating, which gives 144 steps. "Implicit Euler" means each step solves for the *new* temperatures using the *new* conditions, which stays stable even with fairly large steps. One full braking run takes about **1.8 seconds** on a laptop.

### 4.4 The nonlinear part: why Newton's method is needed

**What makes it nonlinear.** Two things:
1. **Radiation** goes like `T⁴`. Double the absolute temperature and the radiated heat goes up sixteen times.
2. **Conductivity changes with temperature**: hot cast iron conducts heat worse than cold iron, so `k` depends on `T`.

So each time step can't be solved in one shot, because the answer depends on itself.

**Newton's method, in plain words.** Guess the temperatures, measure how far off the energy balance is (the "**residual**"), use the slope of that error (the "**tangent matrix**") to make a much better guess, and repeat. Near the answer, the number of correct digits roughly **doubles each iteration** ("quadratic convergence"). This is what our solver actually did on the braking problem:

```
   iteration   residual (how wrong)
       1        2e+02
       2        6e-01
       3        2e-06     ← digits roughly double each time
       4        3e-13     ← done (round-off level)
```

That doubling proves our tangent matrix is right. It is ladder test V8.

**Why it matters for the goal, and what we found.** The handout *requires* a nonlinear component. Ours is type **N2**: the PDE itself is nonlinear and is solved with Newton's method. Doing it properly, with exact gradients and a nonlinear benchmark test, earns the **+5 "fully nonlinear" bonus**.

The handout also asks for the optimized design **with and without** the nonlinearity:

| | Optimized with full physics | Optimized with the linear model |
|---|---|---|
| Mass | 7.20 kg | 6.85 kg |
| Peak face temperature, re-checked with full physics | 695.8 °C ✓ | **722.8 °C ✗ (fails 700 °C)** |

Why the difference: if you pretend hot iron conducts as well as it does at 200 °C, the model under-predicts the peak by about 34 °C. Simplifying radiation goes the other way (it over-predicts by about 9 °C), but that is smaller. So **the temperature-dependent conductivity is what changes the design**, and leaving it out gives a disc that looks fine on paper but would crack.

### 4.5 Scoring a design: a smooth "maximum" (KS aggregation)

**Problem.** "Peak temperature" means the maximum over every node and all 144 time steps. The plain max is not smooth (which point is hottest can jump from one place to another), and optimizers need smooth functions.

**Fix.** The **KS function** (Kreisselmeier–Steinhauser) is a "soft maximum": a smooth formula that always sits slightly *above* the true max. Using it as the constraint is safe (conservative), and it is smooth enough to take gradients of. The report shows the true max next to the KS value.

**Analogy.** Instead of asking "who is the tallest person in the room?" (the answer jumps when someone stands up), you compute a weighted average that is dominated by the tallest few people. It gives almost the same answer, but it changes gradually.

**A tuning decision that mattered.** The KS function has a "sharpness" setting (ρ). With the first setting (0.1 per K), the soft max sat about 47 °C above the true face peak and 65 °C above the hub peak. That would have forced a needlessly heavy disc. At 1.0 per K it sits only 3 to 5 °C above. At the optimum, the face KS value is exactly 700.0 °C while the true peak is 695.8 °C, and the hub KS value is exactly 220.0 °C while the true peak is 215.4 °C.

### 4.6 Verification: proving the solver is right before trusting it

**Why.** A solver written with a coding agent "is wrong until proven otherwise" (handout). We prove it by running problems where we already know the exact answer and showing that our code gets it. The handout requires at least six tests in a table: our number, the exact number, and the relative error. We have ten tests (13 rows in the report table, because V3, V4 and V7 each have two parts), and **all pass** (`python -m verification.ladder`):

| ID | Test, in plain words | Result |
|----|------|--------|
| V1 | A plain wall with air on both sides; the exact answer is a straight-line temperature profile | Error 1e-16 (round-off) |
| V2 | The same wall, rotated and moved to a random spot in 3D. Nothing should change | Identical to 3e-16 |
| V3 | A wall whose conductivity changes with temperature, against the exact formula | Heat flow error 4e-7; error shrinks 4× per mesh halving (order 2.00) |
| V4 | A plate cooling only by glowing (T⁴ radiation), against the exact cooling curve | Error 2e-6; error halves when the time step halves (order 0.99) |
| V5 | "Manufactured solution" on the real curved disc mesh, with every term switched on: we pick the answer first, work out the heating that would produce it, and check the solver recovers it | Error order 1.98 in mesh size (textbook: 2) |
| V6 | The same, but changing in time | Error order 1.00 in time step (textbook: 1) |
| V7 | Energy bookkeeping: heat put in equals the car's lost kinetic energy, and heat in = heat stored + heat lost | Exact; imbalance 2e-16 |
| V8 | Newton on the real braking problem converges quadratically | Order about 2.2 |
| V9 | Mesh refinement of the peak face temperature (1×, 2×, 3× finer) | 683.2, 686.9, 687.6 °C. The optimizer's coarse mesh reads 4.7 °C (0.5 %) below the converged value |
| V10 | Every outside face of the mesh belongs to exactly one boundary type | 184 = 184 |

**Honest caveat to mention in the report:** V9 shows the optimizer's mesh reads the peak about 5 °C low. So the "700 °C" design is really closer to 705 °C on a fine mesh. That belongs in the Limitations section.

**Important rule from the handout:** the expected answers must be derived by hand or from a textbook *before* the agent sees the test, so that the code can't just "agree with itself".

### 4.7 Gradients: "which way should each knob turn?"

**What a gradient is.** For each of the 24 thicknesses, the gradient says: *if I make this thickness slightly bigger, how much does the peak temperature (or the mass) change?* It is a list of 24 sensitivities, one per knob, for each score.

**The slow way (finite differences).** Nudge knob 1, rerun the whole braking simulation, and see what changed. Repeat for all 24 knobs. On our model that takes **87 seconds** for one gradient, and the optimizer needs a gradient every iteration.

**The clever way (the adjoint method).** Run the simulation once forward. Then run **one backward pass**, and out come the sensitivities of both temperature constraints to *all 24 knobs at once*. On our model it takes **2.4 seconds, about 22× faster**. The gap grows with the number of knobs, because the adjoint's cost doesn't depend on how many there are.

**Analogy.** You're at a party and want to know how much each of 24 guests added to the noise. The slow way: ask each guest to leave, one at a time, and measure. The adjoint way: record the room once, then play the recording backward through a model of the room's acoustics, and it tells you every guest's contribution in one pass.

**Why "backward in time".** Our problem changes over time: the temperature during stop 3 depends on the heat left over from stops 1 and 2. So the adjoint starts at the end of the braking run and sweeps backward to the start, asking at each step "how much did this moment's temperatures influence the peak?". That backward sweep is exactly what earns the **+5 "path or time dependence" bonus**.

**Connection to the KKT notes.** The adjoint variable `λ` is the Lagrange multiplier on the heat equation. "Implicit differentiation", "the adjoint method" and "the Lagrange multiplier of the state equation" are three names for the same idea.

**How we proved the gradients are right.**
- **Gradient-check table:** the adjoint answer against the brute-force (central finite difference) answer, for mass, the face constraint and the hub constraint, on both vane and cheek thicknesses. They agree to about **9 digits** (worst case 8e-8 of the gradient size).
- **Taylor test plot:** as the finite-difference step ε shrinks, the mismatch falls like ε² (100× smaller per 10× smaller step) until round-off takes over. A flat line would mean a bug.
- **Cost comparison:** forward 1.8 s, adjoint 2.4 s, finite differences 87 s.

### 4.8 The optimizer (SLSQP)

**What it does.** SLSQP (from `scipy.optimize`) is a standard gradient-based optimizer. At each iteration it takes the mass, the temperature constraints and their gradients, builds a simplified local model of the problem, and jumps to the best design of that model. It repeats until it can't improve without breaking a limit.

**Analogy.** Hiking downhill in fog toward the lowest point of a valley, with fences you may not cross. The gradient tells you which way is downhill; the constraints are the fences. You end up at the lowest spot you can reach, usually pressed against one or two fences.

**What we did, following the handout's rules:**
- **Scaling:** the optimizer works with `x = d / d₀`, mass divided by the starting mass, and constraints written as `KS/limit − 1 ≤ 0`, so every number is around 1. Optimizers behave badly when one number is 0.003 and another is 700.
- **Three starting designs:** the baseline, every thickness at its maximum, and every thickness halfway between its bounds. **All three ended at 7.20 kg**, within 0.0001 kg of each other. That strongly suggests we found the real best design, not a local dip.
- **History saved** for every iteration (design, mass, temperatures) for the convergence plot.

**What the optimized disc looks like, and why it makes sense:**
- **Vanes** at stations 2 to 11 went to the 2 mm minimum. Thin vanes are cheap in weight and the air still flows.
- The **outermost vane** went to its 8 mm maximum: that's where the most heat goes in, and the extra iron conducts it away from the face.
- **Inner cheeks** (stations 2 to 7, inside the pad band) went to the 5 mm minimum. The pad doesn't rub there, so the iron does little work.
- **Right at the hub** (station 1), the optimizer keeps a little extra iron: a 2.4 mm vane and a 6.8 mm cheek. This appeared once the hub limit became active. Our reading (not separately tested) is that this iron acts as a small buffer that soaks up heat before it reaches the hub.
- **Outer cheeks** stay thick (9.3 to 11.0 mm). That iron acts as a heat sponge in the band the pad rubs on, absorbing each stop's heat.

### 4.9 Reading the answer: KKT conditions and multipliers

At the optimum, a constraint is **active** if the design is pressed against its limit, and **inactive** if it is comfortably met. In our design **both temperature limits are active**.

| Constraint | Value at the optimum | Status |
|---|---|---|
| Friction face ≤ 700 °C | KS 700.0 °C (true peak 695.8 °C) | **Active** |
| Hub ≤ 220 °C | KS 220.0 °C (true peak 215.4 °C) | **Active** |
| Thickness bounds | many at their limits (section 4.8) | Active bounds |

Each active constraint has a **Lagrange multiplier**, which has a concrete meaning:

> **The multiplier is how much disc mass we would save if that limit were relaxed by one unit.**

Ours:

| Limit | Multiplier | What it means |
|---|---|---|
| Friction face | **0.0188 kg per °C** (about 19 g per degree) | A pad and disc material that tolerated 710 °C instead of 700 °C would save about **0.19 kg** |
| Hub | **0.0058 kg per °C** (about 6 g per degree) | A grease rated 10 °C higher would save only about **0.06 kg** |

So the face limit is worth about three times as much per degree as the hub limit: it is still the main thing sizing the disc. The hub limit costs a little extra iron near the hub. The rubric asks for exactly this kind of engineering reading.

### 4.10 Bonus: design under uncertainty

A driver won't always brake from exactly 250 km/h. We also designed for three scenarios at once: nominal, **20 km/h faster** (270 → 80 km/h), and **5 % more front brake bias**.

- The nominal 7.20 kg optimum would reach **785 °C** on the face and 235 °C at the hub in the faster scenario, and 746 °C with more front bias. So it fails in both off-nominal cases.
- The **robust design** passes all six constraints (two per scenario) at **9.67 kg**, which is 2.5 kg heavier than the nominal optimum. That is the price of the safety margin.
- The constraints that decide it are the face and hub temperatures in the faster scenario, with multipliers of 0.059 and 0.0046 kg per °C.
- Each scenario's gradients were also checked against finite differences (about 9 digits).

This earns the optional **+5 "design under uncertainty" bonus**.

### 4.11 Results file and website

**`results.json`.** One file, written only by `scripts/run_all.py`, holding every number the report shows and the geometry and temperatures that the 3D viewer draws. **No number is typed into the report by hand**: when the page opens, a script reads `results.json` and fills in every number, table and plot. So if the solver changes, one command (`python scripts/run_all.py`, about 48 minutes) regenerates the whole site. That earns the **Reproducibility 5 points**.

**The website.** A GitHub Pages site built from Markdown (`docs/index.md`, Jekyll), with math rendered by MathJax and a three.js 3D viewer. The viewer meets all 8 handout requirements:

1. The initial and optimized discs side by side, with one shared camera. A dropdown can swap in the linear-model design for the N2 comparison.
2. True 3D geometry: the computed slice is mirrored and copied 48 times into the full disc.
3. Coloured by temperature with a colour bar, on a single-hue orange scale (no rainbow). Over-limit regions are **hatched**, so you can spot them without relying on colour.
4. The temperature field is drawn on the surfaces using the bricks' own blending functions.
5. One animation: play or scrub through the three stops and watch the disc heat up and cool down, with a peak-temperature and speed chart.
6. Hover or tap any face to see its element, radius, current temperature, peak and margin.
7. A metrics panel for each design: mass and every constraint, each marked pass or fail.
8. Loads from `results.json`, works at phone width (tested at 390 px), and loads in about 1 second.

The report itself has the 10 sections the handout lists (Summary, Formulation, Analysis model, Verification, Sensitivities, Optimization, Nonlinear component, Results table, Bonus, Limitations + AI-use statement).

**One thing only you can do:** turn on GitHub Pages in the repository: **Settings → Pages → Deploy from a branch → `main`, folder `/docs`**. Then open the site in a private browser window to check it, and submit both links (repository and site) on Canvas.

---

## 5. Glossary (one line each)

| Term | Meaning |
|---|---|
| **PDE** | Partial differential equation; here, the heat equation |
| **State** `u` / `T` | The physics answer for a given design: the temperatures everywhere, at every time step |
| **Design variables** `d` | The 24 thicknesses the optimizer may change |
| **Functional** | A single number computed from the whole temperature field (e.g. peak temperature) |
| **Finite element / mesh** | Chopping the shape into small bricks so that a computer can solve the PDE |
| **Node** | A brick corner where a temperature is stored |
| **Residual** | How far a guessed solution is from satisfying the equations |
| **Tangent matrix** | The slope of the residual; what Newton's method uses to improve the guess |
| **Newton's method** | Repeatedly fixing a guess using the slope; the digits double each step near the answer |
| **Implicit Euler** | A stable way to step through time |
| **Symmetry reduction** | Simulating one identical slice (1/192) instead of the whole disc |
| **Convection / radiation** | Heat lost to moving air / emitted as glow (∝ T⁴) |
| **k(T)** | Thermal conductivity that depends on temperature (iron conducts worse when hot) |
| **Manufactured solution** | A test where you pick the answer first and work out the problem that produces it |
| **Convergence order** | How fast the error shrinks as the mesh or time step shrinks (2 means 4× less error per halving) |
| **KS function** | A smooth, slightly conservative version of "maximum" |
| **Gradient / sensitivity** | How much a score changes per small change in each design variable |
| **Finite differences** | Getting a gradient by nudging each variable and rerunning: simple but slow |
| **Adjoint method** | Getting the whole gradient with one extra backward pass |
| **Taylor test** | A plot proving the gradients are correct (the error should fall like ε²) |
| **SLSQP** | The optimizer we use (Sequential Least SQuares Programming) |
| **Active constraint** | A limit the optimum is pressed against |
| **Lagrange multiplier** | How much the objective would improve per unit relaxation of a limit |
| **KKT conditions** | The mathematical conditions that hold at a constrained optimum |
| **N2** | The handout's label for "the PDE itself is nonlinear" |
| **results.json** | The single file that all report numbers and the viewer come from |

---

## 6. How the pieces map to the rubric

| Rubric category | Points | Where it comes from in this project | Guide section |
|---|---|---|---|
| Structure and formulation | 10 | Brake disc motivation, 24 variables, mass objective, cited temperature limits | 2 |
| Analysis code and verification | 20 | Our own FE heat solver + the 10-test ladder, all passing | 4.1–4.6 |
| Sensitivity analysis | 20 | Transient adjoint, gradient-check table (about 9 digits), Taylor plot, 20× cost comparison | 4.7 |
| Optimization | 15 | SLSQP, scaling, 3 starts that agree, both active constraints and their multipliers interpreted | 4.8–4.9 |
| Nonlinear component | 15 | k(T) + radiation with Newton; the linear-model design fails at 722.8 °C | 4.4 |
| Website and 3D comparison | 10 | Pages site, all 8 viewer requirements, math renders | 4.11 |
| Reproducibility | 5 | `python scripts/run_all.py` regenerates results.json and so every number | 4.11 |
| Presentation and clarity | 5 | Readable report (and this guide helps you explain it) | all |
| Bonus: fully nonlinear | +5 | N2 with Newton + implicit gradients + nonlinear benchmarks V3, V4 | 4.4, 4.6 |
| Bonus: time dependence | +5 | Adjoint swept backward through the three stops, checked by FD | 4.7 |
| Bonus: uncertainty | +5 | Three braking scenarios, robust 9.67 kg design | 4.10 |

The bonus is capped at +15, and we target all three bonuses.

---

## 7. Things to watch / open risks

- **Structure approval.** The handout says teams proposing their own structure needed instructor approval by Nov 4. A brake disc is close to the menu's "heat sink" (the vanes act like fins), but it is not on the menu. Worth a short email to the instructor.
- **Turn on GitHub Pages** (section 4.11) and check the live site in a private window. Math only renders on the live site, not on github.com.
- **Deadline:** the final site is due Nov 25; the commit at the deadline is graded.
- **Limitations worth stating in the report:** the coarse optimizer mesh reads the peak about 5 °C low (V9); the convection coefficients come from textbook correlations, not measurement; mechanical stress and thermal distortion (coning) are not modelled.

---

## 8. The AI-use statement (required, report section 10)

The handout requires: *"what the coding agent wrote, which parts you derived or checked yourselves, and one bug the verification caught."* Below is a draft built from what the workstreams recorded. **The parts in brackets are yours to fill in honestly**: only you know what you checked yourself, and that is what gets graded.

> **AI-use statement.** The solver, adjoint, optimization scripts and website were written with a coding agent (Claude Code), working in parallel workstreams that coordinated through a shared status file (`ainote.md`). *[Your part, for example: "I chose the structure, the design variables, the objective and the constraint limits and their sources; I reviewed the verification ladder's exact solutions against the textbook; I read every verification and gradient-check result and the optimized design for physical sense."]*
>
> | Component | Written by | How it was checked |
> |---|---|---|
> | Mesh generator (`brakeopt/mesh.py`) | Agent | V2 (rotation invariance), V5 (manufactured solution on the curved mesh), V10 (every boundary face in exactly one set) |
> | Heat solver + Newton (`brakeopt/model.py`, `solver.py`) | Agent | Ladder V1–V9 against exact solutions; Newton quadratic convergence (V8); energy balance to 2e-16 (V7) |
> | Adjoint gradients (`optim/adjoint.py`) | Agent | Adjoint vs central finite differences for every constraint and both variable types (agree to about 9 digits); Taylor test |
> | Optimization (`optim/optimize.py`) | Agent | Three starting designs reach the same 7.20 kg; constraints re-checked at the optimum; linear-model optimum re-checked with full physics |
> | Website and 3D viewer (`docs/`) | Agent | Automated test that every number on the page exists in results.json; checked by hand against the 8 viewer requirements at desktop and phone width |
>
> **Bugs the verification caught:**
> 1. *Sign error in the heat flux (caught by V3).* The agent's first version measured the heat flow leaving the cold end of the test wall with the wrong sign, so V3 failed (−1.29×10⁶ vs +1.29×10⁶ W/m²). Comparing against the exact formula caught it.
> 2. *Newton stopping too early (caught by the gradient check).* On the small stand-in model, the first Newton stopping rule quit after one iteration. The adjoint and finite-difference gradients then disagreed (4×10⁻⁷) and, more tellingly, the Taylor-test curve was **flat** instead of falling like ε², which is the handout's sign of a bug. It was fixed by stopping on the size of the Newton update, and the real solver got an extra Newton iteration.
> 3. *KS too loose (caught by comparing KS with the true maximum).* With the first KS sharpness, the soft maximum overestimated the true peak by 47 to 65 °C, which would have made the optimizer leave in iron it didn't need. Raising ρ brought the gap down to 3 to 5 °C.

---

## 9. Explaining it in 30 seconds

> "We designed a lighter race-car brake disc. Our own finite-element code simulates how heat from the pads spreads through the disc during three hard stops, including the fact that hot iron conducts heat worse and that a hot disc glows heat away. The adjoint method tells us, in one backward pass, how every one of the 24 thicknesses affects the peak temperature, and we checked those gradients against brute force to 9 digits. The optimizer then cut the disc from 9.25 to 7.20 kg, keeping iron mainly in the band the pad rubs on. Two limits stop it: 700 °C on the rubbing surface and 220 °C at the hub, where the wheel-bearing grease would start to run. Each extra degree allowed on the face would save about 19 grams. If you ignore the nonlinear physics, you get a lighter disc that actually reaches 723 °C and cracks."

---

## 10. Change log of this guide

- 2026-11-24, first draft from the project plan, before any code existed.
- 2026-11-24, hub limit changed to 220 °C (grease dropping point). Updated to the rerun: 7.20 kg optimum, both limits active (multipliers 0.0188 and 0.0058 kg/°C), linear-model design 6.85 kg failing at 722.8 °C, robust design 9.67 kg.
- 2026-11-24, checked every number against the final full-run results.json (timings updated: forward 1.8 s, adjoint 2.5 s, FD 85 s).
- 2026-11-24, updated with what was actually built: mesh and scenario details, the 10-test ladder results, gradient-check and cost numbers, the optimized design and its multiplier, the with/without-nonlinearity comparison, the uncertainty bonus, the website features, and a filled-in AI-use statement draft with three caught bugs.
