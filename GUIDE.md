# The beginner's guide to this project

*For someone who only half-follows the lectures and wants to know what is being built, why, and how each piece earns points.*

> **Status of this guide.** Written from the project plan before the code exists. Every "**Status:**" box below says what has actually been built so far. Numbers marked *TBD* get filled in once the solver produces them. The source of truth for progress is [`ainote.md`](ainote.md).

---

## 0. The whole project in one paragraph

We are designing the **brake disc of a race car**: the spinning iron disc that the brake pads squeeze. We want it **as light as possible**, because every kilogram that spins and bounces with the wheel hurts acceleration and grip. But if we take too much iron away, the disc **gets too hot** during hard braking and cracks or loses grip ("brake fade"), and the heat soaks into the hub and **cooks the wheel bearing and brake fluid**. So the question is: *where exactly can we remove iron, and where must we keep it?* With 24 dimensions to tune, guessing by hand is hopeless. We write a program that (1) simulates how heat flows through the disc during a series of hard stops, (2) works out, for every one of the 24 dimensions at once, "if I make this a bit thinner, how much hotter does the disc get?", and (3) hands that information to an optimizer that keeps shaving weight until the temperature limits are just barely met. Then we show the before/after discs in an interactive 3D website.

---

## 1. The object: what a vented brake disc looks like

A race disc is not a solid plate. It is two flat rings ("**cheeks**", the parts the pads rub on) joined by many thin walls ("**vanes**"). Air is flung outward through the gaps between the vanes, like a centrifugal fan, which cools the disc from the inside.

```
   Side cut through the disc (one slice)        Looking at it from the front

      outer edge                                     . - ~ ~ ~ - .
   ┌──────────────┐  ← cheek (pad rubs here)     .'   \  |  /    '.
   │              │                             /  \   \ | /   /   \
   │  air channel │  ← vane between cheeks     |  --  ( hub )  --   |   ← vanes
   │              │                             \  /   / | \   \   /      like
   └──────────────┘  ← cheek (pad rubs here)     '.   /  |  \    .'      spokes
      inner edge → hub / "bell"                     ' - ~ ~ ~ - '
```

**Symmetry trick.** All the vanes are identical, so the disc is the same slice repeated around the circle (like slices of a pizza). We only simulate **one slice** (one "vane pitch sector") and copy it around for display. This is called *cyclic symmetry*, and it makes the computation dozens of times cheaper with no loss of accuracy.

---

## 2. The design problem, written like an engineer

The handout says every problem has the form *"minimize something, subject to limits"*. Ours:

| Piece | In plain words | In our project |
|---|---|---|
| **Design variables** `d` (the knobs we turn) | The 24 numbers the optimizer is allowed to change | Vane thickness at 12 points from inner to outer radius, plus cheek thickness at 12 points from inner to outer radius |
| **Objective** `f` (what we minimize) | The thing we want small | Mass of the disc |
| **Constraint 1** | Don't overheat the rubbing surface | Peak friction-face temperature ≤ about 700 °C (above that grey cast iron cracks and the pads fade; source to be cited in the report) |
| **Constraint 2** | Don't cook the bearing / brake fluid | Temperature where the disc mounts to the hub ≤ a limit (source to be cited) |
| **Constraint 3** | Don't make parts impossibly thin | Every thickness ≥ a minimum you could actually cast and that survives the clamping force |
| **The "state"** `u` | What the physics tells us once the knobs are set | The temperature at every point in the disc at every moment of the braking run |

**Why 12 + 12?** The handout requires at least 20 design variables (so the clever gradient method matters). Letting thickness vary from inner to outer radius lets the optimizer put iron where heat is generated (the outer part, which moves fastest under the pad) and remove it where it does no work.

**The braking scenario.** Several hard stops in a row, e.g. 250 → 80 km/h, with short cool-down gaps between them, like the braking zones of a lap. One stop is easy; repeated stops are what overheat discs, because heat keeps piling up.

---

## 3. The big picture: how the pieces connect

```
          ┌───────────────────────────────────────────────────────────────┐
          │                     THE OPTIMIZATION LOOP                      │
          │                                                                │
  start   │   24 thicknesses ──► [1] BUILD MESH ──► [2] HEAT SIMULATION    │
  design ─┼─►      d                (geometry)        (temperatures over   │
          │       ▲                                     the whole braking   │
          │       │                                     run)               │
          │       │                                         │              │
          │       │                                         ▼              │
          │  [5] OPTIMIZER  ◄── mass, max temps ◄── [3] SCORE THE DESIGN   │
          │   (SLSQP: "move   + their gradients           │                │
          │    the 24 knobs        ▲                      ▼                │
          │    this way")          └──────────── [4] ADJOINT: gradients    │
          │                                      of every score w.r.t.     │
          │                                      all 24 knobs at once      │
          └───────────────────────────────────────────────────────────────┘
                     │ when nothing improves any more
                     ▼
          [6] RESULTS FILE (results.json) ──► [7] WEBSITE with 3D viewer
```

Before any of this runs on the real disc, piece [2] has to pass the **verification ladder** (section 5). The handout's tip says it best: *"A wrong gradient of a wrong model optimizes beautifully."*

The work is split into four workstreams, each in its own thread and branch, coordinated through `ainote.md`:

| Workstream | Builds | Rubric points it mainly earns |
|---|---|---|
| 1. Solver + verification | Pieces [1], [2], [3] and the ladder | Formulation 10, Analysis & verification 20, Nonlinear 15, Fully-nonlinear bonus +5 |
| 2. Adjoint gradients | Piece [4] and the gradient checks | Sensitivity 20, Time-dependence bonus +5 |
| 3. Optimization | Piece [5], the KKT interpretation, with/without nonlinearity | Optimization 15, Nonlinear 15 (the "effect on the design" half) |
| 4. Website | Pieces [6], [7], the report | Website 10, Reproducibility 5, Presentation 5 |

---

## 4. Piece by piece

### 4.1 The physics: the heat equation

**What it is.** One rule governs how temperature moves through solid iron: *heat flows from hot to cold, faster when the iron conducts well and the temperature difference is steep, and a chunk of iron warms up according to how much heat goes in versus how much it can store.* Written as math:

$$\rho c \frac{\partial T}{\partial t} = \nabla\cdot\big(k(T)\,\nabla T\big)$$

Read it as: *(how fast a point warms up) × (how much heat it takes to warm it) = (net heat flowing into that point)*. `ρ` is density, `c` heat capacity, `k` thermal conductivity.

**Where the heat comes from and goes.** This is what the "boundary conditions" describe:

```
   pads rubbing  ──►  heat INTO the cheek faces   (the car's kinetic energy turned to heat;
                                                   a "partition" fraction goes to the disc,
                                                   the rest to the pads)
   air in vanes  ◄──  heat OUT by convection      (h × (T − T_air), h from textbook correlations)
   air on faces  ◄──  heat OUT by convection
   glowing hot   ◄──  heat OUT by radiation       (εσ(T⁴ − T_air⁴): matters a lot at 600+ °C)
```

**Analogy.** Think of the disc as a sponge. The pads pour water (heat) onto its two faces during each stop. Water spreads through the sponge (conduction), and drips out of every surface (convection, radiation). More sponge (more iron) holds more water before it overflows (overheats) but weighs more.

**Why this matters for the goal.** Our constraints are temperatures, so we need the temperature everywhere, at every moment. This equation is the "PDE" in "PDE-constrained optimization".

### 4.2 Finite elements: turning the equation into numbers

**What it is.** A computer can't solve the equation over a smooth shape directly. So we chop the slice of disc into thousands of small bricks ("**hexahedral elements**"), store a temperature at each brick corner ("**node**"), and approximate the temperature inside each brick by blending its 8 corner values. The heat equation then becomes a big system of ordinary equations, one per node:

$$M\,\dot{T} + K(T)\,T = F(t)$$

- `M` ("mass/capacity matrix"): how much heat each node can store.
- `K` ("conductivity/stiffness matrix"): how strongly each node is connected to its neighbours.
- `F`: the heat coming in from the pads.

**The mesh generator is parameterized by the 24 thicknesses.** Change a thickness, and the bricks stretch or squash. That is how the design variables reach the physics.

**Sparse matrices.** Each node only touches its neighbours, so `K` and `M` are almost all zeros. Storing only the non-zeros (`scipy.sparse`) makes solves roughly 10× faster. The handout requires this, and we must write the solver ourselves (no ANSYS or FEniCS for the graded part).

### 4.3 Stepping through time (implicit Euler)

The braking run lasts a while, so we march forward in small time steps Δt: given the temperatures now, compute the temperatures a moment later. "Implicit Euler" means each step solves for the *new* temperatures using the *new* conditions, which is stable even with large steps. Each time step is one solve; the whole braking run is a few hundred of them.

### 4.4 The nonlinear part: why Newton's method is needed

**What makes it nonlinear.** Two things:
1. **Radiation** goes like `T⁴`. Double the (absolute) temperature, sixteen times the radiated heat.
2. **Conductivity changes with temperature**: hot cast iron conducts heat worse than cold, so `k` depends on `T`.

This means we can't solve each time step in one shot: the answer depends on itself.

**Newton's method, in plain words.** Guess the temperatures, see how wrong the energy balance is (the "**residual**"), use the slope of that error (the "**tangent matrix**") to make a much better guess, repeat. Near the answer, the number of correct digits roughly **doubles each iteration** ("quadratic convergence"). That doubling is something we must show in the verification ladder because it proves our tangent matrix is right.

```
   iteration   residual (how wrong)
       1        1e+02
       2        3e-01
       3        2e-05     ← digits roughly double each time
       4        1e-12     ← done
```

**Why it matters for the goal.** The handout *requires* a nonlinear component. Ours is type **N2** (the PDE itself is nonlinear, solved with Newton). Doing N2 properly with exact gradients earns the **+5 "fully nonlinear" bonus**. We must also report the optimized design **with and without** the nonlinearity and explain the difference (for example: ignoring radiation probably makes the disc look hotter than it is, so the optimizer would leave in iron it doesn't need).

### 4.5 Scoring a design: smooth "maximum" (KS aggregation)

**Problem.** "Peak temperature" means the maximum over thousands of nodes and hundreds of time steps. The plain max is not smooth (which point is hottest can jump), and optimizers need smooth functions.

**Fix.** The **KS function** (Kreisselmeier–Steinhauser) is a "soft maximum": a smooth formula that is always slightly *above* the true max. Using it as the constraint is safe (conservative), and it is smooth enough to take gradients. The report must show the true max next to the KS value.

**Analogy.** Instead of asking "who is the tallest person in the room?" (answer jumps when someone stands up), you compute a weighted average that is dominated by the tallest few people. Almost the same answer, but it changes gradually.

### 4.6 Verification: proving the solver is right before trusting it

**Why.** A solver written with a coding agent "is wrong until proven otherwise" (handout). We prove it by running problems where we already know the exact answer and showing our code gets it. The handout requires **at least six tests, in a table**: our number, the exact number, the relative error.

Planned ladder (from the project plan):

| Test | What we check | Exact answer from | Catches |
|---|---|---|---|
| 1D slab | Heat through a flat plate | Textbook formula | Wrong matrices, signs, units |
| Manufactured solution, mesh refinement | Pick a temperature field, compute the heat source that makes it exact, check we recover it; error should shrink at the expected rate as bricks get smaller | We invent it | Wrong discretization or boundary conditions |
| Manufactured solution, time-step refinement | Same, shrinking Δt | We invent it | Wrong time stepping |
| Energy balance | Heat in from pads = heat stored + heat lost to air | Conservation of energy | Lost or double-counted heat |
| Newton convergence | Residual digits double each iteration | Theory | Wrong tangent matrix |
| Nonlinear benchmark | A problem with radiation or k(T) whose answer is known | Textbook / exact | Wrong nonlinear terms (needed for the +5 bonus) |

**Important rule from the handout:** the expected answers must be derived by hand or from a textbook *before* the agent sees the test, so the code can't just "agree with itself".

**Status:** *Not started yet. See ainote.md.*

### 4.7 Gradients: "which way should each knob turn?"

**What a gradient is.** For each of the 24 thicknesses, the gradient says: *if I make this thickness slightly bigger, how much does the peak temperature (or mass) change?* It's a list of 24 sensitivities, one per knob, per score.

**The slow way (finite differences).** Nudge knob 1, rerun the whole braking simulation, see what changed. Repeat for all 24 knobs. That's 24–48 full simulations for *one* gradient, every optimizer iteration. Too slow.

**The clever way (the adjoint method).** Run the simulation once forward. Then run **one extra "backward" calculation** per score, and out come the sensitivities to *all 24 knobs at once*. Cost barely depends on how many knobs there are.

**Analogy.** You're at a party and want to know how much each of 24 guests contributed to the noise level. Slow way: ask each guest to leave one at a time and measure. Adjoint way: record the room once, then play the recording backward through a model of the room's acoustics, and it tells you every guest's contribution in one pass.

**Why "backward in time".** Our problem is transient: the temperature at the end of stop 3 depends on what happened during stops 1 and 2. So the adjoint calculation starts at the end of the braking run and sweeps backward to the start, asking at each step "how much did this moment's temperatures influence the final score?". That backward sweep is exactly what earns the **+5 "path or time dependence" bonus**.

**Connection to the KKT notes.** The adjoint variable `λ` is the Lagrange multiplier on the heat equation. Implicit differentiation, the adjoint method, and "Lagrange multiplier of the state equation" are three names for the same idea.

**How we prove the gradients are right.**
- **Gradient-check table:** adjoint answer vs finite-difference answer, for every constraint and both kinds of variable (vane and cheek thicknesses). They should agree to many digits.
- **Taylor test plot:** shrink the finite-difference step ε; the mismatch should fall like ε² until round-off noise takes over. A flat line means a bug.
- **Cost comparison:** time the adjoint vs finite differences on our model (rubric asks for it).

**Status:** *Not started yet. See ainote.md.*

### 4.8 The optimizer (SLSQP)

**What it does.** SLSQP (from `scipy.optimize`) is a standard gradient-based optimizer. At each iteration it takes the mass, the temperature constraints, and their gradients, builds a simplified local model of the problem, and jumps to the best design of that model. Repeat until it can't improve without breaking a limit.

**Analogy.** Hiking downhill in fog toward the lowest point of a valley, with fences you may not cross. The gradient tells you which way is downhill; the constraints are the fences. You end up at the lowest spot you can reach, usually pressed against one or two fences.

**Practical rules from the handout that we follow:**
- **Scaling:** work with `x = d / d₀` and constraints written as `T/T_limit − 1 ≤ 0`, so every number is around 1. Optimizers behave badly when one number is 0.003 and another is 700.
- **At least two starting designs**, because the answer might depend on where you start (local optima). Report the spread.
- **Save the history** of every iteration (design, mass, temperatures) for the convergence plot and the website animation.

### 4.9 Reading the answer: KKT conditions and multipliers

At the optimum, some constraints are **active** (pressed against their limit) and some are **inactive** (comfortably satisfied). Each active constraint has a **Lagrange multiplier**, which has a concrete meaning:

> **The multiplier = how many grams of disc we would save if that limit were relaxed by one unit.**

So, for example, if the face-temperature multiplier works out to *TBD* g per °C, then a pad compound that tolerates 10 °C more would let us save about 10 × that much iron. The rubric wants exactly this kind of engineering interpretation.

**Status:** *Not started yet. See ainote.md.*

### 4.10 Results file and website

**`results.json`.** One file, written by the solver, holding every number the report shows and the geometry/temperatures the 3D viewer draws. Nothing is copied by hand, so if the solver changes, one command regenerates everything (that is the **Reproducibility 5 points**).

**The website.** A GitHub Pages site built from Markdown (Jekyll), with math rendered by MathJax and a three.js 3D viewer. The viewer must satisfy all 8 handout requirements:

1. Initial and optimized discs side by side with linked cameras.
2. True 3D geometry (the full disc built from the copied slice).
3. Colored by temperature with a colour bar, single-hue scale (no rainbow), over-limit regions marked by more than colour alone.
4. The temperature field drawn on the geometry, interpolated properly.
5. One animation: the disc heating up and cooling through the braking stops.
6. Hover/click to inspect an element and see its temperature.
7. Metrics panel per design: mass and every constraint, each marked pass/fail.
8. Loads from `results.json`, works at phone width, renders in seconds.

The report itself has the 10 sections the handout lists (Summary, Formulation, Analysis model, Verification, Sensitivities, Optimization, Nonlinear component, Results table, Bonus, Limitations + AI-use statement).

**Status:** *Not started yet. See ainote.md.*

---

## 5. Glossary (one line each)

| Term | Meaning |
|---|---|
| **PDE** | Partial differential equation; here, the heat equation |
| **State** `u` / `T` | The physics answer for a given design: temperatures everywhere, every time step |
| **Design variables** `d` | The 24 thicknesses the optimizer may change |
| **Functional** | A single number computed from the whole temperature field (e.g. peak temperature) |
| **Finite element / mesh** | Chopping the shape into small bricks so a computer can solve the PDE |
| **Node** | A brick corner where a temperature is stored |
| **Residual** | How far a guessed solution is from satisfying the equations |
| **Tangent matrix** | The slope of the residual; what Newton's method uses to improve the guess |
| **Newton's method** | Repeatedly fix a guess using the slope; digits double each step near the answer |
| **Implicit Euler** | A stable way to step through time |
| **Cyclic symmetry** | Simulating one identical slice instead of the whole disc |
| **Convection / radiation** | Heat lost to moving air / emitted as glow (∝ T⁴) |
| **KS function** | A smooth, slightly conservative version of "maximum" |
| **Gradient / sensitivity** | How much a score changes per small change in each design variable |
| **Finite differences** | Getting a gradient by nudging each variable and rerunning: simple but slow |
| **Adjoint method** | Getting the whole gradient with one extra backward solve per score |
| **Taylor test** | Plot proving gradients are correct (error should fall like ε²) |
| **SLSQP** | The optimizer we use (Sequential Least SQuares Programming) |
| **Active constraint** | A limit the optimum is pressed against |
| **Lagrange multiplier** | How much the objective would improve per unit relaxation of a limit |
| **KKT conditions** | The mathematical conditions that hold at a constrained optimum |
| **N2** | The handout's label for "the PDE itself is nonlinear" |
| **results.json** | The single file all report numbers and the viewer come from |

---

## 6. How the pieces map to the rubric

| Rubric category | Points | Where it comes from in this project | Section |
|---|---|---|---|
| Structure and formulation | 10 | Brake disc motivation, 24 variables, mass objective, cited temperature limits | 2 |
| Analysis code and verification | 20 | Our own FE heat solver + the ladder table | 4.1–4.6 |
| Sensitivity analysis | 20 | Transient adjoint, gradient-check table, Taylor plot, cost comparison | 4.7 |
| Optimization | 15 | SLSQP, scaling, ≥2 starts, active constraints and multipliers interpreted | 4.8–4.9 |
| Nonlinear component | 15 | Radiation + k(T), Newton; design with vs without | 4.4 |
| Website and 3D comparison | 10 | Pages site, 8 viewer requirements, math renders | 4.10 |
| Reproducibility | 5 | One command regenerates results.json | 4.10 |
| Presentation and clarity | 5 | Readable report (and this guide helps you explain it) | all |
| Bonus: fully nonlinear | +5 | N2 with Newton + implicit gradients + nonlinear benchmark | 4.4, 4.6 |
| Bonus: time dependence | +5 | Adjoint swept backward through the braking stops, checked by FD | 4.7 |
| Bonus: uncertainty (optional) | +5 | Several braking scenarios (initial speed / brake bias) | stretch goal |

Bonus is capped at +15.

---

## 7. Things to watch / open risks

- **Structure approval.** The handout says teams proposing their own structure needed instructor approval by Nov 4. A brake disc is close to the menu's "heat sink" (vanes act like fins), but it is not on the menu. Worth a short email to the instructor.
- **Deadline.** Final site due Nov 25; the commit at the deadline is graded.
- **The AI-use statement is required** (report section 10): what the coding agent wrote, which parts you derived or checked yourself, and **one bug the verification caught**. See section 8 below.

---

## 8. Helping you write the AI-use statement

The handout requires, in section 10 of the report: *"what the coding agent wrote, which parts you derived or checked yourselves, and one bug the verification caught."* It also says to *"record what the agent wrote and how you checked it"* as you go.

A skeleton to fill in as the project progresses (each workstream adds facts to `ainote.md`; this section collects them):

> **AI-use statement.** The solver, adjoint, optimization scripts, and website were written with a coding agent (Claude Code). *[Your part: e.g. "I chose the structure, design variables, objective, and constraint limits and their sources; derived the exact solutions used in the verification ladder by hand before the agent saw the tests; reviewed every verification and gradient-check result."]*
>
> | Component | Written by | How it was checked |
> |---|---|---|
> | Mesh generator | Agent | *TBD* |
> | Heat solver + Newton | Agent | Verification ladder (Section 4) |
> | Adjoint gradients | Agent | Gradient-check table + Taylor test (Section 5) |
> | Optimization | Agent | Constraint values re-checked at the optimum; two starting designs |
> | 3D viewer | Agent | Manual check against the 8 requirements |
>
> **A bug the verification caught:** *TBD: filled in from ainote.md once one turns up.*

Be honest about which parts you checked yourself; that's what is graded.

---

## 9. Change log of this guide

- 2026-11-24: First draft from the project plan, before any code exists.
