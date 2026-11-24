# The beginner's guide to the busbar project

*For someone who only half-follows the lectures and wants to know what was built, why, and how each piece earns points.*

> **Status of this guide (2026-11-24).** Written from the busbar code, its notes in [`ainote.md`](ainote.md) and its full-run `results.json`. The website always reads its numbers from `results.json`, so if this guide and the website ever disagree, **the website is right** and this guide is stale. The earlier brake disc project and its guide are kept in [`alternate_topic/`](alternate_topic/).

---

## 0. The whole project in one paragraph

We designed a **busbar**: the thick copper strip that carries current from one battery module to the next inside an electric car's battery pack. A pack has dozens of them, and copper is heavy and expensive, so we want each one **as light as possible**. But a thinner bar has more electrical resistance, so it **heats up** when the car pulls a big burst of current (full acceleration, or the peak of a fast charge). If it gets too hot, three things go wrong: the plastic insulation around the bar is damaged, heat flows through the bolted ends into the **battery cells**, which age fast when hot, and the extra resistance **wastes energy**. So the question is: *where can we remove copper, and where must we keep it?* We wrote a program that (1) simulates how current flows through the bar and how the heat it makes spreads and escapes during a current burst, (2) works out how every one of 24 dimensions affects the temperatures and the resistance, all at once, and (3) gives that information to an optimizer, which keeps removing copper until a limit is just barely met.

**The headline result:** the bar went from **403 g to 351 g (13 % less copper)**. And the standard 30 × 5 mm bar we started from actually **fails** two of the three requirements, while the optimized one passes all three. If you do the usual hand calculation, which ignores that hot copper conducts electricity worse, you get an even lighter 322 g bar, but it really has **54.3 µΩ of resistance, over the 52 µΩ budget**. So the nonlinear physics is what keeps the design honest.

---

## 1. The object: what a busbar looks like

```
   Battery module A                                  Battery module B
   ┌────────────┐                                    ┌────────────┐
   │   cells    │   bolted                  bolted   │   cells    │
   │  ┌──────┐  │   joint                   joint    │  ┌──────┐  │
   │  │ term │══╪═══╗══════════════════════════╔═════╪══│ term │  │
   │  └──────┘  │   ║    copper busbar, 300 mm ║     │  └──────┘  │
   └────────────┘   ╚══════════════════════════╝     └────────────┘
                    ◄──25 mm──►            ◄──25 mm──►
                    overlap (lap joint)    overlap

   Cross-section of the bar:  ┌─────────────── width w ───────────────┐
                              └───────────────────────────────────────┘  thickness t
```

- The bar is **300 mm** long, and each end overlaps a module terminal by 25 mm and is bolted to it (a "**lap joint**").
- Current enters through one joint, runs along the bar, and leaves through the other.
- **Where the heat comes from:** the copper itself (any resistance turns current into heat, called **Joule heating**), and the bolted joints, which have a small **contact resistance** (two metal faces pressed together never touch perfectly). At about 10 µΩ per joint, the contact resistance is comparable to the resistance of the whole bar.
- **Where the heat goes:** into the air around the bar (natural convection: warm air rises off it), by glowing (radiation), and through the joints into the module terminals, which are held at 40 °C.

**Symmetry trick.** The bar is mirror-symmetric across its width, so we only simulate **half of it** and mirror it for display. That halves the computation with no loss of accuracy.

---

## 2. The design problem, written like an engineer

The handout says every problem has the form *"minimize something, subject to limits"*. Ours:

| Piece | In plain words | In our project |
|---|---|---|
| **Design variables** `d` (the knobs we turn) | The 24 numbers the optimizer may change | The bar's **width** at 12 points along its length (12 to 50 mm) and its **thickness** at the same 12 points (1.5 to 8 mm). Starting design: the standard uniform 30 × 5 mm bar |
| **Objective** `f` (what we minimize) | The thing we want small | Copper mass |
| **Constraint 1** | Don't melt the insulation | Peak bar temperature ≤ 105 °C, the rating of the insulation sleeve (UL 224 105 °C class) |
| **Constraint 2** | Don't cook the battery cells | Peak temperature at the bolted joints ≤ 70 °C, because the joint sits on the cell terminal (Li-ion cell data sheets put the limit around 60 °C for the cells themselves) |
| **Constraint 3** | Don't waste energy | End-to-end resistance at the continuous current ≤ 52 µΩ (an efficiency budget) |
| **Bounds** | Keep it buildable | Width ≥ 12 mm leaves room for the bolt hole; width ≤ 50 mm is the space between terminals; the thickness range is what flat copper bar can be stamped and bent to |
| **The "state"** | What the physics tells us once the knobs are set | The electric voltage and the temperature at every point of the bar, at every moment of the current burst |

**Why 12 + 12?** The handout requires at least 20 design variables so that the clever gradient method matters. Letting width and thickness change along the bar lets the optimizer put copper exactly where heat has to be carried away.

**The current profile (the "load case").** The bar first carries **400 A continuously** (it settles to a steady temperature). Then comes a **60 s boost pulse at 1200 A** (full-power acceleration or a fast-charge peak), followed by **240 s back at 400 A** while it cools. The pack air and the terminals are at 40 °C. The peak temperatures happen at the end of the pulse.

**This is on the handout's menu.** "Busbar or resistive heater" is listed by name, so no instructor approval is needed.

---

## 3. The big picture: how the pieces connect

```
          ┌───────────────────────────────────────────────────────────────┐
          │                     THE OPTIMIZATION LOOP                      │
          │                                                                │
  start   │   24 widths and   ──► [1] BUILD MESH ──► [2] ELECTRO-THERMAL   │
  design ─┼─►  thicknesses d        (geometry)        SIMULATION: voltage  │
          │       ▲                                     and temperature    │
          │       │                                     through the pulse  │
          │       │                                     (~1.6 s)           │
          │       │                                         │              │
          │  [5] OPTIMIZER  ◄── mass, max temps, ◄── [3] SCORE THE DESIGN  │
          │   (SLSQP: "move      resistance              │                 │
          │    the 24 knobs   + their gradients          ▼                 │
          │    this way")          ▲          [4] ADJOINT: gradients of    │
          │                        └───────── every score w.r.t. all 24    │
          │                                   knobs at once (~3.8 s)       │
          └───────────────────────────────────────────────────────────────┘
                     │ when nothing improves any more
                     ▼
          [6] RESULTS FILE (results.json) ──► [7] WEBSITE with 3D viewer
```

Before any of this ran on the real bar, piece [2] had to pass the **verification ladder** (section 4.6). The handout's tip: *"A wrong gradient of a wrong model optimizes beautifully."*

Where each piece lives in the repository:

| File | What's in it | Rubric points it mainly earns |
|---|---|---|
| `busbar/params.py` | Every physical number, with its source | Formulation 10 |
| `busbar/mesh.py`, `model.py`, `solver.py` | Mesh, coupled electro-thermal equations, Newton | Analysis & verification 20, Nonlinear 15, fully-nonlinear bonus +5 |
| `busbar/ladder.py` | The 12 verification tests | Analysis & verification 20 |
| `busbar/problem.py`, `api.py`, `adjoint.py`, `ks.py`, `gradcheck.py`, `optimize.py` | Scores, adjoint, gradient checks, SLSQP (the adjoint, KS, gradient-check and optimizer code was carried over from the brake disc version) | Sensitivity 20, Optimization 15, time-dependence bonus +5, uncertainty bonus +5 |
| `scripts/run_all.py` | The one command that regenerates every number | Reproducibility 5 |
| `docs/` | The website: report (`index.md`), 3D viewer, `results.json` | Website 10, Presentation 5 |

---

## 4. Piece by piece

### 4.1 The physics: two equations that talk to each other

This project has **two** PDEs, and they depend on each other. This is called **two-way coupling**, or "multiphysics".

**Equation 1: how current flows.**

$$\nabla\cdot\big(\sigma(T)\,\nabla\phi\big) = 0$$

`φ` is the voltage at each point, and `σ` is the electrical conductivity of copper. In words: *current flows from high voltage to low voltage, and none is created or lost on the way*. Wherever the bar is thin, the current is squeezed through a smaller area, so it flows faster there and makes more heat.

**Equation 2: how heat moves.**

$$\rho c \frac{\partial T}{\partial t} = \nabla\cdot\big(k(T)\,\nabla T\big) + \sigma(T)\,\lvert\nabla\phi\rvert^2$$

In words: *(how fast a point warms up) × (how much heat it takes to warm it) = (heat conducted in from neighbours) + (heat made right there by the current)*. The last term is the **Joule heating**: it comes straight out of equation 1.

**Why they're coupled both ways.** The current makes heat (equation 1 feeds equation 2). But **hot copper conducts electricity worse**: σ drops about 0.4 % per degree (IEC 60028). So as the bar heats up, its resistance rises, which makes even more heat (equation 2 feeds back into equation 1). That loop is why we can't just solve one equation and then the other: we solve them **together** ("monolithically").

**Where the heat goes (the boundary conditions):**

```
   bolted joints  ──►  heat IN from contact resistance (half into the bar, half into the terminal)
                  ◄──  heat OUT into the 40 °C terminal (a "conductance")
   bar surfaces   ◄──  heat OUT by natural convection: h ∝ ΔT^(1/4) (Holman's formula for still air)
                  ◄──  heat OUT by radiation: εσ(T⁴ − T_air⁴), with ε = 0.85 for a coated bar
```

**Analogy.** Think of the bar as a crowded hallway. People (current) squeeze through narrow spots faster and bump into each other more (heat). The hotter the hallway gets, the more sluggish everyone becomes (resistance rises), which causes more bumping. The air around the hallway and the doors at each end (the joints) let the heat out.

### 4.2 Finite elements: turning the equations into numbers

**What it is.** A computer can't solve the equations over a smooth shape directly. So we chop the half-bar into small bricks ("**hexahedral elements**"), store a voltage and a temperature at each brick corner ("**node**"), and blend those corner values inside each brick. The PDEs then become a big set of ordinary equations, two per node.

**What we actually used:** 132 bricks and 276 nodes. The bricks stretch and squash as the optimizer changes the width and thickness, but the number of bricks and how they connect stay the same. That keeps the gradients smooth.

**Sparse matrices and our own solver.** Each node only touches its neighbours, so the matrices are almost all zeros; we store only the non-zeros (`scipy.sparse`), as the handout requires. The solver is ours, not ANSYS or FEniCS.

### 4.3 Stepping through time (implicit Euler)

The run starts with one very long step that brings the bar to its steady 400 A temperature. Then it marches through the pulse in 2 s steps and through the recovery in 6 s steps. "Implicit Euler" means each step solves for the *new* state using the *new* conditions, which stays stable. One full run takes about **1.6 seconds**.

### 4.4 The nonlinear part: why Newton's method is needed

**What makes it nonlinear.** Four things, all in the model:
1. **σ(T)**: hot copper conducts electricity worse. This is what couples the two equations both ways.
2. **Joule heating** goes like the square of the voltage gradient.
3. **Natural convection** gets stronger as the bar gets hotter (h ∝ ΔT^¼, so the heat loss goes like ΔT^5/4).
4. **Radiation** goes like `T⁴`. (Also a small one: copper's thermal conductivity k drops slightly with temperature.)

So each step can't be solved in one shot: the answer depends on itself.

**Newton's method, in plain words.** Guess the voltages and temperatures, measure how far off the equations are (the "**residual**"), use the slope of that error (the "**tangent matrix**") to make a much better guess, and repeat. Near the answer the number of correct digits roughly **doubles each iteration** ("quadratic convergence"). Our solver shows this on the pulse step (ladder test V9, order 2.1), which proves the tangent matrix is right. We also checked the tangent matrix directly against a numerical derivative (V11).

**What we found by switching the nonlinearity off.** The handout asks for the optimized design **with and without** the nonlinear part. "Without" here means the usual **hand calculation**: copper properties at their 20 °C data-sheet values, a fixed convection coefficient, and simplified radiation.

| | Optimized with full physics | Optimized with the hand-calculation model |
|---|---|---|
| Mass | 351.2 g | 322.1 g |
| Resistance, re-checked with full physics | 50.5 µΩ ✓ | **54.3 µΩ ✗ (over the 52 µΩ budget)** |
| Joint peak, re-checked with full physics | 66.0 °C ✓ | 69.4 °C ✓ |
| Bar peak, re-checked with full physics | 73.7 °C ✓ | 79.6 °C ✓ |

**Why the difference:** at 20 °C copper conducts better than it does at 70 °C, so the hand calculation thinks the bar has less resistance and makes less heat than it really does. We checked which nonlinearity is responsible by switching them off one at a time: **freezing σ at its 20 °C value alone** reproduces almost exactly the hand-calculation answer (50.0 µΩ instead of 54.3 µΩ). So the temperature-dependent conductivity σ(T) is what changes the design, and ignoring it gives a bar that fails its resistance budget.

Doing N2 properly, with Newton's method, exact gradients and a nonlinear benchmark (V4: a coupled σ(T) bar with an exact formula), earns the **+5 "fully nonlinear" bonus**.

### 4.5 Scoring a design: a smooth "maximum" (KS aggregation)

**Problem.** "Peak temperature" means the maximum over every node and every time step. The plain max is not smooth (which point is hottest can jump), and optimizers need smooth functions.

**Fix.** The **KS function** is a "soft maximum": a smooth formula that always sits slightly *above* the true max. Using it as the constraint is safe (conservative), and it is smooth enough to take gradients of. At the optimum, the joint's KS value is exactly 70.0 °C while the true joint peak is 66.0 °C, so the soft max is about 4 °C conservative. The report shows both.

**Analogy.** Instead of asking "who is the tallest person in the room?" (the answer jumps when someone stands up), you compute a weighted average dominated by the tallest few. It gives almost the same answer, but it changes gradually.

The resistance constraint doesn't need KS: it is already one smooth number (voltage drop ÷ current).

### 4.6 Verification: proving the solver is right before trusting it

**Why.** A solver written with a coding agent "is wrong until proven otherwise" (handout). We prove it by running problems where we know the exact answer. The handout requires at least six tests in a table. We have **12 tests (17 rows, since some have two parts), and all pass** (`python -m busbar.ladder`):

| ID | Test, in plain words | Result |
|----|------|--------|
| V1 | Resistance of a plain uniform bar with its two joints, against the textbook formula | Error 2e-16 (round-off) |
| V2 | The same tapered bar, rotated and moved to a random spot in 3D. Nothing should change | Identical to 7e-14 |
| V3 | A bar heated by its own current, cooled at the ends: the exact temperature profile is a parabola | Error 1e-10 |
| V4 | **Nonlinear benchmark:** the two-way coupled σ(T) bar, against its exact cosine-shaped profile | Error 9e-6; error shrinks 4× per mesh halving (order 1.99) |
| V5 | A lump of copper cooling by natural convection (h ∝ ΔT^¼), against the exact cooling curve | Error 3e-7; time-step order 1.00 |
| V6 | A lump cooling only by radiation (T⁴), against the exact curve | Error 4e-7; time-step order 1.00 |
| V7 | "Manufactured solution" in 3D on a tapered mesh, with every nonlinear term on: pick the answer first, work out the sources that produce it, check the solver recovers it | Order 1.96 for temperature, 1.94 for voltage (textbook: 2) |
| V8 | Energy bookkeeping: electrical power in = Joule heat made; heat made = heat stored + heat lost | Exact (2e-16 and 3e-16) |
| V9 | Newton converges quadratically on the pulse step | Order 2.1 |
| V10 | Mesh refinement of the peak temperature (1×, 2×, 3× finer) | Optimizer mesh within 0.04 °C (0.1 %) of the converged value |
| V11 | The Newton tangent matrix against numerical derivatives | Agree to 4e-16 and 2e-9 |
| V12 | Every outside face of the mesh belongs to exactly one boundary type | 232 = 232 |

**Important rule from the handout:** the expected answers must be derived by hand or from a textbook *before* the agent sees the test, so that the code can't just "agree with itself".

### 4.7 Gradients: "which way should each knob turn?"

**What a gradient is.** For each of the 24 widths and thicknesses, the gradient says: *if I make this a bit bigger, how much does the peak joint temperature (or the resistance, or the mass) change?* It is a list of 24 sensitivities per score.

**The slow way (finite differences).** Nudge knob 1, rerun the whole simulation, see what changed, and repeat for all 24 knobs. On our model that takes about **75 seconds** per gradient.

**The clever way (the adjoint method).** Run the simulation once forward. Then run **one backward pass**, and out come the sensitivities of all the scores to *all 24 knobs at once*. On our model it takes **3.8 seconds, about 14× faster**. The gap grows with the number of knobs, because the adjoint's cost doesn't depend on how many there are.

**Analogy.** You want to know how much each of 24 guests at a party added to the noise. The slow way: ask each guest to leave, one at a time, and measure. The adjoint way: record the room once, then play the recording backward through a model of the room's acoustics, and it tells you every guest's contribution in one pass.

**Why "backward in time".** The temperature at the end of the pulse depends on everything that happened before it. So the adjoint starts at the end of the run and sweeps backward to the start, asking at each step "how much did this moment influence the peak?". That backward sweep earns the **+5 "path or time dependence" bonus**.

**Connection to the KKT notes.** The adjoint variable `λ` is the Lagrange multiplier on the PDE. "Implicit differentiation", "the adjoint method" and "the Lagrange multiplier of the state equation" are three names for the same idea.

**How we proved the gradients are right.**
- **Gradient-check table:** 48 checks of adjoint against brute force (central finite differences), covering mass, both temperatures and the resistance, on both widths and thicknesses. The worst disagreement is **4×10⁻¹⁰** of the gradient size: about 10 matching digits.
- **Taylor test plot:** as the finite-difference step ε shrinks, the mismatch falls like ε² until round-off takes over. A flat line would mean a bug.
- **Cost comparison:** forward 1.6 s, adjoint 3.8 s, finite differences 75 s.

### 4.8 The optimizer (SLSQP)

**What it does.** SLSQP (from `scipy.optimize`) is a standard gradient-based optimizer. At each iteration it takes the mass, the constraints and their gradients, builds a simplified local model of the problem, and jumps to the best design of that model. It repeats until it can't improve without breaking a limit.

**Analogy.** Hiking downhill in fog toward the lowest point of a valley, with fences you may not cross. The gradient tells you which way is downhill; the constraints are the fences.

**What we did, following the handout's rules:**
- **Scaling:** the optimizer works with `x = d / d₀`, mass divided by the starting mass, and every constraint written as `value/limit − 1 ≤ 0`, so all numbers are around 1.
- **Three starting designs:** the standard bar, everything at its maximum, and everything halfway. **All three ended at 351.2 g**, essentially identical (32 to 37 iterations each). That strongly suggests it's the real best design, not a local dip.
- **History saved** for every iteration, for the convergence plot.

**What the optimized bar looks like, and why it makes sense:**
- **Every width went to the 50 mm maximum.** For the same amount of copper, a wide, thin bar has more surface for the air to cool, so it runs cooler than a narrow, thick one.
- **The thickness varies along the bar:** about **3.1 mm just inside each joint**, thinning to **2.3 mm in the middle**. Next to the joints, the bar acts as a **cooling fin** for the joint: it carries the contact heat away from the cell terminal into the span, where the air can take it. That copper also stores heat right where the pulse would heat the joint fastest.
- **At the very ends** (inside the bolted overlap) the thickness is at the **1.5 mm minimum**, because the terminal underneath does the cooling there.

### 4.9 Reading the answer: KKT conditions and multipliers

At the optimum, a constraint is **active** if the design is pressed against its limit, and **inactive** if it is comfortably met.

| Constraint | Starting 30 × 5 mm bar | Optimized bar | Status at the optimum |
|---|---|---|---|
| Bar ≤ 105 °C | 72.4 °C ✓ | 73.7 °C ✓ | Inactive |
| Joint ≤ 70 °C | **71.7 °C ✗** | 66.0 °C (KS 70.0 °C) ✓ | **Active** |
| Resistance ≤ 52 µΩ | **54.5 µΩ ✗** | 50.5 µΩ ✓ | Inactive |

Notice that the **lighter bar is also better**: the standard bar fails two requirements, and the optimizer fixes both while removing copper. It does this by moving copper to where it helps (wider, and thicker next to the joints).

**Why only the joint limit is active.** Copper conducts heat so well that the whole bar stays within a few degrees of the joint, so the 70 °C joint limit is always reached long before the 105 °C insulation limit.

Each active constraint has a **Lagrange multiplier**, which has a concrete meaning:

> **The multiplier is how much copper we would save if that limit were relaxed by one unit.**

Ours is **0.0106 kg per °C**, about **11 g per degree**. So if the cells could tolerate a joint at 75 °C instead of 70 °C, each busbar could be about **53 g lighter**. Across the dozens of busbars in a pack, that adds up. The insulation and resistance limits have zero multipliers: relaxing them would save nothing.

### 4.10 Bonus: design under uncertainty

Real packs don't always see exactly the nominal conditions. We also designed for three scenarios at once: nominal, **10 % more peak current** (1320 A), and **pack air and terminals 5 °C warmer**.

- The nominal 351 g optimum would reach **70.6 °C** at the joint with the extra current and **71.3 °C** in the warmer pack. So it fails both off-nominal cases.
- The **robust design** passes all nine constraints (three per scenario) at **420.8 g**, which is 70 g heavier than the nominal optimum. That is the price of the safety margin.
- The constraint that decides it is the joint temperature in the warmer-pack scenario, with a multiplier of 0.017 kg per °C.

This earns the optional **+5 "design under uncertainty" bonus**.

### 4.11 Results file and website

**`results.json`.** One file, written only by `python scripts/run_all.py`, holding every number the report shows and the geometry and fields that the 3D viewer draws. **No number is typed into the report by hand**: when the page opens, a script reads `results.json` and fills in every number, table and plot. So one command (about 25 minutes) regenerates the whole site. That earns the **Reproducibility 5 points**.

**The website.** The same GitHub Pages setup as before (Markdown with Jekyll, MathJax for math, a three.js 3D viewer), live at <https://sgarteiz9090.github.io/Project4/>. The viewer shows the initial and optimized bars side by side, coloured by temperature, with an animation through the current pulse, hover-to-inspect, and a pass/fail metrics panel, and it was checked at desktop and phone width. The report has the 10 sections the handout lists.

---

## 5. Glossary (one line each)

| Term | Meaning |
|---|---|
| **Busbar** | A solid copper strip that connects battery modules and carries the pack current |
| **Lap joint** | Where the bar's end overlaps the terminal and is bolted to it |
| **Contact resistance** | The extra resistance where two metal faces are pressed together |
| **Joule heating** | Heat made by current flowing through a resistance (P = I²R) |
| **σ(T)** | Electrical conductivity that drops as copper heats up |
| **Two-way coupling** | Current makes heat, and heat changes how the current flows |
| **PDE** | Partial differential equation; here, the current equation and the heat equation |
| **State** | The physics answer for a given design: voltage and temperature everywhere, at every time step |
| **Design variables** `d` | The 24 widths and thicknesses the optimizer may change |
| **Finite element / mesh** | Chopping the shape into small bricks so that a computer can solve the PDEs |
| **Residual** | How far a guess is from satisfying the equations |
| **Tangent matrix** | The slope of the residual; what Newton's method uses to improve the guess |
| **Newton's method** | Repeatedly fixing a guess using the slope; the digits double each step near the answer |
| **Implicit Euler** | A stable way to step through time |
| **Natural convection** | Heat carried away by air that rises off a warm surface on its own |
| **Manufactured solution** | A test where you pick the answer first and work out the problem that produces it |
| **Convergence order** | How fast the error shrinks as the mesh or time step shrinks (2 means 4× less error per halving) |
| **KS function** | A smooth, slightly conservative version of "maximum" |
| **Gradient / sensitivity** | How much a score changes per small change in each design variable |
| **Adjoint method** | Getting the whole gradient with one extra backward pass |
| **Taylor test** | A plot proving the gradients are correct (the error should fall like ε²) |
| **SLSQP** | The optimizer we use (Sequential Least SQuares Programming) |
| **Active constraint** | A limit the optimum is pressed against |
| **Lagrange multiplier** | How much the objective would improve per unit relaxation of a limit |
| **N2** | The handout's label for "the PDE itself is nonlinear" |
| **results.json** | The single file that all report numbers and the viewer come from |

---

## 6. How the pieces map to the rubric

| Rubric category | Points | Where it comes from in this project | Guide section |
|---|---|---|---|
| Structure and formulation | 10 | EV busbar motivation, 24 variables, mass objective, three cited limits | 2 |
| Analysis code and verification | 20 | Our own coupled electro-thermal FE solver + the 12-test ladder, all passing | 4.1–4.6 |
| Sensitivity analysis | 20 | Transient adjoint, 48 gradient checks (about 10 digits), Taylor plot, 14× cost comparison | 4.7 |
| Optimization | 15 | SLSQP, scaling, 3 starts that agree, active constraint and its multiplier interpreted | 4.8–4.9 |
| Nonlinear component | 15 | σ(T), Joule, convection, radiation with Newton; the hand-calculation design fails the resistance budget | 4.4 |
| Website and 3D comparison | 10 | Pages site, the 8 viewer requirements, math renders | 4.11 |
| Reproducibility | 5 | `python scripts/run_all.py` regenerates results.json and so every number | 4.11 |
| Presentation and clarity | 5 | Readable report (and this guide helps you explain it) | all |
| Bonus: fully nonlinear | +5 | N2 with Newton + implicit gradients + coupled nonlinear benchmark V4 | 4.4, 4.6 |
| Bonus: time dependence | +5 | Adjoint swept backward through the current pulse, checked by FD | 4.7 |
| Bonus: uncertainty | +5 | Three scenarios, robust 420.8 g design | 4.10 |

The bonus is capped at +15, and we target all three bonuses.

---

## 7. Things to watch / limitations

- **Deadline:** the final site is due Nov 25; the commit at the deadline is graded. Check the live site in a private browser window before submitting both links on Canvas. Math only renders on the live site, not on github.com.
- **Limitations worth stating in the report:** the joints are lumped (a contact resistance plus a conductance into a fixed-temperature terminal), not a full 3D model of the terminal; natural convection uses one textbook formula; radiation ignores neighbouring modules; the bar is straight, with no bends or bolt holes; there is one representative current pulse (the robust design covers three variants); the temperature limits are design values taken from the insulation rating and cell data sheets.

---

## 8. The AI-use statement (required, report section 10)

The handout requires: *"what the coding agent wrote, which parts you derived or checked yourselves, and one bug the verification caught."* Below is a draft built from what was recorded. **The part in brackets is yours to fill in honestly**: only you know what you checked yourself, and that is what gets graded.

> **AI-use statement.** The solver, adjoint, optimization scripts and website were written with a coding agent (Claude Code). The backward adjoint sweep, the KS functional, the gradient-check tools and the SLSQP driver were carried over from our brake disc version of this project; the coupled electro-thermal model, its mesh, the verification ladder and the site are new. *[Your part, for example: "I chose the structure, the design variables, the objective and the limits and their sources; I checked the ladder's exact solutions against the textbook; I read the verification and gradient-check results and the optimized design for physical sense."]*
>
> | Component | Written by | How it was checked |
> |---|---|---|
> | Mesh (`busbar/mesh.py`) | Agent | V2 (rotation invariance), V7 (manufactured solution on a tapered mesh), V12 (every boundary face in exactly one set) |
> | Coupled solver + Newton (`busbar/model.py`, `solver.py`) | Agent | Ladder V1–V10 against exact solutions; Newton quadratic convergence (V9); tangent against numerical derivatives (V11); energy balance to 1e-16 (V8) |
> | Adjoint gradients (`busbar/adjoint.py`, carried over from the brake disc) | Agent | 48 adjoint vs central finite-difference checks across all constraints and both variable types (about 10 digits); Taylor test |
> | Optimization (`busbar/optimize.py`, carried over from the brake disc) | Agent | Three starting designs reach the same 351.2 g; hand-calculation optimum re-checked with full physics |
> | Website and 3D viewer | Agent | Every number read from results.json; checked in a browser at desktop and phone width |
>
> **Bugs the verification caught:**
> 1. *A "correct" answer below absolute zero (caught by V4).* The first version of the coupled nonlinear benchmark used too much current. The solver converged, and it even agreed with the "exact" formula, but the temperatures were **below absolute zero**. At that current the problem has no physical solution: the bar heats faster than the joints can cool it (thermal runaway). Checking whether the answer made physical sense, not just whether it matched, caught it. The benchmark now runs at a current where a real solution exists.
> 2. *A "steady" step that wasn't quite steady (caught by V3).* The 1D Joule-heating test first missed the exact answer by 10⁻⁴, because the "steady" step used a time step of 10⁹ s, which still left a sliver of heat capacity in the equations. Using 10¹⁵ s in the test brought the agreement down to round-off.

---

## 9. Explaining it in 30 seconds

> "We designed a lighter copper busbar for an electric car's battery pack. Our own finite-element code simulates the current flowing through the bar and the heat it makes during a 1200 A boost pulse. Those two pieces of physics feed each other, because hot copper conducts electricity worse. The adjoint method tells us, in one backward pass, how each of the 24 widths and thicknesses affects the temperatures and the resistance, and we checked those gradients against brute force to 10 digits. The optimizer cut the bar from 403 to 351 grams. It made the bar wide and thin, and kept extra copper next to the bolted joints so it works as a cooling fin. The standard bar actually failed two requirements, and the optimized one passes all three. The limit that stops it is 70 °C at the joints, where the battery cells are, and each extra degree allowed would save about 11 grams. If you use the usual 20 °C hand calculation, you get a lighter bar that is actually over its resistance budget."

---

## 10. Change log of this guide

- 2026-11-24, first version for the busbar topic, from the full-run results on the busbar branch.
