---
layout: default
title: Race-car brake disc optimization
---
{% raw %}

<div id="mock-banner" class="mock-banner" hidden></div>

# Lighter race brake discs that still survive hard braking

<p class="lede">PDE-constrained design optimization of a vented cast-iron brake rotor: 24 thickness variables, a transient nonlinear 3D heat-conduction model solved with our own finite-element code, adjoint gradients swept backward in time, and SLSQP.</p>

## 1. Summary

We minimize the mass of a vented race-car brake rotor while keeping the rubbing faces below the grey-iron cracking and fade limit and the hub below the bearing-protection limit through a sequence of <span class="rj" data-param="n_stops" data-d="0"></span> hard stops from <span class="rj" data-param="v_high" data-scale="3.6" data-d="0"></span> to <span class="rj" data-param="v_low" data-scale="3.6" data-d="0"></span> km/h. The optimizer removes <span class="rj" data-calc="mass_saving_kg" data-d="2"></span> kg (<span class="rj" data-calc="mass_saving_pct" data-d="1"></span> %) of iron, going from <span class="rj" data-k="designs.initial.mass_kg" data-d="2"></span> kg to <span class="rj" data-k="designs.optimized.mass_kg" data-d="2"></span> kg. The peak friction-face temperature rises from <span class="rj" data-calc="face_peak_initial" data-d="0"></span> °C to <span class="rj" data-calc="face_peak_optimized" data-d="0"></span> °C against a limit of <span class="rj" data-calc="face_limit" data-d="0"></span> °C, and the peak hub temperature rises from <span class="rj" data-calc="hub_peak_initial" data-d="0"></span> °C to <span class="rj" data-calc="hub_peak_optimized" data-d="0"></span> °C against a limit of <span class="rj" data-calc="hub_limit" data-d="0"></span> °C. Both temperature constraints are active at the optimum: together they stop the optimizer.

<div class="kpis">
  <div class="kpi"><b><span class="rj" data-k="designs.optimized.mass_kg" data-d="2"></span> kg</b><span>optimized rotor mass (initial <span class="rj" data-k="designs.initial.mass_kg" data-d="2"></span> kg)</span></div>
  <div class="kpi"><b>−<span class="rj" data-calc="mass_saving_pct" data-d="1"></span> %</b><span>mass saved</span></div>
  <div class="kpi"><b><span class="rj" data-calc="face_peak_optimized" data-d="0"></span> °C</b><span>peak rubbing-face temperature (limit <span class="rj" data-calc="face_limit" data-d="0"></span> °C)</span></div>
  <div class="kpi"><b><span class="rj" data-calc="hub_peak_optimized" data-d="0"></span> °C</b><span>peak hub temperature (limit <span class="rj" data-calc="hub_limit" data-d="0"></span> °C)</span></div>
  <div class="kpi"><b><span class="rj" data-calc="adjoint_speedup" data-d="0"></span>×</b><span>adjoint speed-up over finite differences</span></div>
</div>

### Interactive 3D comparison

The left disc is the initial design and the right disc the optimized one, coloured by temperature on one shared scale. Press **Play** to watch both discs heat up during each stop and cool while the car accelerates again, or drag the timeline. Hover over (or tap) any face to see its element, surface, temperature and margin to its limit.

<div id="rotor-viewer" class="rv"></div>
<script type="module" src="assets/js/viewer.js"></script>

<div class="rj-plot" data-plot="peaks"></div>

## 2. Structure and formulation

### Why this part matters

A race car's front brake discs turn the car's kinetic energy into heat: a stop from <span class="rj" data-param="v_high" data-scale="3.6" data-d="0"></span> to <span class="rj" data-param="v_low" data-scale="3.6" data-d="0"></span> km/h dumps several hundred kilojoules into each front disc in about three seconds. The disc must store that heat and then shed it, mostly through the radial cooling channels between its two friction rings ("cheeks"), which are joined by vanes that pump air like a centrifugal fan. The disc is also **rotating, unsprung mass**: every gram costs acceleration, braking and grip over bumps. Too little iron and the rubbing faces overheat (pad fade, thermal cracking of grey cast iron) and heat soaks through the hub into the wheel bearing and brake fluid. The question is where iron can be removed and where it must stay.

### Design variables

The rotor has $n_v$ = <span class="rj" data-param="n_vanes" data-d="0"></span> straight radial vanes. At 12 equally spaced radial stations $r_1,\dots,r_{12}$ from the inner radius $r_i$ = <span class="rj" data-param="r_i" data-scale="1000" data-d="0"></span> mm to the outer radius $r_o$ = <span class="rj" data-param="r_o" data-scale="1000" data-d="0"></span> mm we choose

$$
\mathbf d = \big(t_{v,1},\dots,t_{v,12},\; t_{c,1},\dots,t_{c,12}\big),
$$

the vane thickness $t_v$ and the cheek (friction ring) thickness $t_c$, piecewise linear in $r$ between stations. Bounds: <span class="rj" data-param="tv_min" data-scale="1000" data-d="0"></span> mm $\le t_{v,j} \le$ <span class="rj" data-param="tv_max" data-scale="1000" data-d="0"></span> mm (casting limit; the channel must stay open) and <span class="rj" data-param="tc_min" data-scale="1000" data-d="0"></span> mm $\le t_{c,j} \le$ <span class="rj" data-param="tc_max" data-scale="1000" data-d="0"></span> mm (pad-wear allowance and strength).

<div class="rj-table" data-table="design_vars"></div>

### Objective

$$
\min_{\mathbf d}\; m(\mathbf d) = \rho \int_{\Omega(\mathbf d)} \mathrm d V \qquad\text{(full rotor mass, kg).}
$$

### Constraints

Both temperature constraints are functionals of the transient temperature field $T(\mathbf x,t;\mathbf d)$. Each one is aggregated over a surface node set $S$ and the whole braking history with the Kreisselmeier–Steinhauser (KS) function of the normalized nodal temperatures $a_{i,n} = T_{i,n}/T_{\lim}$ (in °C; node $i$, time step $n$):

$$
g = \mathrm{KS}_\rho(a) - 1 \le 0,
\qquad
\mathrm{KS}_\rho(a) = a_{\max} + \frac1\rho \ln \sum_{i\in S}\sum_{n=1}^{N} e^{\rho (a_{i,n} - a_{\max})} \;\ge\; a_{\max},
$$

with $S = \Gamma_{\text{fric}}$ (pad rubbing band) for $g_{\text{face}}$ and $S = \Gamma_{\text{hub}}$ (hub face $r=r_i$) for $g_{\text{hub}}$. We take $\rho = \rho_K T_{\lim}$ with $\rho_K = 1\ \mathrm K^{-1}$, so KS sits within a few kelvin above the true peak; the report shows both.

- $T_{\text{face}}^{\max}$ = <span class="rj" data-calc="face_limit" data-d="0"></span> °C on the pad rubbing band: above about 700 °C grey cast iron rotors suffer heat checking and cracking and pad friction fades (Limpert, *Brake Design and Safety*, 3rd ed., SAE 2011, ch. 3; Day, *Braking of Road Vehicles*, 2014, ch. 3).
- $T_{\text{hub}}^{\max}$ = <span class="rj" data-calc="hub_limit" data-d="0"></span> °C on the hub / bell mounting face, the heat path into the wheel bearing. Wheel-bearing grease of the NLGI GC grade must have a dropping point of at least 220 °C (ASTM D4950; NLGI GC-LB specification, table 10.3), and DOT 4 brake fluid a dry boiling point of at least 230 °C (FMVSS 116, 49 CFR 571.116 S5.1.1). The mounting face must therefore stay at or below the grease dropping point.
- Minimum thicknesses are the lower bounds above.

## 3. Analysis model

**Governing PDE.** Transient nonlinear heat conduction in the iron,

$$
\rho c\,\frac{\partial T}{\partial t} = \nabla\cdot\big(k(T)\,\nabla T\big) \quad\text{in } \Omega(\mathbf d),
\qquad k(T) = k_0\big(1+\beta (T-T_{\text{ref}})\big),
$$

with boundary conditions

$$
-k\,\frac{\partial T}{\partial n} =
\begin{cases}
-\,q_f(t) + h_f(t)(T-T_\infty) + \varepsilon\sigma\,(T^4-T_\infty^4) & \text{friction faces (pad band)}\\
h_c(t)\,(T - T_\infty) & \text{vane channel walls}\\
h_f(t)(T-T_\infty) + \varepsilon\sigma\,(T^4-T_\infty^4) & \text{free cheek face, outer rim}\\
h_{\text{hub}}\,(T - T_{\text{bell}}) & \text{hub face } r=r_i .
\end{cases}
$$

**Heat input.** During each stop the deceleration power $P(t) = m_{\text{car}}\, a\, v(t)$ is split by front bias, between the two front discs, and by the partition coefficient between disc and pad; the disc share is spread over the rubbing band as the flux $q_f(t)$.

**Convection.** Face coefficients $h_f(t)$ from the rotating-disc correlation (Dorfman; Limpert ch. 3), $\mathrm{Nu} = 0.015\,\mathrm{Re}_\omega^{0.8}$; channel coefficients $h_c(t)$ from the vane-channel correlation $\mathrm{Nu} = 0.023\,(1+(d_h/l)^{0.67})\,\mathrm{Re}^{0.8}\mathrm{Pr}^{1/3}$ with air speed from the vane pumping. Both follow the vehicle speed but not the temperature or the design.

**Discretization.** One vane pitch is reduced by two mirror symmetries (vane mid-plane, disc mid-plane) to a quarter sector, meshed with 8-node trilinear hexahedra in true cylindrical geometry (<span class="rj" data-calc="n_elems" data-d="0"></span> elements, <span class="rj" data-calc="n_nodes" data-d="0"></span> nodes). Implicit Euler in time with Newton's method in every step; the residual of step $n$ is

$$
\mathbf R_n = \mathbf M\,\frac{\mathbf T_n-\mathbf T_{n-1}}{\Delta t_n} + \mathbf K(\mathbf T_n)\,\mathbf T_n + \mathbf H_n\,(\mathbf T_n - T_\infty) + \mathbf r_{\text{rad}}(\mathbf T_n) - \mathbf q_n = \mathbf 0 .
$$

**Load case.** <span class="rj" data-param="n_stops" data-d="0"></span> stops from <span class="rj" data-param="v_high" data-scale="3.6" data-d="0"></span> to <span class="rj" data-param="v_low" data-scale="3.6" data-d="0"></span> km/h at <span class="rj" data-param="decel" data-scale="0.10194" data-d="2"></span> g, with <span class="rj" data-param="t_accel" data-d="0"></span> s of re-acceleration between them, starting from a disc at <span class="rj" data-param="T_init" data-offset="-273.15" data-d="0"></span> °C.

**Assumptions.** Axisymmetric (circumferentially averaged) pad heat input, since the disc turns many times per second; straight radial vanes of uniform height; no thermal contact resistance inside the casting; convection coefficients from correlations, not CFD; grey-iron properties as listed below; the pad, caliper and wheel are represented only by the partition coefficient and the hub conductance.

<details><summary>All model parameters (from results.json)</summary>
<div class="rj-table" data-table="params"></div>
</details>

## 4. Verification

Every test is run by `scripts/run_all.py` and reported straight from `results.json`.

<div class="rj-table" data-table="verification"></div>

<figure>
<div class="rj-plot" data-plot="mms"></div>
<figcaption>Manufactured-solution convergence: L2 error against element size (slope 2 expected for trilinear hexes) and against time step (slope 1 for implicit Euler).</figcaption>
</figure>

<figure>
<div class="rj-plot" data-plot="newton"></div>
<figcaption>Newton residual in one braking time step with k(T) and radiation: the number of correct digits roughly doubles every iteration (quadratic convergence).</figcaption>
</figure>

## 5. Sensitivities

**State equation.** Implicit Euler from $t_{n-1}$ to $t_n$ gives one nonlinear residual per step, $\mathbf R_n(\mathbf T_n,\mathbf T_{n-1},\mathbf d) = \mathbf 0$, $n = 1,\dots,N$ (section 3). $\mathbf T_0$ (<span class="rj" data-param="T_init" data-offset="-273.15" data-d="0"></span> °C everywhere) does not depend on $\mathbf d$. Newton's tangent $\mathbf A_n = \partial\mathbf R_n/\partial\mathbf T_n$ contains $k'(T)$ and $4\varepsilon\sigma T^3$; the coupling to the previous step is $\mathbf B_n = \partial\mathbf R_n/\partial\mathbf T_{n-1} = -\mathbf M/\Delta t_n$.

**KS partials.** $\partial g/\partial T_{i,n} = w_{i,n}/T_{\lim}$ with softmax weights $w_{i,n} = e^{\rho(a_{i,n}-\mathrm{KS})}$, $\sum w = 1$: only the hottest few node-instants carry weight. The mass $m(\mathbf d) = \rho_{\text{Fe}}\sum_e \lvert\Omega_e(\mathbf d)\rvert$ is an explicit function of the geometry, so $\mathrm dm/\mathrm d\mathbf d$ needs no adjoint.

**Discrete adjoint, swept backward through the braking history.** With $\mathcal L = g + \sum_{n=1}^{N}\boldsymbol\lambda_n^{\mathsf T}\mathbf R_n$, every $\mathbf R_n = \mathbf 0$ so $\mathrm dg/\mathrm d\mathbf d = \mathrm d\mathcal L/\mathrm d\mathbf d$ for any $\boldsymbol\lambda$. $\mathbf T_n$ appears in $g$, in $\mathbf R_n$ and in $\mathbf R_{n+1}$; choosing $\boldsymbol\lambda$ to remove every $\mathrm d\mathbf T_n/\mathrm d\mathbf d$ gives

$$
\mathbf A_n^{\mathsf T}\boldsymbol\lambda_n = -\Big(\frac{\partial g}{\partial \mathbf T_n}\Big)^{\mathsf T} - \mathbf B_{n+1}^{\mathsf T}\boldsymbol\lambda_{n+1},
\qquad \boldsymbol\lambda_{N+1} = \mathbf 0,\qquad n = N, N-1, \dots, 1,
$$

$$
\frac{\mathrm d g}{\mathrm d \mathbf d} = \frac{\partial g}{\partial \mathbf d} + \sum_{n=1}^{N} \boldsymbol\lambda_n^{\mathsf T}\,\frac{\partial \mathbf R_n}{\partial \mathbf d}.
$$

The sweep runs **backward in time**: the last step is solved first, because the temperature at step $n$ influences the peak through every later step. This is the path-dependent adjoint of the rubric bonus. It uses the stored forward history (small enough that no checkpointing is needed); each step costs one sparse LU of $\mathbf A_n^{\mathsf T}$, shared by both constraints as two right-hand sides. Because $\mathbf A_n$ is the exact Newton tangent at the converged state, this is implicit differentiation of the nonlinear solve (handout implementation iii).

**Shape term.** $\mathbf R_n$ depends on $\mathbf d$ only through the node coordinates $\mathbf X(\mathbf d)$, which the mesh generator moves smoothly at fixed topology. $\sum_n\boldsymbol\lambda_n^{\mathsf T}\partial\mathbf R_n/\partial d_j$ is evaluated by a complex step on the residual with the states held fixed, $\partial\mathbf R_n/\partial d_j = \operatorname{Im}\mathbf R_n(\mathbf T_n,\mathbf T_{n-1},\mathbf d + i h\mathbf e_j)/h$ with $h = 10^{-30}$: exact to round-off, 24 complex residual assemblies per step, and no extra PDE solves.

**How the gradients are checked.** (1) Adjoint vs. central finite differences ($h = 10^{-5} d_j$) for the mass and both KS constraints, for vane and cheek variables, with the full and the linearized model. (2) Taylor test along a random direction $\mathbf v$: $\big\lvert[g(\mathbf x+\epsilon\mathbf v)-g(\mathbf x-\epsilon\mathbf v)]/2\epsilon - \nabla g\cdot\mathbf v\big\rvert$ must fall as $\epsilon^2$ until round-off; a flat curve would mean a wrong gradient. (3) Cost of one forward solve plus one adjoint sweep against $2\times 24$ forward solves for central differences.

### Gradient check (adjoint vs. central finite differences)

<div class="rj-table" data-table="gradient_check"></div>

<figure>
<div class="rj-plot" data-plot="taylor"></div>
<figcaption>Taylor test along a random direction: the remainder of the central difference against the adjoint directional derivative falls as ε² (grey reference) until round-off takes over.</figcaption>
</figure>

### Cost: adjoint vs. finite differences

<div class="rj-table" data-table="cost"></div>

## 6. Optimization

**Algorithm and scaling.** SciPy SLSQP (`ftol` = $10^{-7}$). Variables are scaled $\mathbf x = \mathbf d/\mathbf d_0$, the objective is $f = m/m_0$ and the constraints are $g_k = \mathrm{KS}_k - 1 \le 0$, so everything is of order one, as SLSQP needs. The thickness limits are passed as bounds; the two KS constraints use the adjoint gradients of section 5. Three starts (baseline, all upper bounds, midpoint of the box) show whether the optimum depends on the start; the best feasible result is reported with the spread of the others. Here the final masses of all starts lie within <span class="rj" data-calc="start_spread_g" data-d="1"></span> g of each other.

<div class="rj-table" data-table="starts"></div>

<div class="rj-pair">
<div class="rj-plot" data-plot="history"></div>
<div class="rj-plot" data-plot="history_g"></div>
</div>

<figure>
<div class="rj-plot" data-plot="thickness"></div>
<figcaption>Thickness profiles: solid lines are cheeks, dashed lines vanes.</figcaption>
</figure>

### Active constraints and multipliers

<div class="rj-table" data-table="multipliers"></div>

**KKT read-out.** A constraint is active when $g_k > -10^{-3}$, a bound when the variable sits on it. The multipliers solve $\nabla f + \sum_k \mu_k\nabla g_k + \sum_{\text{bounds}}\nu_i(\pm\mathbf e_i) = \mathbf 0$ with $\mu,\nu\ge 0$ (non-negative least squares), next to SLSQP's own estimates. In scaled form $\mu_k$ is the drop in $m/m_0$ per unit relaxation of $g_k$ (section 11 of the KKT notes), so $\mu_k\,m_0/T_{\lim,k}$ is the **kilograms of iron saved per kelvin of extra allowed temperature**.

## 7. Nonlinear component (N2)

The model has two nonlinear terms: temperature-dependent conductivity $k(T)$ of grey cast iron and surface radiation $\varepsilon\sigma(T^4-T_\infty^4)$. Both enter the residual, the Newton tangent and the adjoint. The table compares the design optimized with the full nonlinear model against the same problem optimized with both terms linearized: radiation linearized about $T_{\text{init}}$ and $k$ frozen at $k(T_{\text{init}})$. That linear optimum is then re-evaluated with the full nonlinear physics, which shows what ignoring the nonlinearity would have cost.

<div class="rj-table" data-table="nonlinear"></div>

<div class="rj-table" data-table="n2"></div>

The design optimized without the nonlinear terms differs in mass by <span class="rj" data-calc="linear_mass_delta_pct" data-d="1"></span> % from the nonlinear optimum.

## 8. Results

<div class="rj-table" data-table="results"></div>

## 9. Bonus sections

- **Fully nonlinear (+5):** section 7; Newton per implicit-Euler step with the exact tangent and the same tangent in the adjoint.
- **Path / time dependence (+5):** the constraints aggregate the whole braking history and the adjoint is swept backward through it (section 5).
- **Design under uncertainty (+5):** three sampled scenarios (nominal, +20 km/h initial speed, +0.05 front brake bias). The robust design satisfies both temperature constraints in every scenario (six constraints, each with its own exact adjoint gradient), i.e. the worst case over the samples.

<div class="rj-table" data-table="robust"></div>

## 10. Limitations and AI-use statement

**Limitations.** Convection coefficients come from correlations, not from a CFD model of the vane channels. The pad heat is treated as circumferentially uniform. Thermal stresses and cracking are represented only through the temperature limit. The nominal design is sized for one braking scenario (the robust design in section 9 covers three). The optimizer's mesh reads the peak face temperature slightly below the mesh-converged (Richardson) value (verification V9 in the table of section 4), so the optimum is slightly unconservative.

**AI-use statement.** The solver, adjoint, optimization scripts and this website were written with a coding agent (Claude Code), working in parallel workstreams that coordinated through a shared status file (`ainote.md`).

<p class="callout"><b>[Sivan: replace this box with what you derived or checked yourself, for example: chose the structure, variables, objective and limits and their sources; checked the ladder's exact solutions; read the verification and gradient-check results and the optimized design for physical sense.]</b></p>

| Component | Written by | How it was checked |
|---|---|---|
| Mesh generator (`brakeopt/mesh.py`) | Agent | Rotation invariance, manufactured solution on the curved mesh, every boundary face in exactly one set |
| Heat solver and Newton (`brakeopt/model.py`, `solver.py`) | Agent | Verification ladder against exact solutions; Newton quadratic convergence; energy balance |
| Adjoint gradients (`optim/adjoint.py`) | Agent | Adjoint vs. central finite differences for every constraint and both variable types; Taylor test |
| Optimization (`optim/optimize.py`) | Agent | Three starting designs reach the same optimum; constraints re-checked at the optimum; linear optimum re-checked with full physics |
| Website and 3D viewer (`docs/`) | Agent | Automated test that every table and plot key exists in `results.json`; checked against the eight viewer requirements at desktop and phone width |

**Bugs the verification caught.**

1. *Sign error in the heat flux (caught by verification V3).* The agent's first version measured the heat flux leaving the cold end of the test wall with the wrong sign, so V3 failed (−1.29×10⁶ vs. +1.29×10⁶ W/m²). Comparing against the exact Kirchhoff-transform flux caught it.
2. *Newton stopping too early (caught by the gradient check).* On the stand-in model, the first Newton stopping rule quit after one iteration. Adjoint and finite-difference gradients disagreed (4×10⁻⁷) and the Taylor-test curve was flat instead of falling as ε². Stopping on the size of the Newton update fixed it.
3. *KS too loose (caught by comparing KS with the true maximum).* The first KS sharpness overestimated the true peak by 47 to 65 °C, which would have left iron in the disc that it did not need. Raising ρ brought the gap down to 3 to 5 °C.

## Reproducibility

One command regenerates `docs/assets/results.json`, and therefore every number, table, plot and the 3D comparison on this page:

```bash
pip install -r requirements.txt
python scripts/run_all.py
```

No number on this page is typed by hand: the page reads them from `results.json` when it loads.

<div class="rj-table" data-table="meta"></div>

{% endraw %}
