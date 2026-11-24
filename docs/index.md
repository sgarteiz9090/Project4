---
layout: default
title: Race-car brake disc optimization
---
{% raw %}

<div id="mock-banner" class="mock-banner" hidden></div>

# Lighter race brake discs that still survive hard braking

<p class="lede">PDE-constrained design optimization of a vented cast-iron brake rotor: 24 thickness variables, a transient nonlinear 3D heat-conduction model solved with our own finite-element code, adjoint gradients swept backward in time, and SLSQP.</p>

## 1. Summary

We minimize the mass of a vented race-car brake rotor while keeping the rubbing faces below the grey-iron cracking and fade limit and the hub below the bearing-protection limit through a sequence of hard stops from 250 to 80 km/h. The optimizer removes <span class="rj" data-calc="mass_saving_kg" data-d="2"></span> kg (<span class="rj" data-calc="mass_saving_pct" data-d="1"></span> %) of iron, going from <span class="rj" data-k="designs.initial.mass_kg" data-d="2"></span> kg to <span class="rj" data-k="designs.optimized.mass_kg" data-d="2"></span> kg. The peak friction-face temperature rises from <span class="rj" data-calc="face_peak_initial" data-d="0"></span> °C to <span class="rj" data-calc="face_peak_optimized" data-d="0"></span> °C against a limit of <span class="rj" data-calc="face_limit" data-d="0"></span> °C, so the face-temperature constraint is the one that stops the optimizer.

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

A race car's front brake discs turn the car's kinetic energy into heat: a stop from 250 to 80 km/h dumps several hundred kilojoules into each front disc in about three seconds. The disc must store that heat and then shed it, mostly through the radial cooling channels between its two friction rings ("cheeks"), which are joined by vanes that pump air like a centrifugal fan. The disc is also **rotating, unsprung mass**: every gram costs acceleration, braking and grip over bumps. Too little iron and the rubbing faces overheat (pad fade, thermal cracking of grey cast iron) and heat soaks through the hub into the wheel bearing and brake fluid. The question is where iron can be removed and where it must stay.

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

Both temperature constraints are functionals of the transient temperature field $T(\mathbf x,t;\mathbf d)$, aggregated over space and the whole braking history with the Kreisselmeier–Steinhauser (KS) function over the nodal values $T_{a,n}$ (node $a$ of the surface set, time step $n$):

$$
g_{\text{face}} = \mathrm{KS}_\rho\big(T_{a,n}\big)_{a\in\Gamma_{\text{fric}},\,n=0..N} \le T_{\text{face}}^{\max},
\qquad
g_{\text{hub}} = \mathrm{KS}_\rho\big(T_{a,n}\big)_{a\in\Gamma_{\text{hub}},\,n=0..N} \le T_{\text{hub}}^{\max},
$$

$$
\mathrm{KS}_\rho(x) = x_{\max} + \frac1\rho \log \sum_i e^{\rho (x_i - x_{\max})} \;\ge\; x_{\max}.
$$

- $T_{\text{face}}^{\max}$ = <span class="rj" data-calc="face_limit" data-d="0"></span> °C on the pad rubbing band: above about 700 °C grey cast iron rotors suffer heat checking and cracking and pad friction fades (Limpert, *Brake Design and Safety*, 3rd ed., SAE 2011, ch. 3; Day, *Braking of Road Vehicles*, 2014, ch. 3).
- $T_{\text{hub}}^{\max}$ = <span class="rj" data-calc="hub_limit" data-d="0"></span> °C on the hub / bell mounting face, protecting wheel bearing grease and brake fluid (source in `brakeopt/params.py`).
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

**Load case.** Three stops 250 → 80 km/h at 1.5 g with re-acceleration between them, starting from a 200 °C disc.

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

**Discrete adjoint, backward in time.** For any functional $J(\mathbf T_0,\dots,\mathbf T_N,\mathbf d)$ and the residuals $\mathbf R_n(\mathbf T_n,\mathbf T_{n-1},\mathbf d) = \mathbf 0$, the Lagrangian $\mathcal L = J + \sum_n \boldsymbol\lambda_n^{\mathsf T}\mathbf R_n$ is stationary in every $\mathbf T_n$ when the adjoint states satisfy, for $n = N, N-1, \dots, 1$,

$$
\Big(\frac{\partial \mathbf R_n}{\partial \mathbf T_n}\Big)^{\mathsf T}\boldsymbol\lambda_n
= -\frac{\partial J}{\partial \mathbf T_n} - \Big(\frac{\partial \mathbf R_{n+1}}{\partial \mathbf T_n}\Big)^{\mathsf T}\boldsymbol\lambda_{n+1},
\qquad \boldsymbol\lambda_{N+1} = \mathbf 0,
\qquad \frac{\partial \mathbf R_{n+1}}{\partial \mathbf T_n} = -\frac{\mathbf M}{\Delta t_{n+1}},
$$

and then

$$
\frac{\mathrm d J}{\mathrm d \mathbf d} = \frac{\partial J}{\partial \mathbf d} + \sum_{n=1}^{N} \boldsymbol\lambda_n^{\mathsf T}\,\frac{\partial \mathbf R_n}{\partial \mathbf d},
\qquad \frac{\partial \mathbf R_n}{\partial \mathbf d} = \frac{\partial \mathbf R_n}{\partial \mathbf X}\,\frac{\partial \mathbf X}{\partial \mathbf d}.
$$

$\partial\mathbf R_n/\partial\mathbf T_n$ is the exact Newton tangent already used in the forward solve, so each adjoint step is one linear solve with the transposed Jacobian. One backward sweep gives the gradient with respect to all 24 variables at once. For the KS constraints, $\partial J/\partial \mathbf T_n$ is the vector of KS weights $e^{\rho(T_{a,n}-T_{\max})}/\sum e^{\rho(\cdot)}$ on the constrained nodes; for the mass, $\boldsymbol\lambda = \mathbf 0$ and only $\partial m/\partial \mathbf d$ remains.

### Gradient check (adjoint vs. central finite differences)

<div class="rj-table" data-table="gradient_check"></div>

<figure>
<div class="rj-plot" data-plot="taylor"></div>
<figcaption>Taylor test along a random direction p: the first-order remainder falls as ε and, with the adjoint gradient included, as ε² until round-off takes over.</figcaption>
</figure>

### Cost: adjoint vs. finite differences

<div class="rj-table" data-table="cost"></div>

## 6. Optimization

**Algorithm and scaling.** <span class="rj" data-k="optimization.algorithm"></span>, with <span class="rj" data-k="optimization.scaling"></span>. Bounds go straight to SLSQP; the two KS temperature constraints use the adjoint gradients above.

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

A multiplier is the mass saved per unit relaxation of its constraint (section 11 of the KKT notes): it prices each limit in grams of iron.

## 7. Nonlinear component (N2)

The model has two nonlinear terms: temperature-dependent conductivity $k(T)$ of grey cast iron and surface radiation $\varepsilon\sigma(T^4-T_\infty^4)$. Both enter the residual, the Newton tangent and the adjoint. The table compares the design optimized with the full nonlinear model against the design optimized with both terms switched off (constant $k$, no radiation), each re-checked with the nonlinear model.

<div class="rj-table" data-table="nonlinear"></div>

The design optimized without the nonlinear terms differs in mass by <span class="rj" data-calc="linear_mass_delta_pct" data-d="1"></span> % from the nonlinear optimum.

## 8. Results

<div class="rj-table" data-table="results"></div>

## 9. Bonus sections

- **Fully nonlinear (+5):** section 7; Newton per implicit-Euler step with the exact tangent and the same tangent in the adjoint.
- **Path / time dependence (+5):** the constraints aggregate the whole braking history and the adjoint is swept backward through it (section 5).

## 10. Limitations and AI-use statement

**Limitations.** Convection coefficients come from correlations, not from a CFD model of the vane channels; the pad heat is treated as circumferentially uniform; thermal stresses and the resulting cracking are represented only through the temperature limit; one braking scenario is designed for.

**AI-use statement.** A coding agent (Claude) wrote most of the code in this repository, split into workstreams coordinated through `ainote.md`. This website, the 3D viewer (`docs/assets/js/viewer.js`) and the report page that fills every number from `results.json` (`docs/assets/js/report.js`) were written by the agent; we checked the viewer against the handout's eight requirements at desktop and phone width. The parts we derived or checked ourselves, and one bug the verification caught, are listed in `ainote.md` and `GUIDE.md`.

## Reproducibility

One command regenerates `docs/assets/results.json`, and therefore every number, table, plot and the 3D comparison on this page:

```bash
pip install -r requirements.txt
python scripts/run_all.py
```

No number on this page is typed by hand: the page reads them from `results.json` when it loads.

<div class="rj-table" data-table="meta"></div>

{% endraw %}
