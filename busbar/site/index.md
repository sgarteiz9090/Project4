---
layout: default
title: EV battery busbar optimization
---
{% raw %}

<div id="mock-banner" class="mock-banner" hidden></div>

# A lighter EV battery busbar that survives a 1200 A boost pulse

<p class="lede">PDE-constrained design optimization of a copper battery-pack busbar: 24 width and thickness variables, a two-way coupled electro-thermal 3D finite-element model (current flow and Joule heating with temperature-dependent conductivity, natural convection and radiation, resistive bolted joints), a transient current pulse, adjoint gradients swept backward in time, and SLSQP.</p>

## 1. Summary

We minimize the copper mass of a 300 mm busbar that links two battery modules, while it carries 400 A continuously and a 60 s boost pulse of 1200 A. Three requirements must hold over the whole history: the bar stays below its insulation rating, the bolted joints (which sit on the cell terminals) stay below the cell limit, and the terminal-to-terminal resistance at the continuous rating stays below the efficiency budget. The optimizer goes from <span class="rj" data-calc="mass_initial_g" data-d="1"></span> g to <span class="rj" data-calc="mass_optimized_g" data-d="1"></span> g (<span class="rj" data-calc="mass_saving_pct" data-d="1"></span> % less copper) and, at the same time, turns an initial design that fails two of the three requirements into one that passes all of them.

<div class="kpis">
  <div class="kpi"><b><span class="rj" data-calc="mass_optimized_g" data-d="0"></span> g</b><span>optimized busbar (initial <span class="rj" data-calc="mass_initial_g" data-d="0"></span> g)</span></div>
  <div class="kpi"><b>−<span class="rj" data-calc="mass_saving_pct" data-d="1"></span> %</b><span>copper saved</span></div>
  <div class="kpi"><b><span class="rj" data-calc="joint_peak_optimized" data-d="1"></span> °C</b><span>peak joint temperature (limit <span class="rj" data-calc="joint_limit" data-d="0"></span> °C)</span></div>
  <div class="kpi"><b><span class="rj" data-calc="R_optimized" data-d="1"></span> µΩ</b><span>resistance (limit <span class="rj" data-calc="R_limit" data-d="0"></span> µΩ)</span></div>
  <div class="kpi"><b><span class="rj" data-calc="adjoint_speedup" data-d="0"></span>×</b><span>adjoint speed-up over finite differences</span></div>
</div>

### Interactive 3D comparison

The left bar is the initial design and the right bar the optimized one, coloured by temperature on one shared scale. Press **Play** to watch both bars heat up during the 1200 A pulse and recover afterwards, or drag the timeline. Hover over (or tap) any face to see its element, temperature, current density and margin to its limit. The thin bars are easier to read with the thickness exaggerated (menu), which is labelled when on.

<div id="rotor-viewer" class="rv"></div>
<script type="module" src="assets/js/viewer.js"></script>

<div class="rj-plot" data-plot="peaks"></div>

## 2. Structure and formulation

### Why this part matters

An electric car's battery pack is built from modules connected in series by copper busbars. At full acceleration or at the peak of a fast charge the pack current reaches 1000 A or more for tens of seconds. Every busbar heats by Joule losses in the copper and, even more locally, in its bolted joints, whose contact resistance is comparable to that of the whole bar. That heat has two bad destinations: the polymer insulation sleeve of the bar, and, through the bolted joint, the module terminal and the cells behind it, which age quickly above about 60 °C. Copper is also heavy and expensive, and a pack has dozens of busbars. The design question is where copper can be removed and where it must stay. The physics is genuinely coupled: copper's resistivity rises 0.39 % per kelvin, so a hotter bar makes more heat (two-way electro-thermal coupling).

### Design variables

At 12 equally spaced stations $x_1,\dots,x_{12}$ along the bar (length $L$ = <span class="rj" data-param="L" data-scale="1000" data-d="0"></span> mm) we choose

$$
\mathbf d = \big(w_1,\dots,w_{12},\; t_1,\dots,t_{12}\big),
$$

the width $w$ and the thickness $t$ of the bar's rectangular cross-section, piecewise linear in $x$ between stations. Bounds: <span class="rj" data-param="w_min" data-scale="1000" data-d="0"></span> mm $\le w_j \le$ <span class="rj" data-param="w_max" data-scale="1000" data-d="0"></span> mm (lower: the bolt hole and the joint; upper: the space between neighbouring module terminals) and <span class="rj" data-param="t_min" data-scale="1000" data-d="1"></span> mm $\le t_j \le$ <span class="rj" data-param="t_max" data-scale="1000" data-d="0"></span> mm (stamping / bending limits of flat copper bar). The initial design is the conventional uniform 30 × 5 mm bar.

<div class="rj-table" data-table="design_vars"></div>

### Objective

$$
\min_{\mathbf d}\; m(\mathbf d) = \rho_{\text{Cu}} \int_{\Omega(\mathbf d)} \mathrm d V \qquad\text{(busbar mass, kg).}
$$

### Constraints

All three are functionals of the coupled transient solution $(\phi, T)(\mathbf x,t;\mathbf d)$. The two temperature constraints are aggregated over a node set $S$ and the whole history with the Kreisselmeier–Steinhauser (KS) function of the normalized nodal temperatures $a_{i,n} = T_{i,n}/T_{\lim}$ (in °C):

$$
g = \mathrm{KS}_\rho(a) - 1 \le 0,
\qquad
\mathrm{KS}_\rho(a) = a_{\max} + \frac1\rho \ln \sum_{i\in S}\sum_{n=1}^{N} e^{\rho (a_{i,n} - a_{\max})} \;\ge\; a_{\max},
$$

with $\rho = \rho_K T_{\lim}$, $\rho_K = 1\ \mathrm K^{-1}$ (KS sits a few kelvin above the true peak; the report shows both).

- $g_{\text{bar}}$: $S$ = every node of the bar, $T_{\lim}$ = <span class="rj" data-calc="bar_limit" data-d="0"></span> °C, the rating of the insulation sleeve (UL 224 105 °C class heat-shrink tubing).
- $g_{\text{joint}}$: $S$ = the nodes of the two bolted joint faces, $T_{\lim}$ = <span class="rj" data-calc="joint_limit" data-d="0"></span> °C. The joint face is clamped to the module terminal; with the terminal itself at 40 °C, this keeps the cell tabs within the maximum operating temperature of typical Li-ion cells (55 to 60 °C at the cell, data sheets) with the 10 K drop of the terminal path.
- $g_R$: resistance between the two terminals at the continuous rating, $R = V_{\text{in}}/I \le R_{\max}$ = <span class="rj" data-calc="R_limit" data-d="0"></span> µΩ (a loss budget of about 8 W per busbar at 400 A, including the two joints).
- Width and thickness bounds as above.

## 3. Analysis model

**Governing PDEs.** Charge conservation and transient heat conduction in the copper, coupled both ways:

$$
\nabla\cdot\big(\sigma(T)\,\nabla \phi\big) = 0,
\qquad
\rho c\,\frac{\partial T}{\partial t} = \nabla\cdot\big(k(T)\,\nabla T\big) + \sigma(T)\,\lvert\nabla\phi\rvert^2
\quad\text{in } \Omega(\mathbf d),
$$

$$
\sigma(T) = \frac{\sigma_0}{1+\alpha (T-T_{\text{ref}})},\qquad k(T) = k_0\big(1+\beta(T-T_{\text{ref}})\big),
$$

with $\sigma_0$ = 5.8×10⁷ S/m and $\alpha$ = 0.00393 K⁻¹ for annealed copper (IEC 60028) and $k_0$ = 398 W/m K.

**Bolted joints.** Each end overlaps a module terminal over $l_c$ = <span class="rj" data-param="l_c" data-scale="1000" data-d="0"></span> mm on its underside. The joint is a contact resistance $\rho_c$ = <span class="rj" data-param="rho_c" data-fmt="sci" data-d="1"></span> Ω m² (about 10 µΩ for a 25 × 30 mm joint; Braunovic et al., *Electrical Contacts*, 2007), which carries the current and generates heat; half of that heat enters the bar. Heat leaves through the terminal with conductance $h_{\text{term}}$ to the 40 °C terminal:

$$
j_n = \frac{\phi - V_t}{\rho_c},\qquad
-k\frac{\partial T}{\partial n} = h_{\text{term}}\,(T - T_{\text{term}}) - \tfrac12\,\rho_c\, j_n^2,
$$

with $V_t = V_{\text{in}}$ at the inlet joint and $V_t = 0$ at the outlet. $V_{\text{in}}$ is an extra unknown, fixed by the total current: $\int_{\text{inlet}} (V_{\text{in}}-\phi)/\rho_c\,\mathrm dA = I(t)$.

**Exposed surfaces.** Electrically insulated; thermally natural convection and radiation to the 40 °C pack air,

$$
-k\frac{\partial T}{\partial n} = h(T)\,(T - T_\infty) + \varepsilon\sigma_{\text{SB}}\big(T^4-T_\infty^4\big),
\qquad h(T) = 1.42\Big(\frac{\lvert T-T_\infty\rvert}{L_c}\Big)^{1/4}
$$

(laminar natural convection in air, Holman, *Heat Transfer*, Table 7-2; $L_c$ = 50 mm, the bar mounted on edge; $\varepsilon$ = 0.85 for the insulated surface). $\lvert\Delta T\rvert$ is smoothed as $(\Delta T^2 + \delta^2)^{1/2}$ with $\delta$ = 0.05 K so the residual stays differentiable at $\Delta T = 0$.

**Load case.** First the steady state at the continuous 400 A (one implicit step with $\Delta t = 10^9$ s, so the starting field depends on the design and is differentiated like every other step), then 60 s at 1200 A ($\Delta t$ = 2 s) and 240 s of recovery at 400 A ($\Delta t$ = 6 s): <span class="rj" data-calc="n_steps" data-d="0"></span> implicit steps.

**Discretization.** Half of the bar (mirror plane $y=0$), 8-node trilinear hexahedra in the true tapered geometry (<span class="rj" data-calc="n_elems" data-d="0"></span> elements, <span class="rj" data-calc="n_nodes" data-d="0"></span> nodes). Each time level has the state $\mathbf u_n = (\boldsymbol\psi_n, \mathbf T_n, V_{\text{in},n})$, <span class="rj" data-calc="n_dof" data-d="0"></span> unknowns (potential in mV, temperature in K, inlet potential), solved **monolithically** by Newton's method with the exact tangent:

$$
\mathbf R_n(\mathbf u_n, \mathbf u_{n-1}) =
\begin{bmatrix}
\mathbf K_\sigma(\mathbf T_n)\,\boldsymbol\psi_n + \mathbf c_\phi(\boldsymbol\psi_n, V_{\text{in},n})\\[2pt]
\mathbf M\,\dfrac{\mathbf T_n-\mathbf T_{n-1}}{\Delta t_n} + \mathbf K_k(\mathbf T_n)\,\mathbf T_n - \mathbf q_J(\boldsymbol\psi_n,\mathbf T_n) + \mathbf r_{\partial}(\mathbf T_n, \boldsymbol\psi_n)\\[2pt]
-\int_{\text{in}} j_n\,\mathrm dA - I_n/2
\end{bmatrix} = \mathbf 0 .
$$

The volume part of the tangent is derived by hand ($\sigma'(T)$, $k'(T)$ and the Joule term $\sigma\lvert\nabla\phi\rvert^2$ give off-diagonal blocks between $\boldsymbol\psi$ and $\mathbf T$); the boundary part is differentiated by complex step over each face's 8 or 9 local unknowns, which is exact to round-off. Verification V11 checks both.

**Assumptions.** Straight bar of rectangular section (no bends); the contact resistance does not depend on temperature or on clamping force; the terminals are a fixed 40 °C with a lumped conductance; natural convection from a simple correlation (no enclosure flow, no view-factor reduction of radiation); the insulation sleeve is thin enough to ignore its thermal resistance; skin effect is negligible (DC traction current).

<details><summary>All model parameters (from results.json)</summary>
<div class="rj-table" data-table="params"></div>
</details>

## 4. Verification

Every test is run by `python -m busbar.run_all` and reported straight from `results.json`. V1 to V3 have closed-form answers the element space reproduces exactly; V4 is a closed-form **nonlinear, two-way coupled** benchmark (in a bar of constant section the current density is fixed by charge conservation, so $\sigma(T)$ makes the Joule source linear in $T$ and $T(x)$ is a cosine); V5 and V6 are exact solutions for a body cooled only by natural convection ($\Delta T(t) = (\Delta T_0^{-1/4} + a t/4)^{-4}$) or only by radiation; V7 is a manufactured solution of the full coupled 3D problem with every nonlinear term on a tapered mesh.

<div class="rj-table" data-table="verification"></div>

<figure>
<div class="rj-plot" data-plot="mms"></div>
<figcaption>Convergence: RMS nodal error against element size (slope 2 expected for trilinear hexes) and against time step (slope 1 for implicit Euler).</figcaption>
</figure>

<figure>
<div class="rj-plot" data-plot="newton"></div>
<figcaption>Newton residual in the first step of the 1200 A pulse: the number of correct digits roughly doubles every iteration (quadratic convergence of the coupled solve).</figcaption>
</figure>

## 5. Sensitivities

**State equation.** One coupled nonlinear residual per implicit step, $\mathbf R_n(\mathbf u_n,\mathbf u_{n-1},\mathbf d) = \mathbf 0$, $n = 1,\dots,N$. $\mathbf u_0$ (40 °C, no current) does not depend on $\mathbf d$; the steady continuous-rating state is step 1, so its dependence on $\mathbf d$ is carried by the same adjoint. The Newton tangent $\mathbf A_n = \partial\mathbf R_n/\partial\mathbf u_n$ is the full coupled Jacobian; the link to the previous step is $\mathbf B_n = \partial\mathbf R_n/\partial\mathbf u_{n-1}$, which is $-\mathbf M/\Delta t_n$ in the temperature block and zero elsewhere (the potential has no time derivative).

**Partials of the functionals.** KS: $\partial g/\partial T_{i,n} = w_{i,n}/T_{\lim}$ with softmax weights $w_{i,n} = e^{\rho(a_{i,n}-\mathrm{KS})}$. Resistance: $g_R = \phi_{\text{scale}} V_{\text{in},1}/(I_1 R_{\max}) - 1$ depends only on one entry of $\mathbf u_1$. The mass is an explicit function of the geometry.

**Discrete adjoint, swept backward through the current history.** With $\mathcal L = g + \sum_{n}\boldsymbol\lambda_n^{\mathsf T}\mathbf R_n$,

$$
\mathbf A_n^{\mathsf T}\boldsymbol\lambda_n = -\Big(\frac{\partial g}{\partial \mathbf u_n}\Big)^{\mathsf T} - \mathbf B_{n+1}^{\mathsf T}\boldsymbol\lambda_{n+1},
\qquad \boldsymbol\lambda_{N+1} = \mathbf 0,\qquad n = N, \dots, 1,
\qquad
\frac{\mathrm d g}{\mathrm d \mathbf d} = \frac{\partial g}{\partial \mathbf d} + \sum_{n=1}^{N} \boldsymbol\lambda_n^{\mathsf T}\,\frac{\partial \mathbf R_n}{\partial \mathbf d}.
$$

The sweep runs **backward in time** from the end of the recovery to the steady state, because the temperature at step $n$ influences the peaks through every later step. Each step costs one sparse LU of $\mathbf A_n^{\mathsf T}$ shared by all three constraints. The adjoint multiplier has an electrical and a thermal part, $\boldsymbol\lambda_n = (\boldsymbol\lambda_{\phi,n}, \boldsymbol\lambda_{T,n}, \lambda_{V,n})$: the electrical part is how a change of the potential field (through Joule heating and through $V_{\text{in}}$) feeds the peak temperature. This is implicit differentiation of the coupled nonlinear solve (handout implementation iii). The backward sweep is the brake-disc project's generic `optim/adjoint.py`, reused unchanged on the coupled state.

**Shape term.** $\mathbf R_n$ depends on $\mathbf d$ only through the node coordinates; $\sum_n\boldsymbol\lambda_n^{\mathsf T}\partial\mathbf R_n/\partial d_j$ is evaluated by complex step on the residual with the states fixed ($h = 10^{-30}$, exact to round-off, 24 complex residual assemblies per step, no extra PDE solves).

### Gradient check (adjoint vs. central finite differences)

<div class="rj-table" data-table="gradient_check"></div>

<figure>
<div class="rj-plot" data-plot="taylor"></div>
<figcaption>Taylor test along a random direction: the error of the central difference against the adjoint directional derivative falls as ε² (grey reference) until round-off takes over.</figcaption>
</figure>

### Cost: adjoint vs. finite differences

<div class="rj-table" data-table="cost"></div>

## 6. Optimization

**Algorithm and scaling.** SciPy SLSQP (`ftol` = $10^{-7}$). Variables $\mathbf x = \mathbf d/\mathbf d_0$, objective $f = m/m_0$, constraints $g_k = \mathrm{KS}_k - 1$ and $g_R = R/R_{\max} - 1$, all of order one. Three starts (the uniform bar, all upper bounds, the midpoint of the box) test whether the optimum depends on the start.

<div class="rj-table" data-table="starts"></div>

<div class="rj-pair">
<div class="rj-plot" data-plot="history"></div>
<div class="rj-plot" data-plot="history_g"></div>
</div>

<figure>
<div class="rj-plot" data-plot="thickness"></div>
<figcaption>Width (dashed) and ten times the thickness (solid) along the bar. The joints occupy the first and last 25 mm.</figcaption>
</figure>

### Active constraints and multipliers

<div class="rj-table" data-table="multipliers"></div>

**KKT read-out.** A constraint is active when $g_k > -10^{-3}$, a bound when the variable sits on it. The multipliers solve $\nabla f + \sum_k \mu_k\nabla g_k + \sum_{\text{bounds}}\nu_i(\pm\mathbf e_i) = \mathbf 0$ with $\mu,\nu\ge 0$ (non-negative least squares). In scaled form $\mu_k$ is the drop in $m/m_0$ per unit relaxation of $g_k$, so $\mu_k m_0/T_{\lim,k}$ is the **kilograms of copper saved per kelvin** of extra allowed temperature, and $\mu_R m_0/R_{\max}$ the kilograms saved per micro-ohm.

**Engineering reading.** The optimizer makes the bar as wide as the space allows and thins it: for the same cross-section a wide, thin bar has more surface for convection and radiation, and copper conducts heat along the bar so well that the cross-section can follow the heat rather than the current. Copper is kept (and added) next to the joints: there the bar acts as a cooling fin for the joint, carrying the contact heat away from the cell terminal into the span, where it can leave by convection. The mid-span is thinner, limited by the resistance budget rather than by temperature.

## 7. Nonlinear component (N2)

The state equation is nonlinear in four ways, all in the residual, the Newton tangent and the adjoint: the temperature-dependent conductivity $\sigma(T)$, which makes the electro-thermal coupling two-way; the Joule source $\sigma\lvert\nabla\phi\rvert^2$ (quadratic in the potential gradient); natural convection $h \propto \Delta T^{1/4}$ (heat loss $\propto \Delta T^{5/4}$); and radiation $\propto T^4$ (plus the weak $k(T)$). We optimize the same problem with a **constant-property model**: $\sigma$, $k$ and $h$ frozen at their values at 70 °C and radiation linearized there, so the potential no longer depends on the temperature. That optimum is then re-evaluated with the full physics.

<div class="rj-table" data-table="nonlinear"></div>

<div class="rj-table" data-table="n2_attr"></div>

The constant-property optimum is <span class="rj" data-calc="linear_mass_delta_pct" data-d="1"></span> % lighter than the full-physics optimum, but re-checked with the full physics its joint reaches <span class="rj" data-calc="joint_peak_linear" data-d="1"></span> °C and its resistance <span class="rj" data-calc="R_linear" data-d="1"></span> µΩ (limits <span class="rj" data-calc="joint_limit" data-d="0"></span> °C and <span class="rj" data-calc="R_limit" data-d="0"></span> µΩ). The second table shows which nonlinearity is responsible.

## 8. Results

<div class="rj-table" data-table="results"></div>

## 9. Bonus sections

- **Fully nonlinear (+5):** section 7; the coupled state equation is solved by Newton's method at every step with the exact tangent, differentiated implicitly, and checked against the closed-form nonlinear coupled benchmark V4 and the exact nonlinear cooling laws V5 and V6.
- **Path / time dependence (+5):** transient heating during the current pulse; both temperature constraints aggregate the whole history and the adjoint is swept backward through it (section 5), checked against finite differences.
- **Design under uncertainty (+5):** three sampled scenarios (nominal, a 10 % higher peak current, 5 K hotter pack air and terminals). The robust design satisfies all nine constraints, each with its own adjoint gradient: the worst case over the samples.

<div class="rj-table" data-table="robust"></div>

## 10. Limitations and AI-use statement

**Limitations.** The joints are lumped (a contact resistance and a conductance to a fixed-temperature terminal), not modelled as the terminal's own 3D conduction. Natural convection uses a single correlation with a fixed length scale, and radiation ignores view factors to neighbouring modules. The bar is straight; real busbars have bends and bolt holes. The current profile is one representative boost pulse (the robust design covers three variants). The temperature limits are design values chosen from the insulation rating and cell data sheets, not from one standard.

**AI-use statement.** The solver, adjoint, optimization scripts and this website were written with a coding agent (Claude Code). The backward adjoint sweep, the KS functional, the gradient-check tools and the SLSQP driver are reused unchanged from our brake-disc version of this project; the coupled electro-thermal model, its mesh, verification ladder and site are new.

<p class="callout"><b>[Sivan: replace this box with what you derived or checked yourself.]</b></p>

**A bug the verification caught.** The first version of the nonlinear coupled benchmark (V4) used too high a current: the solver converged to a temperature profile below absolute zero, and it agreed with the "exact" cosine formula, because at that current the steady problem has no physical solution (thermal runaway: $\tan(\lambda L/2)$ exceeds $h/(k\lambda)$, so the joints cannot remove the heat). Checking the sign of the answer, not just the agreement, caught it; the benchmark now runs in the stable range. Separately, the 1D Joule test (V3) first missed nodal exactness by 10⁻⁴ because the "steady" step used $\Delta t = 10^9$ s, so a sliver of heat capacity remained; raising it to $10^{15}$ s in the test recovered round-off agreement.

## Reproducibility

One command regenerates `busbar/site/assets/results.json`, and therefore every number, table, plot and the 3D comparison on this page:

```bash
pip install -r requirements.txt
python -m busbar.run_all
```

No number on this page is typed by hand: the page reads them from `results.json` when it loads.

<div class="rj-table" data-table="meta"></div>

{% endraw %}
