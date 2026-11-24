# Sensitivities and optimization: derivation (WS2+3)

Text for report sections 5 (Sensitivities) and 6 (Optimization). Every number lives in
`results.json` under `sensitivities` and `optimization`; none is quoted here.

## State equation (one residual per time step)

The semi-discrete heat equation on the sector mesh, with implicit Euler from $t_{n-1}$ to $t_n$, gives one
nonlinear residual per step $n = 1,\dots,N$ (N = 144 for the three-stop sequence):

$$
R_n(T_n, T_{n-1}, d) = \frac{M(d)}{\Delta t_n}(T_n - T_{n-1}) + K(T_n, d)\,T_n + H_n(d)(T_n - T_\text{air})
+ \mathrm{Rad}(T_n, d) - Q_n(d) = 0 .
$$

$T_0$ (200 °C everywhere) does not depend on the design $d$. Newton's method solves each step; its tangent
$A_n = \partial R_n/\partial T_n$ includes $k'(T)$ and $4\varepsilon\sigma T^3$. The coupling to the previous
step is $B_n = \partial R_n/\partial T_{n-1} = -M/\Delta t_n$.

## Constraint functionals over space and time

For a node set $S$ (friction faces, or the hub face), with $a_{i,n} = T_{i,n}/T_\text{lim}$ (°C):

$$
g = \mathrm{KS}_\rho(a) - 1 \le 0,\qquad
\mathrm{KS}_\rho(a) = a_\max + \frac{1}{\rho}\ln\sum_{i\in S}\sum_{n=1}^{N} e^{\rho(a_{i,n}-a_\max)} .
$$

We use $\rho = \rho_K T_\text{lim}$ with $\rho_K = 1\ \mathrm{K}^{-1}$, so KS lies within a few kelvin above the true
peak (the report shows both). Its partials are softmax weights: $\partial g/\partial T_{i,n} = w_{i,n}/T_\text{lim}$,
$w_{i,n} = e^{\rho(a_{i,n}-\mathrm{KS})}$, $\sum w = 1$. Only the hottest few node-instants carry weight.

The objective is the rotor mass $m(d) = \rho_\text{Fe}\sum_e |\Omega_e(d)|$, an explicit function of the geometry,
so $dm/dd$ needs no adjoint.

## Discrete adjoint, swept backward through the braking history

Lagrangian: $\mathcal{L} = g + \sum_{n=1}^{N} \lambda_n^\top R_n$. Because every $R_n = 0$, $dg/dd = d\mathcal L/dd$ for any
$\lambda$. Choose $\lambda$ to kill every $dT_n/dd$ term. $T_n$ appears in $g$, in $R_n$ and in $R_{n+1}$, so

$$
A_n^\top \lambda_n = -\Big(\frac{\partial g}{\partial T_n}\Big)^{\!\top} - B_{n+1}^\top \lambda_{n+1},
\qquad \lambda_{N+1} = 0,\qquad n = N, N-1, \dots, 1,
$$

$$
\frac{dg}{dd} = \frac{\partial g}{\partial d} + \sum_{n=1}^{N} \lambda_n^\top \frac{\partial R_n}{\partial d}.
$$

The sweep runs **backward in time**: the last step is solved first, because the influence of the
temperature at step $n$ on the peak flows forward through all later steps. This is the time-dependent
(path-dependent) adjoint of the rubric bonus. It needs the stored forward history $T_0,\dots,T_N$ (231 nodes ×
145 levels, small, so no checkpointing is needed). Each step costs one sparse LU of $A_n^\top$, shared by both
constraints (two right-hand sides). Because the state equation is nonlinear (radiation and $k(T)$), $A_n$ is the
exact Newton tangent at the converged state: this is implicit differentiation (handout implementation iii).

**Shape term.** $R_n$ depends on $d$ only through the node coordinates $X(d)$, which the mesh generator moves
smoothly (fixed topology): $\partial R_n/\partial d = (\partial R_n/\partial X)(\partial X/\partial d)$. We evaluate
$\sum_n \lambda_n^\top \partial R_n/\partial d_j$ by a complex step in $d_j$ on the residual with the
states fixed: $\partial R_n/\partial d_j = \operatorname{Im} R_n(T_n,T_{n-1},d + i h e_j)/h$, $h = 10^{-30}$. That is exact
to round-off (no subtraction), costs 24 complex residual assemblies per step (vectorized over all steps), and
never re-solves the PDE.

## Verification of the gradients

1. **Finite-difference table.** For each functional (mass, face KS, hub KS), several vane and cheek variables,
   both with full physics and with the linearized model: adjoint vs central FD with $h = 10^{-5} d_j$.
2. **Taylor test.** Along a random direction $v$ in scaled space,
   $\big|[g(x+\epsilon v) - g(x-\epsilon v)]/2\epsilon - \nabla g\cdot v\big|$ must fall as $\epsilon^2$
   until round-off. A flat curve means a wrong gradient, not a bad step.
3. **Cost.** One forward solve plus one adjoint sweep vs $2 \times 24$ forward solves for central FD.

## Optimization setup

- Variables scaled $x = d/d_0$; objective $f = m/m_0$; constraints $g_k = \mathrm{KS}_k - 1 \le 0$; bounds from the
  casting and wear limits. All three are O(1), as SLSQP needs.
- SciPy SLSQP, `ftol = 1e-7`, from three starts (baseline, all upper bounds, midpoint of the box). The best
  feasible result is reported, with the spread of the others.
- **KKT read-out.** At the optimum, active constraints are those with $g_k > -10^{-3}$, active bounds those at
  their limits. Multipliers $\mu$ solve $\nabla f + \sum_k \mu_k \nabla g_k + \sum_{\text{bounds}} \nu_i (\pm e_i) = 0$,
  $\mu,\nu \ge 0$ (NNLS), next to SLSQP's own estimates. In scaled form $\mu_k$ is the fall in $m/m_0$ per unit
  relaxation of $g_k$, so $\mu_k\, m_0 / T_{\text{lim},k}$ is **kilograms saved per kelvin of extra allowed temperature**.
- **N2 comparison.** The same problem with `radiation=False` (radiation linearized about $T_\text{init}$) and
  `kT=False` (k frozen at $k(T_\text{init})$) is optimized, then its optimum is re-evaluated with full physics.
- **Uncertainty (bonus).** Three sampled scenarios (nominal, +20 km/h initial speed, +0.05 front brake bias).
  The robust design satisfies both constraints in every scenario (6 constraints, each with its own exact adjoint
  gradient). That is the worst case over the samples.
