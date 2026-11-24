# Fallback topic: EV battery busbar (handout menu: "Busbar or resistive heater")

Kept on its own branch and in its own folder so it never touches the brake-disc project on `main`.
Use it only if the instructor does not approve the off-menu brake disc.

| | |
|---|---|
| Structure | Copper busbar linking two battery modules, 300 mm, bolted lap joints on the module terminals |
| Design variables | 24: width and thickness at 12 stations (piecewise linear) |
| Objective | Mass |
| Constraints | Peak bar temperature <= 105 C (insulation), peak joint temperature <= 70 C (cell terminal), resistance <= 48 micro-ohm (KS over space and time for the temperatures) |
| PDEs | div(sigma(T) grad phi) = 0 and rho c dT/dt = div(k(T) grad T) + sigma(T) abs(grad phi)^2, solved monolithically |
| Nonlinear (N2, +5) | sigma(T) two-way coupling, Joule source, natural convection h ~ dT^(1/4), radiation T^4, k(T); Newton with exact tangent |
| Time dependence (+5) | 400 A steady, 60 s 1200 A pulse, 240 s recovery; adjoint swept backward through the history |
| Uncertainty (+5) | robust design over 3 scenarios (+10 % peak current, +5 K ambient) |

```bash
python -m busbar.ladder            # verification ladder (V1-V12)
python -m busbar.run_all           # regenerates busbar/site/assets/results.json and figures (~1 h)
python -m pytest -q tests/test_busbar.py
```

Reused unchanged from the brake disc: `brakeopt/fem.py` (hex and quad elements), `optim/adjoint.py`
(backward adjoint sweep), `optim/functionals.py` (KS), `optim/gradcheck.py` (FD table, Taylor test,
cost), `optim/optimize.py` (SLSQP + KKT read-out), and the site's layout, CSS and report engine.

To switch the published site to the busbar: `bash busbar/swap_into_docs.sh` (copies `busbar/site/`
over `docs/`), run the tests, commit, merge to `main`. To preview locally:
`cd busbar/site && jekyll serve`.
