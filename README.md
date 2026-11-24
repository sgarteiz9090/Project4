# Project 4: EV battery busbar, PDE-constrained design optimization

Minimum-mass design of a copper battery-pack busbar that carries 400 A continuously and a 60 s,
1200 A boost pulse. Two-way coupled electro-thermal 3D finite elements (current flow, Joule heating
with temperature-dependent conductivity, resistive bolted joints, natural convection and radiation),
Newton's method in every time step, adjoint gradients swept backward in time, and SLSQP.
Handout menu item: "Busbar or resistive heater".

Report and live 3D comparison: https://sgarteiz9090.github.io/Project4/

- `busbar/`: the Python package (mesh, coupled model, solver, functionals, adjoint, optimizer, verification ladder)
- `scripts/run_all.py`: **one command** that regenerates `docs/assets/results.json` and every number on the site
- `docs/`: the GitHub Pages report site (served from `main:/docs`)
- `GUIDE.md`: plain-language guide to what was built and why
- `ainote.md`: status board
- `alternate_topic/`: the earlier brake-disc version of the project, complete and runnable on its own
  (`cd alternate_topic && python -m pytest`)

```bash
pip install -r requirements.txt
python -m pytest            # fast tests
python -m pytest -m slow    # slow tests (full ladder, adjoint vs FD)
python -m busbar.ladder     # verification ladder only
python scripts/run_all.py   # regenerate results.json (~25 min)
```
