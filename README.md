# Project 4: Race-car brake disc, PDE-constrained design optimization

Minimum-mass design of a vented cast-iron race brake rotor, subject to temperature limits during a
hard-braking sequence. Transient nonlinear 3D heat conduction (k(T) + radiation) solved by our own
finite element code with Newton's method, adjoint gradients swept backward in time, and SLSQP.

- `brakeopt/`: the Python package (mesh, solver, functionals, adjoint, optimizer)
- `verification/`: the verification ladder
- `scripts/run_all.py`: **one command** that regenerates `docs/assets/results.json` and every number on the site
- `docs/`: the GitHub Pages report site
- `ainote.md`: status board shared by the workstreams
- `GUIDE.md`: plain-language guide to what is being built and why

```bash
pip install -r requirements.txt
python -m pytest            # fast tests
python -m pytest -m slow    # slow tests
python scripts/run_all.py   # regenerate results.json
```
