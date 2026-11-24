#!/usr/bin/env bash
# Replace the brake-disc site with the busbar site (only if the instructor rejects the brake disc).
# Run from the repository root on a branch that already contains busbar/ and its generated results:
#     bash busbar/swap_into_docs.sh && python -m pytest -q && git add -A docs tests && git commit -m "Switch site to busbar"
# The brake-disc code (brakeopt/, optim/, scripts/) stays in the repository untouched; busbar reuses
# brakeopt.fem and optim/ (adjoint, KS, gradient checks, SLSQP driver).
set -euo pipefail
test -f busbar/site/assets/results.json || { echo "run: python -m busbar.run_all first"; exit 1; }
rm -rf docs && mkdir docs && cp -R busbar/site/. docs/
# test_site.py checks the brake viewer's n_vanes; the busbar mesh has no vanes
sed -i 's/assert m\["n_vanes"\] > 0/assert m.get("n_vanes", 1) > 0/' tests/test_site.py
echo "docs/ now holds the busbar site. Pages serves main:/docs, so merging this to main publishes it."
