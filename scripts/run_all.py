"""One command that regenerates docs/assets/results.json, the figures and every number on the site.

    python scripts/run_all.py            # full run (~25 min)
    python scripts/run_all.py --quick    # smoke run

Same as `python -m busbar.run_all`.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from busbar.run_all import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
