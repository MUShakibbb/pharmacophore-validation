"""Make `pytest` work from the repo root without an install step.

`python -m pytest` puts the cwd on sys.path; a bare `pytest` does not, so
`from src.gh_score import ...` would fail for anyone who runs the obvious
command. This keeps both working.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
