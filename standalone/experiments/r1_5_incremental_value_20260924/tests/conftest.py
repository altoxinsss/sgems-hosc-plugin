from __future__ import annotations

import sys
from pathlib import Path


EXPERIMENT_DIR = Path(__file__).resolve().parents[1]
if str(EXPERIMENT_DIR) not in sys.path:
    sys.path.insert(0, str(EXPERIMENT_DIR))

# The formal Section 4.4--4.6 orchestrator is intentionally located at the
# experiments root while the R1-5 implementation remains self-contained.
# Include both roots so protocol tests exercise the public entrypoint.
EXPERIMENTS_ROOT = EXPERIMENT_DIR.parent
if str(EXPERIMENTS_ROOT) not in sys.path:
    sys.path.insert(0, str(EXPERIMENTS_ROOT))
