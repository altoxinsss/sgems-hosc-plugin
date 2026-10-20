from __future__ import annotations

import sys
from pathlib import Path


EXPERIMENT_DIRECTORY = Path(__file__).resolve().parents[1]
if str(EXPERIMENT_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(EXPERIMENT_DIRECTORY))

