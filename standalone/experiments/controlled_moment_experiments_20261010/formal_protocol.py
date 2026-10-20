from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType
from typing import Any


EXPERIMENTS_DIRECTORY = Path(__file__).resolve().parents[1]
R1_5_DIRECTORY = EXPERIMENTS_DIRECTORY / "r1_5_incremental_value_20260924"


def _load_r1_5_module(alias: str, filename: str) -> ModuleType:
    """Load a formal R1-5 module by path without colliding with stdlib names."""
    existing = sys.modules.get(alias)
    if existing is not None:
        return existing
    if str(R1_5_DIRECTORY) not in sys.path:
        sys.path.insert(0, str(R1_5_DIRECTORY))
    source_path = R1_5_DIRECTORY / filename
    specification = importlib.util.spec_from_file_location(alias, source_path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"Cannot load formal protocol module {source_path}")
    module = importlib.util.module_from_spec(specification)
    sys.modules[alias] = module
    specification.loader.exec_module(module)
    return module


def local_statistics_from_values(*args: Any, **kwargs: Any) -> Any:
    """Call the one formal local-statistics implementation by explicit path."""
    module = _load_r1_5_module("_r1_5_formal_statistics", "statistics.py")
    return module.local_statistics_from_values(*args, **kwargs)

