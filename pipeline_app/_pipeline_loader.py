"""Loads the three real pipeline scripts as modules, without touching
sys.path or sys.modules under their own names.

`statistics/statistics.py` in particular lives in a folder named
`statistics`, which shadows the stdlib `statistics` module if pipeline_app's
own directory (or above) ever ends up on sys.path — risky, since
pandas/scipy/lifelines may import the real stdlib module internally.
Loading each script via importlib.util.spec_from_file_location under a
private name sidesteps that entirely: nothing is added to sys.path, and
each module is registered in sys.modules only under its own private key.
"""

import importlib.util
from functools import lru_cache
from pathlib import Path

# formatting/, score_patNo/, and statistics/ all live directly under
# pipeline_app/ — this app is self-contained, not a thin layer over
# sibling top-level project folders.
PIPELINE_ROOT = Path(__file__).resolve().parent


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    import sys
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@lru_cache(maxsize=1)
def load_formatting():
    return _load("_pipeline_formatting", PIPELINE_ROOT / "formatting" / "formatting.py")


@lru_cache(maxsize=1)
def load_score_patno():
    return _load("_pipeline_score_patno", PIPELINE_ROOT / "score_patNo" / "score_patno.py")


@lru_cache(maxsize=1)
def load_statistics():
    return _load("_pipeline_statistics", PIPELINE_ROOT / "statistics" / "statistics.py")
