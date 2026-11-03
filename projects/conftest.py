"""
Make the study modules in this folder importable by their names in tests.

pytest runs with ``--import-mode=importlib``, which does not put a test's
folder on :data:`sys.path`; this folder is not an installed package.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
