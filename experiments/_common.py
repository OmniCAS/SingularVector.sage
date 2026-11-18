"""Shared helpers for the experiment scripts (run them from the repository root)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from singular_vector import parse_algebra  # noqa: E402,F401
from singular_vector.export import DATA, Timer, experiment_record, poly_record, write_json  # noqa: E402,F401

# all algebras of the 2026-11 run (delta_s + eps_t): N = 1..6 deltas, orthogonal part of dimension 2..14
ALL_OSP = (["osp(%d|2)" % M for M in range(2, 15)] + ["osp(%d|4)" % M for M in range(2, 13)] +
           ["osp(%d|6)" % M for M in range(2, 11)] + ["osp(%d|8)" % M for M in range(2, 11)] +
           ["osp(%d|10)" % M for M in (2, 4, 5, 6, 7, 8)] + ["osp(4|12)", "osp(6|12)"])
DEFAULT_OSP = ["osp(%d|2)" % M for M in range(2, 11)] + ["osp(%d|4)" % M for M in range(2, 9)] + \
    ["osp(%d|6)" % M for M in range(2, 7)]


def latex_root(g, w):
    return g.weight_name(w, latex=True)
