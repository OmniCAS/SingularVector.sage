"""SingularVector.sage: exact computations of odd singular vectors in Verma modules of Lie superalgebras.

Rewritten in 2026 on top of the SageMath library (passagemath 10.8).  The 2022 tools are kept under
legacy/2022/.  Main objects:

  OSP(M, 2N)                      osp(M|2N) in the conventions of Defs/OSP.lean
  GL(m, n)                        gl(m|n) in the conventions of Defs/GL.lean
  Verma(g, lam), Pencil(g, beta)  Verma modules, the singular system D(lam) v = 0
  symbol_criterion, solve, full_statement, check_at_points   (singular.py)
"""
from .algebra import SuperAlgebra
from .gl import GL, S_formula, c_coeff
from .osp import OSP, paper_matrix
from .singular import (check_at_points, compare_symbols, full_ring, full_statement, hyperplane_ring, solve,
                       symbol_criterion, symbol_derivation)
from .verma import Pencil, Verma, self_check

__version__ = "2026.11"


def parse_algebra(name):
    """'osp(6|2)' -> OSP(6, 2); 'gl(2|3)' -> GL(2, 3)."""
    kind, rest = name.strip().split("(")
    a, b = (int(x) for x in rest.rstrip(")").split("|"))
    return OSP(a, b) if kind == "osp" else GL(a, b)
