"""Fast regression tests (about 10 s): python tests/test_small.py  (or: python -m pytest tests)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from singular_vector import (GL, OSP, Pencil, S_formula, Verma, compare_symbols, full_statement, hyperplane_ring,  # noqa: E402
                             parse_algebra, self_check, solve, symbol_criterion)
from sage.all__sagemath_singular import QQ  # noqa: E402


def test_algebras():
    for name, dim in [("osp(3|2)", 12), ("osp(4|2)", 17), ("osp(6|2)", 30), ("osp(5|4)", 40), ("osp(2|4)", 19),
                      ("gl(2|2)", 16), ("gl(3|1)", 16)]:
        g = parse_algebra(name)
        assert g.dim == dim, name
        assert g.check_closure() and g.check_generation()
        assert self_check(g, 100) == 0, name


def test_rho_and_hyperplane():
    g = OSP(6, 2)
    assert [str(x) for x in g.rho()] == ["-2", "2", "1", "0"]          # the paper's rho for osp(6|2)
    for name in ["osp(6|2)", "osp(7|4)", "osp(2|6)", "gl(3|2)"]:
        g = parse_algebra(name)
        rho = g.rho()
        for k in g.isotropic_odd_positive():
            kind, s, t = g.parse_root(g.weight_name(g.wt[k]))
            S, lam, p = hyperplane_ring(g, kind, s, t)                  # the Lean formula ...
            assert g.pairing([lam[c] + rho[c] for c in range(g.r)], g.wt[k]) == 0   # ... is (lam+rho, beta) = 0


def test_known_numbers():
    # weight space dimensions and numbers of equations of the 2026-11-17 logs (osp(6|2)) and of gl
    g = OSP(6, 2)
    got = {g.weight_name(g.wt[k]): (Pencil(g, g.wt[k]).ncols, Pencil(g, g.wt[k]).nrows) for k in g.isotropic_odd_positive()}
    assert got == {"d1-e1": (1, 1), "d1+e1": (15, 23), "d1-e2": (2, 2), "d1+e2": (8, 13), "d1-e3": (4, 5), "d1+e3": (4, 5)}
    h = GL(3, 3)
    for s in range(1, 4):
        for t in range(1, 4):
            assert Pencil(h, h.beta("-", s, t)).ncols == 2 ** (s + t - 2)
    assert [len(OSP(2 * m, 2).words_of_weight(tuple(-x for x in OSP(2 * m, 2).beta("+", 1, 1)))) for m in range(1, 6)] == [1, 4, 15, 55, 200]


def test_statements():
    for name in ["osp(5|2)", "osp(6|2)", "osp(5|4)", "osp(2|4)"]:
        g = parse_algebra(name)
        for k in g.isotropic_odd_positive():
            P = Pencil(g, g.wt[k])
            assert symbol_criterion(P)[0] and compare_symbols(P) == 0
            kind, s, t = g.parse_root(g.weight_name(g.wt[k]))
            S, lam, p = hyperplane_ring(g, kind, s, t)
            v, res = solve(P, S, lam)
            assert not res
            assert full_statement(g, g.wt[k])["holds"]


def test_example_10_coefficient():
    # the last printed term of Example 10 (coefficient of e_{-eps1-delta} = -f[d1+e1], first term = -(our top word))
    g = OSP(6, 2)
    P = Pencil(g, g.beta("+", 1, 1))
    S, lam, p = hyperplane_ring(g, "+", 1, 1)
    b1, b2, b3 = S.gens()
    v, _ = solve(P, S, lam)
    (j,) = [j for j, w in enumerate(P.words) if len(w) == 1]
    assert v[j] == (b1 - b2 + 1) * (b1 + b2 + 3) * (b1 - b3 + 2) * (b1 + b3 + 2)


def test_formula_7():
    g = GL(2, 2)
    lam = [QQ(3) / 2, QQ(-2), QQ(5), QQ(0)]
    s, t = 2, 2
    p, q = g.coord_of_delta(s), g.coord_of_eps(t)
    pp, f = g.atypical_coordinate("-", s, t)
    lam[pp] = f(lam)
    V = Verma(g, lam)
    Sv = S_formula(V, g, p, q)
    assert Sv and all(not V.act(e, Sv) for e in g.simple)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
