"""Experiment `fbeta-coefficient` (a CONJECTURE, tested, not proved): the coefficient of f[beta] itself.

C1: for beta = delta_s + eps_t in osp(M|2N) and the normalisation 'word in the simple root vectors = 1',
the coefficient of f[beta] is  c * prod_{gamma in Delta_0^+, n = (beta, gamma^vee) > 0} prod_{j<n} ((lam+rho, gamma^vee) - j),
Delta_0^+ = positive roots of g_0 (grading degree 0: no odd simple root in gamma), c = (-1)^(m+1) for M = 2m,
(-1)^(m+1)/2 for M = 2m+1.  For beta = delta_s - eps_t the same product (with another constant) is d_empty of
the paper's formula (7).  Guessed (2026-11-18) from osp(3|2)..osp(12|2) delta_1+eps_1, osp(6|2), osp(7|2)
delta_1+eps_2, osp(5|4), osp(6|4) delta_2+eps_1 and osp(6|4) delta_1+eps_2.

    python experiments/fbeta.py [ALG:root ...]
"""
import argparse

from _common import Timer, experiment_record, parse_algebra
from singular_vector import Pencil, hyperplane_ring, solve
from sage.all__sagemath_singular import QQ
from sage.misc.latex import latex
from singular_vector.export import pretty_factor

GUESSED_FROM = ["osp(%d|2):d1+e1" % M for M in range(3, 13)] + ["osp(6|2):d1+e2", "osp(7|2):d1+e2",
                                                                 "osp(5|4):d2+e1", "osp(6|4):d2+e1", "osp(6|4):d1+e2"]
FRESH = ["osp(13|2):d1+e1", "osp(13|2):d1+e4", "osp(14|2):d1+e1", "osp(14|2):d1+e7", "osp(11|4):d2+e1",
         "osp(11|4):d1+e3", "osp(12|4):d2+e5", "osp(9|6):d3+e1", "osp(9|6):d2+e2", "osp(10|6):d3+e4",
         "osp(8|8):d4+e1", "osp(8|8):d2+e3", "osp(7|8):d4+e2", "osp(9|8):d3+e2", "osp(6|10):d5+e3",
         "osp(5|10):d4+e1", "osp(4|12):d6+e2", "osp(3|8):d4+e1", "osp(2|6):d3+e1", "osp(3|4):d2+e1"]
MINUS = ["osp(6|2):d1-e3", "osp(8|4):d2-e4", "osp(9|4):d2-e3", "osp(7|6):d3-e2", "osp(10|2):d1-e5", "osp(5|8):d4-e2"]


def grading_degree(g, k):
    odd_simple = [i for i, kk in enumerate(g.simple) if g.par[kk]]
    return sum(g.simple_coeffs[k][i] for i in odd_simple)


def predicted(g, S, lam, beta):
    rho = g.rho()
    mu = [lam[c] + rho[c] for c in range(g.r)]
    out, used = S(1), []
    for k in g.pos:
        if g.par[k] or grading_degree(g, k) != 0:
            continue
        w = g.wt[k]
        n = g.coroot_pairing(beta, w)
        if n <= 0:
            continue
        x = S(g.coroot_pairing(mu, w))
        for j in range(int(n)):
            out *= (x - j)
        used.append({"gamma": g.weight_name(w), "gamma_latex": g.weight_name(w, latex=True), "n": int(n),
                     "factor_latex": latex(x)})
    return out, used


def test(spec):
    T = Timer()
    alg, b = spec.split(":")
    g = parse_algebra(alg)
    kind, s, t = g.parse_root(b)
    beta = g.beta(kind, s, t)
    S, lam, p = hyperplane_ring(g, kind, s, t)
    P = Pencil(g, beta)
    v, res = solve(P, S, lam)
    assert not res
    (j,) = [j for j, w in enumerate(P.words) if len(w) == 1]
    pred, used = predicted(g, S, lam, beta)
    K = S.fraction_field()
    ratio = K(v[j]) / K(pred)
    const = S(ratio.numerator()).degree() <= 0 and S(ratio.denominator()).degree() <= 0
    c_rule = QQ((-1) ** (g.M + 1)) / (2 if g.odd else 1)
    return {"case": spec, "kind": kind, "product_matches": bool(const), "constant": str(ratio) if const else None,
            "constant_rule_holds": bool(const and kind == "+" and ratio == c_rule) if kind == "+" else None,
            "roots": used, "degree": int(pred.degree()), "computed_latex": pretty_factor(v[j])[1], "seconds": T()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="*")
    a = ap.parse_args()
    groups = [("guessed-from", GUESSED_FROM), ("fresh", FRESH), ("delta-minus-eps", MINUS)]
    if a.cases:
        groups = [("given", a.cases)]
    rows = []
    for label, specs in groups:
        for spec in specs:
            r = test(spec)
            r["group"] = label
            rows.append(r)
            print("%-16s %-16s product %s constant %s rule %s (%.1f s)" % (label, spec, r["product_matches"], r["constant"],
                                                                          r["constant_rule_holds"], r["seconds"]), flush=True)
    plus = [r for r in rows if r["kind"] == "+"]
    ok = all(r["product_matches"] and r["constant_rule_holds"] for r in plus)
    fresh = [r for r in plus if r["group"] == "fresh"]
    experiment_record(
        "fbeta-coefficient", "Conjecture C1: the coefficient of f[beta] (guess, tested)",
        "Is the coefficient of the root vector f[beta] in the normalised singular vector a product of linear factors, and which?",
        "Compare the computed coefficient with c * prod_{gamma in Delta_0^+, (beta,gamma^vee)>0} prod_{j<(beta,gamma^vee)} "
        "((lam+rho,gamma^vee) - j).",
        "%d delta+eps cases (%d of them computed only after the guess), %d delta-eps cases" % (
            len(plus), len(fresh), len(rows) - len(plus)),
        ("Holds, with the constant rule, in all %d delta+eps cases; for delta-eps the product holds with constant "
         "(-1)^(t-1) (it is d_empty of formula (7))." % len(plus)) if ok else "Fails in some case: see rows.",
        "evidence" if ok else "refuted", "python experiments/fbeta.py", rows,
        {"status": "conjecture (computed evidence, not a proof)"})


if __name__ == "__main__":
    main()
