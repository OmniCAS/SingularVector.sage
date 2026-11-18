"""Experiment `singular-vectors`: the singular vectors themselves, as case files for a website.

For each (algebra, root): the singular vector of weight lam - beta on the hyperplane, normalised so that the
word in the simple root vectors has coefficient 1, written in the PBW basis (even root vectors, then the
one odd root vector), every coefficient factored, with LaTeX.  Writes data/cases/<id>.json and
data/cases/index.json.  Full term lists only up to --max-terms words (default 60); larger cases keep the
summary and the coefficient of f[beta].

    python experiments/singular_cases.py [--max-terms 60] [ALG:root ...]
"""
import argparse

from _common import DATA, Timer, experiment_record, parse_algebra, poly_record, write_json
from singular_vector import Pencil, full_statement, hyperplane_ring, solve, symbol_criterion
from sage.misc.latex import latex

DEFAULT = (["osp(%d|2):d1+e1" % M for M in range(2, 11)] + ["osp(%d|2):d1+e2" % M for M in range(4, 9)] +
           ["osp(6|2):d1+e3", "osp(7|2):d1+e3", "osp(8|2):d1+e3"] +
           ["osp(%d|4):d2+e1" % M for M in range(2, 8)] + ["osp(6|4):d1+e1", "osp(6|4):d2+e2", "osp(6|6):d3+e1"] +
           ["osp(6|2):d1-e3", "osp(8|4):d2-e2"])


def case_id(g, beta_name):
    return "%s-%s" % (g.name.replace("(", "").replace(")", "").replace("|", "-"),
                      beta_name.replace("+", "p").replace("-", "m"))


def record(spec, max_terms):
    T = Timer()
    alg, bname = spec.split(":")
    g = parse_algebra(alg)
    kind, s, t = g.parse_root(bname)
    beta = g.beta(kind, s, t)
    P = Pencil(g, beta)
    sym, blocks = symbol_criterion(P)
    S, lam, p = hyperplane_ring(g, kind, s, t)
    v, res = solve(P, S, lam)
    fs = full_statement(g, beta)
    terms = []
    for j, w in enumerate(P.words):
        c = v.get(j)
        if c is None:
            continue
        odd = [i for i in w if g.negpar[i]]
        terms.append({"word": g.word_name(w), "latex": g.word_name(w, latex=True), "pbw_degree": len(w),
                      "odd_factor": g.neg_name(odd[0]) if odd else None, "coefficient": poly_record(c)})
    deg1 = [tm for tm in terms if tm["pbw_degree"] == 1]
    rec = {
        "id": case_id(g, g.weight_name(beta)), "algebra": g.name, "algebra_latex": "\\mathfrak{osp}(%d|%d)" % (g.dimV0, 2 * g.N),
        "beta": g.weight_name(beta), "beta_latex": g.weight_name(beta, latex=True),
        "conventions": "generators of Defs/OSP.lean; lam = (a_N..a_1 | b_1..b_M) = lam(T a a); f_gamma = root vector "
                       "of weight -gamma; PBW order: even root vectors (by height), then the odd one",
        "rho": [str(x) for x in g.rho()], "rho_latex": [latex(x) for x in g.rho()],
        "simple_roots_latex": [g.weight_name(g.wt[k], latex=True) for k in g.simple],
        "hyperplane": {"text": "a%d = %s" % (s, lam[p]), "latex": "a_{%d} = %s" % (s, latex(lam[p]))},
        "free_variables": list(S.variable_names()),
        "pbw_order": [{"name": g.neg_name(i), "latex": g.neg_name(i, latex=True), "odd": bool(g.negpar[i])}
                      for i in range(len(g.neg))],
        "weight_space_dim": P.ncols, "equations": P.nrows, "top_degree": max(P.coldeg),
        "symbol_blocks": {str(d): list(x) for d, x in blocks.items()},
        "statements": {"symbol_criterion": sym, "exactly_one_on_H": not res, "every_lambda": fs["holds"]},
        "normalisation": "coefficient of the word in the simple root vectors = 1",
        "n_terms": len(terms), "n_linear_factor_coefficients": sum(tm["coefficient"]["linear_factors_only"] for tm in terms),
        "f_beta_coefficient": deg1[0]["coefficient"] if deg1 else None,
        "terms": terms if len(terms) <= max_terms else None,
        "terms_omitted": len(terms) > max_terms, "seconds": T(),
    }
    write_json(DATA / "cases" / ("%s.json" % rec["id"]), rec)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="*")
    ap.add_argument("--max-terms", type=int, default=60)
    a = ap.parse_args()
    index = []
    for spec in a.cases or DEFAULT:
        rec = record(spec, a.max_terms)
        index.append({k: rec[k] for k in ("id", "algebra", "algebra_latex", "beta", "beta_latex", "weight_space_dim",
                                          "top_degree", "n_terms", "n_linear_factor_coefficients", "statements",
                                          "terms_omitted")})
        print("%-22s dim %4d terms %4d linear-factor coefficients %4d (%.1f s)" % (
            rec["id"], rec["weight_space_dim"], rec["n_terms"], rec["n_linear_factor_coefficients"], rec["seconds"]), flush=True)
    write_json(DATA / "cases" / "index.json", index)
    experiment_record(
        "singular-vectors", "The singular vectors, written out (case files)",
        "What does the singular vector of weight lam - beta look like, term by term?",
        "Exact solution on the hyperplane, normalised by the word in the simple root vectors; coefficients factored.",
        "%d (algebra, root) cases; full term lists up to %d terms" % (len(index), a.max_terms),
        "Every PBW word of the weight space occurs with a non-zero coefficient in every case: %s." % all(
            x["n_terms"] == x["weight_space_dim"] for x in index),
        "confirmed", "python experiments/singular_cases.py", index)


if __name__ == "__main__":
    main()
