"""Experiment `example-10`: Example 10 of Liu-Luo-Wang (osp(6|2), beta = delta_1 + eps_1, a = b1 + 4).

The paper's root vectors are built from the matrices printed in Section 2.4 (osp.paper_matrix) and identified
with the generators of Defs/OSP.lean.  The singular vector (unique on the hyperplane) is written in the 15
printed monomials and compared with the 15 printed coefficients.  Also: which printed monomials have the
wrong weight, and every way of repairing them by changing one factor.

    python experiments/example10.py
"""
from itertools import product

from _common import experiment_record
from singular_vector import OSP, Pencil, Verma, hyperplane_ring, solve
from sage.all__sagemath_singular import QQ, matrix
from sage.misc.latex import latex

g = OSP(6, 2)
S, lam, p = hyperplane_ring(g, "+", 1, 1)          # free b1, b2, b3; a1 = b1 + 4
b1, b2, b3 = S.gens()
A = (b1 - b2) * (b1 + b2 + 3)
# (printed coefficient, printed monomial); names e_{...} of the paper: 'e3-e2' = e_{eps3-eps2}, 'e1-d' = e_{eps1-delta}
PRINTED = [
    (S(1), ["e3-e2", "-e3-e2", "e2-e1", "e2-e1", "e1-d"]),
    (-(b2 + b3 + 2), ["e3-e2", "-e3-e1", "e2-e1", "e1-d"]),
    (-(b2 - b3 + 2), ["-e3-e2", "e3-e2", "e2-e1", "e1-d"]),
    (-(2 * b1 + 4), ["e3-e2", "-e3-e2", "e2-e1", "e2-d"]),
    (-(A + (b2 - b3 + 2) * (b2 + b3 + 2)), ["-e2-e1", "e2-e1", "e1-d"]),
    (-A, ["-e3-e1", "-e3-e1", "e1-d"]),
    ((b1 + b2 + 3) * (b1 + b3 + 2), ["e3-e2", "e2-e1", "-e3-d"]),
    ((b1 + b2 + 3) * (b1 - b3 + 2), ["-e3-e2", "e2-e1", "e3-d"]),
    ((b1 - b2 + 1) * (b2 + b3 + 2) - A, ["e3-e2", "-e3-e1", "e2-d"]),
    ((b1 - b2 + 1) * (b2 - b3 + 2) - A, ["-e3-e2", "e3-e1", "e2-d"]),
    ((b1 - b2) ** 2 * (b1 + b2 + 3) + (b1 - b2 + 1) * (b2 - b3 + 2) * (b2 + b3 + 2), ["-e2-e1", "e2-d"]),
    ((b1 + b2 + 3) * (b1 - b3 + 2) * (b1 + b3 + 2), ["e2-e1", "-e2-d"]),
    (A * (b1 + b3 + 2), ["e3-e1", "-e3-d"]),
    (A * (b1 - b3 + 2), ["-e3-e1", "e3-d"]),
    ((b1 - b2 + 1) * (b1 + b2 + 3) * (b1 - b3 + 2) * (b1 + b3 + 2), ["-e1-d"]),
]
NAMES = ["e3-e2", "e2-e1", "e3-e1", "-e3-e2", "-e3-e1", "-e2-e1", "e1-d", "e2-d", "e3-d", "-e1-d", "-e2-d", "-e3-d"]


def gamma_name(paper):
    """'e1-d' -> '-d1+e1' (the root in the d/e notation of osp.paper_matrix)."""
    return paper.replace("-d", "-d1") if paper.endswith("-d") else paper


FAC = {nm: g.paper_root_vector(gamma_name(nm)) for nm in NAMES}


def tex_name(nm):
    s = nm.replace("-d", "-\\delta")
    for k in "123":
        s = s.replace("e%s" % k, "\\epsilon_{%s}" % k)
    return "e_{%s}" % s


def main():
    beta = g.beta("+", 1, 1)
    P = Pencil(g, beta)
    v, res = solve(P, S, lam)
    assert not res
    V = Verma(g, [S(x) for x in lam])
    target = tuple(-x for x in beta)

    def weight(names):
        return tuple(sum(g.wt[FAC[nm][0]][c] for nm in names) for c in range(g.r))

    def vec_of(names):
        sc = QQ(1)
        for nm in names:
            sc *= FAC[nm][1]
        u = V.apply_product([FAC[nm][0] for nm in names])
        return [sc * QQ(u.get(w, 0)) for w in P.words]
    sing = [v.get(j, S(0)) for j in range(P.ncols)]

    def coords(terms):
        Mx = matrix(QQ, [vec_of(nm) for _, nm in terms]).transpose()
        if Mx.rank() < len(terms):
            return None
        Mi = Mx.inverse()
        x = [sum(Mi[i, j] * sing[j] for j in range(P.ncols)) for i in range(len(terms))]
        return [S(xi / x[0]) for xi in x]
    wrong = [i + 1 for i, (_, nm) in enumerate(PRINTED) if weight(nm) != target]
    cands = {}
    for t in wrong:
        names = PRINTED[t - 1][1]
        cands[t] = [alt for pos, new in product(range(len(names)), NAMES)
                    for alt in [names[:pos] + [new] + names[pos + 1:]] if weight(alt) == target and alt != names]
    working = []
    for choice in product(*[cands[t] for t in wrong]):
        terms = list(PRINTED)
        for t, alt in zip(wrong, choice):
            terms[t - 1] = (terms[t - 1][0], alt)
        x = coords(terms)
        if x is not None and all(xi == cf for xi, (cf, _) in zip(x, terms)):
            working.append([" ".join(c) for c in choice])
    terms = list(PRINTED)
    terms[2] = (terms[2][0], ["-e3-e2", "e3-e1", "e2-e1", "e1-d"])
    terms[5] = (terms[5][0], ["-e3-e1", "e3-e1", "e1-d"])
    x = coords(terms)
    rows = []
    for i, ((cf, nm), xi) in enumerate(zip(terms, x)):
        rows.append({"term": i + 1, "monomial": " ".join(nm), "monomial_latex": " ".join(tex_name(n) for n in nm),
                     "printed_monomial": " ".join(PRINTED[i][1]), "monomial_misprinted": (i + 1) in wrong,
                     "printed_coefficient_latex": latex(cf), "computed_coefficient_latex": latex(xi),
                     "equal": xi == cf})
    ok = all(r["equal"] for r in rows) and len(wrong) == 2
    experiment_record(
        "example-10", "Example 10 recomputed (osp(6|2), beta = delta_1 + eps_1)",
        "Are the 15 printed terms of Example 10 right?",
        "Root vectors from the matrices printed in Section 2.4; the unique singular vector on a = b1 + 4 written in the "
        "15 printed monomials; exhaustive one-factor repairs of the monomials whose weight is wrong.",
        "weight space of dimension %d; %s one-factor repairs tried" % (P.ncols, {str(k): len(c) for k, c in cands.items()}),
        ("Monomials %s are misprinted (wrong weight); after the repair (unique up to the order of two commuting "
         "factors) all 15 printed coefficients are correct." % wrong) if ok else "Mismatch: see rows.",
        "confirmed" if ok else "refuted", "python experiments/example10.py", rows,
        {"repairs_that_work": working, "scalars_paper_vs_generators": {nm: [str(g.labels[k]), str(c)] for nm, (k, c) in FAC.items()}})
    print("misprinted monomials:", wrong, "| repairs that work:", working, "| all 15 coefficients equal:", ok)


if __name__ == "__main__":
    main()
