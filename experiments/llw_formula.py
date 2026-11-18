"""Experiment `llw-formula`: the closed formula S_{-beta} of Liu-Luo-Wang, (3)-(7), against the computation.

gl(m|n), beta = delta_s - eps_t (Theorem 1), and osp(M|2N), beta = delta_s - eps_t (Proposition 9: the same
formula with the generators T a b).  For random rational points of the hyperplane: S_{-beta} v is non-zero,
killed by the simple root vectors, and spans the (one-dimensional) solution space; the weight space has
dimension 2^(s+t-2) for gl.  Control: S_{-beta} v at a point off the hyperplane is not singular.

    python experiments/llw_formula.py [--max-rank 6]
"""
import argparse
import random

from _common import Timer, experiment_record
from singular_vector import GL, OSP, Pencil, S_formula, Verma
from sage.all__sagemath_singular import QQ, matrix

rng = random.Random(1903046)


def check(g, s, t, chain=None):
    T = Timer()
    kind = "-"
    beta = g.beta(kind, s, t)
    p, q = g.coord_of_delta(s), g.coord_of_eps(t)
    P = Pencil(g, beta)
    out = {"algebra": g.name, "beta": g.weight_name(beta), "dim_weight_space": P.ncols}
    good = True
    for _ in range(2):
        lam = [QQ(rng.randint(-20, 20)) / rng.randint(1, 5) for _ in range(g.r)]
        pp, f = g.atypical_coordinate(kind, s, t)
        lam[pp] = QQ(f(lam))
        V = Verma(g, lam)
        Sv = S_formula(V, g, p, q, chain=chain)
        killed = all(not V.act(e, Sv) for e in g.simple)
        Mx = P.at(lam)
        ker = Mx.right_kernel()
        vec = [Sv.get(w, 0) for w in P.words]
        good &= bool(bool(Sv) and killed and ker.dimension() == 1 and
                     matrix(QQ, [vec, list(ker.basis()[0])]).rank() == 1)
    rho = g.rho()
    while True:                                   # a point genuinely off the hyperplane
        lam2 = [QQ(rng.randint(-200, 200)) / rng.randint(1, 13) for _ in range(g.r)]
        if g.pairing([lam2[c] + rho[c] for c in range(g.r)], beta) != 0:
            break
    V2 = Verma(g, lam2)
    S2 = S_formula(V2, g, p, q, chain=chain)
    control = bool(S2) and not all(not V2.act(e, S2) for e in g.simple)
    if not control:
        out["control_point"] = [str(x) for x in lam2]
    out.update({"formula_is_the_singular_vector": bool(good), "control_off_hyperplane_rejected": control,
                "seconds": T()})
    if isinstance(g, GL):
        out["dim_equals_2^(s+t-2)"] = P.ncols == 2 ** (s + t - 2)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-rank", type=int, default=6)
    a = ap.parse_args()
    rows = []
    for tot in range(2, a.max_rank + 1):
        for m in range(1, tot):
            g = GL(m, tot - m)
            for s in range(1, m + 1):
                for t in range(1, tot - m + 1):
                    rows.append(check(g, s, t))
    for (M, two_N) in [(4, 2), (5, 2), (6, 2), (7, 2), (8, 2), (4, 4), (5, 4), (6, 4), (6, 6), (3, 4), (2, 4)]:
        g = OSP(M, two_N)
        chain = lambda a, b, g=g: (g.index[("T", a, b)], QQ(1))
        for s in range(1, g.N + 1):
            for t in range(1, g.M + 1):
                rows.append(check(g, s, t, chain))
    for r in rows:
        print("%-9s %-6s dim %3d formula %s control %s" % (r["algebra"], r["beta"], r["dim_weight_space"],
                                                          r["formula_is_the_singular_vector"], r["control_off_hyperplane_rejected"]))
    ok = all(r["formula_is_the_singular_vector"] and r["control_off_hyperplane_rejected"] and r.get("dim_equals_2^(s+t-2)", True)
             for r in rows)
    ngl = sum(1 for r in rows if r["algebra"].startswith("gl"))
    experiment_record(
        "llw-formula", "The formula S_{-beta} of Theorem 1 and Proposition 9, recomputed",
        "Does formula (7) give the singular vector for gl(m|n) (Theorem 1) and, with the generators T a b, for "
        "beta = delta_s - eps_t in osp (Proposition 9)?",
        "Apply S_{-beta} to the highest weight vector at random rational atypical weights; compare with the kernel of the "
        "singular system; control at a non-atypical weight.",
        "gl(m|n) with m + n <= %d: %d cases; osp: %d cases" % (a.max_rank, ngl, len(rows) - ngl),
        "In every case S_{-beta} v spans the solution space (one-dimensional); the controls fail as they should." if ok
        else "Mismatch: see rows.",
        "confirmed" if ok else "refuted", "python experiments/llw_formula.py", rows,
        {"note": "points only (random rational atypical weights); Theorem 1 itself is proved in Lean"})


if __name__ == "__main__":
    main()
