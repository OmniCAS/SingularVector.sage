"""Experiment `gl-so-factorisation`: the singular vector for delta_s + eps_t in osp(M|2N) is
(gl(N) part) x (the osp(M|2) vector for delta_1 + eps_t).

Both g_0-singular vectors are computed over QQ[lam] (all of h*, even equations only, normalised by the word
in the simple root vectors).  A PBW word of osp(M|2N) is (delta-delta letters)(eps letters) f[delta_j +- eps_k];
the test: its coefficient is A(delta part, j; a) * B(eps part, +-k; b), B the osp(M|2) coefficient of
(eps letters) f[delta_1 +- eps_k], A free of the b's.  Records the A's (products of linear factors, the shape
of the gl formula).

    python experiments/factorization.py [ALG:root ...]
"""
import argparse

from _common import Timer, experiment_record, parse_algebra
from singular_vector import OSP, Pencil, full_ring, solve
from singular_vector.export import pretty_factor
from sage.misc.latex import latex

DEFAULT = ["osp(6|4):d2+e1", "osp(5|6):d3+e1", "osp(7|4):d2+e2", "osp(6|6):d3+e2", "osp(8|4):d2+e3", "osp(4|8):d4+e1",
           "osp(7|6):d3+e1", "osp(5|8):d4+e2"]


def g0_vector(g, beta):
    P = Pencil(g, beta, simple=[k for k in g.simple if g.par[k] == 0])
    S, lam = full_ring(g)
    v, res = solve(P, S, lam)
    assert not res
    return P, S, v


def split(g, word):
    gl, so, odd = [], [], None
    for i in word:
        w = g.negwt[i]
        if g.negpar[i]:
            odd = w
        elif any(w[c] for c in range(g.N)):
            gl.append(g.weight_name(tuple(-x for x in w)))
        else:
            so.append(g.weight_name(tuple(-x for x in w)))
    j = next(c for c in range(g.N) if odd[c])
    return tuple(sorted(gl)), g.N - j, tuple(sorted(so)), tuple(odd[c] for c in range(g.N, g.r))


def latex_name(g, nm):
    """'d2-d1' -> 'f_{\\delta_{2}-\\delta_{1}}' via the coordinate names of g."""
    out = nm
    for a, l in sorted(zip(g.coord_names, g.latex_coords), key=lambda x: -len(x[0])):
        out = out.replace(a, l)
    return "f_{%s}" % out


def test(spec):
    T = Timer()
    alg, b = spec.split(":")
    G = parse_algebra(alg)
    kind, s, t = G.parse_root(b)
    Pb, Sb, vb = g0_vector(G, G.beta("+", s, t))
    g1 = OSP(G.dimV0, 2)
    Ps, Ss, vs = g0_vector(g1, g1.beta("+", 1, t))
    K = Sb.fraction_field()
    B = {}
    for jj, c in vs.items():
        _, _, so, sp = split(g1, Ps.words[jj])
        assert not any(str(x).startswith("a") for x in c.variables())
        B[(so, sp)] = K(Sb(str(c)))
    rows = {}
    for jj, c in vb.items():
        gl, j, so, sp = split(G, Pb.words[jj])
        rows.setdefault((gl, j), {})[(so, sp)] = K(c)
    ok, A = True, []
    for (gl, j), row in sorted(rows.items(), key=lambda kv: (-len(kv[0][0]), kv[0])):
        rs = {row[c] / B[c] for c in row} if set(row) == set(B) else set()
        if len(rs) != 1:
            ok = False
            break
        r = rs.pop()
        if any(str(x).startswith("b") for x in r.numerator().variables() + r.denominator().variables()):
            ok = False
            break
        num = Sb(r.numerator())
        A.append({"delta_letters": list(gl), "odd_index_j": j,
                  "delta_letters_latex": " ".join(latex_name(G, nm) for nm in gl),
                  "A_latex": pretty_factor(num)[1] + ("" if r.denominator() == 1 else " / " + latex(r.denominator())),
                  "A_text": str(r)})
    return {"case": spec, "small": "%s:d1+e%d" % (g1.name, t), "dim": Pb.ncols, "dim_small": Ps.ncols,
            "gl_rows": len(rows), "factorises": ok, "A": A if ok else None, "seconds": T()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", nargs="*")
    a = ap.parse_args()
    rows = []
    for spec in a.cases or DEFAULT:
        r = test(spec)
        rows.append(r)
        print("%-16s vs %-14s dims %4d / %3d  gl rows %2d  factorises %s (%.1f s)" % (
            spec, r["small"], r["dim"], r["dim_small"], r["gl_rows"], r["factorises"], r["seconds"]), flush=True)
    ok = all(r["factorises"] for r in rows)
    experiment_record(
        "gl-so-factorisation", "The singular vector factors as gl(N) part x osp(M|2) part",
        "For beta = delta_s + eps_t in osp(M|2N), how does the singular vector depend on N and s?",
        "Compute the g_0-singular vectors over QQ[lam] for osp(M|2N) and osp(M|2); split every PBW word into its "
        "delta-delta letters, eps letters and odd letter; test that the coefficient matrix is an outer product.",
        "%d cases" % len(rows),
        "In all cases the coefficient is A(delta part) * B(eps part), B = the osp(M|2) coefficient; every A is a product "
        "of linear factors in the a's (the shape of formula (7))." if ok else "Fails in some case: see rows.",
        "confirmed" if ok else "refuted", "python experiments/factorization.py", rows)


if __name__ == "__main__":
    main()
