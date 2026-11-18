"""Singular vectors of weight lam - beta: the symbol criterion, the exact solution, the full statement.

For a Pencil P (rows: simple root vectors e, PBW words u; columns: PBW words of weight -beta):

* `symbol_criterion(P)`: the degree-preserving blocks of D (independent of lam) are injective below
  the top degree and the top degree has a single word (the monomial in the simple root vectors).
  Then for EVERY lam the space of vectors of weight lam - beta killed by those e has dimension <= 1.
* `solve(P, S, lam)`: D(lam) v = 0 over a polynomial ring S (lam a list of affine elements of S),
  degree by degree from the top, inverting only the constant symbol blocks, v(top word) = 1.
  Returns (v, residual); an empty residual means v(lam) is a solution at every point.
* `full_statement(g, beta)`: for EVERY lam, dim (singular vectors of weight lam - beta) = 1 if
  (lam + rho, beta) = 0, else 0 -- via the even simple roots (symbol criterion + solution over all of
  QQ[lam]) and the odd simple roots (each e_odd . v0 divisible by (lam + rho, beta), quotient with a
  non-zero constant term).
"""
import random

from sage.all__sagemath_singular import QQ, PolynomialRing, matrix

from .verma import Pencil, Verma


def symbol_blocks(P):
    blocks = {}
    for (i, j), c in P.ent.items():
        if P.rowdeg[i] == P.coldeg[j]:
            assert c.degree() <= 0, "a symbol entry depends on lam"
            blocks.setdefault(P.coldeg[j], {})[(i, j)] = QQ(c)
        else:
            assert P.rowdeg[i] < P.coldeg[j], "D does not lower the PBW degree"
    return blocks


def symbol_criterion(P, keep=None):
    """(holds, {degree: (columns, rows, rank)}).  `keep(e)` restricts the simple root vectors."""
    top = max(P.coldeg)
    blocks = symbol_blocks(P)
    ok = sum(1 for d in P.coldeg if d == top) == 1
    details = {}
    for d in sorted(set(P.coldeg)):
        cols = [j for j in range(P.ncols) if P.coldeg[j] == d]
        rows = [i for i in range(P.nrows) if P.rowdeg[i] == d and (keep is None or keep(P.rowkeys[i][0]))]
        ci = {j: k for k, j in enumerate(cols)}
        ri = {i: k for k, i in enumerate(rows)}
        S = matrix(QQ, len(rows), len(cols), sparse=True)
        for (i, j), c in blocks.get(d, {}).items():
            if i in ri:
                S[ri[i], ci[j]] = c
        rk = int(S.rank()) if rows else 0
        details[d] = (len(cols), len(rows), rk)
        kernel = len(cols) - rk
        if (d < top and kernel) or (d == top and kernel != 1):
            ok = False
    return ok, details


def symbol_derivation(g, e, word):
    """D_e(word) in the supersymmetric algebra S(n^-), D_e(f) = the n^- component of [e, f]."""
    out, sign = {}, 1
    for i, j in enumerate(word):
        for c, cc in g.bracket(e, g.neg[j]).items():
            if g.kind[c] != "-":
                continue
            letters = list(word[:i]) + [g.negpos[c]] + list(word[i + 1:])
            s, k = sign, i
            while k > 0 and letters[k - 1] > letters[k]:
                if g.negpar[letters[k - 1]] and g.negpar[letters[k]]:
                    s = -s
                letters[k - 1], letters[k] = letters[k], letters[k - 1]
                k -= 1
            while k + 1 < len(letters) and letters[k + 1] < letters[k]:
                if g.negpar[letters[k + 1]] and g.negpar[letters[k]]:
                    s = -s
                letters[k + 1], letters[k] = letters[k], letters[k + 1]
                k += 1
            if any(letters[q] == letters[q + 1] and g.negpar[letters[q]] for q in range(len(letters) - 1)):
                continue
            out[tuple(letters)] = out.get(tuple(letters), 0) + s * cc
        if g.par[e] and g.negpar[j]:
            sign = -sign
    return {k: x for k, x in out.items() if x}


def compare_symbols(P):
    """Mismatches between the degree-preserving part of e . word (PBW) and D_e(word) (S(n^-))."""
    V = Verma(P.g, list(P.L.gens()))
    bad = 0
    for e in P.simple:
        for w in P.words:
            lead = {u: c for u, c in V.act_word(e, w).items() if len(u) == len(w)}
            der = symbol_derivation(P.g, e, w)
            if set(lead) != set(der) or any(lead[u] != der[u] for u in der):
                bad += 1
    return bad


def hyperplane_ring(g, kind, s, t):
    """Polynomial ring of the free coordinates of the hyperplane, and lam as a list in it."""
    p, f = g.atypical_coordinate(kind, s, t)
    free = [c for c in range(g.r) if c != p]
    S = PolynomialRing(QQ, [g.var_names[c] for c in free], order="degrevlex")
    lam = [None] * g.r
    for c, x in zip(free, S.gens()):
        lam[c] = x
    lam[p] = S(f(lam))
    return S, lam, p


def full_ring(g):
    S = PolynomialRing(QQ, g.var_names)
    return S, list(S.gens())


def solve(P, S, lam):
    """D(lam) v = 0 over S, degree by degree from the top; v(top word) = 1.  Returns (v, residual)."""
    top = max(P.coldeg)
    (jtop,) = [j for j in range(P.ncols) if P.coldeg[j] == top]
    Dl = {}
    for (i, j), c in P.ent.items():
        val = S(c.constant_coefficient())
        for li, x in zip(lam, P.L.gens()):
            co = c.monomial_coefficient(x)
            if co:
                val += QQ(co) * li
        if val != 0:
            Dl[(i, j)] = val
    bycol = {}
    for (i, j), c in Dl.items():
        bycol.setdefault(j, []).append((i, c))
    rows_of, cols_of = {}, {}
    for i in range(P.nrows):
        rows_of.setdefault(P.rowdeg[i], []).append(i)
    for j in range(P.ncols):
        cols_of.setdefault(P.coldeg[j], []).append(j)
    v, acc = {jtop: S(1)}, {}

    def push(j):
        for i, c in bycol.get(j, ()):
            acc[i] = acc.get(i, S(0)) + c * v[j]
    push(jtop)
    residual = []
    for d in range(top - 1, -1, -1):
        rows, cols = rows_of.get(d, []), cols_of.get(d, [])
        if cols:
            rhs = [-acc.get(i, S(0)) for i in rows]
            ri = {i: k for k, i in enumerate(rows)}
            ci = {j: k for k, j in enumerate(cols)}
            Sig = matrix(QQ, len(rows), len(cols), sparse=True)
            for j in cols:
                for i, c in bycol.get(j, ()):
                    if i in ri:
                        Sig[ri[i], ci[j]] = QQ(c)
            Sig = Sig.dense_matrix()
            piv = Sig.transpose().pivots()
            assert len(piv) == len(cols), "symbol block not injective (symbol criterion fails)"
            Ainv = Sig.matrix_from_rows(piv).inverse()
            b = [rhs[k] for k in piv]
            for jj, j in enumerate(cols):
                val = S(0)
                for kk in range(len(cols)):
                    if Ainv[jj, kk]:
                        val += Ainv[jj, kk] * b[kk]
                if val != 0:
                    v[j] = val
                    push(j)
        for i in rows:
            if acc.get(i, S(0)) != 0:
                residual.append((P.rowkeys[i], acc[i]))
    return v, residual


def full_statement(g, beta):
    """The 'every lam' statement; returns a dict of the checks (see the module docstring)."""
    even = [k for k in g.simple if g.par[k] == 0]
    odds = [k for k in g.simple if g.par[k] == 1]
    P = Pencil(g, beta, simple=even)
    sym_ok, _ = symbol_criterion(P)
    S, lam = full_ring(g)
    v0, res = solve(P, S, lam) if sym_ok else ({}, ["symbol criterion fails"])
    V = Verma(g, lam)
    vec = {P.words[j]: c for j, c in v0.items()}
    rho = g.rho()
    form = S(g.pairing([lam[c] + rho[c] for c in range(g.r)], beta))
    evs = [c for k in odds for c in V.act(k, vec).values()]
    divisible = all((c % form) == 0 for c in evs)
    quot = [c // form for c in evs] if divisible else []
    const = any(q != 0 and q.degree() <= 0 for q in quot)
    return {"even_symbol_criterion": sym_ok, "even_system_consistent_for_all_lam": not res,
            "odd_action_divisible_by_atypicality": divisible, "quotient_has_constant_term": const,
            "atypicality_form": str(form),
            "holds": bool(sym_ok and not res and divisible and const),
            "v0_uses": sorted({str(x) for c in v0.values() for x in c.variables()})}


def check_at_points(P, S, lam, v, kind, s, t, npts=3, rng=None):
    """At random rational points of the hyperplane: v(lam) != 0, D(lam) v(lam) = 0, nullity 1;
    D(lam) is rebuilt from a fresh numeric Verma module.  Off the hyperplane: nullity 0."""
    rng = rng or random.Random(46830)
    g = P.g
    ok = True
    for _ in range(npts):
        x = {var: QQ(rng.randint(-40, 40)) / rng.randint(1, 9) for var in S.gens()}
        lamx = [QQ(l.subs(x)) for l in lam]
        vec = [QQ(v[j].subs(x)) if j in v else QQ(0) for j in range(P.ncols)]
        Vx = Verma(g, lamx)
        rows = {}
        for e in g.simple:
            for jw, w in enumerate(P.words):
                for u, c in Vx.act_word(e, w).items():
                    rows.setdefault((e, u), {})[jw] = c
        Mx = matrix(QQ, len(rows), P.ncols, sparse=True)
        for i, (_, rw) in enumerate(sorted(rows.items())):
            for j, c in rw.items():
                Mx[i, j] = c
        ok &= bool((Mx * matrix(QQ, vec).transpose()).is_zero() and any(vec) and P.ncols - Mx.rank() == 1)
    rho, off = g.rho(), set()
    while len(off) < npts:                    # points genuinely off the hyperplane
        y = [QQ(rng.randint(-60, 60)) / rng.randint(1, 11) for _ in range(g.r)]
        if g.pairing([y[c] + rho[c] for c in range(g.r)], P.beta) != 0:
            off.add((tuple(y), P.nullity(y)))
    return ok, sorted({n for _, n in off})
