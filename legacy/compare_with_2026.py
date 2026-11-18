"""Run the 2022 tools (legacy/2022/src, unchanged) on small cases and compare with the 2026 package.

The old .sage files are preparsed and executed here, unchanged, in a private namespace.  Its classes:
SpO(m, n) = spo(2m|2n) (m = rank of sp = number of deltas, n = rank of o = number of eps; only even
orthogonal part, n >= 2), structure constants from matrices (AlgebraBSC.lie_sup_SC), super PBW
straightening (SPBWElement.tuple2dict), Verma projection (to_verma).  The repository has no driver for the
singular vector equations (it exported them to MATLAB via rex.py); the driver below is new:

  1. candidate words of weight -beta from SpO.positive_dec (the repository's 'P function'),
  2. e_i * word for the simple root vectors e_i, straightened with SPBWElement, projected with to_verma
     at a symbolic weight on the hyperplane,
  3. exact kernel over QQ(b...) (no MATLAB),
  4. the solution, written in the old root vectors, is pushed into the new engine (each old matrix is
     +-1 times a generator of Defs/OSP.lean after the index map below) and compared with the new
     singular vector (must be proportional).

Index map old -> new (V = C^{2m|2n} in gl(2m+2n); new: u(a) = a, v(a) = r + a, r = m + n):
  old k < m (delta_{k+1}) -> u(k);  m <= k < 2m -> v(k-m);  2m <= k < 2m+n (eps_{k-2m+1}) -> u(k-m);
  2m+n <= k -> v(k-m-n).  Old delta_{i} is new delta_{m+1-i} (the old odd simple root is delta_m - eps_1).

    python legacy/compare_with_2026.py [--points] [m n s t ...]   (SpO(m, n); s, t in the NEW notation)

Needs the symbolic ring: pip install passagemath-symbolics==10.8.14.
"""
import json
import sys
import time
from pathlib import Path

from sage.all__sagemath_symbolics import *  # noqa: F401,F403  (the old code expects the Sage globals)
from sage.repl.preparse import preparse_file
from sage.symbolic.expression import Expression

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from singular_vector import OSP, Pencil, Verma, hyperplane_ring, solve  # noqa: E402

OLD = ROOT / "legacy" / "2022" / "src"


def load_old():
    ns = dict(globals())
    ns["Expression"] = Expression
    for f in ["Lie.sage", "MatSp.sage", "AlgebraBSC.sage", "pbw.sage", "pbws.sage", "spo.sage"]:
        exec(compile(preparse_file(open(OLD / f).read()), str(OLD / f), "exec"), ns)
    return ns


def run(m, n, s, t, mode="symbolic"):
    out = {"old": "SpO(%d,%d) = spo(%d|%d)" % (m, n, 2 * m, 2 * n), "new": "osp(%d|%d)" % (2 * n, 2 * m),
           "beta_new_notation": "d%d+e%d" % (s, t)}
    ns = load_old()
    t0 = time.time()
    L = ns["SpO"](m, n)
    P = len(L.positive_matrixs)
    r = m + n
    mats = list(L.negative_mats) + list(L.h_matrixs) + list(L.positive_matrixs)
    odd = list(L.odd_mats)
    N = ns["AlgebraBSC"].lie_sup_SC(mats, odd)
    syms = tuple(var("f%d" % i) for i in range(P)) + tuple(var("h%d" % i) for i in range(r)) + \
        tuple(var("e%d" % i) for i in range(P))
    odds = [syms[i] for i, M_ in enumerate(mats) if M_ in odd]
    SP = lambda d: ns["SPBWElement"](d, N, syms, odds)
    out["t_old_algebra"] = round(time.time() - t0, 1)
    # beta: new delta_s = old delta_{m+1-s}
    d, e = L.basis[:m], L.basis[m:]
    beta_vec = d[m - s] + e[t - 1]
    beta_root = L.vector2root(beta_vec)
    decs = L.positive_dec(beta_root)
    words = set()
    for co in decs:
        w = []
        ok = True
        for i, c in enumerate(co):
            if c and L.positive_matrixs[i] in odd and c > 1:
                ok = False
            w += [syms[i]] * int(c)
        if ok:
            words.add(tuple(w))
    words = sorted(words, key=lambda w: (-len(w), [syms.index(x) for x in w]))
    out["old_candidate_words"] = len(words)
    # weight on the hyperplane, in NEW coordinates: coordinate c of the new engine; old delta_{k+1} <-> new coordinate k
    g = OSP(2 * n, 2 * m)
    Snew, lam_new, p = hyperplane_ring(g, "+", s, t)
    bnames = [str(x) for x in Snew.gens()]
    bvars = {nm: var(nm) for nm in bnames}
    lam_sr = [SR(str(c)) for c in lam_new]                       # symbolic copies
    weight = sum(lam_sr[k] * d[k] for k in range(m)) + sum(lam_sr[m + j] * e[j] for j in range(n))
    # the system with the old engine
    t1 = time.time()
    rows = {}
    for i in range(r):
        ei = SP({(syms[P + r + i],): 1})
        for jw, w in enumerate(words):
            res = (ei * SP({w: 1})).to_verma(weight, L, r)
            for key, val in res.element.items():
                rows.setdefault((i, key), {})[jw] = val
    out["old_equations"] = len(rows)
    out["t_old_system"] = round(time.time() - t1, 1)
    K = Snew.fraction_field()
    rowlist = sorted(rows.items(), key=str)
    entries = [(ri, jw, Snew(str(SR(val).expand()))) for ri, (key, row) in enumerate(rowlist) for jw, val in row.items()]
    if mode == "symbolic":
        A = matrix(K, len(rows), len(words), sparse=True)
        for ri, jw, val in entries:
            A[ri, jw] = K(val)
        ker = A.right_kernel().basis()
        out["old_kernel_dim_on_H"] = len(ker)
    out["t_old_solve"] = round(time.time() - t1, 1)
    # push into the new engine
    def new_index(k):
        if k < m:
            return k
        if k < 2 * m:
            return r + (k - m)
        if k < 2 * m + n:
            return k - m
        return r + (k - m - n)
    def to_new(Mx):
        X = {}
        for (i, j), c in Mx.dict().items():
            X[(new_index(i), new_index(j))] = QQ(c)
        for kk in range(g.r, g.dim):
            B = g.mats[kk]
            if set(B) == set(X):
                (i0, j0) = next(iter(B))
                c = X[(i0, j0)] / B[(i0, j0)]
                if all(X[q] == c * B[q] for q in B):
                    return kk, c
        raise ValueError("old matrix is not a multiple of a new generator")
    fmap = [to_new(Mx) for Mx in L.negative_mats]
    out["old_root_vectors_are_pm_new_generators"] = sorted({str(c) for _, c in fmap})
    Pn = Pencil(g, g.beta("+", s, t))
    v_new, res = solve(Pn, Snew, lam_new)
    def mapped(kv, Ring, lamv):
        Vx = Verma(g, [Ring(x) for x in lamv])
        vec = {}
        for jw, w in enumerate(words):
            if kv[jw] == 0:
                continue
            idxs = [fmap[syms.index(x)][0] for x in w]
            scal = QQ(1)
            for x in w:
                scal *= fmap[syms.index(x)][1]
            for ww, c in Vx.apply_product(idxs).items():
                vec[ww] = vec.get(ww, Ring(0)) + Ring(c) * scal * kv[jw]
        return vec

    def proportional(vec, newvec, Ring):
        keys = set(vec) | set(newvec)
        ratios = {vec.get(k, Ring(0)) / newvec[k] for k in keys if newvec.get(k, 0) != 0}
        ok = len(ratios) == 1 and all((vec.get(k, Ring(0)) == 0) == (newvec.get(k, Ring(0)) == 0) for k in keys)
        return ok, (ratios.pop() if len(ratios) == 1 else None)
    if mode == "symbolic":
        vec = mapped(ker[0], K, lam_new)
        newvec = {Pn.words[j]: K(c) for j, c in v_new.items()}
        ok, ratio = proportional(vec, newvec, K)
        out["old_solution_proportional_to_new"] = ok
        out["ratio"] = str(ratio)
    else:
        import random
        rng = random.Random(2022)
        checks = []
        for _ in range(3):
            x = {v_: QQ(rng.randint(-40, 40)) / rng.randint(1, 9) for v_ in Snew.gens()}
            A = matrix(QQ, len(rows), len(words), sparse=True)
            for ri, jw, val in entries:
                A[ri, jw] = QQ(val.subs(x))
            kerx = A.right_kernel().basis()
            lamx = [QQ(c.subs(x)) for c in lam_new]
            if len(kerx) != 1:
                checks.append(("kernel dim", len(kerx)))
                continue
            vec = mapped(kerx[0], QQ, lamx)
            newvec = {Pn.words[j]: QQ(c.subs(x)) for j, c in v_new.items()}
            checks.append(proportional(vec, newvec, QQ))
        out["old_kernel_dim_on_H"] = "1 at 3 random points" if all(c[0] is True for c in checks) else str(checks)
        out["old_solution_proportional_to_new"] = all(c[0] is True for c in checks)
        out["ratio"] = sorted({str(c[1]) for c in checks})
    out["new_dim_weight_space"] = Pn.ncols
    out["t_total"] = round(time.time() - t0, 1)
    return out


if __name__ == "__main__":
    args = sys.argv[1:]
    mode = "symbolic"
    if "--points" in args:
        mode = "points"
        args.remove("--points")
    a = [int(x) for x in args] or [1, 3, 1, 1]
    for i in range(0, len(a), 4):
        res = run(*a[i:i + 4], mode=mode)
        res["mode"] = mode
        print(json.dumps(res, default=str), flush=True)
        with open(ROOT / "data" / "experiments" / "legacy-comparison.jsonl", "a") as fh:
            fh.write(json.dumps(res, default=str) + "\n")
