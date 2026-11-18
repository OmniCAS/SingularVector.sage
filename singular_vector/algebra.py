"""Matrix Lie superalgebras with a triangular decomposition, and U(n^-) in a PBW basis.

A `SuperAlgebra` is given by homogeneous matrices: r diagonal Cartan matrices H_0..H_{r-1} and root
vectors (one per root; every root space here is one-dimensional).  Brackets are matrix super
commutators, decomposed in the given basis; nothing about roots or structure constants is typed in.
Weights are the eigenvalues of ad(H_c); a root is positive iff its first non-zero coordinate is
positive (the coordinate order of the builders makes this the distinguished Borel).

PBW order of the negative root vectors: even ones first, then odd ones; inside each parity by
height, then by weight.  A word is a non-decreasing tuple of positions in that order; an odd position
occurs at most once (x x = 1/2 [x, x] is rewritten).
"""
from sage.all__sagemath_singular import QQ, ZZ, matrix, vector


def mat_mul(A, B):
    out = {}
    rows = {}
    for (j, k), b in B.items():
        rows.setdefault(j, []).append((k, b))
    for (i, j), a in A.items():
        for k, b in rows.get(j, ()):
            out[(i, k)] = out.get((i, k), 0) + a * b
    return {k: x for k, x in out.items() if x}


class SuperAlgebra:
    """cartan: list of r diagonal matrices {(i, i): c}; roots: list of (label, matrix, parity);
    parV: parities of the basis of V (to check homogeneity); signs: the bilinear form on weights is
    (x, y) = sum_c signs[c] x_c y_c; coord_names: names of the weight coordinates."""

    def __init__(self, name, parV, cartan, roots, signs, coord_names, latex_coords=None):
        self.name = name
        self.parV = parV
        self.r = len(cartan)
        self.signs = list(signs)
        self.coord_names = list(coord_names)
        self.latex_coords = list(latex_coords or coord_names)
        self.labels = [("H", c) for c in range(self.r)] + [lab for lab, _, _ in roots]
        self.mats = [dict(H) for H in cartan] + [dict(X) for _, X, _ in roots]
        self.par = [0] * self.r + [p for _, _, p in roots]
        self.dim = len(self.mats)
        self.index = {lab: i for i, lab in enumerate(self.labels)}
        for k in range(self.r, self.dim):
            assert all((parV[i] + parV[j]) % 2 == self.par[k] for (i, j) in self.mats[k]), self.labels[k]
        D = len(parV)
        diag = [[self.mats[c].get((i, i), 0) for i in range(D)] for c in range(self.r)]
        self.wt = [(0,) * self.r] * self.r
        for k in range(self.r, self.dim):
            ws = {tuple(diag[c][i] - diag[c][j] for c in range(self.r)) for (i, j) in self.mats[k]}
            assert len(ws) == 1, ("not a weight vector", self.labels[k])
            self.wt.append(tuple(int(x) for x in ws.pop()))
        self.byweight = {}
        for k in range(self.r, self.dim):
            assert self.wt[k] not in self.byweight, "root space of dimension > 1"
            self.byweight[self.wt[k]] = k
        self._check_independent(D)
        ispos = lambda w: next(x for x in w if x != 0) > 0
        self.pos = [k for k in range(self.r, self.dim) if ispos(self.wt[k])]
        self.negs = [k for k in range(self.r, self.dim) if not ispos(self.wt[k])]
        self._br = {}
        self._ml = {}
        self._simple_and_heights()
        self._pbw_order()

    # -- structure ----------------------------------------------------------------------------
    def _check_independent(self, D):
        rows = []
        for X in self.mats:
            v = [0] * (D * D)
            for (i, j), c in X.items():
                v[i * D + j] = c
            rows.append(v)
        assert matrix(QQ, rows).rank() == self.dim, "the matrices are linearly dependent"

    def bracket(self, x, y):
        """[x, y] of basis elements as {basis index: QQ} (matrix super commutator)."""
        key = (x, y)
        if key in self._br:
            return self._br[key]
        s = -1 if (self.par[x] and self.par[y]) else 1
        C = dict(mat_mul(self.mats[x], self.mats[y]))
        for k, c in mat_mul(self.mats[y], self.mats[x]).items():
            C[k] = C.get(k, 0) - s * c
        C = {k: c for k, c in C.items() if c}
        res = {}
        if C:
            w = tuple(a + b for a, b in zip(self.wt[x], self.wt[y]))
            if any(w):
                z = self.byweight[w]
                assert self.par[z] == (self.par[x] + self.par[y]) % 2
                (i, j), zc = next(iter(self.mats[z].items()))
                c = C.get((i, j), 0) / zc
                assert {k: c * v for k, v in self.mats[z].items()} == C, ("bracket leaves the span", x, y)
                res = {z: c}
            else:
                res = self._cartan_coordinates(C)
        self._br[key] = res
        return res

    def _cartan_coordinates(self, C):
        if not hasattr(self, "_Hsolve"):
            idx = sorted({i for c in range(self.r) for (i, _) in self.mats[c]})
            A = matrix(QQ, [[self.mats[c].get((i, i), 0) for c in range(self.r)] for i in idx])
            self._Hsolve = (idx, A)
        idx, A = self._Hsolve
        assert all(i == j for (i, j) in C), "bracket of weight 0 is not diagonal"
        b = vector(QQ, [C.get((i, i), 0) for i in idx])
        cs = A.solve_right(b)
        assert A * cs == b, "bracket leaves the span of the Cartan"
        return {c: QQ(x) for c, x in enumerate(cs) if x}

    def check_closure(self):
        for x in range(self.dim):
            for y in range(self.dim):
                self.bracket(x, y)
        return True

    def check_generation(self):
        """The simple root vectors generate n^+."""
        reached, frontier = set(self.simple), list(self.simple)
        while frontier:
            new = []
            for x in frontier:
                for e in self.simple:
                    for z, c in self.bracket(e, x).items():
                        if c and z not in reached and z in self.pos:
                            reached.add(z)
                            new.append(z)
            frontier = new
        return reached == set(self.pos)

    # -- roots --------------------------------------------------------------------------------
    def _simple_and_heights(self):
        P = {self.wt[k] for k in self.pos}
        self.simple = [k for k in self.pos
                       if not any(tuple(a - b for a, b in zip(self.wt[k], p)) in P for p in P if p != self.wt[k])]
        self._S = matrix(QQ, [list(self.wt[k]) for k in self.simple])
        assert self._S.rank() == len(self.simple), "simple roots are linearly dependent"
        self.height, self.simple_coeffs = {}, {}
        for k in range(self.r, self.dim):
            co = self.simple_coordinates(self.wt[k])
            assert all(c in ZZ for c in co)
            sgn = 1 if k in self.pos else -1
            assert all(sgn * c >= 0 for c in co), ("not a combination of simple roots", self.labels[k])
            self.simple_coeffs[k] = tuple(int(sgn * c) for c in co)
            self.height[k] = int(sum(sgn * c for c in co))

    def simple_coordinates(self, w):
        """The coefficients of w in the simple roots (w must lie in their span; gl has a centre)."""
        co = self._S.solve_left(vector(QQ, list(w)))
        assert co * self._S == vector(QQ, list(w)), "not in the span of the simple roots"
        return co

    def _pbw_order(self):
        key = lambda k: (self.par[k], self.height[k], tuple(-x for x in self.wt[k]))
        self.neg = sorted(self.negs, key=key)
        self.negpos = {k: i for i, k in enumerate(self.neg)}
        self.negwt = [self.wt[k] for k in self.neg]
        self.negpar = [self.par[k] for k in self.neg]
        self.kind = {}
        for k in range(self.r):
            self.kind[k] = "h"
        for k in self.pos:
            self.kind[k] = "+"
        for k in self.neg:
            self.kind[k] = "-"

    def pairing(self, x, y):
        return sum(s * a * b for s, a, b in zip(self.signs, x, y))

    def rho(self):
        out = [QQ(0)] * self.r
        for k in self.pos:
            sgn = -1 if self.par[k] else 1
            for c, x in enumerate(self.wt[k]):
                out[c] += QQ(sgn * x) / 2
        return out

    def coroot_pairing(self, mu, gamma):
        """(mu, gamma^vee): 2(mu, gamma)/(gamma, gamma), or (mu, gamma) for isotropic gamma."""
        nn = self.pairing(gamma, gamma)
        val = self.pairing(mu, gamma)
        return val * (QQ(2) / nn) if nn != 0 else val

    def isotropic_odd_positive(self):
        return [k for k in self.pos if self.par[k] and self.pairing(self.wt[k], self.wt[k]) == 0]

    # -- names --------------------------------------------------------------------------------
    def weight_name(self, w, latex=False):
        names = self.latex_coords if latex else self.coord_names
        parts = []
        for c, x in enumerate(w):
            if x:
                parts.append(("+" if x > 0 else "-") + ("" if abs(x) == 1 else str(abs(x))) + names[c])
        return "".join(parts).lstrip("+") or "0"

    def neg_name(self, i, latex=False):
        gamma = tuple(-x for x in self.negwt[i])
        return ("f_{%s}" if latex else "f[%s]") % self.weight_name(gamma, latex)

    def word_name(self, word, latex=False):
        out, k = [], 0
        while k < len(word):
            j = k
            while j < len(word) and word[j] == word[k]:
                j += 1
            e = j - k
            nm = self.neg_name(word[k], latex)
            out.append(nm + (("^{%d}" if latex else "^%d") % e if e > 1 else ""))
            k = j
        return (" " if not latex else " ").join(out) if out else "1"

    # -- U(n^-) -------------------------------------------------------------------------------
    def words_of_weight(self, target):
        """All PBW words of weight `target` (minus a sum of positive roots), longest first."""
        n = len(self.neg)
        if not hasattr(self, "_negco"):
            self._negco = [tuple(-c for c in self.simple_coeffs[k]) for k in self.neg]
        co_t = self.simple_coordinates(target)
        assert all(c in ZZ and c <= 0 for c in co_t), "target is not minus a sum of positive roots"
        memo = {}

        def rec(i, rem):
            key = (i, rem)
            if key in memo:
                return memo[key]
            if not any(rem):
                res = [()]
            elif i == n:
                res = []
            else:
                res = list(rec(i + 1, rem))
                cur, k = rem, 0
                while True:
                    cur = tuple(a - b for a, b in zip(cur, self._negco[i]))
                    k += 1
                    if any(a > 0 for a in cur):
                        break
                    res.extend((i,) * k + tail for tail in rec(i + 1, cur))
                    if self.negpar[i]:
                        break
            memo[key] = res
            return res
        out = rec(0, tuple(int(c) for c in co_t))
        out.sort(key=lambda w: (-len(w), w))
        return out

    def word_weight(self, word):
        w = [0] * self.r
        for i in word:
            for c, x in enumerate(self.negwt[i]):
                w[c] += x
        return tuple(w)

    def mul_left(self, i, word):
        """f_i * word in U(n^-) as {word: QQ} (i a PBW position).  Memoised."""
        key = (i, word)
        res = self._ml.get(key)
        if res is not None:
            return res
        if not word or i < word[0]:
            res = {(i,) + word: QQ(1)}
        else:
            j, rest = word[0], word[1:]
            if i == j:
                if not self.negpar[i]:
                    res = {(i,) + word: QQ(1)}
                else:
                    res = {}
                    a = self.neg[i]
                    for c, cc in self.bracket(a, a).items():            # f f = 1/2 [f, f]
                        for w, cw in self.mul_left(self.negpos[c], rest).items():
                            res[w] = res.get(w, 0) + cc * cw / 2
            else:
                s = -1 if (self.negpar[i] and self.negpar[j]) else 1
                res = {}
                for w, cw in self.mul_left(i, rest).items():
                    for w2, c2 in self.mul_left(j, w).items():
                        res[w2] = res.get(w2, 0) + s * cw * c2
                for c, cc in self.bracket(self.neg[i], self.neg[j]).items():
                    for w, cw in self.mul_left(self.negpos[c], rest).items():
                        res[w] = res.get(w, 0) + cc * cw
            res = {w: x for w, x in res.items() if x}
        self._ml[key] = res
        return res
