"""Verma modules M(lam) = U(n^-) v and the singular system D(lam) v = 0 as an affine pencil."""
import random

from sage.all__sagemath_singular import QQ, PolynomialRing, matrix


class Verma:
    """M(lam); `lam` is a list of r elements of a commutative ring R containing QQ."""

    def __init__(self, g, lam):
        self.g = g
        self.lam = list(lam)
        self.R = lam[0].parent() if hasattr(lam[0], "parent") else QQ
        self._act = {}

    def act_word(self, x, word):
        """x . (word v) for a basis index x of g, as {word: coefficient}."""
        key = (x, word)
        res = self._act.get(key)
        if res is not None:
            return res
        g = self.g
        kind = g.kind[x]
        if kind == "-":
            res = dict(g.mul_left(g.negpos[x], word))
        elif kind == "h":
            s = self.lam[x] + g.word_weight(word)[x]
            res = {word: s} if s != 0 else {}
        else:
            res = {}
            if word:
                j, rest = word[0], word[1:]
                s = -1 if (g.par[x] and g.negpar[j]) else 1
                for w, cw in self.act_word(x, rest).items():          # s f_j (x . rest)
                    for w2, c2 in g.mul_left(j, w).items():
                        res[w2] = res.get(w2, 0) + s * c2 * cw
                for c, cc in g.bracket(x, g.neg[j]).items():          # [x, f_j] . rest
                    for w, cw in self.act_word(c, rest).items():
                        res[w] = res.get(w, 0) + cc * cw
                res = {w: c for w, c in res.items() if c != 0}
        self._act[key] = res
        return res

    def act(self, x, vec):
        out = {}
        for word, c in vec.items():
            for w, d in self.act_word(x, word).items():
                out[w] = out.get(w, 0) + c * d
        return {w: c for w, c in out.items() if c != 0}

    def apply_product(self, factors, vec=None):
        """x_1 ... x_k . vec (rightmost first); default vec = the highest weight vector."""
        v = {(): self.R(1)} if vec is None else vec
        for x in reversed(factors):
            v = self.act(x, v)
        return v


def self_check(g, samples=200, rng=None):
    """Number of failures of x.(y.w) - (-1)^{|x||y|} y.(x.w) = [x,y].w on random words."""
    rng = rng or random.Random(1903)
    lam = [QQ(rng.randint(-9, 9)) / rng.randint(1, 4) for _ in range(g.r)]
    V = Verma(g, lam)
    n, bad = len(g.neg), 0
    for _ in range(samples):
        word = tuple(sorted(rng.randrange(n) for _ in range(rng.randint(0, 3))))
        if any(word[k] == word[k + 1] and g.negpar[word[k]] for k in range(len(word) - 1)):
            continue
        x, y = rng.randrange(g.dim), rng.randrange(g.dim)
        lhs = V.act(x, V.act_word(y, word))
        s = -1 if (g.par[x] and g.par[y]) else 1
        for w, c in V.act(y, V.act_word(x, word)).items():
            lhs[w] = lhs.get(w, 0) - s * c
        rhs = {}
        for z, cz in g.bracket(x, y).items():
            for w, c in V.act_word(z, word).items():
                rhs[w] = rhs.get(w, 0) + cz * c
        if any(lhs.get(w, 0) != rhs.get(w, 0) for w in set(lhs) | set(rhs)):
            bad += 1
    return bad


class Pencil:
    """The system 'the simple root vectors in `simple` kill v' on the PBW basis of weight lam - beta,
    computed once with symbolic lam: D(lam) = D0 + sum_c lam_c D_c.
    Rows (e, u): e a simple root vector, u a PBW word; columns: the PBW words of weight -beta."""

    def __init__(self, g, beta, simple=None):
        self.g, self.beta = g, tuple(beta)
        self.L = PolynomialRing(QQ, ["L%d" % c for c in range(g.r)])
        V = Verma(g, list(self.L.gens()))
        self.words = g.words_of_weight(tuple(-x for x in beta))
        self.col = {w: i for i, w in enumerate(self.words)}
        self.simple = list(g.simple if simple is None else simple)
        rows, ent = {}, {}
        for e in self.simple:
            for w in self.words:
                for u, c in V.act_word(e, w).items():
                    rows.setdefault((e, u), len(rows))
                    ent[(rows[(e, u)], self.col[w])] = self.L(c)
        self.rowkeys = sorted(rows, key=lambda k: (-len(k[1]), k[0], k[1]))
        perm = {rows[k]: i for i, k in enumerate(self.rowkeys)}
        self.ent = {(perm[i], j): c for (i, j), c in ent.items()}
        self.nrows, self.ncols = len(self.rowkeys), len(self.words)
        self.rowdeg = [len(k[1]) for k in self.rowkeys]
        self.coldeg = [len(w) for w in self.words]
        self.D0, self.Dc = {}, [dict() for _ in range(g.r)]
        for k, c in self.ent.items():
            assert c.degree() <= 1, "D(lam) is not affine in lam"
            if c.constant_coefficient():
                self.D0[k] = QQ(c.constant_coefficient())
            for i, x in enumerate(self.L.gens()):
                if c.monomial_coefficient(x):
                    self.Dc[i][k] = QQ(c.monomial_coefficient(x))

    def at(self, lam):
        """D(lam) as a dense matrix over QQ."""
        Mx = matrix(QQ, self.nrows, self.ncols, sparse=True)
        for (i, j), c in self.D0.items():
            Mx[i, j] += c
        for li, Di in zip(lam, self.Dc):
            if li:
                for (i, j), c in Di.items():
                    Mx[i, j] += li * c
        return Mx.dense_matrix()

    def nullity(self, lam, rows=None):
        Mx = self.at(lam)
        if rows is not None:
            if not rows:
                return self.ncols
            Mx = Mx.matrix_from_rows(rows)
        return int(self.ncols - Mx.rank())
