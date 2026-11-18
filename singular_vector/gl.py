"""gl(m|n) with the standard Borel, in the conventions of Defs/GL.lean (OmniCAS/odd-singular-vector).

Positions 0..m+n-1; p < m is even and is the paper's index (m-p)-bar (delta_{m-p}); position m+j-1 is
the index j (eps_j).  Root vectors E[p, q] (p != q), Cartan H_p = E[p, p]; (delta, delta) = 1,
(eps, eps) = -1 (the paper's normalisation).  The formula S_{-beta} of Theorem 1, (3)-(7), is
`S_formula`; it is checked against the computed singular vector in experiments/theorem1.py.
"""
from itertools import combinations

from sage.all__sagemath_singular import QQ

from .algebra import SuperAlgebra


class GL(SuperAlgebra):
    def __init__(self, m, n):
        self.m, self.n = m, n
        N = m + n
        parV = [0] * m + [1] * n
        cartan = [{(p, p): QQ(1)} for p in range(N)]
        roots = [(("E", p, q), {(p, q): QQ(1)}, (parV[p] + parV[q]) % 2)
                 for p in range(N) for q in range(N) if p != q]
        names = ["d%d" % (m - p) if p < m else "e%d" % (p - m + 1) for p in range(N)]
        latex = ["\\delta_{%d}" % (m - p) if p < m else "\\epsilon_{%d}" % (p - m + 1) for p in range(N)]
        signs = [1] * m + [-1] * n
        super().__init__("gl(%d|%d)" % (m, n), parV, cartan, roots, signs, names, latex)
        self.var_names = ["a%d" % (m - p) if p < m else "b%d" % (p - m + 1) for p in range(N)]
        self.N, self.M = m, n          # 'delta' coordinates first, as for OSP

    def coord_of_delta(self, i):
        return self.m - i

    def coord_of_eps(self, k):
        return self.m + k - 1

    def beta(self, kind, s, t):
        assert kind == "-", "the odd positive roots of gl(m|n) are delta_s - eps_t"
        w = [0] * self.r
        w[self.coord_of_delta(s)] = 1
        w[self.coord_of_eps(t)] = -1
        return tuple(w)

    def parse_root(self, name):
        d, e = name.split("-")
        return "-", int(d[1:]), int(e[1:])

    def atypical_coordinate(self, kind, s, t):
        """Atypical of Defs/GL.lean: lam_p + lam_q + (m - p) - (q - m + 1) = 0, solved for lam_p."""
        p, q = self.coord_of_delta(s), self.coord_of_eps(t)
        return p, (lambda lam: -lam[q] - (self.m - p) + (q - self.m + 1))

    def E(self, p, q):
        return self.index[("E", p, q)]


def c_coeff(lam, m, p, q, r):
    """c_r of (3)-(4) in positions (cCoeff of Defs/GL.lean): r < m: lam_p - lam_r + (r - p - 1);
    r >= m: lam_q - lam_r - (q - r)."""
    if r < m:
        return lam[p] - lam[r] + (r - p - 1)
    return lam[q] - lam[r] - (q - r)


def S_formula(V, g, p, q, c=c_coeff, chain=None):
    """S_{-beta} v = sum_J d_J E_J v of Theorem 1, E_J = E[q, j_k] ... E[j_1, p] (J increasing),
    d_J = prod_{r in (p, q) not in J} c_r.  `chain(a, b)` gives (basis index, scalar) of the root vector
    taking position b to position a (default: E[a, b] of gl); for osp the same formula with T a b."""
    chain = chain or (lambda a, b: (g.E(a, b), QQ(1)))
    lam = V.lam
    inner = range(p + 1, q)
    out = {}
    for k in range(len(inner) + 1):
        for J in combinations(inner, k):
            d = V.R(1)
            for r in inner:
                if r not in J:
                    d *= c(lam, g.N, p, q, r)
            if d == 0:
                continue
            lst = [q] + sorted(J, reverse=True) + [p]
            vec = {(): V.R(1)}
            for i in reversed(range(len(lst) - 1)):
                x, sc = chain(lst[i], lst[i + 1])
                vec = {w: sc * c_ for w, c_ in V.act(x, vec).items()}
            for w, c_ in vec.items():
                out[w] = out.get(w, 0) + d * c_
    return {w: c_ for w, c_ in out.items() if c_ != 0}
