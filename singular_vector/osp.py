"""osp(2M|2N) and osp(2M+1|2N) in the conventions of the Lean file Defs/OSP.lean.

(Transcribed in OmniCAS/odd-singular-vector, lean/scripts/check_osp.py; same names and signs.)

  positions a = 0..N+M-1:  a < N  <->  delta_{N-a},   a >= N  <->  eps_{a-N+1}
  V has basis u(a) = a, v(a) = r + a (r = N + M) and z = 2r (only for osp(2M+1|2N));
  the delta-vectors u(a), v(a) (a < N) are odd, the others even;
  T a b = E[u a, u b] + sigT(a,b) E[v b, v a],  P a b = E[u a, v b] + tauP(a,b) E[u b, v a],
  Q a b = E[v a, u b] + tauQ(a,b) E[v b, u a],  Rp a = E[u a, z] + rhoP(a) E[z, v a],
  Rm a = E[v a, z] + rhoM(a) E[z, u a];
  Cartan H_a = T a a, weight coordinates lam_a = lam(H_a); (delta, delta) = -1, (eps, eps) = 1.

Root vectors are the generators themselves (T a b for a != b; P a b and Q a b for a <= b when
non-zero; Rp a, Rm a).  The weight lam = (a_N, ..., a_1 | b_1, ..., b_M).

The root vectors printed in Section 2.4 of Liu-Luo-Wang (paper notation, with the sign of
e_{delta_i - delta_j} corrected) are recovered by `paper_root_vector`.
"""
from sage.all__sagemath_singular import QQ

from .algebra import SuperAlgebra


def sigT(N, a, b):
    return 1 if (a < N and b >= N) else -1


def tauP(N, a, b):
    return -1 if (a >= N and b >= N) else 1


def tauQ(N, a, b):
    return 1 if (a < N and b < N) else -1


def rhoP(N, a):
    return 1 if a < N else -1


def rhoM(N, a):
    return -1


def generators(odd, N, M):
    """{generator label: matrix} for all generators of Defs/OSP.lean (zero matrices included)."""
    r = N + M
    u = lambda a: a
    v = lambda a: r + a
    z = 2 * r

    def E(*entries):
        m = {}
        for (i, j, c) in entries:
            if c:
                m[(i, j)] = m.get((i, j), 0) + c
        return {k: QQ(x) for k, x in m.items() if x}
    gens = {}
    for a in range(r):
        for b in range(r):
            gens[("T", a, b)] = E((u(a), u(b), 1), (v(b), v(a), sigT(N, a, b)))
            gens[("P", a, b)] = E((u(a), v(b), 1), (u(b), v(a), tauP(N, a, b)))
            gens[("Q", a, b)] = E((v(a), u(b), 1), (v(b), u(a), tauQ(N, a, b)))
        if odd:
            gens[("Rp", a)] = E((u(a), z, 1), (z, v(a), rhoP(N, a)))
            gens[("Rm", a)] = E((v(a), z, 1), (z, u(a), rhoM(N, a)))
    return gens


def gpar(N, g):
    if g[0] in ("T", "P", "Q"):
        return int((g[1] < N) != (g[2] < N))
    return int(g[1] < N)


class OSP(SuperAlgebra):
    def __init__(self, M_dim, two_N):
        """osp(M_dim | two_N), e.g. OSP(6, 2) = osp(6|2)."""
        odd, M, N = M_dim % 2 == 1, M_dim // 2, two_N // 2
        assert M >= 1 and N >= 1
        self.odd, self.N, self.M = odd, N, M
        self.dimV0 = M_dim
        r = N + M
        gens = generators(odd, N, M)
        self.gens = gens
        parV = [1 if a < N else 0 for a in range(r)] * 2 + ([0] if odd else [])
        cartan = [gens[("T", a, a)] for a in range(r)]
        roots = []
        for g, X in gens.items():
            if not X or (g[0] == "T" and g[1] == g[2]):
                continue
            if g[0] in ("P", "Q") and g[1] > g[2]:
                continue
            roots.append((g, X, gpar(N, g)))
        names = ["d%d" % (N - c) if c < N else "e%d" % (c - N + 1) for c in range(r)]
        latex = ["\\delta_{%d}" % (N - c) if c < N else "\\epsilon_{%d}" % (c - N + 1) for c in range(r)]
        signs = [-1 if c < N else 1 for c in range(r)]
        super().__init__("osp(%d|%d)" % (M_dim, two_N), parV, cartan, roots, signs, names, latex)
        expect = M_dim * (M_dim - 1) // 2 + N * (2 * N + 1) + 2 * N * M_dim
        assert self.dim == expect, (self.dim, expect)
        for k in range(self.r, self.dim):     # IsPos of Defs/OSP.lean = lexicographic positivity
            g = self.labels[k]
            assert ((g[0] == "T" and g[1] < g[2]) or g[0] in ("P", "Rp")) == (k in self.pos), g
        self.var_names = ["a%d" % (N - c) if c < N else "b%d" % (c - N + 1) for c in range(r)]

    # -- roots in the paper's notation -------------------------------------------------------
    def coord_of_delta(self, i):
        return self.N - i

    def coord_of_eps(self, k):
        return self.N + k - 1

    def beta(self, kind, s, t):
        """delta_s + eps_t (kind '+') or delta_s - eps_t (kind '-')."""
        w = [0] * self.r
        w[self.coord_of_delta(s)] = 1
        w[self.coord_of_eps(t)] = 1 if kind == "+" else -1
        return tuple(w)

    def parse_root(self, name):
        """'d1+e2' -> ('+', 1, 2)."""
        kind = "+" if "+" in name else "-"
        d, e = name.replace("+", "-").split("-")
        return kind, int(d[1:]), int(e[1:])

    def atypical_coordinate(self, kind, s, t):
        """The hyperplane (lam + rho, beta) = 0 of Defs/OSP.lean, solved for the coordinate of
        delta_s: returns (p, f) with lam[p] = f(lam).  '+': a_s - b_t + s + t - dim V_0 = 0;
        '-': a_s + b_t + s - t = 0."""
        p, q = self.coord_of_delta(s), self.coord_of_eps(t)
        if kind == "+":
            return p, (lambda lam: lam[q] - s - t + self.dimV0)
        return p, (lambda lam: -lam[q] - s + t)

    # -- the paper's root vectors (Section 2.4), for osp(2M|2N) ------------------------------------
    def paper_root_vector(self, gamma_name):
        """(basis index, scalar) with paper's e_gamma = scalar * generator, for the root vectors of
        Section 2.4 (osp(2m|2n), m = M, n = N; indices i-bar -> u(N-i), (i+n)-bar -> v(N-i),
        k -> u(N+k-1), k+m -> v(N+k-1)).  gamma_name like '-d1+e2', 'e3-e2', '-e2-e3', '2d1'."""
        assert not self.odd, "Section 2.4 prints the matrices for osp(2m|2n) only"
        X = paper_matrix(self.N, self.M, gamma_name)
        for k in range(self.r, self.dim):
            B = self.mats[k]
            if set(B) == set(X):
                (i, j) = next(iter(B))
                c = X[(i, j)] / B[(i, j)]
                if all(X[key] == c * B[key] for key in B):
                    return k, c
        raise ValueError("not a multiple of a generator: %s" % gamma_name)


def paper_matrix(N, M, name):
    """The matrix of e_gamma printed in Section 2.4 (e_{delta_i - delta_j} with the corrected sign),
    in the basis u, v of Defs/OSP.lean.  Paper indices: i-bar = delta_i-vector, (i+n)-bar, k, k+m."""
    n, m, r = N, M, N + M
    B = lambda i: N - i                      # i-bar
    Bn = lambda i: r + N - i                 # (i+n)-bar
    U = lambda k: N + k - 1                  # k
    Um = lambda k: r + N + k - 1             # k+m

    def lin(*terms):
        out = {}
        for c, i, j in terms:
            out[(i, j)] = out.get((i, j), 0) + QQ(c)
        return {k: x for k, x in out.items() if x}
    terms = parse_signed(name)
    ds = [(s, i) for s, x, i in terms if x == "d"]
    es = [(s, k) for s, x, k in terms if x == "e"]
    if len(ds) == 1 and not es and abs(ds[0][0]) == 2:     # +-2 delta_i
        s, i = ds[0]
        return lin((1, B(i), Bn(i))) if s > 0 else lin((1, Bn(i), B(i)))
    if len(ds) == 2:                                        # delta_i +- delta_j
        (s1, i), (s2, j) = sorted(ds, key=lambda t: -t[0])
        if s1 > 0 and s2 > 0:
            return lin((1, B(i), Bn(j)), (1, B(j), Bn(i)))
        if s1 < 0 and s2 < 0:
            return lin((1, Bn(j), B(i)), (1, Bn(i), B(j)))
        return lin((1, B(i), B(j)), (-1, Bn(j), Bn(i)))     # delta_i - delta_j (corrected sign)
    if len(es) == 2:                                        # eps_k +- eps_l
        (s1, k), (s2, l) = sorted(es, key=lambda t: -t[0])
        if s1 > 0 and s2 < 0:
            return lin((1, U(k), U(l)), (-1, Um(l), Um(k)))
        k, l = sorted((k, l))
        if s1 > 0:
            return lin((1, U(k), Um(l)), (-1, U(l), Um(k)))
        return lin((1, Um(l), U(k)), (-1, Um(k), U(l)))
    (sd, i), (se, k) = ds[0], es[0]                         # +-delta_i +- eps_k
    if sd > 0 and se > 0:
        return lin((1, U(k), Bn(i)), (1, B(i), Um(k)))
    if sd < 0 and se < 0:
        return lin((1, Um(k), B(i)), (-1, Bn(i), U(k)))
    if sd > 0 and se < 0:
        return lin((1, Um(k), Bn(i)), (1, B(i), U(k)))
    return lin((1, U(k), B(i)), (-1, Bn(i), Um(k)))


def parse_signed(name):
    """'-d1+e2' -> [(-1, 'd', 1), (1, 'e', 2)];  '2d1' -> [(2, 'd', 1)]."""
    out, i, sign = [], 0, 1
    while i < len(name):
        ch = name[i]
        if ch in "+-":
            sign = 1 if ch == "+" else -1
            i += 1
            continue
        mult = 1
        if ch.isdigit():
            mult = int(ch)
            i += 1
            ch = name[i]
        j = i + 1
        while j < len(name) and name[j].isdigit():
            j += 1
        out.append((sign * mult, ch, int(name[i + 1:j])))
        sign = 1
        i = j
    return out
