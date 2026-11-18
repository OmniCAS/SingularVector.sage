"""Experiment `conjecture-8`: Conjecture 8 of Liu-Luo-Wang for osp, beta = delta_s + eps_t.

For every isotropic odd positive root beta of the chosen kind and every algebra:
  (a) symbol criterion (all simple roots, and the even ones alone)  => dim <= 1 for every lam;
  (b) exact polynomial solution on the hyperplane H                  => dim >= 1 on H;
  (c) the full statement: for EVERY lam, dim = 1 iff (lam + rho, beta) = 0, else 0;
  (d) checks at random rational points on and off H (fresh numeric Verma module).
Writes data/experiments/conjecture-8.json (or conjecture-8-minus.json for --kind -).

    python experiments/conjecture8.py [--all] [--kind +|-] [ALG ...]
"""
import argparse
import json
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from _common import ALL_OSP, DEFAULT_OSP, Timer, experiment_record, parse_algebra
from singular_vector import Pencil, check_at_points, full_statement, hyperplane_ring, solve, symbol_criterion


def run_root(g, k):
    T = Timer()
    w = g.wt[k]
    nm = g.weight_name(w)
    kind, s, t = g.parse_root(nm)
    P = Pencil(g, w)
    sym, _ = symbol_criterion(P)
    sym_even, _ = symbol_criterion(P, keep=lambda e: g.par[e] == 0)
    S, lam, p = hyperplane_ring(g, kind, s, t)
    v, res = solve(P, S, lam) if sym else ({}, ["-"])
    fs = full_statement(g, w)
    pts, off = check_at_points(P, S, lam, v, kind, s, t) if sym and not res else (False, None)
    return {"algebra": g.name, "beta": nm, "beta_latex": g.weight_name(w, latex=True),
            "hyperplane": "a%d = %s" % (s, lam[p]), "dim_weight_space": P.ncols, "equations": P.nrows,
            "top_degree": max(P.coldeg), "symbol_criterion": sym, "even_simple_roots_suffice": sym_even,
            "exactly_one_on_H": bool(sym and not res), "every_lambda": fs["holds"],
            "checked_at_points": pts, "nullity_off_H": off, "solution_degree": max((c.degree() for c in v.values()), default=None),
            "seconds": T()}


def worker(name, kind, out):
    g = parse_algebra(name)
    assert g.check_generation()
    rows = []
    for k in g.isotropic_odd_positive():
        if g.parse_root(g.weight_name(g.wt[k]))[0] != kind:
            continue
        rows.append(run_root(g, k))
    Path(out).write_text(json.dumps(rows, default=str))


def run_isolated(name, kind):
    """One algebra per process (long single processes have crashed inside libsingular's allocator)."""
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as fh:
        out = fh.name
    proc = subprocess.run([sys.executable, __file__, "--worker", name, "--kind", kind, "--out", out],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        return [{"algebra": name, "error": "worker exited with %d" % proc.returncode}]
    rows = json.loads(Path(out).read_text())
    Path(out).unlink()
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("algebras", nargs="*")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--kind", default="+", choices=["+", "-"])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--worker", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    if a.worker:
        return worker(a.worker, a.kind, a.out)
    names = a.algebras or (ALL_OSP if a.all else DEFAULT_OSP)
    rows = []
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        for name, part in zip(names, pool.map(lambda nm: run_isolated(nm, a.kind), names)):
            for row in part:
                rows.append(row)
                if "error" in row:
                    print("%-10s ERROR %s" % (name, row["error"]), flush=True)
                else:
                    print("%-10s %-6s dim %5d  every lam: %s  (%.1f s)" % (row["algebra"], row["beta"], row["dim_weight_space"],
                                                                         row["every_lambda"], row["seconds"]), flush=True)
    ok = all("error" not in r and r["exactly_one_on_H"] and r["every_lambda"] and r["checked_at_points"] for r in rows)
    biggest = max((r for r in rows if "error" not in r), key=lambda r: r["dim_weight_space"])
    sign = "+" if a.kind == "+" else "-"
    path = experiment_record(
        "conjecture-8" + ("" if a.kind == "+" else "-minus"),
        "Conjecture 8 for osp, beta = delta_s %s eps_t: exact computation" % sign,
        "Is the space of singular vectors of weight lam - beta in M(lam) one-dimensional for every beta-atypical lam?",
        "Symbol criterion (lam-independent ranks over QQ), exact solution over the polynomial ring of the hyperplane "
        "(constant pivots only), and the even-part reduction with divisibility by (lam + rho, beta); checks at random points.",
        "%d algebras, %d (algebra, root) cases; largest weight space %d (%s, %s)" % (
            len(names), len(rows), biggest["dim_weight_space"], biggest["algebra"], biggest["beta"]),
        "All cases: for every lam, dim = 1 iff (lam + rho, beta) = 0, else 0." if ok else "NOT all cases hold: see rows.",
        "confirmed" if ok else "refuted",
        "python experiments/conjecture8.py" + (" --all" if a.all else "") + (" --kind -" if a.kind == "-" else ""),
        rows)
    print("wrote", path)


if __name__ == "__main__":
    main()
