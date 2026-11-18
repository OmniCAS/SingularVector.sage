"""Machine-readable output: polynomials with LaTeX, case records, experiment records.

Experiment records follow the 'exploration log' entries of the sister site (question, method, scale,
result, outcome, reproduce) so that a site can show them directly; case records carry the data of one
(algebra, root) with LaTeX strings ready for KaTeX.
"""
import json
import platform
import time
from datetime import datetime, timezone
from pathlib import Path

from sage.misc.latex import latex

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


def software():
    import sage.version
    return {"sage_library": "passagemath " + sage.version.version, "python": platform.python_version()}


def pretty_factor(c):
    """(text, latex) of a factorisation with every factor normalised to a positive leading coefficient."""
    if c == 0:
        return "0", "0"
    if c.degree() <= 0:
        return str(c), latex(c)
    F = c.factor()
    unit = F.unit()
    facs = []
    for f, e in F:
        if f.lc() < 0:
            f, unit = -f, unit * (-1) ** e
        facs.append((f, e))
    facs.sort(key=lambda fe: (fe[0].degree(), str(fe[0])))
    def part(f, e, tex):
        s = latex(f) if tex else str(f)
        if len(f.monomials()) > 1 and (len(facs) > 1 or e > 1 or unit != 1):
            s = "(%s)" % s
        return s + (("^{%d}" % e if tex else "^%d" % e) if e > 1 else "")
    body_t = "*".join(part(f, e, False) for f, e in facs)
    body_l = " ".join(part(f, e, True) for f, e in facs)
    if unit == 1:
        return body_t, body_l
    if unit == -1:
        return "-" + body_t, "-" + body_l
    return "%s*%s" % (unit, body_t), "%s\\, %s" % (latex(unit), body_l)


def poly_record(c):
    ft, fl = pretty_factor(c)
    linear = c.degree() <= 0 or all(f.degree() == 1 for f, _ in c.factor())
    return {"text": str(c), "factored": ft, "latex": fl, "degree": int(c.degree()) if c != 0 else -1,
            "linear_factors_only": bool(linear)}


def jdefault(o):
    try:
        return int(o)
    except Exception:
        return str(o)


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=1, ensure_ascii=False, default=jdefault) + "\n")
    return path


def experiment_record(id, title, question, method, scale, result, outcome, reproduce, rows, extra=None):
    """outcome: 'confirmed' (proved by exact computation in the stated range), 'evidence', 'refuted',
    'partial' (the vocabulary of the sister site's exploration log)."""
    rec = {"id": id, "title": title, "question": question, "method": method, "scale": scale,
           "result": result, "outcome": outcome, "reproduce": reproduce,
           "generated_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
           "software": software(), "rows": rows}
    if extra:
        rec.update(extra)
    return write_json(DATA / "experiments" / ("%s.json" % id), rec)


class Timer:
    def __init__(self):
        self.t = time.time()

    def __call__(self):
        return round(time.time() - self.t, 2)
