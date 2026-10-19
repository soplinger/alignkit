#!/usr/bin/env python3
"""alignkit: offline pairwise sequence alignment from FASTA files.

Global, local and semi-global alignment with affine gaps (Gotoh), exact
Levenshtein distance with CIGAR strings, reverse-complement search, pairwise
identity matrices, FASTA linting with line-numbered errors and CI identity
gates. Standard library only; no network access.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

NEG = -(10 ** 12)
MODES = ("global", "local", "semiglobal")
DEFAULT_MAX_CELLS = 2_000_000
DNA_BASES = set("ACGTN")
IUPAC_CODES = set("RYKMSWBDHV")
_COMPLEMENT = {
    "A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
    "R": "Y", "Y": "R", "K": "M", "M": "K", "B": "V", "V": "B",
    "D": "H", "H": "D", "S": "S", "W": "W",
}


class AlignError(Exception):
    """A problem with the request itself: bad input, bad options or too large."""


@dataclass(frozen=True)
class Scoring:
    match: int = 2
    mismatch: int = -1
    gap_open: int = -4
    gap_extend: int = -1

    def check(self) -> None:
        if self.gap_open > 0 or self.gap_extend > 0:
            raise AlignError("gap penalties must be zero or negative")
        if self.match <= self.mismatch:
            raise AlignError("match score must be higher than the mismatch score")

    def score(self, x: str, y: str) -> int:
        return self.match if x == y else self.mismatch


@dataclass(frozen=True)
class Record:
    name: str
    seq: str
    line: int


@dataclass(frozen=True)
class Issue:
    line: int
    severity: str
    message: str


@dataclass
class Alignment:
    mode: str
    score: int
    aligned_a: str
    aligned_b: str
    a_start: int
    a_end: int
    b_start: int
    b_end: int
    strand: str = "+"

    def stats(self) -> dict:
        aa, ab = self.aligned_a, self.aligned_b
        columns = len(aa)
        matches = sum(1 for x, y in zip(aa, ab) if x != "-" and x == y)
        mismatches = sum(1 for x, y in zip(aa, ab) if x != "-" and y != "-" and x != y)
        gap_columns = sum(1 for x, y in zip(aa, ab) if x == "-" or y == "-")
        return {
            "columns": columns,
            "matches": matches,
            "mismatches": mismatches,
            "gap_columns": gap_columns,
            "gap_opens": _gap_runs(aa) + _gap_runs(ab),
            "identity": Fraction(matches, columns) if columns else Fraction(0),
        }


def _gap_runs(row: str) -> int:
    return sum(1 for k, ch in enumerate(row) if ch == "-" and (k == 0 or row[k - 1] != "-"))


def parse_fasta(text: str, alphabet: str = "any") -> tuple[list[Record], list[Issue]]:
    """Parse FASTA text. Returns (valid records, issues) with 1-based line numbers."""
    if alphabet not in ("any", "dna"):
        raise AlignError(f"unknown alphabet {alphabet!r}")
    issues: list[Issue] = []
    blocks: list[dict] = []
    current = None
    for line_no, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            current = {"name": line[1:].strip(), "line": line_no, "lines": []}
            blocks.append(current)
        elif current is None:
            issues.append(Issue(line_no, "error", "sequence data before the first '>' header"))
        else:
            current["lines"].append((line_no, line))

    records: list[Record] = []
    seen: dict[str, int] = {}
    for block in blocks:
        name, header_line = block["name"], block["line"]
        if not name:
            issues.append(Issue(header_line, "error", "header has no record name"))
            continue
        if name in seen:
            issues.append(Issue(header_line, "error",
                                f"duplicate record name {name!r} (first at line {seen[name]})"))
            continue
        seen[name] = header_line
        residues: list[str] = []
        lowercase = False
        for line_no, text_line in block["lines"]:
            for col, ch in enumerate(text_line, start=1):
                if ch.islower():
                    lowercase = True
                up = ch.upper()
                if not (ch.isascii() and ch.isalpha()):
                    issues.append(Issue(line_no, "error", f"invalid character {ch!r} at column {col}"))
                    continue
                if alphabet == "dna" and up not in DNA_BASES:
                    if up in IUPAC_CODES:
                        issues.append(Issue(line_no, "warning",
                                            f"ambiguity code {ch!r} at column {col}"))
                    else:
                        issues.append(Issue(line_no, "error",
                                            f"{ch!r} at column {col} is not a DNA base"))
                        continue
                residues.append(up)
        if lowercase:
            issues.append(Issue(header_line, "warning",
                                f"record {name!r} has lowercase residues (normalised to upper case)"))
        if not residues:
            issues.append(Issue(header_line, "error", f"record {name!r} has no sequence"))
            continue
        records.append(Record(name, "".join(residues), header_line))
    issues.sort(key=lambda issue: issue.line)
    return records, issues


def reverse_complement(seq: str) -> str:
    try:
        return "".join(_COMPLEMENT[c] for c in reversed(seq.upper()))
    except KeyError as exc:
        raise AlignError(f"cannot reverse-complement non-DNA symbol {exc.args[0]!r}") from None


def check_cells(n: int, m: int, max_cells: int) -> None:
    if (n + 1) * (m + 1) > max_cells:
        raise AlignError(f"alignment needs {(n + 1) * (m + 1)} cells, above the limit of "
                         f"{max_cells} (raise --max-cells)")


def _choose(value: int, options: list[tuple[str, int]]) -> str:
    for state, candidate in options:
        if candidate == value:
            return state
    raise RuntimeError("traceback found no predecessor")  # pragma: no cover


def _trace(a, b, M, X, Y, i, j, state, scoring, local):
    go, ge = scoring.gap_open, scoring.gap_extend
    out_a: list[str] = []
    out_b: list[str] = []
    while True:
        if local and state == "M" and M[i][j] <= 0:
            break
        if i == 0 and j == 0:
            break
        if i == 0 or j == 0:
            while i > 0:
                out_a.append(a[i - 1])
                out_b.append("-")
                i -= 1
            while j > 0:
                out_a.append("-")
                out_b.append(b[j - 1])
                j -= 1
            break
        if state == "M":
            out_a.append(a[i - 1])
            out_b.append(b[j - 1])
            prev = M[i][j] - scoring.score(a[i - 1], b[j - 1])
            i, j = i - 1, j - 1
            state = _choose(prev, [("M", M[i][j]), ("X", X[i][j]), ("Y", Y[i][j])])
        elif state == "X":
            out_a.append(a[i - 1])
            out_b.append("-")
            state = _choose(X[i][j], [("M", M[i - 1][j] + go), ("X", X[i - 1][j] + ge),
                                      ("Y", Y[i - 1][j] + go)])
            i -= 1
        else:
            out_a.append("-")
            out_b.append(b[j - 1])
            state = _choose(Y[i][j], [("M", M[i][j - 1] + go), ("Y", Y[i][j - 1] + ge),
                                      ("X", X[i][j - 1] + go)])
            j -= 1
    return "".join(reversed(out_a)), "".join(reversed(out_b)), i, j


def align(a: str, b: str, mode: str = "global", scoring: Scoring | None = None,
          max_cells: int = DEFAULT_MAX_CELLS) -> Alignment:
    """Align two sequences with affine gaps (Gotoh).

    global: end-to-end. local: Smith-Waterman best segment (score 0 if none).
    semiglobal: end gaps are free, so overhangs cost nothing.
    """
    if mode not in MODES:
        raise AlignError(f"unknown mode {mode!r}")
    scoring = scoring or Scoring()
    scoring.check()
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        raise AlignError("cannot align an empty sequence")
    check_cells(n, m, max_cells)
    go, ge = scoring.gap_open, scoring.gap_extend
    local = mode == "local"
    free = mode == "semiglobal"

    M = [[NEG] * (m + 1) for _ in range(n + 1)]
    X = [[NEG] * (m + 1) for _ in range(n + 1)]
    Y = [[NEG] * (m + 1) for _ in range(n + 1)]
    if local:
        for i in range(n + 1):
            M[i][0] = 0
        for j in range(m + 1):
            M[0][j] = 0
    else:
        M[0][0] = 0
        for i in range(1, n + 1):
            X[i][0] = 0 if free else go + (i - 1) * ge
        for j in range(1, m + 1):
            Y[0][j] = 0 if free else go + (j - 1) * ge

    for i in range(1, n + 1):
        ai = a[i - 1]
        for j in range(1, m + 1):
            diag = max(M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1])
            step = diag + scoring.score(ai, b[j - 1])
            M[i][j] = max(0, step) if local else step
            X[i][j] = max(M[i - 1][j] + go, X[i - 1][j] + ge, Y[i - 1][j] + go)
            Y[i][j] = max(M[i][j - 1] + go, Y[i][j - 1] + ge, X[i][j - 1] + go)

    if local:
        best, ei, ej = 0, 0, 0
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                if M[i][j] > best:
                    best, ei, ej = M[i][j], i, j
        if best <= 0:
            return Alignment(mode, 0, "", "", 0, 0, 0, 0)
        score, state = best, "M"
    elif free:
        cands = [(M[n][m], "M", n, m)]
        cands += [(X[n][j], "X", n, j) for j in range(m + 1)]
        cands += [(Y[i][m], "Y", i, m) for i in range(n + 1)]
        score, state, ei, ej = max(cands, key=lambda c: c[0])
    else:
        cands = [(M[n][m], "M"), (X[n][m], "X"), (Y[n][m], "Y")]
        score, state = max(cands, key=lambda c: c[0])
        ei, ej = n, m

    aa, ab, si, sj = _trace(a, b, M, X, Y, ei, ej, state, scoring, local)
    a_end, b_end = ei, ej
    if free:
        tail_a, tail_b = a[ei:], b[ej:]
        aa += tail_a + "-" * len(tail_b)
        ab += "-" * len(tail_a) + tail_b
        a_end, b_end = n, m
    return Alignment(mode, score, aa, ab, si, a_end, sj, b_end)


def align_strands(a: str, b: str, mode: str, scoring: Scoring | None, both: bool,
                  max_cells: int = DEFAULT_MAX_CELLS) -> Alignment:
    """Align forward, and also against b's reverse complement when both is set."""
    forward = align(a, b, mode, scoring, max_cells)
    if not both:
        return forward
    reverse = align(a, reverse_complement(b), mode, scoring, max_cells)
    if reverse.score > forward.score:
        reverse.strand = "-"
        return reverse
    return forward


def levenshtein(a: str, b: str, max_cells: int = DEFAULT_MAX_CELLS) -> tuple[int, str, str]:
    """Unit-cost edit distance with one optimal alignment (prefers diagonal moves)."""
    n, m = len(a), len(b)
    check_cells(n, m, max_cells)
    D = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        D[i][0] = i
    for j in range(m + 1):
        D[0][j] = j
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            cost = 0 if a[i - 1] == b[j - 1] else 1
            D[i][j] = min(D[i - 1][j - 1] + cost, D[i - 1][j] + 1, D[i][j - 1] + 1)
    out_a: list[str] = []
    out_b: list[str] = []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and D[i][j] == D[i - 1][j - 1] + (0 if a[i - 1] == b[j - 1] else 1):
            out_a.append(a[i - 1])
            out_b.append(b[j - 1])
            i, j = i - 1, j - 1
        elif i > 0 and D[i][j] == D[i - 1][j] + 1:
            out_a.append(a[i - 1])
            out_b.append("-")
            i -= 1
        else:
            out_a.append("-")
            out_b.append(b[j - 1])
            j -= 1
    return D[n][m], "".join(reversed(out_a)), "".join(reversed(out_b))


def cigar(aligned_a: str, aligned_b: str) -> str:
    """Run-length ops: = match, X mismatch, I base only in b, D base only in a."""
    runs: list[list] = []
    for x, y in zip(aligned_a, aligned_b):
        if x == "-":
            op = "I"
        elif y == "-":
            op = "D"
        elif x == y:
            op = "="
        else:
            op = "X"
        if runs and runs[-1][1] == op:
            runs[-1][0] += 1
        else:
            runs.append([1, op])
    return "".join(f"{count}{op}" for count, op in runs)


def pairwise(records: list[Record], mode: str = "global", scoring: Scoring | None = None,
             max_cells: int = DEFAULT_MAX_CELLS) -> list[tuple[Record, Record, Alignment]]:
    rows = []
    for i, ra in enumerate(records):
        for rb in records[i + 1:]:
            rows.append((ra, rb, align(ra.seq, rb.seq, mode, scoring, max_cells)))
    return rows


def render_blocks(al: Alignment, a_name: str, b_name: str, width: int = 60) -> str:
    if not al.aligned_a:
        return "(no aligned residues)"
    label = max(len(a_name), len(b_name))
    pad = label + 9
    a_pos, b_pos = al.a_start, al.b_start
    blocks = []
    for k in range(0, len(al.aligned_a), width):
        ca = al.aligned_a[k:k + width]
        cb = al.aligned_b[k:k + width]
        na = len(ca) - ca.count("-")
        nb = len(cb) - cb.count("-")
        mid = "".join(_bar(x, y) for x, y in zip(ca, cb))
        sa = a_pos + 1 if na else a_pos
        sb = b_pos + 1 if nb else b_pos
        blocks.append("\n".join([
            f"{a_name:<{label}} {sa:>7} {ca} {a_pos + na:<7}".rstrip(),
            f"{'':<{pad}}{mid}".rstrip(),
            f"{b_name:<{label}} {sb:>7} {cb} {b_pos + nb:<7}".rstrip(),
        ]))
        a_pos += na
        b_pos += nb
    return "\n\n".join(blocks)


def _bar(x: str, y: str) -> str:
    if x == "-" or y == "-":
        return " "
    return "|" if x == y else "."


def _span(start: int, end: int) -> str:
    return f"{start + 1}-{end}" if end > start else "none"


def render_report(al: Alignment, a_name: str, b_name: str, scoring: Scoring, width: int) -> str:
    st = al.stats()
    b_label = b_name + (" (reverse complement)" if al.strand == "-" else "")
    head = [
        f"alignkit: {al.mode} alignment, strand {al.strand}",
        f"  a: {a_name}  residues {_span(al.a_start, al.a_end)}",
        f"  b: {b_label}  residues {_span(al.b_start, al.b_end)}",
        f"  scoring: match {scoring.match}, mismatch {scoring.mismatch}, "
        f"gap open {scoring.gap_open}, gap extend {scoring.gap_extend}",
        f"  score {al.score} | columns {st['columns']} | matches {st['matches']} | "
        f"mismatches {st['mismatches']} | gap columns {st['gap_columns']} | "
        f"gap opens {st['gap_opens']}",
        f"  identity {float(st['identity']) * 100:.1f}% ({st['matches']}/{st['columns']}) | "
        f"cigar {cigar(al.aligned_a, al.aligned_b) or '-'}",
    ]
    return "\n".join(head) + "\n\n" + render_blocks(al, a_name, b_label, width)


def alignment_dict(al: Alignment, a_name: str, b_name: str) -> dict:
    st = al.stats()
    return {
        "a": a_name,
        "b": b_name,
        "mode": al.mode,
        "strand": al.strand,
        "score": al.score,
        "a_range": [al.a_start, al.a_end],
        "b_range": [al.b_start, al.b_end],
        "columns": st["columns"],
        "matches": st["matches"],
        "mismatches": st["mismatches"],
        "gap_columns": st["gap_columns"],
        "gap_opens": st["gap_opens"],
        "identity": round(float(st["identity"]), 4),
        "cigar": cigar(al.aligned_a, al.aligned_b),
        "aligned_a": al.aligned_a,
        "aligned_b": al.aligned_b,
    }


# ---------------------------------------------------------------- CLI

def _load_records(path: str, alphabet: str = "any") -> list[Record]:
    text = Path(path).read_text(encoding="utf-8")
    records, issues = parse_fasta(text, alphabet)
    errors = [issue for issue in issues if issue.severity == "error"]
    if errors:
        for issue in errors:
            print(f"{path}:{issue.line}: error: {issue.message}", file=sys.stderr)
        raise AlignError(f"{path} has {len(errors)} error(s); run 'alignkit lint {path}'")
    if not records:
        raise AlignError(f"{path} has no records")
    return records


def _first_record(path: str) -> Record:
    return _load_records(path)[0]


def _scoring(args) -> Scoring:
    return Scoring(args.match, args.mismatch, args.gap_open, args.gap_extend)


def cmd_align(args) -> int:
    scoring = _scoring(args)
    ra, rb = _first_record(args.a), _first_record(args.b)
    al = align_strands(ra.seq, rb.seq, args.mode, scoring, args.strand == "both", args.max_cells)
    if args.json:
        print(json.dumps(alignment_dict(al, ra.name, rb.name), indent=2))
    else:
        print(render_report(al, ra.name, rb.name, scoring, args.width))
    if args.min_identity is not None and al.stats()["identity"] < Fraction(args.min_identity):
        print(f"alignkit: identity {float(al.stats()['identity']) * 100:.1f}% is below "
              f"--min-identity {args.min_identity}", file=sys.stderr)
        return 1
    return 0


def cmd_pairs(args) -> int:
    records = _load_records(args.file)
    if len(records) < 2:
        raise AlignError(f"{args.file} needs at least two records for pairs")
    scoring = _scoring(args)
    rows = pairwise(records, args.mode, scoring, args.max_cells)
    gate = Fraction(args.min_identity) if args.min_identity is not None else None
    failing = []
    labels = [f"{ra.name} vs {rb.name}" for ra, rb, _ in rows]
    width = max(len(label) for label in labels)
    if args.json:
        out = []
        for label, (ra, rb, al) in zip(labels, rows):
            st = al.stats()
            out.append({"a": ra.name, "b": rb.name, "score": al.score,
                        "columns": st["columns"], "identity": round(float(st["identity"]), 4)})
        print(json.dumps(out, indent=2))
    else:
        print(f"{'pair':<{width}}  identity  score  columns")
        for label, (ra, rb, al) in zip(labels, rows):
            st = al.stats()
            print(f"{label:<{width}}  {float(st['identity']) * 100:7.1f}%  "
                  f"{al.score:>5}  {st['columns']:>7}")
    for label, (_, _, al) in zip(labels, rows):
        if gate is not None and al.stats()["identity"] < gate:
            failing.append(label)
    for label in failing:
        print(f"alignkit: FAIL {label}: identity below {args.min_identity}", file=sys.stderr)
    return 1 if failing else 0


def cmd_distance(args) -> int:
    ra, rb = _first_record(args.a), _first_record(args.b)
    dist, aa, ab = levenshtein(ra.seq, rb.seq, args.max_cells)
    ops = cigar(aa, ab)
    hamming = (sum(x != y for x, y in zip(ra.seq, rb.seq))
               if len(ra.seq) == len(rb.seq) else None)
    if args.json:
        print(json.dumps({"a": ra.name, "b": rb.name, "edit_distance": dist,
                          "hamming": hamming, "cigar": ops,
                          "aligned_a": aa, "aligned_b": ab}, indent=2))
        return 0
    print(f"edit distance {dist} | hamming {hamming if hamming is not None else 'n/a'} | "
          f"cigar {ops or '-'}")
    print()
    print(render_blocks(Alignment("edit", -dist, aa, ab, 0, len(ra.seq), 0, len(rb.seq)),
                        ra.name, rb.name, args.width))
    return 0


def cmd_lint(args) -> int:
    text = Path(args.file).read_text(encoding="utf-8")
    records, issues = parse_fasta(text, args.alphabet)
    for issue in issues:
        print(f"{args.file}:{issue.line}: {issue.severity}: {issue.message}")
    errors = sum(1 for issue in issues if issue.severity == "error")
    warnings = len(issues) - errors
    print(f"{args.file}: {errors} error(s), {warnings} warning(s), "
          f"{len(records)} valid record(s)")
    return 1 if errors or (args.strict and warnings) else 0


def cmd_revcomp(args) -> int:
    lines = []
    for record in _load_records(args.file):
        rc = reverse_complement(record.seq)
        lines.append(f">{record.name} reverse complement")
        lines.extend(rc[k:k + 60] for k in range(0, len(rc), 60))
    print("\n".join(lines))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="alignkit",
                                     description="Offline pairwise sequence alignment.")
    sub = parser.add_subparsers(dest="command", required=True)

    def scoring_args(p):
        g = p.add_argument_group("scoring")
        g.add_argument("--mode", choices=MODES, default="global")
        g.add_argument("--match", type=int, default=2)
        g.add_argument("--mismatch", type=int, default=-1)
        g.add_argument("--gap-open", type=int, default=-4)
        g.add_argument("--gap-extend", type=int, default=-1)
        g.add_argument("--max-cells", type=int, default=DEFAULT_MAX_CELLS)

    p = sub.add_parser("align", help="align the first record of two FASTA files")
    p.add_argument("a")
    p.add_argument("b")
    scoring_args(p)
    p.add_argument("--strand", choices=("forward", "both"), default="forward",
                   help="'both' also tries b's reverse complement (DNA)")
    p.add_argument("--min-identity", help="exit 1 if identity is below this fraction, e.g. 0.9")
    p.add_argument("--width", type=int, default=60)
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_align)

    p = sub.add_parser("pairs", help="pairwise identity table for every record in a file")
    p.add_argument("file")
    scoring_args(p)
    p.add_argument("--min-identity", help="exit 1 if any pair is below this fraction")
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_pairs)

    p = sub.add_parser("distance", help="Levenshtein edit distance, Hamming and CIGAR")
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("--max-cells", type=int, default=DEFAULT_MAX_CELLS)
    p.add_argument("--width", type=int, default=60)
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_distance)

    p = sub.add_parser("lint", help="check a FASTA file with line-numbered findings")
    p.add_argument("file")
    p.add_argument("--alphabet", choices=("any", "dna"), default="any")
    p.add_argument("--strict", action="store_true", help="warnings also fail")
    p.set_defaults(func=cmd_lint)

    p = sub.add_parser("revcomp", help="print the reverse complement of every DNA record")
    p.add_argument("file")
    p.set_defaults(func=cmd_revcomp)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (AlignError, OSError, ValueError) as exc:
        print(f"alignkit: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
