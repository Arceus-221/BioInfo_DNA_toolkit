"""Task 2: nucleotide composition, derived metrics and sliding windows.

Convention (stated in every report): only the unambiguous bases A, T, G, C
are used for composition metrics; L = A+T+G+C is the denominator.  Ambiguous
symbols are counted separately.  Undefined values are reported as NA (None).
"""
import math

from validation import AMBIGUOUS

CONVENTION = ("Convention: GC/AT content, skews and entropy use ONLY unambiguous A,T,G,C; "
              "L = A+T+G+C is the denominator (not total length). Ambiguous symbols are "
              "counted separately; undefined values are reported as NA.")


def fmt(value, digits=4):
    """Format a number, or 'NA' for None."""
    return "NA" if value is None else f"{value:.{digits}f}"


def shannon_entropy(counts: dict, total: int):
    """H = -sum(p*log2 p) over bases with p>0; None when total == 0."""
    if total == 0:
        return None
    h = 0.0
    for n in counts.values():
        if n > 0:
            p = n / total
            h -= p * math.log2(p)
    return h + 0.0


def compute_stats(seq: str) -> dict:
    counts = {"A": 0, "T": 0, "G": 0, "C": 0}
    ambiguous = invalid = 0
    for ch in seq:
        if ch in counts:
            counts[ch] += 1
        elif ch in AMBIGUOUS:
            ambiguous += 1
        else:
            invalid += 1
    a, t, g, c = counts["A"], counts["T"], counts["G"], counts["C"]
    L = a + t + g + c
    total = len(seq)

    def pct(x):
        return None if L == 0 else x / L * 100

    return {
        "length": total, "A": a, "T": t, "G": g, "C": c, "L": L,
        "A_pct": pct(a), "T_pct": pct(t), "G_pct": pct(g), "C_pct": pct(c),
        "gc_content": pct(g + c),
        "at_content": pct(a + t),
        "gc_skew": None if (g + c) == 0 else (g - c) / (g + c),
        "at_skew": None if (a + t) == 0 else (a - t) / (a + t),
        "entropy": shannon_entropy(counts, L),
        "ambiguous_count": ambiguous,
        "ambiguous_pct": None if total == 0 else ambiguous / total * 100,
        "invalid_count": invalid,
    }


def sliding_windows(seq: str, window: int, step: int, include_partial: bool = False) -> list:
    """Window metrics with 1-based inclusive coordinates.

    Only complete windows are used by default.  With include_partial=True a
    final shorter window covering the tail is appended (documented option).
    """
    if window < 1 or step < 1:
        raise ValueError("Window size and step size must be positive integers.")
    n = len(seq)
    if n == 0:
        raise ValueError("Cannot run sliding-window analysis on an empty sequence.")
    if window > n:
        raise ValueError(f"Window size ({window}) is larger than the sequence length ({n}).")
    rows = []
    start = 0
    last_end = 0
    while start + window <= n:
        rows.append(_window_row(seq, len(rows) + 1, start, start + window))
        last_end = start + window
        start += step
    if include_partial and last_end < n:
        rows.append(_window_row(seq, len(rows) + 1, start, n))
    return rows


def _window_row(seq, idx, start0, end0):
    s = compute_stats(seq[start0:end0])
    return {
        "window": idx, "start": start0 + 1, "end": end0,
        "gc_content": s["gc_content"], "at_content": s["at_content"],
        "gc_skew": s["gc_skew"], "at_skew": s["at_skew"],
        "entropy": s["entropy"], "L": s["L"], "ambiguous_count": s["ambiguous_count"],
    }


def format_summary(records_stats: list) -> str:
    """records_stats: list of (id, stats dict)."""
    lines = ["TASK 2 - NUCLEOTIDE COMPOSITION SUMMARY", "=" * 60, CONVENTION, ""]
    for rid, s in records_stats:
        lines += [
            f"Sequence ID: {rid}",
            f"  Total length          : {s['length']}",
            f"  Counts A/T/G/C        : {s['A']}/{s['T']}/{s['G']}/{s['C']}   (L = {s['L']})",
            f"  Percentages A/T/G/C   : {fmt(s['A_pct'], 2)}/{fmt(s['T_pct'], 2)}/"
            f"{fmt(s['G_pct'], 2)}/{fmt(s['C_pct'], 2)} %",
            f"  GC content            : {fmt(s['gc_content'], 2)} %",
            f"  AT content            : {fmt(s['at_content'], 2)} %",
            f"  GC skew               : {fmt(s['gc_skew'])}",
            f"  AT skew               : {fmt(s['at_skew'])}",
            f"  Shannon entropy       : {fmt(s['entropy'])} bits",
            f"  Ambiguous symbols     : {s['ambiguous_count']} ({fmt(s['ambiguous_pct'], 2)} % of total length)",
            "-" * 60,
        ]
    return "\n".join(lines)


def summary_rows(records_stats: list) -> list:
    rows = []
    for rid, s in records_stats:
        row = {"sequence_id": rid}
        row.update({k: ("NA" if v is None else (round(v, 6) if isinstance(v, float) else v))
                    for k, v in s.items()})
        rows.append(row)
    return rows


def window_csv_rows(rid: str, windows: list) -> list:
    out = []
    for w in windows:
        row = {"sequence_id": rid}
        row.update({k: ("NA" if v is None else (round(v, 6) if isinstance(v, float) else v))
                    for k, v in w.items()})
        out.append(row)
    return out
