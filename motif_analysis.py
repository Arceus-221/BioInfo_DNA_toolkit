"""Task 4: pattern matching, palindromes, tandem repeats and benchmarking.

All positions returned are 1-based (start of the match).
"""
import random
import statistics
import time

from transformation import reverse_complement
from validation import IUPAC


def naive_search(text: str, pattern: str, overlapping: bool = True) -> list:
    n, m = len(text), len(pattern)
    result = []
    if m == 0 or m > n:
        return result
    i = 0
    last = n - m
    while i <= last:
        j = 0
        while j < m and text[i + j] == pattern[j]:
            j += 1
        if j == m:
            result.append(i + 1)
            i += 1 if overlapping else m
        else:
            i += 1
    return result


def build_lps(pattern: str) -> list:
    """lps[i] = length of longest proper prefix of pattern[:i+1] that is also a suffix."""
    m = len(pattern)
    lps = [0] * m
    length = 0
    i = 1
    while i < m:
        if pattern[i] == pattern[length]:
            length += 1
            lps[i] = length
            i += 1
        elif length:
            length = lps[length - 1]
        else:
            lps[i] = 0
            i += 1
    return lps


def kmp_search(text: str, pattern: str, overlapping: bool = True) -> list:
    n, m = len(text), len(pattern)
    result = []
    if m == 0 or m > n:
        return result
    lps = build_lps(pattern)
    i = j = 0
    while i < n:
        if text[i] == pattern[j]:
            i += 1
            j += 1
            if j == m:
                result.append(i - m + 1)
                j = lps[j - 1] if overlapping else 0
        elif j:
            j = lps[j - 1]
        else:
            i += 1
    return result


def clean_motif(motif: str) -> str:
    motif = "".join(motif.split()).upper()
    if not motif:
        raise ValueError("Motif is empty.")
    bad = sorted({c for c in motif if c not in IUPAC})
    if bad:
        raise ValueError(f"Motif '{motif}' contains invalid symbol(s): {', '.join(bad)}")
    return motif


def is_rc_palindrome(s: str) -> bool:
    return len(s) > 0 and len(s) % 2 == 0 and s == reverse_complement(s)


def find_palindromes(seq: str, min_len: int = 4, max_len: int = 12) -> list:
    """All substrings (even length) equal to their reverse complement.

    Returns (start_1based, end_1based, substring).  O(n * (max-min) * L).
    """
    found = []
    n = len(seq)
    for length in range(min_len + (min_len % 2), max_len + 1, 2):
        for i in range(0, n - length + 1):
            sub = seq[i:i + length]
            try:
                if sub == reverse_complement(sub):
                    found.append((i + 1, i + length, sub))
            except ValueError:
                continue
    found.sort()
    return found


def _is_primitive(unit: str) -> bool:
    k = len(unit)
    return not any(k % d == 0 and unit[:d] * (k // d) == unit for d in range(1, k))


def find_tandem_repeats(seq: str, min_unit: int = 1, max_unit: int = 6, min_copies: int = 3) -> list:
    """Scan for maximal runs of adjacent copies of a repeat unit.

    Returns dicts (start, end, unit, copies) with 1-based inclusive coordinates.
    The scan is left-to-right; after a run is reported the scan resumes after it
    (runs that begin inside a reported run are not reported again).
    """
    n = len(seq)
    runs = []
    i = 0
    while i < n:
        best = None
        for k in range(min_unit, max_unit + 1):
            if i + k * min_copies > n:
                break
            unit = seq[i:i + k]
            if not _is_primitive(unit):
                continue
            copies, j = 1, i + k
            while seq[j:j + k] == unit:
                copies += 1
                j += k
            if copies >= min_copies and (best is None or copies * k > best[0] * len(best[1])):
                best = (copies, unit)
        if best:
            copies, unit = best
            runs.append({"start": i + 1, "end": i + copies * len(unit), "unit": unit, "copies": copies})
            i += copies * len(unit)
        else:
            i += 1
    return runs


def tandem_runs_for_motif(seq: str, motif: str, min_copies: int = 2) -> list:
    """Runs of the user-given motif repeated back-to-back (uses KMP)."""
    starts = kmp_search(seq, motif, overlapping=True)
    m = len(motif)
    runs, k = [], 0
    while k < len(starts):
        copies = 1
        while k + copies < len(starts) and starts[k + copies] == starts[k] + copies * m:
            copies += 1
        if copies >= min_copies:
            runs.append({"start": starts[k], "end": starts[k] + copies * m - 1, "unit": motif, "copies": copies})
        k += copies
    return runs


def _time_ms(func, text, pattern, reps):
    times = []
    hits = 0
    for _ in range(reps):
        t0 = time.perf_counter()
        hits = len(func(text, pattern, True))
        times.append((time.perf_counter() - t0) * 1000)
    return statistics.median(times), hits


def benchmark(sizes=(1000, 100000, 1000000), reps=3, seed=42, motif_len=8) -> list:
    """Time naive vs KMP (median of `reps`, search only - no I/O or plotting).

    Two scenarios: 'random' DNA and an adversarial 'repeat' case
    (text = A*n, pattern = A*(m-1)+T) which forces naive toward O(n*m).
    """
    rng = random.Random(seed)
    rows = []
    for scenario in ("random", "repeat"):
        for n in sizes:
            if scenario == "random":
                text = "".join(rng.choice("ACGT") for _ in range(n))
                pattern = text[n // 2:n // 2 + motif_len]
            else:
                if n > 100000:
                    continue
                text = "A" * n
                pattern = "A" * 19 + "T"
            naive_ms, naive_hits = _time_ms(naive_search, text, pattern, reps)
            kmp_ms, kmp_hits = _time_ms(kmp_search, text, pattern, reps)
            rows.append({"scenario": scenario, "size": n, "motif_length": len(pattern),
                         "naive_ms": round(naive_ms, 3), "kmp_ms": round(kmp_ms, 3),
                         "matches": kmp_hits, "results_agree": naive_hits == kmp_hits})
    return rows


def format_benchmark(rows: list) -> str:
    lines = ["BENCHMARK - NAIVE vs KMP (median of repeated runs, time.perf_counter, search only)",
             f"{'scenario':<9}{'size':>10}{'motif':>7}{'naive (ms)':>13}{'KMP (ms)':>12}{'matches':>9}{'agree':>7}"]
    for r in rows:
        lines.append(f"{r['scenario']:<9}{r['size']:>10}{r['motif_length']:>7}{r['naive_ms']:>13.3f}"
                     f"{r['kmp_ms']:>12.3f}{r['matches']:>9}{str(r['results_agree']):>7}")
    lines += ["",
              "Discussion: theory says naive is O(nm) and KMP is O(n+m). On random DNA mismatches occur",
              "after ~1.3 comparisons, so naive is effectively O(n) and the two are close; KMP's extra",
              "bookkeeping can even make it slower in pure Python. On the repetitive adversarial input",
              "naive pays ~m comparisons per position while KMP stays linear, so the gap widens with m.",
              "Other effects: interpreter overhead, caching, timer resolution and machine load."]
    return "\n".join(lines)


def analyse_record(rid: str, seq: str, motifs: list, overlapping: bool,
                   pal_min=4, pal_max=12, tr_max_unit=6, tr_min_copies=3):
    """Returns (text, match_rows)."""
    mode = "overlapping" if overlapping else "non-overlapping"
    lines = [f"Sequence ID: {rid}   (length {len(seq)}; mode: {mode}; 1-based start positions)"]
    rows = []
    for motif in motifs:
        nv = naive_search(seq, motif, overlapping)
        kp = kmp_search(seq, motif, overlapping)
        agree = "agree" if nv == kp else "MISMATCH"
        shown = ", ".join(map(str, kp[:40])) + (" ..." if len(kp) > 40 else "")
        pal = " [reverse-complement palindrome]" if is_rc_palindrome(motif) else ""
        lines.append(f"  Motif {motif}{pal}: {len(kp)} match(es); naive vs KMP {agree}")
        if kp:
            lines.append(f"    positions: {shown}")
        for p in kp:
            rows.append({"sequence_id": rid, "motif": motif, "start": p, "end": p + len(motif) - 1,
                         "mode": mode})
        runs = tandem_runs_for_motif(seq, motif, 2)
        for r in runs:
            lines.append(f"    motif tandem run: {motif} x{r['copies']} at {r['start']}-{r['end']}")
    pals = find_palindromes(seq, pal_min, pal_max)
    lines.append(f"  Reverse-complement palindromes (even length {pal_min}-{pal_max}): {len(pals)}")
    for s, e, sub in pals[:25]:
        lines.append(f"    {s}-{e}: {sub}")
    if len(pals) > 25:
        lines.append(f"    ... {len(pals) - 25} more (see CSV)")
    trs = find_tandem_repeats(seq, 1, tr_max_unit, tr_min_copies)
    lines.append(f"  Tandem repeats (unit 1-{tr_max_unit}, >= {tr_min_copies} copies): {len(trs)}")
    for r in trs[:25]:
        lines.append(f"    {r['start']}-{r['end']}: unit {r['unit']} x{r['copies']}")
    if len(trs) > 25:
        lines.append(f"    ... {len(trs) - 25} more (see CSV)")
    lines.append("-" * 60)
    return "\n".join(lines), rows, pals, trs
