"""Task 3: DNA sequence transformations (all implemented from scratch)."""

COMPLEMENT = {
    "A": "T", "T": "A", "G": "C", "C": "G",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N",
}


def complement(seq: str) -> str:
    try:
        return "".join(COMPLEMENT[ch] for ch in seq)
    except KeyError as exc:
        raise ValueError(f"Cannot complement unsupported symbol {exc.args[0]!r}.") from None


def reverse(seq: str) -> str:
    return seq[::-1]


def reverse_complement(seq: str) -> str:
    return complement(reverse(seq))


def dna_to_rna(seq: str) -> str:
    """Coding-strand representation: T -> U."""
    return seq.replace("T", "U")


def rna_to_dna(seq: str) -> str:
    return seq.replace("U", "T")


def subsequence(seq: str, start: int, end: int) -> str:
    """1-based inclusive [start, end] -> Python slice seq[start-1:end]."""
    if start < 1:
        raise ValueError(f"Invalid start {start}: must be >= 1.")
    if end > len(seq):
        raise ValueError(f"Invalid end {end}: sequence length is {len(seq)}.")
    if start > end:
        raise ValueError(f"Invalid range: start ({start}) is greater than end ({end}).")
    return seq[start - 1:end]


def transformation_report(rid: str, seq: str, start=None, end=None) -> str:
    comp = complement(seq)
    rev = reverse(seq)
    rc = reverse_complement(seq)
    rna = dna_to_rna(seq)
    back = rna_to_dna(rna)
    lines = [f"Sequence ID: {rid}",
             f"  Original (5'->3')   : {seq}",
             f"  Complement          : {comp}",
             f"  Reverse             : {rev}",
             f"  Reverse complement  : {rc}",
             f"  RNA (T->U)          : {rna}",
             f"  DNA recovered (U->T): {back}"]
    if start is not None and end is not None:
        try:
            lines.append(f"  Subsequence {start}-{end} (1-based inclusive): {subsequence(seq, start, end)}")
        except ValueError as exc:
            lines.append(f"  Subsequence {start}-{end}: ERROR - {exc}")
    checks = [
        ("reverse-complement twice returns original", reverse_complement(rc) == seq),
        ("complement preserves length", len(comp) == len(seq)),
        ("reverse-complement preserves length", len(rc) == len(seq)),
        ("RNA -> DNA round trip", back == seq),
    ]
    lines.append("  Self-checks: " + "; ".join(f"{n}: {'PASS' if ok else 'FAIL'}" for n, ok in checks))
    lines.append("-" * 60)
    return "\n".join(lines)
