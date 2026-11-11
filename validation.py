"""Task 1 (part 2): validation against the IUPAC DNA alphabet."""

UNAMBIGUOUS = frozenset("ATGC")
AMBIGUOUS = frozenset("RYSWKMBDHVN")
IUPAC = UNAMBIGUOUS | AMBIGUOUS

IUPAC_MEANING = {
    "R": "A/G", "Y": "C/T", "S": "G/C", "W": "A/T", "K": "G/T", "M": "A/C",
    "B": "C/G/T", "D": "A/G/T", "H": "A/C/T", "V": "A/C/G", "N": "A/C/G/T",
}


def validate_sequence(seq: str) -> dict:
    """Classify every symbol as unambiguous, ambiguous or invalid.

    Positions are 1-based.  Time complexity O(n).
    """
    ambiguous_counts = {}
    ambiguous_positions = []
    invalid = []
    for pos, ch in enumerate(seq, start=1):
        if ch in UNAMBIGUOUS:
            continue
        if ch in AMBIGUOUS:
            ambiguous_counts[ch] = ambiguous_counts.get(ch, 0) + 1
            ambiguous_positions.append((ch, pos))
        else:
            invalid.append((ch, pos))

    if len(seq) == 0:
        status = "Invalid (empty sequence)"
    elif invalid:
        status = "Invalid"
    elif ambiguous_counts:
        status = "Valid with ambiguity"
    else:
        status = "Valid"
    return {
        "length": len(seq),
        "status": status,
        "ambiguous_count": sum(ambiguous_counts.values()),
        "ambiguous_counts": ambiguous_counts,
        "ambiguous_positions": ambiguous_positions,
        "invalid": invalid,
    }


def validate_record(rec) -> dict:
    result = validate_sequence(rec.sequence)
    result["id"] = rec.id
    result["description"] = rec.description
    return result


def is_usable(result: dict) -> bool:
    """True when the record may be used for downstream analysis."""
    return not result["status"].startswith("Invalid")


def format_validation_report(results: list) -> str:
    lines = ["TASK 1 - SEQUENCE VALIDATION REPORT", "=" * 60,
             "Alphabet: 15 IUPAC DNA symbols (A,T,G,C + R,Y,S,W,K,M,B,D,H,V,N).",
             "Ambiguous symbols are VALID; characters outside the alphabet are INVALID.",
             "Positions are 1-based.", ""]
    for r in results:
        lines.append(f"Sequence ID : {r['id']}")
        if r["description"]:
            lines.append(f"Description : {r['description']}")
        lines.append(f"Length      : {r['length']}")
        lines.append(f"Status      : {r['status']}")
        if r["ambiguous_count"]:
            kinds = ", ".join(f"{c}x{n}" for c, n in sorted(r["ambiguous_counts"].items()))
            where = ", ".join(f"{c}@{p}" for c, p in r["ambiguous_positions"][:30])
            more = " ..." if len(r["ambiguous_positions"]) > 30 else ""
            lines.append(f"Ambiguous   : {r['ambiguous_count']} symbol(s) [{kinds}] at {where}{more}")
        else:
            lines.append("Ambiguous   : 0")
        if r["invalid"]:
            bad = ", ".join(f"'{c}' at position {p}" for c, p in r["invalid"][:30])
            more = " ..." if len(r["invalid"]) > 30 else ""
            lines.append(f"Invalid     : {len(r['invalid'])} character(s): {bad}{more}")
        else:
            lines.append("Invalid     : none")
        lines.append("-" * 60)
    return "\n".join(lines)


def validation_rows(results: list) -> list:
    """Rows for CSV export."""
    rows = []
    for r in results:
        rows.append({
            "sequence_id": r["id"],
            "length": r["length"],
            "status": r["status"],
            "ambiguous_count": r["ambiguous_count"],
            "ambiguous_symbols": ";".join(f"{c}:{n}" for c, n in sorted(r["ambiguous_counts"].items())),
            "invalid_count": len(r["invalid"]),
            "invalid_characters": ";".join(f"{c}@{p}" for c, p in r["invalid"]),
        })
    return rows
