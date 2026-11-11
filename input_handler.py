"""Task 1 (part 1): sequence input.

Reads a manually entered sequence, a (multi-record) FASTA file, or a plain-text
sequence file.  Whitespace is removed and letters are upper-cased, but invalid
characters are NEVER silently dropped - validation.py reports them.
"""
import os
import re
from dataclasses import dataclass


@dataclass
class SeqRecord:
    id: str
    description: str
    sequence: str


def normalize(raw: str) -> str:
    """Remove all whitespace (spaces, tabs, line breaks) and uppercase."""
    return re.sub(r"\s+", "", raw).upper()


def parse_fasta_text(text: str, default_id: str = "seq") -> list:
    """Parse FASTA (or header-less plain text) into SeqRecord objects.

    * A line starting with '>' starts a new record; the first word is the ID,
      the rest is the description.
    * Sequence lines are concatenated per record.
    * Text before any header (plain-text input) becomes a record named
      ``default_id``.
    * Lines starting with ';' are treated as comments and ignored.
    """
    records = []
    current = None

    def flush():
        if current is not None:
            records.append(SeqRecord(current[0], current[1], normalize("".join(current[2]))))

    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(";"):
            continue
        if stripped.startswith(">"):
            flush()
            parts = stripped[1:].strip().split(None, 1)
            rid = parts[0] if parts else f"{default_id}_{len(records) + 1}"
            desc = parts[1] if len(parts) > 1 else ""
            current = [rid, desc, []]
        else:
            if current is None:
                current = [default_id, "", []]
            current[2].append(stripped)
    flush()
    return records


def load_file(path: str) -> list:
    """Load FASTA / plain-text file.  Raises clear errors for bad paths."""
    path = path.strip().strip('"').strip("'")
    if not path:
        raise ValueError("No file path was given.")
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: '{path}'")
    if not os.path.isfile(path):
        raise ValueError(f"'{path}' is not a file.")
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            text = fh.read()
    except OSError as exc:
        raise ValueError(f"Could not read '{path}': {exc}") from exc
    if not text.strip():
        raise ValueError(f"File '{path}' is empty.")
    stem = os.path.splitext(os.path.basename(path))[0]
    records = parse_fasta_text(text, default_id=stem)
    if not records:
        raise ValueError(f"No sequence records found in '{path}'.")
    return records


def from_manual(text: str, index: int = 1) -> list:
    """Create a record from manually typed text (FASTA syntax also accepted)."""
    if not text.strip():
        raise ValueError("Empty input: no sequence was entered.")
    return parse_fasta_text(text, default_id=f"manual_{index}")
