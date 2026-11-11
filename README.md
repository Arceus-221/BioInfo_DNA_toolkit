# DNA Sequence Analysis and Computational Genomics Toolkit (Assignment 1A)

Menu-driven Python toolkit covering Tasks 1-4: input/validation, composition and
sliding windows, sequence transformations, and motif detection (naive vs KMP).
All core algorithms are implemented from scratch (no Biopython needed).

## Installation
Python 3.10+ and:
```
pip install numpy pandas matplotlib
```
(Only `matplotlib` is actually imported; the others/Biopython are optional.)

## Running
```
cd dna_toolkit
python main.py                       # interactive menu
python main.py --run-all             # complete analysis on data/sample_sequences.fasta
python main.py --input my.fasta --run-all --window 50 --step 25 --motifs GAATTC,ATG --coords 3 12
python main.py --run-all --quick     # shorter benchmark (1k and 10k bases)
python -m unittest discover -s tests -v
```
Menu: 1 Load, 2 Validate, 3 Composition + sliding window, 4 Transformations,
5 Motifs, 6 Complete analysis, 7 Export combined report, 0 Exit.

## Input
* FASTA (one or many records; `>id description` then sequence lines), or plain text
  (header-less; id = file name), or a sequence typed at the prompt.
* Whitespace is removed and letters upper-cased. Invalid characters are kept and reported.
* Alphabet: the 15 IUPAC DNA symbols. Ambiguous symbols are valid but are excluded
  from composition statistics (denominator L = A+T+G+C) and reported separately.
* Records that fail validation are skipped for Tasks 2-4 with a message.
* Coordinates in all reports are 1-based inclusive.

## Output (`output/`)
| Folder | Contents |
|---|---|
| `reports/` | `task1_validation_report.txt`, `task2_composition_report.txt`, `task3_transformation_report.txt`, `task4_motif_report.txt`, `combined_report.txt` |
| `csv/` | validation, summary statistics, sliding windows, motif matches, palindromes, tandem repeats, benchmark |
| `plots/` | per sequence: GC content, GC skew, AT skew, entropy vs position; plus `benchmark_naive_vs_kmp.png` |

## Modules
`input_handler.py` (parse/load), `validation.py`, `composition.py`, `transformation.py`,
`motif_analysis.py`, `plotting.py`, `main.py`; tests in `tests/`; data in `data/`.

## Algorithms and complexity
| Item | Time |
|---|---|
| Validation, composition, entropy | O(n) |
| Sliding windows (size w, step s) | O(n/s * w) |
| Complement / reverse / reverse complement / RNA | O(n) |
| Naive search | O(nm) worst case |
| KMP (LPS build + search) | O(n + m) |
| Palindrome scan (lengths 4-12) | O(n * k) |
| Tandem-repeat scan (unit <= u) | O(n * u * copies) |

Design notes: undefined metrics (L=0, G+C=0, A+T=0) print `NA`; sliding windows use complete
windows only unless a partial final window is requested; the tandem-repeat scan is left-to-right
and reports primitive units only (e.g. `ATATAT` -> `AT` x3, default minimum 3 copies); the
palindrome scan reports every nested even-length palindrome.

## Limitations
Motif search is literal (IUPAC codes in a motif are not expanded). Tandem repeats are limited to
unit length <= 6 by default. Benchmark timings depend on the machine. Assignment parts that
appear only in the lab manual (restriction sites, ORFs, Hamming/edit distance, alignment) are
outside the four tasks of the assignment sheet and are not implemented.
