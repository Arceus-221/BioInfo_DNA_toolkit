"""DNA Sequence Analysis Toolkit - menu-driven application (Assignment 1A).

Interactive:      python main.py
Non-interactive:  python main.py --input data/sample_sequences.fasta --run-all
"""
import argparse
import csv
import os
import sys
from datetime import datetime

import composition as comp
import input_handler as ih
import motif_analysis as ma
import plotting
import transformation as tr
import validation as val

BASE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_MOTIFS = ["GAATTC", "ATG", "TATA", "GGATCC"]

MENU = """
==================== DNA SEQUENCE ANALYSIS TOOLKIT ====================
 1. Load sequence(s)
 2. Validate sequence(s)
 3. Nucleotide composition and sliding-window analysis
 4. Sequence transformations
 5. Motif and pattern detection
 6. Run complete analysis
 7. Export reports
 0. Exit
======================================================================="""


class Session:
    def __init__(self, outdir):
        self.records = []
        self.validations = {}
        self.outdir = outdir
        self.sections = {}
        for sub in ("reports", "csv", "plots"):
            os.makedirs(os.path.join(outdir, sub), exist_ok=True)

    def path(self, sub, name):
        return os.path.join(self.outdir, sub, name)

    def write_text(self, name, text):
        p = self.path("reports", name)
        with open(p, "w", encoding="utf-8") as fh:
            fh.write(text + "\n")
        return p

    def write_csv(self, name, rows, fields=None):
        p = self.path("csv", name)
        if not rows:
            return None
        fields = fields or list(rows[0].keys())
        with open(p, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        return p

    def add_records(self, records, replace=True):
        if replace:
            self.records = []
            self.validations = {}
        seen = {r.id for r in self.records}
        for r in records:
            base, k = r.id, 2
            while r.id in seen:
                r.id = f"{base}_{k}"
                k += 1
            seen.add(r.id)
            self.records.append(r)
            self.validations[r.id] = val.validate_record(r)

    def usable(self, task):
        """Return records that passed validation; explain which were skipped."""
        if not self.records:
            print("  No sequences loaded. Use option 1 first.")
            return []
        good = []
        for r in self.records:
            if val.is_usable(self.validations[r.id]):
                good.append(r)
            else:
                print(f"  [skipped for {task}] {r.id}: {self.validations[r.id]['status']} "
                      f"(fix the invalid characters first)")
        if not good:
            print("  No valid sequences available for this task.")
        return good


def task_load_file(s, path):
    recs = ih.load_file(path)
    s.add_records(recs)
    print(f"  Loaded {len(recs)} record(s) from '{path}'.")


def task_validate(s):
    if not s.records:
        print("  No sequences loaded. Use option 1 first.")
        return
    results = [s.validations[r.id] for r in s.records]
    text = val.format_validation_report(results)
    print(text)
    s.sections["Task 1 - Validation"] = text
    s.write_text("task1_validation_report.txt", text)
    s.write_csv("task1_validation.csv", val.validation_rows(results))
    print("  Saved task1_validation_report.txt and task1_validation.csv")


def task_composition(s, window, step, include_partial=False):
    recs = s.usable("composition")
    if not recs:
        return
    all_stats, all_rows, plots = [], [], []
    notes = []
    for r in recs:
        all_stats.append((r.id, comp.compute_stats(r.sequence)))
        w, st = window, step
        n = len(r.sequence)
        if 0 < n < window:
            w = max(1, n // 4)
            st = max(1, w // 2)
            notes.append(f"  Note: {r.id} ({n} bp) is shorter than window {window}; "
                         f"used window={w}, step={st} for this sequence.")
        try:
            wins = comp.sliding_windows(r.sequence, w, st, include_partial)
        except ValueError as exc:
            notes.append(f"  Sliding window for {r.id}: {exc}")
            continue
        all_rows += comp.window_csv_rows(r.id, wins)
        plots += plotting.plot_window_metrics(r.id, wins, w, st, s.path("plots", ""))
    text = comp.format_summary(all_stats)
    text += (f"\n\nSliding window: size={window}, step={step}, "
             f"{'partial final window included' if include_partial else 'complete windows only'}; "
             "coordinates 1-based inclusive.")
    if notes:
        text += "\n" + "\n".join(notes)
    text += f"\nPlots written: {len(plots)} PNG file(s) in output/plots/"
    print(text)
    s.sections["Task 2 - Composition"] = text
    s.write_text("task2_composition_report.txt", text)
    s.write_csv("task2_summary_statistics.csv", comp.summary_rows(all_stats))
    s.write_csv("task2_sliding_windows.csv", all_rows)
    print("  Saved task2 report + CSVs.")


def task_transform(s, start=None, end=None):
    recs = s.usable("transformation")
    if not recs:
        return
    parts = ["TASK 3 - DNA SEQUENCE TRANSFORMATION REPORT", "=" * 60]
    for r in recs:
        parts.append(tr.transformation_report(r.id, r.sequence, start, end))
    text = "\n".join(parts)
    print(text if sum(len(r.sequence) for r in recs) < 5000 else text[:3000] + "\n  ... (full text in file)")
    s.sections["Task 3 - Transformation"] = text
    s.write_text("task3_transformation_report.txt", text)
    print("  Saved task3_transformation_report.txt")


def task_motifs(s, motifs, overlapping, run_benchmark=True, sizes=(1000, 100000, 1000000)):
    recs = s.usable("motif detection")
    if not recs:
        return
    parts = ["TASK 4 - PATTERN MATCHING AND MOTIF DETECTION REPORT", "=" * 60,
             f"Motifs: {', '.join(motifs)}", ""]
    match_rows, pal_rows, tr_rows = [], [], []
    for r in recs:
        text, rows, pals, trs = ma.analyse_record(r.id, r.sequence, motifs, overlapping)
        parts.append(text)
        match_rows += rows
        pal_rows += [{"sequence_id": r.id, "start": a, "end": b, "palindrome": c} for a, b, c in pals]
        tr_rows += [{"sequence_id": r.id, **t} for t in trs]
    if run_benchmark:
        print("  Running benchmark (this can take a few seconds)...")
        bench = ma.benchmark(sizes=sizes)
        parts += ["", ma.format_benchmark(bench)]
        s.write_csv("task4_benchmark.csv", bench)
        plot = plotting.plot_benchmark(bench, s.path("plots", ""))
        parts.append(f"Benchmark plot: {plot}")
    text = "\n".join(parts)
    print(text)
    s.sections["Task 4 - Motifs"] = text
    s.write_text("task4_motif_report.txt", text)
    s.write_csv("task4_motif_matches.csv", match_rows, ["sequence_id", "motif", "start", "end", "mode"])
    s.write_csv("task4_palindromes.csv", pal_rows, ["sequence_id", "start", "end", "palindrome"])
    s.write_csv("task4_tandem_repeats.csv", tr_rows, ["sequence_id", "start", "end", "unit", "copies"])
    print("  Saved task4 report + CSVs.")


def export_combined(s):
    if not s.sections:
        print("  Nothing to export yet - run at least one analysis first.")
        return
    head = [f"DNA TOOLKIT COMBINED REPORT - generated {datetime.now():%Y-%m-%d %H:%M:%S}",
            f"Sequences analysed: {', '.join(r.id for r in s.records)}", ""]
    body = []
    for name in sorted(s.sections):
        body += [f"##### {name} #####", s.sections[name], ""]
    p = s.write_text("combined_report.txt", "\n".join(head + body))
    print(f"  Combined report saved: {p}")


def run_all(s, window, step, motifs, overlapping, start, end, sizes):
    task_validate(s)
    task_composition(s, window, step)
    task_transform(s, start, end)
    task_motifs(s, motifs, overlapping, True, sizes)
    export_combined(s)


def ask(prompt, default=None):
    suffix = f" [{default}]" if default is not None else ""
    ans = input(f"{prompt}{suffix}: ").strip()
    return ans if ans else (str(default) if default is not None else "")


def ask_int(prompt, default=None):
    while True:
        raw = ask(prompt, default)
        try:
            return int(raw)
        except ValueError:
            print(f"  '{raw}' is not a whole number; try again.")


def ask_yes(prompt, default="y"):
    return ask(prompt + " (y/n)", default).lower().startswith("y")


def menu_load(s):
    print("  a) Type a sequence   b) FASTA / plain-text file")
    choice = ask("Choose", "b").lower()
    try:
        if choice.startswith("a"):
            print("  Enter sequence (or FASTA text on one line; blank to cancel):")
            raw = input("  > ")
            recs = ih.from_manual(raw, index=len(s.records) + 1)
            s.add_records(recs, replace=not ask_yes("Add to already loaded sequences?", "n"))
            print(f"  Loaded {len(recs)} record(s).")
        else:
            path = ask("File path", os.path.join(BASE, "data", "sample_sequences.fasta"))
            task_load_file(s, path)
    except (ValueError, OSError) as exc:
        print(f"  ERROR: {exc}")
        return
    for r in s.records:
        v = s.validations[r.id]
        print(f"   - {r.id}: {v['length']} symbols, {v['status']}")


def interactive(s):
    window, step = 100, 50
    while True:
        print(MENU)
        try:
            choice = ask("Select option")
            if choice == "0":
                print("Goodbye.")
                return
            elif choice == "1":
                menu_load(s)
            elif choice == "2":
                task_validate(s)
            elif choice == "3":
                window = ask_int("Window size", window)
                step = ask_int("Step size", step)
                partial = ask_yes("Include a final partial window?", "n")
                task_composition(s, window, step, partial)
            elif choice == "4":
                use = ask_yes("Extract a subsequence too?", "y")
                st = en = None
                if use:
                    st = ask_int("Start (1-based)", 1)
                    en = ask_int("End (1-based, inclusive)", 10)
                task_transform(s, st, en)
            elif choice == "5":
                raw = ask("Motif(s), comma-separated", ",".join(DEFAULT_MOTIFS))
                try:
                    motifs = [ma.clean_motif(m) for m in raw.split(",") if m.strip()]
                except ValueError as exc:
                    print(f"  ERROR: {exc}")
                    continue
                if not motifs:
                    print("  ERROR: no motif given.")
                    continue
                overlapping = ask_yes("Report overlapping matches?", "y")
                bench = ask_yes("Run benchmark (1k, 100k, 1M bases)?", "y")
                task_motifs(s, motifs, overlapping, bench)
            elif choice == "6":
                if not s.records:
                    print("  No sequences loaded - loading the sample file.")
                    task_load_file(s, os.path.join(BASE, "data", "sample_sequences.fasta"))
                print(f"  Using defaults: window=100, step=50, motifs={DEFAULT_MOTIFS}, overlapping, subsequence 2-10.")
                run_all(s, 100, 50, DEFAULT_MOTIFS, True, 2, 10, (1000, 100000, 1000000))
            elif choice == "7":
                export_combined(s)
            else:
                print("  Unknown option; enter a number from the menu.")
        except (ValueError, OSError) as exc:
            print(f"  ERROR: {exc}")
        except (EOFError, KeyboardInterrupt):
            print("\nExiting.")
            return


def main():
    ap = argparse.ArgumentParser(description="DNA sequence analysis toolkit")
    ap.add_argument("--input", help="FASTA or plain-text file to load")
    ap.add_argument("--run-all", action="store_true", help="run the complete analysis non-interactively")
    ap.add_argument("--window", type=int, default=100)
    ap.add_argument("--step", type=int, default=50)
    ap.add_argument("--motifs", default=",".join(DEFAULT_MOTIFS))
    ap.add_argument("--non-overlapping", action="store_true")
    ap.add_argument("--coords", nargs=2, type=int, metavar=("START", "END"), default=[2, 10])
    ap.add_argument("--outdir", default=os.path.join(BASE, "output"))
    ap.add_argument("--quick", action="store_true", help="benchmark only 1k and 10k bases")
    a = ap.parse_args()

    s = Session(a.outdir)
    try:
        if a.input:
            task_load_file(s, a.input)
        if a.run_all:
            if not s.records:
                task_load_file(s, os.path.join(BASE, "data", "sample_sequences.fasta"))
            motifs = [ma.clean_motif(m) for m in a.motifs.split(",") if m.strip()]
            sizes = (1000, 10000) if a.quick else (1000, 100000, 1000000)
            run_all(s, a.window, a.step, motifs, not a.non_overlapping, a.coords[0], a.coords[1], sizes)
            return 0
    except (ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    interactive(s)
    return 0


if __name__ == "__main__":
    sys.exit(main())
