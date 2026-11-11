"""Unit tests (run:  python -m unittest discover -s tests -v   from dna_toolkit/)."""
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import composition as comp
import input_handler as ih
import motif_analysis as ma
import transformation as tr
import validation as val


class TestValidation(unittest.TestCase):
    def status(self, s):
        return val.validate_sequence(ih.normalize(s))["status"]

    def test_cases_from_manual(self):
        self.assertEqual(self.status("ATGCGT"), "Valid")
        self.assertEqual(self.status("atgcgt"), "Valid")
        self.assertEqual(self.status("ATGNCGT"), "Valid with ambiguity")
        self.assertEqual(self.status("ATGRY"), "Valid with ambiguity")
        self.assertEqual(self.status("ATGXZ"), "Invalid")
        self.assertTrue(self.status("").startswith("Invalid"))

    def test_sample_record(self):
        recs = ih.parse_fasta_text(">seq_A\nATGcnXTA\n")
        self.assertEqual(recs[0].sequence, "ATGCNXTA")
        r = val.validate_record(recs[0])
        self.assertEqual(r["length"], 8)
        self.assertEqual(r["invalid"], [("X", 6)])
        self.assertEqual(r["ambiguous_positions"], [("N", 5)])
        self.assertEqual(r["status"], "Invalid")

    def test_fasta_multi_record_and_whitespace(self):
        recs = ih.parse_fasta_text(">a desc here\nAT GC\n\tGG\n>b\nTT\n")
        self.assertEqual([r.id for r in recs], ["a", "b"])
        self.assertEqual(recs[0].description, "desc here")
        self.assertEqual(recs[0].sequence, "ATGCGG")

    def test_missing_file(self):
        with self.assertRaises(FileNotFoundError):
            ih.load_file("definitely_missing_file.fasta")

    def test_empty_manual(self):
        with self.assertRaises(ValueError):
            ih.from_manual("   ")


class TestComposition(unittest.TestCase):
    def test_worked_example_ambiguous(self):
        s = comp.compute_stats("ATGCRYN")
        self.assertEqual((s["A"], s["T"], s["G"], s["C"], s["L"]), (1, 1, 1, 1, 4))
        self.assertEqual(s["ambiguous_count"], 3)
        self.assertAlmostEqual(s["gc_content"], 50)
        self.assertAlmostEqual(s["gc_skew"], 0)
        self.assertAlmostEqual(s["entropy"], 2)
        self.assertAlmostEqual(s["ambiguous_pct"], 3 / 7 * 100)

    def test_example_12(self):
        s = comp.compute_stats("ATGCGCGCATAT")
        self.assertEqual((s["A"], s["T"], s["G"], s["C"]), (3, 3, 3, 3))
        self.assertAlmostEqual(s["entropy"], 2)

    def test_na_cases(self):
        s = comp.compute_stats("NNNN")
        self.assertIsNone(s["gc_content"])
        self.assertIsNone(s["entropy"])
        self.assertIsNone(comp.compute_stats("AAAA")["gc_skew"])
        self.assertIsNone(comp.compute_stats("GGGG")["at_skew"])
        self.assertEqual(comp.compute_stats("AAAA")["entropy"], 0)

    def test_windows(self):
        w = comp.sliding_windows("ATGCGCGCATAT", 6, 6)
        self.assertEqual(len(w), 2)
        self.assertAlmostEqual(w[0]["gc_content"], 200 / 3)
        self.assertAlmostEqual(w[1]["gc_content"], 100 / 3)
        self.assertEqual((w[0]["start"], w[0]["end"], w[1]["start"], w[1]["end"]), (1, 6, 7, 12))

    def test_partial_window_and_errors(self):
        self.assertEqual(len(comp.sliding_windows("ATGCGCGCATAT", 5, 5)), 2)
        self.assertEqual(len(comp.sliding_windows("ATGCGCGCATAT", 5, 5, True)), 3)
        with self.assertRaises(ValueError):
            comp.sliding_windows("ATG", 10, 1)
        with self.assertRaises(ValueError):
            comp.sliding_windows("ATG", 0, 1)


class TestTransformation(unittest.TestCase):
    def test_manual_example(self):
        d = "ATGCCG"
        self.assertEqual(tr.complement(d), "TACGGC")
        self.assertEqual(tr.reverse(d), "GCCGTA")
        self.assertEqual(tr.reverse_complement(d), "CGGCAT")
        self.assertEqual(tr.dna_to_rna(d), "AUGCCG")
        self.assertEqual(tr.rna_to_dna("AUGCCG"), d)
        self.assertEqual(tr.subsequence(d, 2, 5), "TGCC")

    def test_iupac_and_involution(self):
        s = "ATGCRYSWKMBDHVN"
        self.assertEqual(tr.complement(s), "TACGYRSWMKVHDBN")
        self.assertEqual(tr.reverse_complement(tr.reverse_complement(s)), s)
        self.assertEqual(len(tr.reverse_complement(s)), len(s))

    def test_bad_coordinates(self):
        for a, b in ((0, 3), (2, 99), (5, 2)):
            with self.assertRaises(ValueError):
                tr.subsequence("ATGCCG", a, b)


class TestMotifs(unittest.TestCase):
    def test_overlap(self):
        for f in (ma.naive_search, ma.kmp_search):
            self.assertEqual(f("AAAAA", "AAA", True), [1, 2, 3])
            self.assertEqual(f("AAAAA", "AAA", False), [1])
            self.assertEqual(f("ACGT", "TTTTT"), [])

    def test_lps(self):
        self.assertEqual(ma.build_lps("AABAACAABAA"), [0, 1, 0, 1, 2, 0, 1, 2, 3, 4, 5])

    def test_algorithms_agree_randomly(self):
        import random
        rng = random.Random(1)
        for _ in range(200):
            t = "".join(rng.choice("AC") for _ in range(rng.randint(0, 40)))
            p = "".join(rng.choice("AC") for _ in range(rng.randint(1, 5)))
            for ov in (True, False):
                self.assertEqual(ma.naive_search(t, p, ov), ma.kmp_search(t, p, ov))

    def test_palindrome(self):
        self.assertTrue(ma.is_rc_palindrome("GAATTC"))
        self.assertFalse(ma.is_rc_palindrome("GAATTG"))
        self.assertIn((3, 8, "GAATTC"), ma.find_palindromes("AAGAATTCAA", 6, 6))

    def test_tandem(self):
        r = ma.find_tandem_repeats("ATATAT", 1, 6, 3)
        self.assertEqual(len(r), 1)
        self.assertEqual((r[0]["unit"], r[0]["copies"], r[0]["start"], r[0]["end"]), ("AT", 3, 1, 6))
        r = ma.tandem_runs_for_motif("GCATCATCATGG", "CAT", 2)
        self.assertEqual((r[0]["start"], r[0]["copies"]), (2, 3))

    def test_clean_motif(self):
        self.assertEqual(ma.clean_motif(" gaa ttc"), "GAATTC")
        with self.assertRaises(ValueError):
            ma.clean_motif("GAXTC")

    def test_benchmark_runs(self):
        rows = ma.benchmark(sizes=(500,), reps=1)
        self.assertTrue(all(r["results_agree"] for r in rows))


if __name__ == "__main__":
    unittest.main()
