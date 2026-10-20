"""Automated tests for the FASTA module."""

import json
import tempfile
import unittest
from pathlib import Path

from bioeasy.fasta import (
    FastaFormatError,
    calculate_fasta_stats,
    parse_fasta,
)


class TestFastaStats(unittest.TestCase):

    def setUp(self):
        self.temp_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_directory.cleanup)
        self.root = Path(self.temp_directory.name)
        self.fasta_path = self.root / "sample.fasta"

        self.fasta_path.write_text(
            ">seq1 Example sequence\n"
            "ACGT\n"
            ">seq2\n"
            "GG\n",
            encoding="utf-8",
        )

    def test_parse_fasta_records(self):
        records = list(parse_fasta(self.fasta_path))

        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].identifier, "seq1")
        self.assertEqual(records[0].sequence, "ACGT")
        self.assertEqual(records[1].sequence, "GG")

    def test_calculate_statistics(self):
        stats = calculate_fasta_stats(self.fasta_path)

        self.assertEqual(stats.sequence_count, 2)
        self.assertEqual(stats.total_length, 6)
        self.assertEqual(stats.longest_sequence, 4)
        self.assertEqual(stats.n50, 4)
        self.assertEqual(stats.gc_content_percent, 66.67)
        self.assertEqual(stats.canonical_nucleotides, 6)

    def test_reject_sequence_before_header(self):
        bad_file = self.root / "invalid.fasta"
        bad_file.write_text("ACGT\n", encoding="utf-8")

        with self.assertRaises(FastaFormatError):
            list(parse_fasta(bad_file))

    def test_reject_empty_file(self):
        empty_file = self.root / "empty.fasta"
        empty_file.write_text("", encoding="utf-8")

        with self.assertRaises(FastaFormatError):
            calculate_fasta_stats(empty_file)

    def test_reject_empty_sequence_record(self):
        bad_file = self.root / "empty_record.fasta"
        bad_file.write_text(
            ">seq1\n>seq2\nACGT\n",
            encoding="utf-8",
        )

        with self.assertRaises(FastaFormatError):
            list(parse_fasta(bad_file))

    def test_report_is_json_serializable(self):
        stats = calculate_fasta_stats(self.fasta_path)
        encoded = json.dumps(stats.to_dict())
        decoded = json.loads(encoded)

        self.assertEqual(decoded["sequence_count"], 2)


if __name__ == "__main__":
    unittest.main()
