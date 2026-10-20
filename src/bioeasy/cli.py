"""Command-line interface for BioEasy."""

import argparse
import json
import sys
from pathlib import Path

from bioeasy import __version__
from bioeasy.fasta import FastaFormatError, calculate_fasta_stats


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        prog="bioeasy",
        description=(
            "Run common bioinformatics tasks with simple commands."
        ),
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    commands = parser.add_subparsers(
        dest="command",
        required=True,
    )

    fasta_parser = commands.add_parser(
        "fasta-stats",
        help="Calculate statistics from a FASTA file.",
        description=(
            "Calculate sequence count, total length, longest sequence, "
            "N50 and GC content."
        ),
    )

    fasta_parser.add_argument(
        "input",
        type=Path,
        help="Path to the input FASTA file.",
    )

    fasta_parser.add_argument(
        "--output",
        type=Path,
        help="Optional path for a JSON report.",
    )

    return parser


def _print_fasta_report(stats: dict) -> None:
    """Print a readable FASTA statistics report."""
    print()
    print("BioEasy | FASTA Statistics")
    print("=" * 30)
    print(f"Input file:             {stats['input_file']}")
    print(f"Number of sequences:    {stats['sequence_count']:,}")
    print(f"Total length:           {stats['total_length']:,}")
    print(f"Longest sequence:       {stats['longest_sequence']:,}")
    print(f"N50:                    {stats['n50']:,}")

    gc = stats["gc_content_percent"]

    if gc is None:
        print("GC content (%):         N/A")
        print("Reason:                 No canonical nucleotides found")
    else:
        print(f"GC content (%):         {gc:.2f}")

    print(
        "Canonical nucleotides:  "
        f"{stats['canonical_nucleotides']:,}"
    )
    print("=" * 30)
    print()


def main(argv: list[str] | None = None) -> int:
    """Run BioEasy and return a process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "fasta-stats":
        try:
            result = calculate_fasta_stats(args.input)
            report = result.to_dict()

            _print_fasta_report(report)

            if args.output is not None:
                output_path = args.output
                output_path.parent.mkdir(
                    parents=True,
                    exist_ok=True,
                )

                output_path.write_text(
                    json.dumps(report, indent=2) + "\n",
                    encoding="utf-8",
                )

                print(f"JSON report saved to: {output_path}")

            return 0

        except (OSError, FastaFormatError, ValueError) as exc:
            print(f"bioeasy: error: {exc}", file=sys.stderr)
            return 2

    parser.error("Unknown command.")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
