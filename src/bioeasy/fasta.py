"""FASTA parsing and sequence-statistics calculations."""

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterator


class FastaFormatError(ValueError):
    """Raised when a FASTA file is empty or incorrectly formatted."""


@dataclass(frozen=True)
class FastaRecord:
    """A single FASTA record."""

    identifier: str
    description: str
    sequence: str


@dataclass(frozen=True)
class FastaStats:
    """Statistics calculated from a FASTA file."""

    input_file: str
    sequence_count: int
    total_length: int
    longest_sequence: int
    n50: int
    gc_content_percent: float | None
    canonical_nucleotides: int

    def to_dict(self) -> dict:
        """Return a JSON-compatible representation."""
        return asdict(self)


def _build_record(header: str, sequence_parts: list[str]) -> FastaRecord:
    """Convert a FASTA header and sequence lines into a record."""
    fields = header.split(maxsplit=1)

    return FastaRecord(
        identifier=fields[0],
        description=fields[1] if len(fields) > 1 else "",
        sequence="".join(sequence_parts),
    )


def parse_fasta(file_path: str | Path) -> Iterator[FastaRecord]:
    """Parse a FASTA file one record at a time.

    Supports multiline sequences, blank lines, and UTF-8 text files.
    Does not restrict the sequence alphabet, so it can also read
    protein FASTA files.
    """
    path = Path(file_path)

    if not path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    header = None
    sequence_parts: list[str] = []

    try:
        with path.open("r", encoding="utf-8-sig") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                line = raw_line.strip()

                if not line:
                    continue

                if line.startswith(">"):
                    if header is not None:
                        if not sequence_parts:
                            raise FastaFormatError(
                                f"Empty sequence record before line "
                                f"{line_number} in {path}"
                            )

                        yield _build_record(header, sequence_parts)

                    header = line[1:].strip()

                    if not header:
                        raise FastaFormatError(
                            f"Empty FASTA header on line {line_number}"
                        )

                    sequence_parts = []
                    continue

                if header is None:
                    raise FastaFormatError(
                        f"Sequence data before the first FASTA header "
                        f"on line {line_number} in {path}"
                    )

                # Remove whitespace inside sequence lines and normalize case.
                cleaned = "".join(line.split()).upper()

                if cleaned:
                    sequence_parts.append(cleaned)

    except UnicodeDecodeError as exc:
        raise FastaFormatError(
            f"Could not decode {path} as UTF-8 text."
        ) from exc

    if header is not None:
        if not sequence_parts:
            raise FastaFormatError(
                f"Empty sequence record at the end of {path}"
            )

        yield _build_record(header, sequence_parts)


def calculate_fasta_stats(file_path: str | Path) -> FastaStats:
    """Calculate common sequence statistics for a FASTA file."""
    path = Path(file_path)
    records = list(parse_fasta(path))

    if not records:
        raise FastaFormatError(f"No FASTA records found in {path}")

    lengths = sorted(
        (len(record.sequence) for record in records),
        reverse=True,
    )

    total_length = sum(lengths)

    # N50: the length at which cumulative sequence length reaches
    # at least half the total length.
    cumulative_length = 0
    n50 = 0

    for length in lengths:
        cumulative_length += length

        if cumulative_length >= total_length / 2:
            n50 = length
            break

    gc_count = 0
    canonical_count = 0

    for record in records:
        sequence = record.sequence

        gc_count += sequence.count("G") + sequence.count("C")
        canonical_count += sum(
            sequence.count(base) for base in "ACGTU"
        )

    gc_percent = (
        round(100.0 * gc_count / canonical_count, 2)
        if canonical_count
        else None
    )

    return FastaStats(
        input_file=str(path),
        sequence_count=len(records),
        total_length=total_length,
        longest_sequence=max(lengths),
        n50=n50,
        gc_content_percent=gc_percent,
        canonical_nucleotides=canonical_count,
    )
