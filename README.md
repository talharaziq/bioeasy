# BioEasy

BioEasy is an open-source command-line toolkit that makes repetitive
bioinformatics tasks easier to run and reproduce.

## Current features

- FASTA parsing
- Sequence count and total sequence length
- Longest sequence and N50
- GC content calculated over canonical A, C, G, T and U characters
- Human-readable reports and optional JSON output
- Input validation and automated tests

## Installation

Python 3.10 or newer is required.

```bash
python -m pip install -e .
```

## Usage

```bash
bioeasy --help
bioeasy fasta-stats examples/sample.fasta
bioeasy fasta-stats examples/sample.fasta --output reports/sample_stats.json
```

You can also run BioEasy as a Python module:

```bash
python -m bioeasy fasta-stats examples/sample.fasta
```

## Development tests

```bash
python -m unittest discover -s tests -v
```

## GC content and sequence statistics

GC content is calculated using only canonical A, C, G, T and U
characters. Other characters contribute to sequence length but not
to the GC-content denominator. The N50 calculation uses all sequence
lengths in the input file.

BioEasy is a research utility. Its output should be checked before
being used in scientific publications or clinical settings.

## License

MIT. See LICENSE.
