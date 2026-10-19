# alignkit

**Offline pairwise sequence alignment from FASTA files, with exact scores, CIGAR strings and CI identity gates.**

alignkit aligns DNA or protein sequences with global, local and semi-global dynamic programming, using affine gap penalties. It reports the exact score, the aligned columns, the CIGAR string and identity, and it lints FASTA files with line-numbered errors. It is a single Python file that uses only the standard library. It makes no network calls and needs no accounts or API keys.

## Features

- **Three alignment modes:** global (Needleman-Wunsch), local (Smith-Waterman) and semi-global, where leading and trailing overhangs are free.
- **Affine gaps (Gotoh):** a gap of length *k* costs `gap_open + (k-1) * gap_extend`, so one long gap beats several short ones.
- **Reverse-complement search:** `--strand both` also tries the reverse complement of the second sequence and keeps the better hit, with IUPAC codes supported.
- **Exact Levenshtein distance:** unit-cost edit distance with one optimal alignment, plus Hamming distance when the lengths match.
- **CIGAR strings** for every alignment (`=` match, `X` mismatch, `I` base only in the second sequence, `D` base only in the first).
- **Readable alignment blocks:** three-row output with residue positions on each side and a match bar (`|` match, `.` mismatch).
- **Full statistics:** score, columns, matches, mismatches, gap columns, gap opens and identity as an exact fraction.
- **FASTA linting** with `file:line` messages for invalid characters (with column), duplicate or empty names, empty records, sequence before the first header, and lowercase soft-masking. A `--alphabet dna` mode checks bases and flags IUPAC ambiguity codes.
- **Pairwise identity table** for every pair in a file, with `--min-identity` as a CI gate that exits 1 if any pair falls short.
- **Identity gate for single alignments** with `align --min-identity`.
- **Reverse-complement FASTA output** through `revcomp`.
- **JSON output** for `align`, `pairs` and `distance`, for scripts and CI.
- **Safety limit:** `--max-cells` stops runaway memory use on very long sequences before the matrices are allocated.
- **Deterministic results:** ties are broken in a fixed order, so repeated runs give the same alignment.

## Quick start

```sh
python3 alignkit.py align examples/pair.fasta examples/panel.fasta
```

This aligns the first record of each file, `wildtype exon 3` against `isolate_01`, using the default scoring of match +2, mismatch -1, gap open -4 and gap extend -1.

Other useful commands:

```sh
python3 alignkit.py align examples/pair.fasta examples/panel.fasta --mode local --json
python3 alignkit.py pairs examples/panel.fasta --min-identity 0.8
python3 alignkit.py distance examples/pair.fasta examples/panel.fasta
python3 alignkit.py lint examples/broken.fasta
```

## Commands

| Command | Purpose | Exit codes |
| --- | --- | --- |
| `align A B` | Align the first record of each file. Options: `--mode`, `--match`, `--mismatch`, `--gap-open`, `--gap-extend`, `--strand`, `--min-identity`, `--width`, `--json`, `--max-cells`. | 0 ok, 1 identity below gate, 2 bad input or options |
| `pairs FILE` | Identity table for every pair of records, with the same scoring options and `--min-identity`. | 0 ok, 1 a pair is below the gate, 2 bad input |
| `distance A B` | Levenshtein distance, Hamming distance, CIGAR and alignment blocks. | 0 ok, 2 bad input |
| `lint FILE` | Line-numbered FASTA findings. Options: `--alphabet any\|dna`, `--strict`. | 0 clean, 1 errors (or warnings with `--strict`), 2 unreadable file |
| `revcomp FILE` | Print the reverse complement of every DNA record. | 0 ok, 2 bad input |

Each input file must be valid FASTA. `align` and `distance` use the first record in each file, and `pairs` uses all of them.

## How to run

Needs Python 3.9 or newer. There is nothing to install:

```sh
python3 alignkit.py --help
python3 -m alignkit align examples/pair.fasta examples/panel.fasta
```

## How to test

The tests use pytest. Install it once, then run the suite from the project root:

```sh
python3 -m pip install pytest
python3 -m pytest
```

The tests cover FASTA parsing and line-numbered errors, identical and embedded-motif alignments, semi-global free end gaps, affine gap runs, a randomized check that rescores every global and local alignment column by column, Levenshtein distances and CIGAR strings, reverse complements, strand selection, scoring and memory validation, identity gates and lint exit codes.

## Layout

```
alignkit.py            the whole tool: parser, DP core, CLI
examples/pair.fasta    two sequences, for pairwise commands
examples/panel.fasta   four clean sequences, for pairs and lint
examples/broken.fasta  deliberate FASTA errors, for lint
tests/                 pytest suite
```

## Limits

- Scoring is match/mismatch with affine gaps. Substitution matrices such as BLOSUM are not built in.
- Only the first record of each file is aligned by `align` and `distance`.
- The dynamic programme keeps three matrices of size (n+1)×(m+1) in memory. The default limit is 2,000,000 cells; raise it with `--max-cells` if you have the memory.
