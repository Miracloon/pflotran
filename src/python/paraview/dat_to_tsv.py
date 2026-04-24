"""Convert a PFLOTRAN .dat file to a tab-delimited .tsv file.

The .dat format has:
  - Line 1: comma-separated, quoted column headers
  - Remaining lines: whitespace-delimited numeric data

Usage:
    python dat_to_tsv.py <input.dat> [output.tsv]

If output path is omitted, the output is written alongside the input file
with the same name but a .tsv extension.
"""

import sys
import csv
from pathlib import Path


def convert(input_path: Path, output_path: Path) -> None:
    with open(input_path, newline='') as fh:
        lines = fh.readlines()

    if not lines:
        print(f"ERROR: {input_path} is empty.")
        sys.exit(1)

    # Parse quoted, comma-separated header (skipinitialspace handles leading whitespace
    # before the first quoted field, preventing stray quotes in the output)
    header = next(csv.reader([lines[0]], skipinitialspace=True))
    header = [col.strip() for col in header]

    with open(output_path, 'w', newline='') as fh:
        fh.write('\t'.join(header) + '\n')
        for lineno, line in enumerate(lines[1:], start=2):
            line = line.strip()
            if not line:
                continue
            fields = line.split()
            fh.write('\t'.join(fields) + '\n')

    print(f"Converted {input_path} -> {output_path}  ({len(lines) - 1} data rows)")


def main():
    if len(sys.argv) < 2:
        print("Usage: dat_to_tsv.py <input.dat> [output.tsv]")
        sys.exit(1)

    input_path = Path(sys.argv[1])
    if not input_path.exists():
        print(f"ERROR: File not found: {input_path}")
        sys.exit(1)

    output_path = Path(sys.argv[2]) if len(sys.argv) >= 3 else input_path.with_suffix('.tsv')
    convert(input_path, output_path)


if __name__ == '__main__':
    main()
