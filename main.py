#!/usr/bin/env python3
"""
Recode aligned DNA sequences into a 0–3 STANDARD NEXUS matrix
according to the procedure in To_John_Franco.txt:

Part 1 (global mode):
    - Count A, C, G, T over the entire matrix (ignoring '?').
    - Assign weights 1, 2, 4, 8 so that the most frequent nucleotide
      gets weight 1, next 2, next 4, least frequent 8 (minimizes sum).
    - Map 1 -> 0, 2 -> 1, 4 -> 2, 8 -> 3 to get the final symbols 0–3.
    - Use that same mapping for all sites.

Part 2 (column mode):
    - Do the *same* procedure independently for each column:
        * count A, C, G, T in that column (ignoring '?'),
        * assign 1,2,4,8 by frequency as above,
        * map to 0–3 column-specifically.
    - Thus each column can have its own code mapping.

Usage examples (from the shell):

    python recode_alignment.py input.txt --mode global      > output_global.nex
    python recode_alignment.py input.txt --mode column      > output_column.nex
    python recode_alignment.py input.txt --mode global --add-outgroup > with_outgroup.nex
"""

import argparse
from collections import Counter, defaultdict
from typing import Dict, List, Tuple


NUCS = ("A", "C", "G", "T")
WEIGHTS = [1, 2, 4, 8]
WEIGHT_TO_CODE = {1: 0, 2: 1, 4: 2, 8: 3}


def parse_alignment(path: str) -> Tuple[Dict[str, str], int, int]:
    """
    Parse a simple alignment text file.

    Expected format:
        <ntax> <nchar>
        <taxon_name><whitespace><sequence>
        ...

    - First non-empty, non-comment line: two integers.
    - Next <ntax> non-empty, non-comment lines: taxon + sequence.
      Taxon name is taken as the first field; sequence as the last field.
    """
    alignment: Dict[str, str] = {}

    with open(path, "r", encoding="utf-8") as f:
        # Read first non-empty, non-comment line for ntax, nchar
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError("First non-empty line must contain NTAX and NCHAR.")
            ntax_declared = int(parts[0])
            nchar_declared = int(parts[1])
            break
        else:
            raise ValueError("Empty file or no valid NTAX/NCHAR line found.")

        # Read taxon lines
        for line in f:
            line = line.rstrip("\n")
            if not line.strip() or line.strip().startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            taxon = parts[0]
            seq = parts[-1].strip()
            alignment[taxon] = seq

    ntax_actual = len(alignment)
    if ntax_actual != ntax_declared:
        raise ValueError(
            f"Declared NTAX={ntax_declared}, but read {ntax_actual} sequences."
        )

    for taxon, seq in alignment.items():
        if len(seq) != nchar_declared:
            raise ValueError(
                f"Sequence length mismatch for {taxon}: "
                f"expected {nchar_declared}, got {len(seq)}."
            )

    return alignment, ntax_declared, nchar_declared


def compute_global_mapping(alignment: Dict[str, str]) -> Dict[str, int]:
    """
    Count A/C/G/T over the entire alignment (ignoring '?' and '-')
    and produce a nucleotide -> {0,1,2,3} mapping using:

        - sort by decreasing frequency (ties -> alphabetical),
        - assign weights [1,2,4,8] in that order,
        - map weights to 0–3 via WEIGHT_TO_CODE.
    """
    total_counts = Counter()
    for seq in alignment.values():
        for ch in seq:
            ch = ch.upper()
            if ch in NUCS:
                total_counts[ch] += 1

    # If some nucleotides are absent, keep them with count 0 for determinism
    for n in NUCS:
        total_counts.setdefault(n, 0)

    # Sort by frequency desc, then alphabetically
    sorted_nucs = sorted(NUCS, key=lambda n: (-total_counts[n], n))

    mapping: Dict[str, int] = {}
    for nuc, weight in zip(sorted_nucs, WEIGHTS):
        mapping[nuc] = WEIGHT_TO_CODE[weight]

    return mapping


def compute_column_mappings(
    alignment: Dict[str, str], nchar: int
) -> List[Dict[str, int]]:
    """
    For each column j:
        - count A/C/G/T in that column (ignoring '?' and '-'),
        - sort by decreasing frequency (ties alphabetically),
        - assign weights 1,2,4,8 and map to 0–3.

    Returns list of length nchar, where each element is a dict nuc->code.
    """
    # Initialize per-column counters
    col_counts: List[Counter] = [Counter() for _ in range(nchar)]

    for seq in alignment.values():
        for j, ch in enumerate(seq):
            ch = ch.upper()
            if ch in NUCS:
                col_counts[j][ch] += 1

    col_mappings: List[Dict[str, int]] = []

    for j in range(nchar):
        counts = col_counts[j]

        # Ensure all four nucleotides are present (with count 0 if truly absent)
        for n in NUCS:
            counts.setdefault(n, 0)

        sorted_nucs = sorted(NUCS, key=lambda n: (-counts[n], n))

        mapping: Dict[str, int] = {}
        for nuc, weight in zip(sorted_nucs, WEIGHTS):
            mapping[nuc] = WEIGHT_TO_CODE[weight]

        col_mappings.append(mapping)

    return col_mappings


def encode_alignment_global(
    alignment: Dict[str, str], mapping: Dict[str, int]
) -> Dict[str, str]:
    """
    Apply a single global mapping nuc->code to every position.
    '?' and '-' are left as-is.
    """
    encoded: Dict[str, str] = {}

    for taxon, seq in alignment.items():
        out_chars: List[str] = []
        for ch in seq:
            uch = ch.upper()
            if uch in mapping:
                out_chars.append(str(mapping[uch]))
            elif ch in ("?", "-"):
                out_chars.append(ch)
            else:
                # Treat unknown/ambiguous as missing
                out_chars.append("?")
        encoded[taxon] = "".join(out_chars)

    return encoded


def encode_alignment_columnwise(
    alignment: Dict[str, str], col_mappings: List[Dict[str, int]]
) -> Dict[str, str]:
    """
    Apply column-specific mappings.
    col_mappings[j] is a dict nuc->code for column j.
    '?' and '-' are left as-is.
    """
    encoded: Dict[str, str] = {}
    nchar = len(col_mappings)

    for taxon, seq in alignment.items():
        if len(seq) != nchar:
            raise ValueError(
                f"Column mapping length mismatch for {taxon}: "
                f"expected {nchar}, got {len(seq)}"
            )

        out_chars: List[str] = []
        for j, ch in enumerate(seq):
            uch = ch.upper()
            if uch in col_mappings[j]:
                out_chars.append(str(col_mappings[j][uch]))
            elif ch in ("?", "-"):
                out_chars.append(ch)
            else:
                out_chars.append("?")
        encoded[taxon] = "".join(out_chars)

    return encoded


def format_nexus_matrix(
    encoded: Dict[str, str],
    add_outgroup: bool = False,
    outgroup_name: str = "Out",
) -> str:
    """
    Return the NEXUS MATRIX block as a string.

    If add_outgroup is True, append an 'Out' taxon coded as all zeros.
    (Useful if you want to replicate the example output exactly.)
    """
    taxa = list(encoded.keys())
    nchar = len(next(iter(encoded.values())))
    ntax = len(taxa)

    lines: List[str] = []
    lines.append("#NEXUS")
    lines.append("BEGIN DATA;")
    lines.append(f"    DIMENSIONS NTAX={ntax + (1 if add_outgroup else 0)} NCHAR={nchar};")
    lines.append(
        '    FORMAT DATATYPE=STANDARD MISSING=? GAP=- SYMBOLS="0123";'
    )
    lines.append("    MATRIX")

    # Determine width for taxon name column
    max_name_len = max(len(t) for t in taxa + ([outgroup_name] if add_outgroup else []))
    name_field = max_name_len + 2  # a bit of spacing

    for t in taxa:
        seq = encoded[t]
        lines.append(f"    {t.ljust(name_field)}{seq}")

    if add_outgroup:
        lines.append(
            f"    {outgroup_name.ljust(name_field)}" + ("0" * nchar)
        )

    lines.append("    ;")
    lines.append("END;")
    lines.append("")
    lines.append("BEGIN ASSUMPTIONS;")
    lines.append(f"\tTYPESET * UNTITLED   =  ord:  1- {nchar};")
    lines.append("END;")

    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recode nucleotide alignment into NEXUS 0–3 matrix "
        "using a global or per-column frequency-based mapping."
    )
    parser.add_argument("input", help="Path to input alignment text file.")
    parser.add_argument(
        "--mode",
        choices=["global", "column"],
        default="global",
        help=(
            "Recoding mode:\n"
            "  global  – use a single mapping for all sites (Part 1).\n"
            "  column  – compute mapping independently for each column (Part 2)."
        ),
    )
    parser.add_argument(
        "--add-outgroup",
        action="store_true",
        help="If set, add an 'Out' taxon with all zeros as in the example.",
    )

    args = parser.parse_args()

    alignment, ntax, nchar = parse_alignment(args.input)

    if args.mode == "global":
        mapping = compute_global_mapping(alignment)
        encoded = encode_alignment_global(alignment, mapping)
    else:  # column-wise
        col_mappings = compute_column_mappings(alignment, nchar)
        encoded = encode_alignment_columnwise(alignment, col_mappings)

    nexus_text = format_nexus_matrix(encoded, add_outgroup=args.add_outgroup)
    print(nexus_text)


if __name__ == "__main__":
    main()
