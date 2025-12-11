#!/usr/bin/env python3
"""
FRECode: Frequency-based Recoder for DNA Alignments

Implements the procedure described in the To_John_Franco.txt note:

- Global mode:
    * Count A, C, G, T over the entire alignment (ignoring '?' and '-').
    * Rank nucleotides by decreasing frequency (ties -> alphabetical).
    * Assign weights 1, 2, 4, 8 in that order.
    * Map weights to codes via: 1->0, 2->1, 4->2, 8->3.
    * Apply the same nuc->code mapping to all sites.

- Column mode:
    * Do the same ranking/weighting independently for each column.
    * Each column has its own nuc->code mapping based on column frequencies.

- Both mode:
    * Output both the global and column-wise recodings in one NEXUS file:
      first the global DATA block, then the column-wise DATA block.

Input format:

    <NTAX> <NCHAR>
    <TaxonName> <Sequence>
    <TaxonName> <Sequence>
    ...

Example:

    5 12
    Tax1   ACGTACGTACGT
    Tax2   A?GTACCTACGA
    Tax3   TCGTACGTACGA
    Tax4   ACGTACGTTCGT
    Tax5   GCGTACGTACGA
"""

import argparse
from collections import Counter
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

    - First non-empty, non-comment line: two integers (ntax, nchar).
    - Next ntax non-empty, non-comment lines: taxon + sequence.
      Taxon name is taken as the first field; sequence as the last field.
    """
    alignment: Dict[str, str] = {}

    with open(path, "r", encoding="utf-8") as f:
        # Get ntax, nchar from the first non-empty, non-comment line
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(
                    "First non-empty line must contain NTAX and NCHAR."
                )
            ntax_declared = int(parts[0])
            nchar_declared = int(parts[1])
            break
        else:
            raise ValueError("Empty file or no valid NTAX/NCHAR line found.")

        # Read taxon lines
        for line in f:
            raw = line.rstrip("\n")
            if not raw.strip() or raw.strip().startswith("#"):
                continue
            parts = raw.split()
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

    # Ensure all four nucleotides are present so the mapping is deterministic
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
    col_counts: List[Counter] = [Counter() for _ in range(nchar)]

    for seq in alignment.values():
        for j, ch in enumerate(seq):
            ch = ch.upper()
            if ch in NUCS:
                col_counts[j][ch] += 1

    col_mappings: List[Dict[str, int]] = []

    for j in range(nchar):
        counts = col_counts[j]

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
    '?' and '-' are left as-is; unknown symbols become '?'.
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
                out_chars.append("?")
        encoded[taxon] = "".join(out_chars)

    return encoded


def encode_alignment_columnwise(
    alignment: Dict[str, str], col_mappings: List[Dict[str, int]]
) -> Dict[str, str]:
    """
    Apply column-specific mappings.
    col_mappings[j] is a dict nuc->code for column j.
    '?' and '-' are left as-is; unknown symbols become '?'.
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
    include_header: bool = True,
) -> str:
    """
    Return a NEXUS DATA block (plus ASSUMPTIONS) as a string.

    If add_outgroup is True, append an 'Out' taxon coded as all zeros.
    If include_header is True, add the '#NEXUS' line at the top.
    """
    taxa = list(encoded.keys())
    nchar = len(next(iter(encoded.values())))
    ntax = len(taxa)

    lines: List[str] = []
    if include_header:
        lines.append("#NEXUS")
    lines.append("BEGIN DATA;")
    lines.append(
        f"    DIMENSIONS NTAX={ntax + (1 if add_outgroup else 0)} "
        f"NCHAR={nchar};"
    )
    lines.append(
        '    FORMAT DATATYPE=STANDARD MISSING=? GAP=- SYMBOLS="0123";'
    )
    lines.append("    MATRIX")

    max_name_len = max(len(t) for t in taxa + ([outgroup_name] if add_outgroup else []))
    name_field = max_name_len + 2

    for t in taxa:
        seq = encoded[t]
        lines.append(f"    {t.ljust(name_field)}{seq}")

    if add_outgroup:
        lines.append(f"    {outgroup_name.ljust(name_field)}" + ("0" * nchar))

    lines.append("    ;")
    lines.append("END;")
    lines.append("")
    lines.append("BEGIN ASSUMPTIONS;")
    lines.append(f"\tTYPESET * UNTITLED   =  ord:  1- {nchar};")
    lines.append("END;")

    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "FRECode: recode nucleotide alignments into 0–3 NEXUS matrices "
            "using global or column-wise frequency-based mappings."
        )
    )
    parser.add_argument("input", help="Path to input alignment text file.")
    parser.add_argument(
        "--mode",
        choices=["global", "column", "both"],
        default="global",
        help=(
            "Recoding mode:\n"
            "  global  – use a single mapping for all sites.\n"
            "  column  – compute mapping independently for each column.\n"
            "  both    – output both global and column-wise recodings."
        ),
    )
    parser.add_argument(
        "--add-outgroup",
        action="store_true",
        help="If set, add an 'Out' taxon with all zeros.",
    )

    args = parser.parse_args()

    alignment, ntax, nchar = parse_alignment(args.input)

    if args.mode == "global":
        mapping = compute_global_mapping(alignment)
        encoded = encode_alignment_global(alignment, mapping)
        nexus_text = format_nexus_matrix(
            encoded,
            add_outgroup=args.add_outgroup,
            include_header=True,
        )
        print(nexus_text)

    elif args.mode == "column":
        col_mappings = compute_column_mappings(alignment, nchar)
        encoded = encode_alignment_columnwise(alignment, col_mappings)
        nexus_text = format_nexus_matrix(
            encoded,
            add_outgroup=args.add_outgroup,
            include_header=True,
        )
        print(nexus_text)

    elif args.mode == "both":
        # Global recoding
        mapping = compute_global_mapping(alignment)
        encoded_global = encode_alignment_global(alignment, mapping)
        nexus_global = format_nexus_matrix(
            encoded_global,
            add_outgroup=args.add_outgroup,
            include_header=True,
        )

        # Column-wise recoding
        col_mappings = compute_column_mappings(alignment, nchar)
        encoded_col = encode_alignment_columnwise(alignment, col_mappings)
        # No second #NEXUS header
        nexus_col = format_nexus_matrix(
            encoded_col,
            add_outgroup=args.add_outgroup,
            include_header=False,
        )

        # Output both blocks, separated by NEXUS comments
        print("[FRECode GLOBAL recoding]")
        print(nexus_global)
        print("")
        print("[FRECode COLUMN-WISE recoding]")
        print(nexus_col)


if __name__ == "__main__":
    main()
