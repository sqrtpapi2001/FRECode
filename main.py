#!/usr/bin/env python3
"""
FRECode: Frequency-based Recoder for DNA Alignments
Now supports:
  - PHYLIP (.phy) input (sequential or interleaved; wrapped lines supported)
  - Simple TXT input (ntax nchar + taxon seq)
  - Recoding modes: global | column | both
  - Output: NEXUS
  - Optional ZIP packaging of output files

Examples
--------
# Column-wise recoding from PHYLIP, zip output:
python frecode.py Test_1.phy --mode column --zip FRECode_output.zip

# Global recoding from PHYLIP:
python frecode.py Test_1.phy --mode global > out.nex

# Both recodings, write to files and zip them:
python frecode.py Test_1.phy --mode both --out-prefix Test_1 --zip Test_1_frecode.zip

# Simple TXT format:
python frecode.py alignment.txt --mode column --zip out.zip
"""

from __future__ import annotations

import argparse
import zipfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional


NUCS = ("A", "C", "G", "T")
WEIGHTS = [1, 2, 4, 8]
WEIGHT_TO_CODE = {1: 0, 2: 1, 4: 2, 8: 3}


@dataclass
class Alignment:
    taxa: List[str]              # ordered taxa
    seqs: Dict[str, str]         # taxon -> sequence
    nchar: int


def _clean_seq(s: str) -> str:
    return "".join(ch for ch in s.replace(" ", "").replace("\t", "") if ch)


def read_txt_alignment(path: Path) -> Alignment:
    """
    Read the simple TXT format:

        <NTAX> <NCHAR>
        <TaxonName> <Sequence>
        ...

    - First non-empty, non-comment line gives ntax and nchar.
    - Next ntax non-empty, non-comment lines provide taxon and sequence.
      First token = taxon, last token = sequence.
    """
    lines: List[str] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            lines.append(line.rstrip("\n"))

    if not lines:
        raise ValueError("Empty alignment file.")

    ntax, nchar = map(int, lines[0].split()[:2])
    taxa: List[str] = []
    seqs: Dict[str, str] = {}

    for line in lines[1:]:
        if len(taxa) >= ntax:
            break
        parts = line.split()
        if len(parts) < 2:
            continue
        taxon = parts[0]
        seq = parts[-1].strip()
        seqs[taxon] = seq
        taxa.append(taxon)

    if len(taxa) != ntax:
        raise ValueError(f"Declared NTAX={ntax}, but read {len(taxa)} sequences.")

    for t in taxa:
        if len(seqs[t]) != nchar:
            raise ValueError(f"Sequence length mismatch for {t}: expected {nchar}, got {len(seqs[t])}.")

    return Alignment(taxa=taxa, seqs=seqs, nchar=nchar)


def read_phylip(path: Path) -> Alignment:
    """
    Robust PHYLIP reader:
      - Handles sequential and interleaved PHYLIP
      - Handles wrapped lines
      - Assumes taxon name is the first field (often fixed-width 10 chars, but we accept relaxed)
      - For interleaved: after first block, subsequent blocks may omit taxon names

    Strategy:
      1) Read header NTAX NCHAR
      2) Parse first block of up to NTAX taxa lines with names and initial sequence chunks
      3) Then keep reading sequence chunks (interleaved or sequential-wrapped) until each taxon reaches NCHAR
    """
    raw_lines: List[str] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                raw_lines.append(line.rstrip("\n"))

    if not raw_lines:
        raise ValueError("Empty PHYLIP file.")

    header = raw_lines[0].split()
    if len(header) < 2:
        raise ValueError("PHYLIP header must begin with: <NTAX> <NCHAR>")

    ntax, nchar = int(header[0]), int(header[1])
    body = raw_lines[1:]

    taxa: List[str] = []
    seqs: Dict[str, str] = {}

    i = 0

    # ---- Parse first block: expect NTAX taxa lines with names + sequence chunk
    while i < len(body) and len(taxa) < ntax:
        line = body[i]
        i += 1
        if not line.strip():
            continue

        # Common PHYLIP: name is first 10 chars; relaxed PHYLIP: name is first token
        # We try fixed-width first; if that yields empty, fall back to token.
        name_fw = line[:10].strip()
        rest_fw = line[10:]
        if name_fw:
            name = name_fw
            chunk = _clean_seq(rest_fw)
        else:
            parts = line.split()
            if len(parts) < 2:
                continue
            name = parts[0]
            chunk = _clean_seq("".join(parts[1:]))

        taxa.append(name)
        seqs[name] = chunk

    if len(taxa) != ntax:
        raise ValueError(f"PHYLIP header says NTAX={ntax}, but could not parse {ntax} taxon lines.")

    # ---- Now collect remaining sequence chunks until each taxon reaches nchar
    def done() -> bool:
        return all(len(seqs[t]) >= nchar for t in taxa)

    # Decide whether remaining looks interleaved (blocks of NTAX lines) or sequential-wrapped.
    # We'll handle either by two passes:
    # - If a line begins with a known taxon name (or fixed-width name), treat as named line.
    # - Else, treat as unnamed chunk line for taxa in order (interleaved blocks).

    while i < len(body) and not done():
        # Skip blank separators
        if not body[i].strip():
            i += 1
            continue

        # Peek next non-empty line
        line = body[i]

        # Try interpret as named line (interleaved blocks that repeat taxon names)
        name_fw = line[:10].strip()
        parts = line.split()
        named = False
        name = None
        chunk = ""

        if name_fw in seqs:
            named = True
            name = name_fw
            chunk = _clean_seq(line[10:])
        elif parts and parts[0] in seqs:
            named = True
            name = parts[0]
            chunk = _clean_seq("".join(parts[1:]))

        if named:
            # Consume named lines as they appear
            i += 1
            seqs[name] += chunk
            continue

        # Otherwise, assume an unnamed interleaved block: next NTAX lines are chunks (possibly with whitespace)
        for t in taxa:
            # advance to next non-empty
            while i < len(body) and not body[i].strip():
                i += 1
            if i >= len(body):
                break
            chunk_line = body[i]
            i += 1
            seqs[t] += _clean_seq(chunk_line)

    # Trim sequences to nchar and validate
    for t in taxa:
        seqs[t] = seqs[t][:nchar]
        if len(seqs[t]) != nchar:
            raise ValueError(
                f"PHYLIP parse incomplete for {t}: expected {nchar} characters, got {len(seqs[t])}."
            )

    return Alignment(taxa=taxa, seqs=seqs, nchar=nchar)


def read_alignment_auto(path: Path) -> Alignment:
    ext = path.suffix.lower()
    if ext in (".phy", ".phylip"):
        return read_phylip(path)
    # fallback: try txt format
    return read_txt_alignment(path)


def compute_global_mapping(aln: Alignment) -> Dict[str, int]:
    total = Counter()
    for seq in aln.seqs.values():
        for ch in seq:
            u = ch.upper()
            if u in NUCS:
                total[u] += 1
    for n in NUCS:
        total.setdefault(n, 0)
    order = sorted(NUCS, key=lambda n: (-total[n], n))
    mapping: Dict[str, int] = {}
    for nuc, w in zip(order, WEIGHTS):
        mapping[nuc] = WEIGHT_TO_CODE[w]
    return mapping


def compute_column_mappings(aln: Alignment) -> List[Dict[str, int]]:
    col_counts = [Counter() for _ in range(aln.nchar)]
    for seq in aln.seqs.values():
        for j, ch in enumerate(seq):
            u = ch.upper()
            if u in NUCS:
                col_counts[j][u] += 1

    col_maps: List[Dict[str, int]] = []
    for c in col_counts:
        for n in NUCS:
            c.setdefault(n, 0)
        order = sorted(NUCS, key=lambda n: (-c[n], n))
        mapping: Dict[str, int] = {}
        for nuc, w in zip(order, WEIGHTS):
            mapping[nuc] = WEIGHT_TO_CODE[w]
        col_maps.append(mapping)
    return col_maps


def encode_global(aln: Alignment, mapping: Dict[str, int]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for t in aln.taxa:
        seq = aln.seqs[t]
        enc = []
        for ch in seq:
            u = ch.upper()
            if u in mapping:
                enc.append(str(mapping[u]))
            elif ch in ("?", "-"):
                enc.append(ch)
            else:
                enc.append("?")
        out[t] = "".join(enc)
    return out


def encode_columnwise(aln: Alignment, col_maps: List[Dict[str, int]]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for t in aln.taxa:
        seq = aln.seqs[t]
        enc = []
        for j, ch in enumerate(seq):
            u = ch.upper()
            if u in col_maps[j]:
                enc.append(str(col_maps[j][u]))
            elif ch in ("?", "-"):
                enc.append(ch)
            else:
                enc.append("?")
        out[t] = "".join(enc)
    return out


def format_nexus(
    taxa: List[str],
    encoded: Dict[str, str],
    nchar: int,
    title_comment: Optional[str] = None,
    include_header: bool = True,
    add_outgroup: bool = False,
    outgroup_name: str = "Out",
) -> str:
    ntax = len(taxa) + (1 if add_outgroup else 0)
    pad = max([len(t) for t in taxa] + ([len(outgroup_name)] if add_outgroup else [0])) + 2

    lines: List[str] = []
    if include_header:
        lines.append("#NEXUS")
    if title_comment:
        lines.append(f"[{title_comment}]")

    lines.append("BEGIN DATA;")
    lines.append(f"    DIMENSIONS NTAX={ntax} NCHAR={nchar};")
    lines.append('    FORMAT DATATYPE=STANDARD MISSING=? GAP=- SYMBOLS="0123";')
    lines.append("    MATRIX")

    for t in taxa:
        lines.append(f"    {t.ljust(pad)}{encoded[t]}")
    if add_outgroup:
        lines.append(f"    {outgroup_name.ljust(pad)}" + ("0" * nchar))

    lines.append("    ;")
    lines.append("END;")
    lines.append("")
    lines.append("BEGIN ASSUMPTIONS;")
    lines.append(f"    TYPESET * UNTITLED = ord: 1-{nchar};")
    lines.append("END;")

    return "\n".join(lines)


def write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def main() -> None:
    p = argparse.ArgumentParser(description="FRECode: frequency-based recoder (PHYLIP/TXT -> NEXUS).")
    p.add_argument("input", help="Input alignment (.phy/.phylip or simple txt format).")
    p.add_argument("--mode", choices=["global", "column", "both"], default="column",
                   help="Recoding mode: global | column | both (default: column).")
    p.add_argument("--add-outgroup", action="store_true", help="Append an all-zero outgroup named 'Out'.")
    p.add_argument("--out-prefix", default=None,
                   help="Output file prefix (default: derived from input filename).")
    p.add_argument("--zip", dest="zip_path", default=None,
                   help="If set, write outputs to files and package into this ZIP.")
    args = p.parse_args()

    in_path = Path(args.input)
    if not in_path.exists():
        raise SystemExit(f"Input file not found: {in_path}")

    aln = read_alignment_auto(in_path)

    prefix = args.out_prefix or in_path.stem
    out_files: List[Path] = []

    def emit(name_suffix: str, nexus_text: str) -> None:
        nonlocal out_files
        if args.zip_path:
            out_path = Path(f"{prefix}_{name_suffix}.nex")
            # write next to current working dir
            write_text(out_path, nexus_text)
            out_files.append(out_path)
        else:
            # stdout
            print(nexus_text)

    if args.mode == "global":
        gmap = compute_global_mapping(aln)
        enc = encode_global(aln, gmap)
        nexus = format_nexus(
            taxa=aln.taxa,
            encoded=enc,
            nchar=aln.nchar,
            title_comment="FRECode GLOBAL recoding",
            include_header=True,
            add_outgroup=args.add_outgroup,
        )
        emit("global", nexus)

    elif args.mode == "column":
        cmaps = compute_column_mappings(aln)
        enc = encode_columnwise(aln, cmaps)
        nexus = format_nexus(
            taxa=aln.taxa,
            encoded=enc,
            nchar=aln.nchar,
            title_comment="FRECode COLUMN-WISE recoding",
            include_header=True,
            add_outgroup=args.add_outgroup,
        )
        emit("column", nexus)

    else:  # both
        gmap = compute_global_mapping(aln)
        enc_g = encode_global(aln, gmap)
        nexus_g = format_nexus(
            taxa=aln.taxa,
            encoded=enc_g,
            nchar=aln.nchar,
            title_comment="FRECode GLOBAL recoding",
            include_header=True,
            add_outgroup=args.add_outgroup,
        )

        cmaps = compute_column_mappings(aln)
        enc_c = encode_columnwise(aln, cmaps)
        nexus_c = format_nexus(
            taxa=aln.taxa,
            encoded=enc_c,
            nchar=aln.nchar,
            title_comment="FRECode COLUMN-WISE recoding",
            include_header=False,  # no second #NEXUS
            add_outgroup=args.add_outgroup,
        )

        if args.zip_path:
            emit("global", nexus_g)
            emit("column", nexus_c)
        else:
            # stdout: print both blocks in one stream
            print(nexus_g)
            print("")
            print(nexus_c)

    # ZIP packaging (if requested)
    if args.zip_path:
        zip_path = Path(args.zip_path)
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
            for f in out_files:
                z.write(f, arcname=f.name)

        # Optional: clean up loose files after zipping
        # Comment out the next two lines if you want to keep the .nex files.
        for f in out_files:
            try:
                f.unlink()
            except OSError:
                pass

        print(f"Wrote ZIP: {zip_path.resolve()}")


if __name__ == "__main__":
    main()
