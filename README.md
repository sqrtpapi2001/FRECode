# FRECode  
*A frequency-based DNA alignment recoder for phylogenetic analysis*

FRECode implements a flexible nucleotide-recoding method in which each aligned
site (or the entire alignment) is transformed into a 4-state numerical matrix
(`0,1,2,3`) based on relative nucleotide frequencies.  

This produces a compact, model-agnostic representation suitable for downstream
phylogenetic inference or exploratory analyses.

FRECode supports two major modes:

1. **Global mode**  
   The four nucleotides A/C/G/T are ranked once using their **total frequencies
   over the entire alignment**; the same mapping is used for all sites.

2. **Column-wise mode**  
   Each alignment column independently ranks A/C/G/T using **per-site
   frequencies**, generating a site-specific mapping.

Both modes follow a consistent weight → code scheme designed to minimize a
simple cost function over a four-symbol encoding.

---

## Features

- Simple text alignment input (taxon name + sequence)
- NEXUS formatted output (STANDARD, symbols `0–3`)
- Two recoding strategies:
  - Global frequency mapping
  - Column-wise (site-specific) frequency mapping
- Handles missing data (`?`) and gaps (`-`) without alteration
- Deterministic even under ties (alphabetical tie-break)
- Optional synthetic outgroup (all-zero character vector)

---

## Input Format

FRECode expects a plain-text file with this structure:

```text
<NTAX> <NCHAR>
<TaxonName> <Sequence>
<TaxonName> <Sequence>
```

where

```text
NTAX = number of taxa (sequences)

NCHAR = alignment length (characters per sequence)

Taxon names must be single tokens (no spaces)

All sequences must be exactly NCHAR long

Valid sequence characters: A C G T ? -
```

## Example
```text
5 12
Tax1   ACGTACGTACGT
Tax2   A?GTACCTACGA
Tax3   TCGTACGTACGA
Tax4   ACGTACGTTCGT
Tax5   GCGTACGTACGA
```

Usage

Basic usage:

```python frecode.py input.txt --mode global```

```python frecode.py input.txt --mode column```


Redirect output to a file:

```python frecode.py input.txt --mode column > recoded.nex```


Add an all-zero synthetic outgroup (named Out by default):

```python frecode.py input.txt --mode global --add-outgroup > recoded_with_outgroup.nex```

Command-line options:

``input (positional)``` – path to the alignment .txt file

```--mode {global,column}``` – recoding mode (default: global)

```--add-outgroup``` – append an all-zero outgroup taxon to the matrix

Recoding Method
1. Frequency counting

Depending on the mode:

Global mode:
FRECode counts A/C/G/T across the entire alignment.

Column-wise mode:
FRECode counts A/C/G/T within each column independently.

Missing data (?) and gaps (-) are ignored in the counts.

2. Ranking and weights

For each recoding context (global or per-column):

Sort nucleotides by:

frequency (descending)

alphabetically as a tie-breaker

Assign weights:

most frequent     → 1
2nd most frequent → 2
3rd most frequent → 4
least frequent    → 8


Map weights to codes:

1 → 0
2 → 1
4 → 2
8 → 3


Example (for a single column):

Counts:  A=3, T=2, C=0, G=0
Order:   A, T, C, G
Weights: A=1, T=2, C=4, G=8
Codes:   A=0, T=1, C=2, G=3

3. Handling missing data

Any ? or - in the input is copied unchanged to the output.

Any unknown/other symbol is recoded as ?.

Output Format

FRECode writes a NEXUS file to standard output.

Example (column-wise mode):

#NEXUS
BEGIN DATA;
    DIMENSIONS NTAX=5 NCHAR=12;
    FORMAT DATATYPE=STANDARD MISSING=? GAP=- SYMBOLS="0123";
    MATRIX
    Tax1   203120312031
    Tax2   2?3120032130
    Tax3   303120312130
    Tax4   203120313031
    Tax5   103120312130
    ;
END;

BEGIN ASSUMPTIONS;
    TYPESET * UNTITLED   =  ord:  1- 12;
END;


This is directly usable in:

PAUP*

MrBayes

RevBayes

Mesquite

TNT (via simple conversion)

Other tools that accept 4-state STANDARD data

Roadmap

Planned or easy-to-add features:

FASTA / PHYLIP input formats

Handling of ambiguous IUPAC bases (R, Y, S, W, K, M, N)

Optional export of:

global mapping table (A/C/G/T → 0/1/2/3)

per-column mapping tables

Alignment diagnostics:

base frequency plots

GC content by site or taxon

Integration into workflows (Snakemake, Nextflow)

##Citation

If FRECode helps your research, please consider citing:

FRECode: A Frequency-Based Nucleotide Recoder for Four-State Phylogenetic Matrices
<Your Name>, <Year>, GitHub repository: https://github.com/yourname/frecode


(You can add a DOI via Zenodo for formal citation.)

##License

FRECode is released under the MIT License.
You are free to use, modify, and distribute it in academic or commercial projects.

##Contributing

Contributions, bug reports, and feature requests are welcome.

Open an issue for questions or ideas

Submit a pull request for code changes

Suggest new recoding strategies or input formats

Happy recoding! 🧬
