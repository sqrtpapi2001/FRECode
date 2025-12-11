# FRECode  
*A frequency-based DNA alignment recoder for phylogenetic analysis*

FRECode implements a flexible nucleotide-recoding method in which each aligned
site (or the entire alignment) is converted into a 4-state numerical matrix
(0,1,2,3) based on relative nucleotide frequencies. This method reduces
compositional noise, highlights informative site patterns, and produces compact
matrices usable by major phylogenetic software.

FRECode supports three recoding modes:

1. GLOBAL mode  
   Rank A, C, G, T once for the entire alignment based on overall frequencies.

2. COLUMN-WISE mode  
   Rank nucleotides independently for each alignment column.

3. BOTH mode  
   Output both recodings within a single file: first global, then column-wise.

--------------------------------------------------------------------

## Features

- Simple plain-text alignment input.
- Global, column-wise, and dual recoding modes.
- NEXUS output compatible with PAUP*, MrBayes, RevBayes, Mesquite, TNT.
- Missing data (?) and gaps (-) preserved.
- Deterministic alphabetical tie-breaking for equal frequencies.
- Optional automatically generated all-zero outgroup.
- Pure Python (no dependencies).

--------------------------------------------------------------------

## Input Format

The input must follow this structure:

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

Rules:

- Taxon names must have no spaces.
- All sequences must be exactly NCHAR characters long.
- Allowed characters: A, C, G, T, ?, -

--------------------------------------------------------------------

## Installation

Clone the repository:

    git clone https://github.com/yourname/frecode.git
    cd frecode

Optional: make executable:

    chmod +x frecode.py

--------------------------------------------------------------------

## Usage

GLOBAL recoding:

    python frecode.py input.txt --mode global

COLUMN-WISE recoding:

    python frecode.py input.txt --mode column

Output BOTH recodings (two NEXUS blocks):

    python frecode.py input.txt --mode both

Add a synthetic all-zero outgroup:

    python frecode.py input.txt --mode global --add-outgroup

Redirect output:

    python frecode.py input.txt --mode column > recoded.nex

--------------------------------------------------------------------

## Recoding Algorithm

1. Frequency Counting  
   - In global mode: count A, C, G, T across the entire alignment.  
   - In column-wise mode: count within each site separately.  
   Missing data (?) and gaps (-) do not contribute to counts.

2. Ranking  
   Nucleotides are sorted by:
   - decreasing frequency  
   - alphabetical order for ties  

3. Weight Assignment  
   Frequencies map to weights in descending rank:
       1, 2, 4, 8

4. Code Conversion  
   Weights are transformed into final states:
       1 → 0  
       2 → 1  
       4 → 2  
       8 → 3

5. Missing Data  
   ? and - are preserved.  
   Unknown symbols become ?.

--------------------------------------------------------------------

## Output Format

FRECode outputs a valid NEXUS matrix.

A typical DATA block will look like:

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
        TYPESET * UNTITLED = ord: 1-12;
    END;

In BOTH mode, the output contains two DATA blocks in sequence:
first the GLOBAL recoding, then the COLUMN-WISE recoding.

--------------------------------------------------------------------

## Roadmap

Possible future enhancements:

- FASTA and PHYLIP input support
- Exporting frequency tables and mapping tables
- Handling IUPAC ambiguity codes (R, Y, S, W, K, M, N)
- Visualization tools for per-site composition
- Snakemake / Nextflow workflow integration
- JSON/CSV export for ML pipelines

--------------------------------------------------------------------

## Citation

If FRECode contributes to your work, you may cite:

    FRECode: A Frequency-Based Recoder for Four-State Phylogenetic Matrices
    <Your Name>, <Year>
    GitHub: https://github.com/yourname/frecode

--------------------------------------------------------------------

## License

MIT License – free for academic and commercial use.

--------------------------------------------------------------------

## Contributing

Pull requests and feature suggestions are welcome.
Open an Issue for bugs, ideas, or improvements.

