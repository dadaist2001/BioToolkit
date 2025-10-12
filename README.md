# BioToolkit

This is a small bioinformatics toolkit for working with DNA and RNA sequences be dadaist2001.

## What it can do

- DNA tools: check if a sequence is DNA, transcribe DNA to RNA, reverse, complement, reverse complement.
- RNA tools: check if a sequence is RNA, reverse transcribe RNA to DNA, reverse, complement, reverse complement.
- Filter FASTQ sequences by GC content, length, or quality..
- Bioinformatics file utilities (HW5):
  - `convert_multiline_fasta_to_oneline(input_fasta, output_fasta=None)` – converts multi-line FASTA sequences to single-line format.
  - `parse_blast_output(input_file, output_file)` – extracts top hits from BLAST reports and saves them sorted.
  - `select_genes_from_gbk_to_fasta` – extracts neighboring genes from GenBank files into FASTA format (monster-function).

## How to use it

Clone the repository:

```bash
git clone https://github.com/YOUR_USERNAME/BioToolkit.git
cd BioToolkit
```

## For Python

from BioToolkit import run_dna_rna_tools
from BioToolkit import filter_fastq
from BioToolkit import convert_multiline_fasta_to_oneline
from BioToolkit import parse_blast_output
from BioToolkit import select_genes_from_gbk_to_fasta