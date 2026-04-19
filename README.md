# BioToolkit

This is a small bioinformatics toolkit for working with DNA and RNA sequences by dadaist2001.

## What it can do

The toolkit now uses an object-oriented approach for working with biological sequences:

- **DNASequence**
  - transcription to RNA
  - reverse, complement, reverse complement

- **RNASequence**
  - reverse, complement, reverse complement

- **AminoAcidSequence**
  - validation of protein sequences
  - molecular weight calculation

- **Bioinformatics file utilities**:
  - `convert_multiline_fasta_to_oneline(input_fasta, output_fasta=None)` – converts multi-line FASTA sequences to single-line format.
  - `parse_blast_output(input_file, output_file)` – extracts top hits from BLAST reports and saves them sorted.
  - `select_genes_from_gbk_to_fasta` – extracts neighboring genes from GenBank files into FASTA format (monster-function).

- **FASTQ filtering (Biopython)**
FASTQ filtering is implemented using the Biopython library:

  - filtering by sequence length
  - filtering by GC content

## How to use it

Clone the repository:

```bash
git clone https://github.com/dadaist2001/BioToolkit.git
cd BioToolkit

