def convert_multiline_fasta_to_oneline(input_fasta, output_fasta=None):
    """
    A function that converts a FASTA file where sequences
    are split across multiple lines
    into a format where each sequence is on a single line.

    Args:
        input_fasta (str): path to the input file.
        output_fasta (str, optional): path for the converted file.

    Returns:
        str: path to the converted FASTA file.
    """

    if output_fasta is None:
        output_fasta = "converted_" + input_fasta.split("/")[-1]

    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:

        header = None
        seq = ""

        for line in infile:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    outfile.write(f"{header}\n{seq}\n")

                header = line
                seq = ""
            else:
                seq += line

        if header is not None:
            outfile.write(f"{header}\n{seq}\n")


def parse_blast_output(input_file, output_file):
    """
    A function for extracting the first name for QUERY
    from a BLAST report and saves all found names.

    Args:
        input_file (str): path to the BLAST file.
        output_file (str): path to the file will be saved.

    Returns:
        None
    """

    results = []

    with open(input_file, "r") as f:
        lines = f.readlines()

    for i in range(len(lines)):
        line = lines[i].strip()

        if "Sequences producing significant alignments" in line:
            if i + 1 < len(lines):
                next_line = lines[i + 1].strip()
                if next_line:
                    results.append(next_line)
    results.sort()

    with open(output_file, "w") as out:
        for res in results:
            out.write(res + "\n")


def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes: "str | list[str]",
    n_before: int = 1,
    n_after: int = 1,
    output_fasta: str = "selected_genes.fasta",
) -> None:
    """
    The function for extracting protein sequences
    for genes near specified genes of interest.
    """

    with open(input_gbk, "r") as f:
        all_lines = f.readlines()

    all_genes = []
    gene = ""
    translation = ""

    for line in all_lines:
        line = line.rstrip()
        if line.startswith("/gene="):
            gene = line.split('"')[1]
        elif line.startswith("/translation="):
            translation = line.split('"')[1]
            all_genes.append({"gene_name": gene, "protein": translation})
            gene = ""
            translation = ""

    if isinstance(genes, str):
        genes_to_find = [genes]
    else:
        genes_to_find = genes

    genes_to_write = []
    index = 0
    while index < len(all_genes):
        g = all_genes[index]
        if g["gene_name"] in genes_to_find:
            start_index = max(0, index - n_before)
            end_index = min(len(all_genes), index + n_after + 1)
            for k in range(start_index, end_index):
                if k != index:
                    genes_to_write.append(all_genes[k])
        index += 1

    with open(output_fasta, "w") as out_file:
        for item in genes_to_write:
            out_file.write(f">{item['gene_name']}\n{item['protein']}\n")
