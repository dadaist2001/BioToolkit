def read_fastq(input_fastq: str):
    """
    A function that reads a FASTQ file sequence by sequence.

    Yields: tuple(str, str, str): (sequence_name, sequence, quality)
    """
    with open(input_fastq, "r") as f:
        while True:
            name = f.readline().rstrip()
            if not name:
                break
            seq = f.readline().rstrip()
            f.readline()  # plus line
            qual = f.readline().rstrip()
            yield name[1:], seq, qual  # remove '@' from name

def write_fastq(sequence_data: tuple, output_fastq: str, overwrite=False):
    """
    A function that writes a single sequence to the output FASTQ file.

    Args:
    sequence_data: tuple(str, str, str) — (name, sequence, quality)
    output_fastq: str — path to the output file
    overwrite: bool — if True, clears the file before writing
    """
    mode = "w" if overwrite else "a"
    name, seq, qual = sequence_data
    with open(output_fastq, mode) as f:
        f.write(f"@{name}\n{seq}\n+\n{qual}\n")
