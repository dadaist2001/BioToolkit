from tools.utils import gc_content, mean_quality

from tools.rna_dna_tools import (
    is_dna,
    is_rna,
    transcribe,
    reverse_transcribe,
    reverse,
    complement,
    reverse_complement,
)
from tools.read_tools import read_fastq, write_fastq


def run_dna_rna_tools(*args):
    """
    The main function.

    Args: The last argument is the procedure name.

     Returns: All previous arguments are sequences.
    """
    if len(args) < 2:
        raise ValueError("You need to pass the sequence and procedure.")

    *sequences, procedure = args
    results = []

    for seq in sequences:
        if procedure == "is_rna":
            results.append(is_rna(seq))
        elif procedure == "is_dna":
            results.append(is_dna(seq))
        elif procedure == "transcribe":
            results.append(transcribe(seq))
        elif procedure == "reverse_transcribe":
            results.append(reverse_transcribe(seq))
        elif procedure == "reverse":
            results.append(reverse(seq))
        elif procedure == "complement":
            results.append(complement(seq))
        elif procedure == "reverse_complement":
            results.append(reverse_complement(seq))
        else:
            raise ValueError(f"There is no such {procedure} in the function")

    return results[0] if len(results) == 1 else results


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold=0,
) -> None:
    """
    Function for Filtering sequences from a FASTQ file
    and writes only those passing
    GC content, length, and quality thresholds to an output FASTQ file.

    Args:
        input_fastq (str): path to the input file.
        output_fastq (str): path to save the filtered sequences.
        gc_bounds (tuple or float): GC content range.
        length_bounds (tuple or int): length range.
        quality_threshold (int): min average quality.

    Returns:
        None
    """

    if isinstance(gc_bounds, (int, float)):
        gc_min, gc_max = 0, float(gc_bounds)
    else:
        gc_min, gc_max = gc_bounds

    if isinstance(length_bounds, (int, float)):
        len_min, len_max = 0, int(length_bounds)
    else:
        len_min, len_max = length_bounds

    mode = "w" if overwrite else "a"

    with open(output_fastq, mode) as out:
        for name, seq, qual in read_fastq(input_fastq):
            if (
                gc_min <= gc_content(seq) <= gc_max
                and len_min <= len(seq) <= len_max
                and mean_quality(qual) >= quality_threshold
            ):
                out.write(f"@{name}\n{seq}\n+\n{qual}\n")
