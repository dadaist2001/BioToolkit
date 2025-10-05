def is_nucleic_acid(seq: str) -> bool:
    """
    A function that checks whether DNA and RNA have mixed.

    Args: seq (str): Input string to check.

    Returns: bool: True if the sequence is DNA or RNA.

    """
    return is_dna(seq) or is_rna(seq)

def gc_content(seq: str) -> float:
    """
    A function for calculating the GC content of a nucleic acid sequence.

    Args: seq (str): DNA or RNA sequence.

    Returns: float: GC content in percentage.
    """
    if not seq:
        return 0.0
    gc_count = sum(1 for base in seq.upper() if base in "GC")
    return (gc_count / len(seq)) * 100

def mean_quality(quality: str) -> float:
    """
    A function for calculating the average quality of a FASTQ read.

    Args: qual (str): String of ASCII characters to Phred33 quality.

    Returns: float: Average quality score.
    """
    if not quality:
        return 0.0
    score = [ord(c) - 33 for c in quality]  
    return sum(score) / len(score)
