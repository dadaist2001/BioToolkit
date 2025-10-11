def is_dna(seq: str) -> bool:
    """
    A function for checking whether a sequence is a DNA sequence

    Args: seq (str): Input string to check.

    Returns: bool: True if the sequence contains only A, T, G, C.

    """
    return all(base in "aAtTgGcC" for base in seq)


def is_rna(seq: str) -> bool:
    """
    A function for checking whether a sequence is a RNA sequence

    Args: seq (str): Input string to check.

    Returns: bool: True if the sequence contains only A, U, G, C.

    """
    return all(base in "aAuUgGcC" for base in seq)


def transcribe(seq: str) -> str:
    """
    A function that transcribes DNA into RNA.

    Args: seq (str): Input string to be transcribed.

    Returns: Complementary DNA sequence of RNA.

    """
    if not is_dna(seq):
        raise ValueError(
            "This sequence cannot be transcribed because it is not DNA."
            "For retroviruses, use reverse_transcribe."
        )
    return seq.replace("T", "U").replace("t", "u")


def reverse_transcribe(seq: str) -> str:
    """
    A function that reverse-transcribes RNA into DNA.

    Args: seq (str): Input string to be reverse transcribed.

    Returns: str: DNA sequence where U are replaced by T.

    """
    if not is_rna(seq):
        raise ValueError(
            "This sequence cannot be transcribed because it is not RNA."
            "For DNA use transcribe"
        )
    return seq.replace("U", "T").replace("u", "t")


def reverse(seq: str) -> str:
    """
    A function that reverse sequence.

    Args: seq (str): Input string to be reversed.

    Returns: str: Sequence written in reverse order.

    """
    return seq[::-1]


def complement(seq: str) -> str:
    """
    Return the complementary DNA or RNA sequence.

    Args: seq (str): DNA or RNA sequence.

    Returns: str: Complementary sequence.

    """
    dna_table = {
        "A": "T", "T": "A", "G": "C", "C": "G",
        "a": "t", "t": "a", "g": "c", "c": "g",
    }
    rna_table = {
        "A": "U", "U": "A", "G": "C", "C": "G",
        "a": "u", "u": "a", "g": "c", "c": "g",
    }

    if is_dna(seq):
        return "".join(dna_table[base] for base in seq)
    elif is_rna(seq):
        return "".join(rna_table[base] for base in seq)
    else:
        raise ValueError("Sequence should be DNA or RNA.")


def reverse_complement(seq: str) -> str:
    """
    A function that reverse complement.

    Args: seq (str): Input string to be reversed.

    Returns: str: Complement sequence written in reverse order.

    """
    return reverse(complement(seq))
