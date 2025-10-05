def is_rna(seq: str) -> bool:
    """
    A function for checking whether a sequence is a RNA sequence

    Args: seq (str): Input string to check.

    Returns: bool: True if the sequence contains only A, U, G, C.

    """
    return all(base in "aAuUgGcC" for base in seq)

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
    Return the complementary RNA sequence.

    Args: seq (str): RNA sequence.

    Returns: str: Complementary sequence.

    """
    rna = {
        "A": "U",
        "U": "A",
        "G": "C",
        "C": "G",
        "a": "u",
        "u": "a",
        "g": "c",
        "c": "g",
    }

    if is_rna(seq):
        return "".join(rna[base] for base in seq)
    else:
        raise ValueError("Sequence should be RNA.")


def reverse_complement(seq: str) -> str:
    """
    A function that reverse complement.

    Args: seq (str): Input string to be reversed.

    Returns: str: Complement sequence written in reverse order.

    """
    if not is_rna(seq):
        raise ValueError("Sequence should b DNA or RNA.")

    return reverse(complement(seq))
