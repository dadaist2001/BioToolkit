from tools.utils import gc_content, mean_quality

from tools.rna_dna_tools import (
    is_dna, is_rna, transcribe, reverse_transcribe,
    reverse, complement, reverse_complement
)

def run_dna_rna_tools(*args):
    """
    The main function.

    Args: The last argument is the procedure name (for example, 'transcribe').

     Returns: All previous arguments are sequences.
    """
    if len(args) < 2:
        raise ValueError("You need to pass the sequence and procedure.")

    *sequences, procedure = args
    results = []

    for seq in sequences:
        if procedure == "is_nucleic_acid":
            results.append(is_nucleic_acid(seq))
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
            raise ValueError(f"There is no such procedure in the function: {procedure}")

    return results[0] if len(results) == 1 else results

def filter_fastq(seqs: dict, 
                 gc_bounds=(0, 100), 
                 length_bounds=(0, 2**32), 
                 quality_threshold=0) -> dict:
    """
    A function for filtering FASTQ sequences.
    
    Args: seqs (dict): Dictionary with sequences and their quality.
    gc_bounds: GC content range (in percentage).
    length_bounds: Range of read lengths.
    quality_threshold: Minimum average quality value.

    Returns: dict: Only sequences that passed two filters.
    """

    filtered = {}

    if isinstance(gc_bounds, (int, float)):
        gc_min, gc_max = 0, float(gc_bounds)
    else:
        gc_min, gc_max = gc_bounds
        
    if isinstance(length_bounds, (int, float)):
        len_min, len_max = 0, int(length_bounds)
    else:
        len_min, len_max = length_bounds

    for name, (seq, qual) in seqs.items():
        seq_gc = gc_content(seq)
        seq_len = len(seq)
        seq_qual = mean_quality(qual)

        if (gc_min <= seq_gc <= gc_max and
            len_min <= seq_len <= len_max and
            seq_qual >= quality_threshold):
            filtered[name] = (seq, qual)

    return filtered
