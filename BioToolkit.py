from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import logging

# simple logger
logging.basicConfig(filename="tool.log", level=logging.INFO)
logger = logging.getLogger(__name__)

class BiologicalSequence(ABC):
    """
    Abstract base class for biological sequences.
    """

    def __init__(self, sequence: str):
        """
        Initialize biological sequence.
        """
        self.sequence = sequence.upper()
        self._check_alphabet()

    def __len__(self):
        """
        Return length of the sequence.
        """
        return len(self.sequence)

    def __getitem__(self, item):
        """
        Return element or slice of the sequence.
        """
        result = self.sequence[item]
        if isinstance(item, slice):
            return self.__class__(result)
        return result

    def __str__(self):
        """
        Return human readable string representation.
        """
        return f"{self.__class__.__name__}: {self.sequence}"

    @abstractmethod
    def _check_alphabet(self):
        """
        Check that the sequence contains only valid characters.
        """
        pass


class NucleicAcidSequence(BiologicalSequence):
    """
    Abstract class for nucleic acid sequences.
    """

    alphabet = None
    complement_map = None

    def __init__(self, sequence: str):
        if self.__class__ == NucleicAcidSequence:
            raise NotImplementedError("Cannot instantiate abstract class")
        super().__init__(sequence)

    def _check_alphabet(self):
        """
        Validate nucleic acid alphabet.
        """
        if self.alphabet is None:
            raise NotImplementedError("Alphabet not defined")

        if not all(base in self.alphabet for base in self.sequence):
            raise ValueError("Invalid alphabet for this sequence type")

    def complement(self):
        """
        Return complementary sequence.
        """
        return self.__class__(
            "".join(self.complement_map[b] for b in self.sequence)
        )

    def reverse(self):
        """
        Return reversed sequence.
        """
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self):
        """
        Return reverse complementary sequence.
        """
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    """
    Class representing DNA sequence.
    """

    alphabet = {"A", "T", "G", "C"}
    complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}

    def transcribe(self):
        """
        Transcribe DNA into RNA.

        Returns:
            RNASequence: transcribed RNA sequence
        """
        return RNASequence(self.sequence.replace("T", "U"))


class RNASequence(NucleicAcidSequence):
    """
    Class representing RNA sequence.
    """

    alphabet = {"A", "U", "G", "C"}
    complement_map = {"A": "U", "U": "A", "G": "C", "C": "G"}


class AminoAcidSequence(BiologicalSequence):
    """
    Class representing amino acid (protein) sequence.
    """

    alphabet = set("ACDEFGHIKLMNPQRSTVWY")

    weights = {
        "A": 89, "C": 121, "D": 133, "E": 147,
        "F": 165, "G": 75,  "H": 155, "I": 131,
        "K": 146, "L": 131, "M": 149, "N": 132,
        "P": 115, "Q": 146, "R": 174, "S": 105,
        "T": 119, "V": 117, "W": 204, "Y": 181,
    }

    def _check_alphabet(self):
        """
        Validate amino acid alphabet.
        """
        if not all(aa in self.alphabet for aa in self.sequence):
            raise ValueError("Invalid amino acid sequence")

    def molecular_weight(self):
        """
        Calculate molecular weight of the protein sequence.
        """
        return sum(self.weights[aa] for aa in self.sequence)


def filter_fastq(
    input_fastq,
    output_fastq,
    length_bounds=(0, float("inf")),
    gc_bounds=(0, 100),
):
    """
    Filter FastQ file by sequence length and GC content.
    """

    logger.info("Start filtering...")

    try:
        records = SeqIO.parse(input_fastq, "fastq")
    except Exception:
        logger.error("Cannot read file")
        raise

    filtered = []

    min_len, max_len = length_bounds
    min_gc, max_gc = gc_bounds

    for record in records:
        seq_len = len(record)
        gc = gc_fraction(record.seq)

        if min_len <= seq_len <= max_len and min_gc <= gc <= max_gc:
            filtered.append(record)

    SeqIO.write(filtered, output_fastq, "fastq")

    logger.info(f"Done filtering!: {output_fastq}")

    return len(filtered)