import pytest
from BioToolkit import filter_fastq
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os


# helper - without thшs thing, the code was falling
def make_record(seq, id="1"):
    return SeqRecord(
        Seq(seq),
        id=id,
        letter_annotations={"phred_quality": [40] * len(seq)}
   )


# create file
def test_file_created(tmp_path):
    inp = tmp_path / "in.fastq"
    out = tmp_path / "out.fastq"

    SeqIO.write([make_record("ATGC")], inp, "fastq")

    filter_fastq(str(inp), str(out))

    assert os.path.exists(out)


# empty file
def test_empty(tmp_path):
    inp = tmp_path / "in.fastq"
    out = tmp_path / "out.fastq"

    SeqIO.write([], inp, "fastq")

    filter_fastq(str(inp), str(out))

    assert len(list(SeqIO.parse(out, "fastq"))) == 0


# file error
def test_wrong():
    with pytest.raises(Exception):
        filter_fastq("no.fastq", "out.fastq")


# one record
def test_record(tmp_path):
    inp = tmp_path / "in.fastq"
    out = tmp_path / "out.fastq"

    SeqIO.write([make_record("ATGC")], inp, "fastq")

    filter_fastq(str(inp), str(out))

    assert len(list(SeqIO.parse(out, "fastq"))) == 1


# empty end
def test_not_empty(tmp_path):
    inp = tmp_path / "a.fastq"
    out = tmp_path / "b.fastq"

    SeqIO.write([make_record("ATGC", "x")], inp, "fastq")

    filter_fastq(str(inp), str(out))

    result = list(SeqIO.parse(out, "fastq"))

    assert result != []


# filtering GC
def test_gc(tmp_path):
    inp = tmp_path / "a.fastq"
    out = tmp_path / "b.fastq"

    SeqIO.write([make_record("GGGG")], inp, "fastq")

    filter_fastq(str(inp), str(out), gc_bounds=(0.75, 1))

    assert len(list(SeqIO.parse(out, "fastq"))) == 1


# filtering by len
def test_length(tmp_path):
    inp = tmp_path / "a.fastq"
    out = tmp_path / "b.fastq"

    SeqIO.write([
        make_record("A"),
        make_record("AAAAA")
    ], inp, "fastq")

    filter_fastq(str(inp), str(out), length_bounds=(5, 10))

    result = list(SeqIO.parse(out, "fastq"))

    assert len(result) == 1