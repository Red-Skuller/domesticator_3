import copy
import os
from typing import Generator, List, Dict
import warnings
from pathlib import Path

from Bio import SeqRecord, SeqIO, BiopythonParserWarning
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature

from dnachisel import reverse_translate


def load_vector_record(vector_filepath: Path) -> SeqRecord.SeqRecord:
    """Return Biopython SeqRecord.SeqRecords.

    Essentially a wrapper around Biopython SeqIO machinery to read a single vector file.

    Args:
        vector_filepath: Path to the vector file.

    Returns:
        A SeqRecord object representing the vector.
    """
    records = list(SeqIO.parse(vector_filepath, "genbank"))

    if len(records) != 1:
        raise RuntimeError(
            "correctly formatted vector files only have one record -- the sequence of the vector"
        )
    record = records[0]

    return record

def load_gb_as_naive(filepaths: List[str]) -> List[SeqRecord.SeqRecord]:
    """Load GenBank files as naive records.

    Args:
        filepaths: List of paths to GenBank files.

    Returns:
        List of SeqRecord objects.
    """
    records = []
    for filepath in filepaths:
        records.append(SeqIO.read(filepath, "genbank"))
    return records


def load_insert_seq(sequence: str) -> SeqRecord.SeqRecord:
    """Returns a Biopython SeqRecord holding DNA encoding the protein sequence.

    Args:
        sequence: Amino acid sequence.

    Returns:
        A SeqRecord object representing the insert DNA.
    """

    new_dna_seq = Seq(reverse_translate(sequence))
    assert sequence == new_dna_seq.translate()
    new_record = SeqRecord.SeqRecord(
        seq=new_dna_seq,
        id="insert",
        name="insert",
        description="synthetic gene",
        annotations={"molecule_type": "DNA", "chain": "A"},
    )
    return new_record


def load_inserts(
    filenames: List[str], increasing_chain_fasta: bool
) -> Generator[List[SeqRecord.SeqRecord], None, None]:
    """Yields lists of Biopython SeqRecords.

    A generator function which returns lists of SeqRecords holding DNA
    encoding the protein sequences taken from the input files.

    Args:
        filenames: List of file paths to load.
        increasing_chain_fasta: Whether to use increasing chain letters for FASTA records.

    Yields:
        Lists of SeqRecord objects.
    """
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for filename in filenames:
        assert os.path.isfile(filename)
        basename = os.path.basename(filename)
        name, ext = os.path.splitext(basename)
        assert ext in [".pdb", ".fasta"]

        if ext == ".fasta":
            records = []
            for chain_num, record in enumerate(SeqIO.parse(filename, "fasta")):
                if not increasing_chain_fasta:
                    chain_num = 0
                orig_aa_seq = record.seq.upper()
                new_dna_seq = Seq(reverse_translate(record.seq.upper()))
                assert orig_aa_seq == new_dna_seq.translate()
                new_record = SeqRecord.SeqRecord(
                    seq=new_dna_seq,
                    id=record.id,
                    name=record.id,
                    description=record.description,
                    annotations={"molecule_type": "DNA", "chain": alphabet[chain_num]},
                )
                records.append(new_record)
            yield records
        else:
            records = []
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=PDBConstructionWarning)
                warnings.filterwarnings("ignore", category=BiopythonParserWarning)
                for record in SeqIO.parse(filename, "pdb-atom"):
                    orig_aa_seq = record.seq.upper()
                    new_dna_seq = Seq(reverse_translate(record.seq.upper()))
                    assert orig_aa_seq == new_dna_seq.translate()
                    new_id = name + "_" + record.annotations["chain"]
                    new_name = name + "_" + record.annotations["chain"]
                    records.append(
                        SeqRecord.SeqRecord(
                            seq=new_dna_seq,
                            id=new_id,
                            name=new_name,
                            description=record.description,
                            annotations={
                                "molecule_type": "DNA",
                                "chain": record.annotations["chain"],
                            },
                        )
                    )
                yield records


def get_insert_locations(record: SeqRecord.SeqRecord) -> Dict[str, FeatureLocation]:
    """Returns a dict of chain ID to feature locations.

    This function looks through the record for features of type "misc_feature"
    and name matching the form "!insert(?)" where "?" is a chain letter.

    Args:
        record: The SeqRecord to search.

    Returns:
        A dictionary mapping chain letters to FeatureLocation objects.
    """

    location_dict = {}

    for feature in record.features:
        if feature.type == "misc_feature" and feature.qualifiers["label"][0].startswith(
                "!insert("
        ):
            chain_letter = feature.qualifiers["label"][0].split("(")[1].split(")")[0]
            location = feature.location
            location_dict[chain_letter] = feature.location
    return location_dict


def replace_sequence_in_record(
    record: SeqRecord.SeqRecord, location: FeatureLocation, insert: SeqRecord.SeqRecord
) -> SeqRecord.SeqRecord:
    """Returns a modified SeqRecord with the insert sequence.

    Replaces the sequence at the specified location with the insert sequence
    and adjusts all feature locations accordingly.

    Args:
        record: The original SeqRecord.
        location: The location to replace.
        insert: The SeqRecord to insert.

    Returns:
        The modified SeqRecord.
    """

    # I don't know if this is right. What does the strand number mean? --rdkibler 210320
    if location.strand >= 0:
        adjusted_seq = (
                record.seq[: location.start] + insert.seq + record.seq[location.end:]
        )
    else:
        adjusted_seq = (
                record.seq[: location.start]
                + insert.reverse_complement().seq
                + record.seq[location.end:]
        )

    record.seq = adjusted_seq

    seq_diff = len(insert) - len(location)
    orig_start = location.start
    orig_end = location.end

    processed_features = []

    # add a feature for the insert
    processed_features.append(
        SeqFeature(
            location=FeatureLocation(
                location.start, location.end + seq_diff, strand=location.strand
            ),
            type="protein",
            qualifiers={"label": [insert.name]},
        )
    )

    for feat in record.features:

        f_loc = feat.location

        loc_list = []

        for subloc in f_loc.parts:

            assert subloc.start <= subloc.end

            # type 1: where the start and end are contained within the original location
            # -> do not add it to the processed_features list because I have no idea where it should go
            if (
                    subloc.start > location.start
                    and subloc.start < location.end
                    and subloc.end > location.start
                    and subloc.end < location.end
            ):
                continue

            # type 1b: where the start and end are the same which will happen a lot for storing constraints and objectives
            elif subloc.start == location.start and subloc.end == location.end:
                new_loc = FeatureLocation(
                    location.start, location.end + seq_diff, strand=subloc.strand
                )

            # type 2a: where they span the location but start or end at the start or end of the location
            elif subloc.start == location.start and subloc.end > location.end:
                new_loc = FeatureLocation(
                    location.start, subloc.end + seq_diff, strand=subloc.strand
                )
            elif subloc.start < location.start and subloc.end == location.end:
                new_loc = FeatureLocation(
                    subloc.start, location.end + seq_diff, strand=subloc.strand
                )

            # type 2b: where they start or end inside the location but don't fully span the location
            # Here I assume that the total length of the annotation is important, so adjust the annotation to have the correct
            # length anchored outside of the insert, unless it'd extend through the insert
            elif subloc.start >= location.start and subloc.start <= location.end:
                # we already caught the case where it's fully within, so I know the end of the subloc extends after the end of the insert
                new_loc = FeatureLocation(
                    max(subloc.end + seq_diff - len(subloc), location.start),
                    subloc.end + seq_diff,
                    strand=subloc.strand,
                )
                assert len(new_loc) == len(subloc)

            elif subloc.end >= location.start and subloc.end <= location.end:
                new_loc = FeatureLocation(
                    subloc.start,
                    min(subloc.end, location.end + seq_diff),
                    strand=subloc.strand,
                )
                assert len(new_loc) == len(subloc)

            # type 3: where they span the location
            # -> keep the leftmost point same and add diff to rightmost. do not split
            elif (
                    location.start >= subloc.start
                    and location.start <= subloc.end
                    and location.end >= subloc.start
                    and location.end <= subloc.end
            ):
                new_loc = FeatureLocation(
                    subloc.start, subloc.end + seq_diff, strand=subloc.strand
                )

            # type 4: where they start and end before location
            # -> add it to list unchanged
            elif subloc.start <= location.start and subloc.end <= location.start:
                new_loc = subloc

            # type 5: where they start and end after location
            # -> add diff to whole location
            elif subloc.start >= location.end and subloc.end >= location.end:
                new_loc = subloc + seq_diff

            loc_list.append(new_loc)

        # if the list is empty, it means that all the sublocs were contained within the insert
        if len(loc_list) > 0:
            feat.location = sum(loc_list)
            processed_features.append(feat)

    record.features = processed_features

    return record


def make_naive_vector_records(
    base_vector_record: SeqRecord.SeqRecord,
    protein_filepaths: List[str],
    increasing_chain_fasta: bool = False,
    do_not_append_vector_name: bool = False,
) -> List[SeqRecord.SeqRecord]:
    """Returns a list of SeqRecords with randomly reverse-translated inserts.

    Args:
        base_vector_record: The base vector SeqRecord.
        protein_filepaths: Paths to protein sequence files.
        increasing_chain_fasta: Whether to use increasing chain letters for FASTA.
        do_not_append_vector_name: Whether to omit the vector name from the record name.

    Returns:
        List of generated SeqRecord objects.
    """

    output_records = []

    insert_locations = get_insert_locations(base_vector_record)
    # print(insert_locations)
    for inserts in load_inserts(protein_filepaths, increasing_chain_fasta):

        if len(insert_locations) == 1:
            for insert in inserts:
                intermediate_vector_record = copy.deepcopy(base_vector_record)
                intermediate_vector_record = replace_sequence_in_record(
                    intermediate_vector_record, insert_locations["A"], insert
                )
                vec_name = intermediate_vector_record.name
                insert_name = insert.name
                if not do_not_append_vector_name:
                    intermediate_vector_record.name = f"{insert_name}__{vec_name}"
                else:
                    intermediate_vector_record.name = f"{insert_name}"
                output_records.append(intermediate_vector_record)
        else:
            intermediate_vector_record = copy.deepcopy(base_vector_record)
            for insert in inserts:
                intermediate_insert_locations = get_insert_locations(
                    intermediate_vector_record
                )
                intermediate_vector_record = replace_sequence_in_record(
                    intermediate_vector_record,
                    intermediate_insert_locations[insert.annotations["chain"]],
                    insert,
                )

            vec_name = intermediate_vector_record.name
            insert_name = insert.name[:-2]  # cuts off _A or whatever chain ID it is
            if not do_not_append_vector_name:
                intermediate_vector_record.name = f"{insert_name}__{vec_name}"
            else:
                intermediate_vector_record.name = f"{insert_name}"
            output_records.append(intermediate_vector_record)

    return output_records


def make_naive_vector_record_by_seq(
    base_vector_record: SeqRecord.SeqRecord, amino_acid_sequence: str
) -> SeqRecord.SeqRecord:
    """Returns a SeqRecord with the protein sequence inserted into the vector.

    Args:
        base_vector_record: The base vector SeqRecord.
        amino_acid_sequence: The amino acid sequence to insert.

    Returns:
        The generated SeqRecord.
    """

    insert_locations = get_insert_locations(base_vector_record)
    # print(insert_locations)
    insert = load_insert_seq(amino_acid_sequence)

    intermediate_vector_record = copy.deepcopy(base_vector_record)
    intermediate_vector_record = replace_sequence_in_record(
        intermediate_vector_record, insert_locations["A"], insert
    )
    vec_name = intermediate_vector_record.name
    insert_name = insert.name
    intermediate_vector_record.name = f"{insert_name}__{vec_name}"
    return intermediate_vector_record
