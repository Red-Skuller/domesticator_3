import logging
import re
from collections import Counter
from typing import Tuple, Optional, Callable
from pathlib import Path

import dnachisel as dc
from dnachisel import Specification, SpecEvaluation, reverse_translate, DEFAULT_SPECIFICATIONS_DICT
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

logger = logging.getLogger(__name__)


class MinimizeNumKmers(Specification):
    """Minimize a no-kmers score."""
    best_possible_score = 0

    def __init__(self, k=9, location=None, boost=1.0):
        self.location = location
        self.k = k
        self.boost = boost

    def initialize_on_problem(self, problem, role=None):
        return self._copy_with_full_span_if_no_location(problem)

    def evaluate(self, problem):
        sequence = self.location.extract_sequence(problem.sequence)
        all_kmers = [sequence[i: i + self.k] for i in range(len(sequence) - self.k)]
        number_of_non_unique_kmers = sum(
            count for kmer, count in Counter(all_kmers).items() if count > 1
        )
        score = -(float(self.k) * number_of_non_unique_kmers) / len(sequence)
        return SpecEvaluation(
            self, problem, score=score, locations=[self.location],
            message=f"Score: {score:.2f} ({number_of_non_unique_kmers} non-unique {self.k}-mers)"
        )

    def label_parameters(self): return [("k", str(self.k))]

    def short_label(self): return f"Avoid {self.k}mers {self.boost}"

    def __str__(self): return "MinimizeNum%dmers" % self.k


def _get_codon_table(template_record: SeqRecord) -> str:
    for feature in template_record.features:
        labels = feature.qualifiers.get("label", []) + feature.qualifiers.get("note", [])
        if any("!INSERT" in str(item).upper() for item in labels):
            for item in labels:
                match = re.search(r'@EnforceTranslation\s*\(\s*genetic_table\s*=\s*["\']([^"\']+)["\']\s*\)', str(item))
                if match: return match.group(1).strip().lower()
            break
    return "Standard"


def _insert_into_template(template_record: SeqRecord, insert_dna: str) -> SeqRecord:
    insert_feature = None
    for feature in template_record.features:
        labels = feature.qualifiers.get("label", []) + feature.qualifiers.get("note", [])
        if any("INSERT" in label.upper() for label in labels):
            insert_feature = feature
            break

    if not insert_feature:
        raise ValueError("Template vector must contain a feature with label or note 'INSERT'.")

    start, end = int(insert_feature.location.start), int(insert_feature.location.end)
    old_length, new_length = end - start, len(insert_dna)

    left_flank = template_record[:start]
    right_flank = template_record[end:]
    insert_record = SeqRecord(Seq(insert_dna), id="target_insert", annotations={"molecule_type": "DNA"})

    for feature in template_record.features:
        f_start, f_end = int(feature.location.start), int(feature.location.end)
        if f_start >= start and f_end <= end:
            rel_start, rel_end = f_start - start, f_end - start
            new_start = 0 if (rel_start == 0 and rel_end == old_length) else rel_start
            new_end = new_length if (rel_start == 0 and rel_end == old_length) else min(rel_end, new_length)

            insert_record.features.append(SeqFeature(
                FeatureLocation(new_start, new_end, strand=feature.location.strand),
                type=feature.type, qualifiers=feature.qualifiers
            ))

    merged_record = left_flank + insert_record + right_flank
    merged_record.annotations["topology"] = template_record.annotations.get("topology", "linear")
    return merged_record


def optimize_sequence(
        protein_sequence: Optional[str],
        template_path: Path,
        evaluator: Optional[Callable[[str], Tuple[bool, Optional[float]]]] = None
) -> Tuple[Optional[str], str, bool, Optional[float], Optional[SeqRecord]]:
    records = list(SeqIO.parse(template_path, "genbank"))
    if len(records) != 1: raise ValueError("Template must contain exactly one record.")
    template_record = records[0]

    # Structural branch based on amino acid sequence availability
    if protein_sequence and protein_sequence.strip():
        naive_dna = reverse_translate(protein_sequence, table=_get_codon_table(template_record))
        merged_record = _insert_into_template(template_record, naive_dna)
    else:
        naive_dna = None
        merged_record = template_record

    custom_specs = dict(DEFAULT_SPECIFICATIONS_DICT)
    custom_specs["MinimizeNumKmers"] = MinimizeNumKmers

    trial, max_tries, current_max_iters = 0, 10, 1000
    problem = None

    while trial < max_tries:
        try:
            problem = dc.DnaOptimizationProblem.from_record(merged_record, specifications_dict=custom_specs)
            problem.resolve_constraints()
            problem.optimize()
            problem.resolve_constraints(final_check=True)
            break
        except dc.NoSolutionError:
            current_max_iters += 1000
            trial += 1
            if trial >= max_tries:
                raise dc.NoSolutionError("Failed to converge on a valid optimization solution.")

    if protein_sequence and protein_sequence.strip():
        insert_start, insert_end = 0, 0
        for feature in problem.record.features:
            labels = feature.qualifiers.get("label", []) + feature.qualifiers.get("note", [])
            if any("INSERT" in label.upper() for label in labels):
                insert_start, insert_end = int(feature.location.start), int(feature.location.end)
                break
        optimized_output_dna = problem.sequence[insert_start:insert_end]
    else:
        optimized_output_dna = problem.sequence

    is_accepted, final_score = True, None
    if evaluator:
        is_accepted, final_score = evaluator(problem.sequence)

    return naive_dna, optimized_output_dna, is_accepted, final_score, problem.record