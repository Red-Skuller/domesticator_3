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
    # 1. Identify all features containing the target INSERT annotation
    insert_features = []
    for feature in template_record.features:
        labels = (
            feature.qualifiers.get("label", [])
            + feature.qualifiers.get("note", [])
            + feature.qualifiers.get("locus_tag", [])
        )
        if any("!INSERT" in str(item).upper() for item in labels):
            insert_features.append(feature)

    if not insert_features:
        logger.debug("No INSERT annotation discovered. Defaulting to 'Standard'.")
        return "Standard"

    # 2. Define regex to capture the 'genetic_table' parameter value
    regex = r'@EnforceTranslation\s*\((?=[^)]*?\bgenetic_table\s*=\s*(?P<quote>["\']?)(?P<value>[^,"\')\s]+)(?P=quote))[^)]*\)'

    # 3. Scan for features containing @EnforceTranslation that overlap with the INSERT features
    for insert_feat in insert_features:
        for feature in template_record.features:
            # Check for physical sequence overlap: max(start1, start2) < min(end1, end2)
            if max(feature.location.start, insert_feat.location.start) < min(feature.location.end, insert_feat.location.end):
                labels = (
                    feature.qualifiers.get("label", [])
                    + feature.qualifiers.get("note", [])
                    + feature.qualifiers.get("locus_tag", [])
                )
                for item in labels:
                    match = re.search(regex, str(item))
                    if match:
                        table_name = match.group("value").strip()
                        logger.debug("Extracted genetic table rule: '%s' from overlapping feature qualifiers.", table_name)
                        return table_name

    logger.debug("No custom genetic table definition discovered in overlapping intervals. Defaulting to 'Standard'.")
    return "Standard"


def _insert_into_template(template_record: SeqRecord, insert_dna: str) -> SeqRecord:
    insert_feature = None
    for feature in template_record.features:
        labels = feature.qualifiers.get("label", []) + feature.qualifiers.get("note", [])
        if any("INSERT" in label.upper() for label in labels):
            insert_feature = feature
            break

    if not insert_feature:
        logger.error("Template structural invalidation: Failed to discover an 'INSERT' tag indicator within sequence features.")
        raise ValueError("Template vector must contain a feature with label or note 'INSERT'.")

    start, end = int(insert_feature.location.start), int(insert_feature.location.end)
    old_length, new_length = end - start, len(insert_dna)
    delta = new_length - old_length
    logger.debug("Modifying sequence template target interval locus boundaries: [%d, %d] (Original Len: %d, Insert Len: %d)", start, end, old_length, new_length)

    left_flank = template_record[:start]
    right_flank = template_record[end:]
    insert_record = SeqRecord(Seq(insert_dna), id="target_insert", annotations={"molecule_type": "DNA"})

    def fully_inside(f_start, f_end):
        return f_start >= start and f_end <= end

    def fully_left(f_end):
        return f_end <= start

    def fully_right(f_start):
        return f_start >= end

    for feature in template_record.features:
        f_start, f_end = int(feature.location.start), int(feature.location.end)
        if fully_inside(f_start, f_end):
            rel_start, rel_end = f_start - start, f_end - start
            new_start = 0 if (rel_start == 0 and rel_end == old_length) else rel_start
            new_end = new_length if (rel_start == 0 and rel_end == old_length) else min(rel_end, new_length)

            insert_record.features.append(SeqFeature(
                FeatureLocation(new_start, new_end, strand=feature.location.strand),
                type=feature.type, qualifiers=feature.qualifiers
            ))

    merged_record = left_flank + insert_record + right_flank
    merged_record.annotations["topology"] = template_record.annotations.get("topology", "linear")

    # Re-attach features that straddle the INSERT window boundary, which are
    # otherwise dropped: Biopython's slicing only keeps features fully inside
    # a slice, and the loop above only keeps features fully inside the window.
    left_len = len(left_flank)

    def map_pos(p):
        if p <= start:
            return p
        if p >= end:
            return p + delta
        return left_len + (p - start)

    for feature in template_record.features:
        f_start, f_end = int(feature.location.start), int(feature.location.end)
        if fully_inside(f_start, f_end) or fully_left(f_end) or fully_right(f_start):
            continue  # already handled above, via left_flank, or via right_flank

        new_start, new_end = map_pos(f_start), map_pos(f_end)
        logger.debug("Re-attaching boundary-spanning feature %s: [%d, %d] -> [%d, %d]",
                     feature.type, f_start, f_end, new_start, new_end)
        merged_record.features.append(SeqFeature(
            FeatureLocation(new_start, new_end, strand=feature.location.strand),
            type=feature.type, qualifiers=feature.qualifiers
        ))
    return merged_record


def optimize_sequence(
        protein_sequence: Optional[str],
        template_path: Path,
        evaluator: Optional[Callable[[str], Tuple[bool, Optional[float]]]] = None
) -> Tuple[Optional[str], str, bool, Optional[float], Optional[SeqRecord]]:
    logger.debug("Loading single record vector structural profile from path: %s", template_path)
    records = list(SeqIO.parse(template_path, "genbank"))
    if len(records) != 1:
        logger.error("Template structure constraint error: Expected 1 record profile, discovered %d.", len(records))
        raise ValueError("Template must contain exactly one record.")
    template_record = records[0]

    if protein_sequence and protein_sequence.strip():
        logger.debug("Protein sequence input confirmed. Initiating algorithmic reverse-translation mappings.")
        naive_dna = reverse_translate(protein_sequence, table=_get_codon_table(template_record))
        merged_record = _insert_into_template(template_record, naive_dna)
    else:
        logger.info("Empty or null protein sequence entry. Directing sequence execution strategy onto base template.")
        naive_dna = None
        merged_record = template_record

    custom_specs = dict(DEFAULT_SPECIFICATIONS_DICT)
    custom_specs["MinimizeNumKmers"] = MinimizeNumKmers

    trial, max_tries, current_max_iters = 0, 10, 1000
    problem = None

    while trial < max_tries:
        logger.debug("Beginning sequence optimization cycle loop attempt %d/%d (Iteration Limit: %d)", trial + 1, max_tries, current_max_iters)
        try:
            problem = dc.DnaOptimizationProblem.from_record(merged_record, specifications_dict=custom_specs)
            logger.debug("Evaluating heuristic local constraint resolution criteria...")
            problem.resolve_constraints()
            logger.debug("Executing core dnachisel local optimization operations...")
            problem.optimize()
            #problem.optimize_with_report(target="/home/lukah/Downloads/report.zip")
            logger.debug("Executing final structural global constraint consistency verification checks...")
            problem.resolve_constraints(final_check=True)
            logger.debug("Optimization problem constraints converged successfully on attempt loop %d.", trial + 1)
            break
        except dc.NoSolutionError as nse:
            logger.warning("Optimization pass iteration %d failed to converge under current problem parameters: %s", trial + 1, str(nse))
            current_max_iters += 1000
            trial += 1
            if trial >= max_tries:
                logger.error(
                    "Critical convergence breakdown: Complete search space exhausted without locating acceptable structural solutions.")
                raise RuntimeError("Failed to converge on a valid optimization solution.") from nse
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
        logger.info("Invoking vendor API sequence optimization score validation checks.")
        try:
            is_accepted, final_score = evaluator(problem.sequence)
            logger.info("Vendor evaluation assessment processing completed. Acceptance Status: %s, Assigned Metrics: %s", is_accepted, final_score)
        except Exception:
            logger.exception("Critical communication or validation exception thrown during external vendor complexity evaluation processing.")
            raise

    return naive_dna, optimized_output_dna, is_accepted, final_score, problem.record