import logging
import re
from collections import Counter
from typing import Tuple, Optional, Callable
from pathlib import Path

import dnachisel as dc
from dnachisel import (
    Specification,
    SpecEvaluation,
    reverse_translate,
    DEFAULT_SPECIFICATIONS_DICT,
)
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from domestica.specs.base import get_all_specifications

logger = logging.getLogger(__name__)


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
            if max(feature.location.start, insert_feat.location.start) < min(
                feature.location.end, insert_feat.location.end
            ):
                labels = (
                    feature.qualifiers.get("label", [])
                    + feature.qualifiers.get("note", [])
                    + feature.qualifiers.get("locus_tag", [])
                )
                for item in labels:
                    match = re.search(regex, str(item))
                    if match:
                        table_name = match.group("value").strip()
                        logger.debug(
                            "Extracted genetic table rule: '%s' from overlapping feature qualifiers.",
                            table_name,
                        )
                        return table_name

    logger.debug(
        "No custom genetic table definition discovered in overlapping intervals. Defaulting to 'Standard'."
    )
    return "Standard"


def _insert_into_template(template_record: SeqRecord, insert_dna: str) -> SeqRecord:
    insert_feature = None
    for feature in template_record.features:
        labels = feature.qualifiers.get("label", []) + feature.qualifiers.get(
            "note", []
        )
        if any("INSERT" in label.upper() for label in labels):
            insert_feature = feature
            break

    if not insert_feature:
        logger.error(
            "Template structural invalidation: Failed to discover an 'INSERT' tag indicator within sequence features."
        )
        raise ValueError(
            "Template vector must contain a feature with label or note 'INSERT'."
        )

    start, end = int(insert_feature.location.start), int(insert_feature.location.end)
    old_length, new_length = end - start, len(insert_dna)
    delta = new_length - old_length
    logger.debug(
        "Modifying sequence template target interval locus boundaries: [%d, %d] (Original Len: %d, Insert Len: %d)",
        start,
        end,
        old_length,
        new_length,
    )

    left_flank = template_record[:start]
    right_flank = template_record[end:]

    # Explicitly clear features to prevent Biopython boundary-slicing omission bugs
    left_flank.features = []
    right_flank.features = []

    insert_record = SeqRecord(
        Seq(insert_dna), id="target_insert", annotations={"molecule_type": "DNA"}
    )
    merged_record = left_flank + insert_record + right_flank

    # Transfer original annotations and metadata
    merged_record.id = template_record.id
    merged_record.name = template_record.name
    merged_record.description = template_record.description
    merged_record.annotations = template_record.annotations.copy()

    merged_record.annotations["topology"] = template_record.annotations.get(
        "topology", "linear"
    )  # TODO make it so topology is inherited from template i.e. implement circular templates and outputs
    left_len = len(left_flank)

    def map_pos(p):
        if p <= start:
            return p
        if p >= end:
            return p + delta
        return left_len + (p - start)

    # Re-attach all features manually using exact relative coordinates
    for feature in template_record.features:
        f_start, f_end = int(feature.location.start), int(feature.location.end)

        if f_start >= start and f_end <= end:
            # Feature is fully contained within the replaced insert region
            rel_start, rel_end = f_start - start, f_end - start
            if rel_start == 0 and rel_end == old_length:
                new_start = left_len
                new_end = left_len + new_length
            else:
                new_start = left_len + rel_start
                new_end = left_len + min(rel_end, new_length)
        else:
            # Feature spans boundaries, is fully left, or fully right
            new_start = map_pos(f_start)
            new_end = map_pos(f_end)

        merged_record.features.append(
            SeqFeature(
                FeatureLocation(new_start, new_end, strand=feature.location.strand),
                type=feature.type,
                qualifiers=feature.qualifiers,
            )
        )

    return merged_record


def optimize_sequence(
    protein_sequence: Optional[str],
    template_path: Path,
    min_length: int = 300,
    evaluator: Optional[Callable[[str], Tuple[bool, Optional[float]]]] = None,
) -> Tuple[Optional[str], str, bool, Optional[float], Optional[SeqRecord]]:
    logger.debug(
        "Loading single record vector structural profile from path: %s", template_path
    )
    records = list(SeqIO.parse(template_path, "genbank"))
    if len(records) != 1:
        logger.error(
            "Template structure constraint error: Expected 1 record profile, discovered %d.",
            len(records),
        )
        raise ValueError("Template must contain exactly one record.")
    template_record = records[0]

    if protein_sequence and protein_sequence.strip():
        logger.debug(
            "Protein sequence input confirmed. Initiating algorithmic reverse-translation mappings."
        )
        naive_dna = reverse_translate(
            protein_sequence, table=_get_codon_table(template_record)
        )
        merged_record = _insert_into_template(template_record, naive_dna)
    else:
        logger.info(
            "Empty or null protein sequence entry. Directing sequence execution strategy onto base template."
        )
        naive_dna = None
        merged_record = template_record

    current_length = len(merged_record.seq)
    if current_length < min_length:
        num_pad = min_length - current_length
        logger.info(
            "Sequence length (%d bp) is below %d bp. Appending %d bp of ACGT padding.",
            current_length,
            min_length,
            num_pad,
        )

        # Pad with a neutral repeating sequence instead of "A"s to prevent massive overlapping
        # AvoidPattern(AAAAAA) breaches which trigger a localization bug in dnachisel
        padding_seq = ("ATGG" * (num_pad // 4 + 1))[:num_pad]
        merged_record.seq = merged_record.seq + padding_seq
        padding_feature = SeqFeature(
            FeatureLocation(current_length, min_length, strand=1),
            type="misc_feature",
            qualifiers={
                "note": [
                    "@AvoidHairpins() & @AvoidPattern(AAAAAA) & @AvoidPattern(CCCCCC) & "
                    "@AvoidPattern(GAGACC) & @AvoidPattern(GGGGGG) & @AvoidPattern(TTTTTT) & "
                    "@EnforceGCContent(25-80%/50bp) & @EnforceGCContent(40-65%) & ~MinimizeNumKmers(8, boost=10)"
                ],
                "label": ["padding_as"],
            },
        )
        merged_record.features.append(padding_feature)

    custom_specs = dict(DEFAULT_SPECIFICATIONS_DICT)
    custom_specs.update(get_all_specifications())

    trial, max_tries, current_max_iters = 0, 10, 1000
    problem = None

    while trial < max_tries:
        logger.debug(
            "Beginning sequence optimization cycle loop attempt %d/%d (Iteration Limit: %d)",
            trial + 1,
            max_tries,
            current_max_iters,
        )
        try:
            problem = dc.DnaOptimizationProblem.from_record(
                merged_record, specifications_dict=custom_specs
            )
            logger.debug("Evaluating heuristic local constraint resolution criteria...")
            problem.resolve_constraints()
            logger.debug("Executing core dnachisel local optimization operations...")
            problem.optimize()
            # REMOVED hardcoded debug paths here
            logger.debug(
                "Executing final structural global constraint consistency verification checks..."
            )
            logger.info(problem.constraints_text_summary())
            logger.info(problem.objectives_text_summary())
            logger.debug(
                "Optimization problem constraints converged successfully on attempt loop %d.",
                trial + 1,
            )
            break
        except dc.NoSolutionError as nse:
            logger.warning(
                "Optimization pass iteration %d failed to converge: %s",
                trial + 1,
                str(nse),
            )
            current_max_iters += 1000
            trial += 1
            if trial >= max_tries:
                logger.error(
                    "Critical convergence breakdown: Search space exhausted without locating acceptable structures."
                )
                raise RuntimeError(
                    "Failed to converge on a valid optimization solution."
                ) from nse

    merged_record.seq = Seq(problem.sequence)
    optimized_output_dna = str(merged_record.seq)

    is_accepted, final_score = True, None
    if evaluator:
        logger.info(
            "Invoking vendor API sequence optimization score validation checks."
        )
        try:
            is_accepted, final_score = evaluator(problem.sequence)
        except Exception:
            logger.exception("Critical communication or validation exception.")
            raise

    return naive_dna, optimized_output_dna, is_accepted, final_score, merged_record
