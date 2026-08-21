import pytest
import dnachisel as dc
from domestica.specs.base import get_all_specifications
from domestica.specs.minimize_kmers import MinimizeNumKmers


def test_specification_registration() -> None:
    """Verify that the custom specification is correctly registered."""
    specs = get_all_specifications()
    assert "MinimizeNumKmers" in specs
    assert specs["MinimizeNumKmers"] == MinimizeNumKmers


def test_minimize_num_kmers_scoring() -> None:
    """Verify that MinimizeNumKmers assigns worse scores to sequences with repeated k-mers."""
    k_val = 4

    # Sequence with high repetition of 4-mers
    seq_repeated = "ATCGATCGATCGATCG"
    problem_rep = dc.DnaOptimizationProblem(sequence=seq_repeated)

    spec_rep = MinimizeNumKmers(k=k_val)
    spec_rep = spec_rep.initialize_on_problem(problem_rep)
    eval_rep = spec_rep.evaluate(problem_rep)

    # Sequence with unique 4-mers
    seq_unique = "ATCGAGCTTAGCGACT"
    problem_uniq = dc.DnaOptimizationProblem(sequence=seq_unique)

    spec_uniq = MinimizeNumKmers(k=k_val)
    spec_uniq = spec_uniq.initialize_on_problem(problem_uniq)
    eval_uniq = spec_uniq.evaluate(problem_uniq)

    # Assertions
    assert eval_rep.score < eval_uniq.score
    assert eval_rep.best_possible_score == 0
    assert eval_uniq.best_possible_score == 0


def test_specification_labels() -> None:
    """Verify parameter labeling and string representations."""
    spec = MinimizeNumKmers(k=9, boost=2.0)
    assert spec.short_label() == "Avoid 9mers 2.0"
    assert spec.label_parameters() == [("k", "9")]
    assert str(spec) == "MinimizeNum9mers"