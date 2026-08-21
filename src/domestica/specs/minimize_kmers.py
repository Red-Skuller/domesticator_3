from collections import Counter
from dnachisel import Specification, SpecEvaluation
from domestica.specs.base import register_specification


@register_specification("MinimizeNumKmers")
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
        all_kmers = [sequence[i : i + self.k] for i in range(len(sequence) - self.k)]
        number_of_non_unique_kmers = sum(
            count for kmer, count in Counter(all_kmers).items() if count > 1
        )
        score = -(float(self.k) * number_of_non_unique_kmers) / len(sequence)
        return SpecEvaluation(
            self,
            problem,
            score=score,
            locations=[self.location],
            message=f"Score: {score:.2f} ({number_of_non_unique_kmers} non-unique {self.k}-mers)",
        )

    def label_parameters(self):
        return [("k", str(self.k))]

    def short_label(self):
        return f"Avoid {self.k}mers {self.boost}"

    def __str__(self):
        return "MinimizeNum%dmers" % self.k
