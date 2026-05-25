from abc import ABC, abstractmethod
from typing import List, Tuple, Dict, Any
from Bio.SeqRecord import SeqRecord

class SequenceOptimizer(ABC):
    @abstractmethod
    def optimize(self, record: SeqRecord, **kwargs) -> Any:
        pass

class ComplexityEvaluator(ABC):
    @abstractmethod
    def evaluate(self, sequence: str, **kwargs) -> Tuple[float, List[Dict[str, Any]]]:
        pass
