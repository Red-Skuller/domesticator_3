import importlib
import pkgutil
from abc import ABC, abstractmethod
from typing import Tuple, Optional, Dict, Type
from pathlib import Path

_REGISTRY: Dict[str, Type["ComplexityEvaluator"]] = {}

class ComplexityEvaluator(ABC):
    def __init__(self, product: str):
        self.product = product

    @abstractmethod
    def evaluate(self, sequence: str) -> Tuple[bool, Optional[float]]:
        pass

def register_vendor(vendor_id: str):
    def decorator(cls: Type[ComplexityEvaluator]) -> Type[ComplexityEvaluator]:
        _REGISTRY[vendor_id.lower()] = cls
        return cls
    return decorator

def load_vendors() -> None:
    package_dir = Path(__file__).resolve().parent
    for _, module_name, _ in pkgutil.iter_modules([str(package_dir)]):
        importlib.import_module(f"domestica.vendors.{module_name}")

def get_evaluator(vendor_id: str, product: str) -> ComplexityEvaluator:
    load_vendors()
    vendor_lower = vendor_id.lower()
    if vendor_lower not in _REGISTRY:
        raise ValueError(f"Vendor '{vendor_id}' unsupported. Available: {list(_REGISTRY.keys())}")
    return _REGISTRY[vendor_lower](product=product)