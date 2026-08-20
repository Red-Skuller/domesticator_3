import importlib
import pkgutil
from typing import Dict, Type
from pathlib import Path
from dnachisel import Specification

_SPEC_REGISTRY: Dict[str, Type[Specification]] = {}

def register_specification(spec_name: str):
    def decorator(cls: Type[Specification]) -> Type[Specification]:
        _SPEC_REGISTRY[spec_name] = cls
        return cls
    return decorator

def load_specifications() -> None:
    package_dir = Path(__file__).resolve().parent
    for _, module_name, _ in pkgutil.iter_modules([str(package_dir)]):
        importlib.import_module(f"domestica.specs.{module_name}")

def get_all_specifications() -> Dict[str, Type[Specification]]:
    load_specifications()
    return _SPEC_REGISTRY.copy()