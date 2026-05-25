from Bio.SeqUtils.ProtParam import ProteinAnalysis

PKA = {
    "Cterm": 3.55, "Nterm": 7.50,
    "D": 3.90, "E": 4.07, "C": 8.50, "Y": 10.10,
    "H": 6.04, "K": 10.54, "R": 12.48,
}


def net_charge(seq: str, pH: float = 7.4) -> float:
    """Calculate the net charge of a protein sequence at a given pH.

    Args:
        seq: The amino acid sequence.
        pH: The pH value for the calculation.

    Returns:
        The calculated net charge.
    """
    if not seq: return 0.0

    nterm = (10 ** (PKA["Nterm"] - pH)) / (1 + 10 ** (PKA["Nterm"] - pH))
    cterm = - (1 / (1 + 10 ** (PKA["Cterm"] - pH)))

    counts = {aa: seq.count(aa) for aa in "DECYHKR"}

    neg = (
            counts["D"] * (1 / (1 + 10 ** (PKA["D"] - pH))) +
            counts["E"] * (1 / (1 + 10 ** (PKA["E"] - pH))) +
            counts["C"] * (1 / (1 + 10 ** (PKA["C"] - pH))) +
            counts["Y"] * (1 / (1 + 10 ** (PKA["Y"] - pH)))
    )
    pos = (
            counts["H"] * ((10 ** (PKA["H"] - pH)) / (1 + 10 ** (PKA["H"] - pH))) +
            counts["K"] * ((10 ** (PKA["K"] - pH)) / (1 + 10 ** (PKA["K"] - pH))) +
            counts["R"] * ((10 ** (PKA["R"] - pH)) / (1 + 10 ** (PKA["R"] - pH)))
    )
    return nterm + cterm + pos - neg


# Extensible metric registry.
# Easily add new metrics by mapping a command-line key to a callable function.
METRICS_REGISTRY = {
    "mw": lambda seq, ph: ProteinAnalysis(seq).molecular_weight(),
    "pi": lambda seq, ph: ProteinAnalysis(seq).isoelectric_point(),
    "gravy": lambda seq, ph: ProteinAnalysis(seq).gravy(),
    "net_charge": lambda seq, ph: net_charge(seq, ph),
    "length": lambda seq, ph: len(seq)
}


from typing import Callable, List, Dict, Any

def register_metric(name: str, func: Callable[[str, float], Any]):
    """Programmatically register a custom metric function.

    Args:
        name: The name of the metric.
        func: The function to calculate the metric, taking (sequence, ph) and returning any value.
    """
    METRICS_REGISTRY[name] = func


def calculate_selected_params(sequence: str, requested_params: List[str], ph: float = 7.4) -> Dict[str, Any]:
    """Calculates only the requested protein parameters by iterating over the registry.

    Args:
        sequence: The amino acid sequence.
        requested_params: List of parameters to calculate.
        ph: The pH value for relevant calculations.

    Returns:
        A dictionary mapping parameter names to their calculated values.
    """
    if not sequence or not requested_params:
        return {}

    results = {}
    calc_all = "all" in requested_params

    for param, func in METRICS_REGISTRY.items():
        if calc_all or param in requested_params:
            try:
                results[param] = func(sequence, ph)
            except Exception:
                results[param] = None  # Gracefully handle calculation errors

    return results