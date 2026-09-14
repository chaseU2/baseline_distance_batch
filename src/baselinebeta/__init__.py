"""baselinebeta: distance-to-individual-baseline analysis for longitudinal beta diversity."""

from .analysis import BaselineDistanceAnalysis, BaselineResults
from .distances import load_distance_matrix
from .statistics import cliffs_delta
from .subsets import apply_subset, generate_subsets

__all__ = [
    "BaselineDistanceAnalysis",
    "BaselineResults",
    "load_distance_matrix",
    "cliffs_delta",
    "apply_subset",
    "generate_subsets",
]

__version__ = "0.1.0"
