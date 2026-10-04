"""SPECTRA for cobrapy.

SPECTRA reconstructs metabolic networks from multi-omics data at a range of
biological scales: minimal reactomes, context-specific models, gap-filled
reconstructions, minimal microbiomes, microbial community models and
multi-tissue models. What changes between them is the universal model, the
evidence supplied and the objective chosen, not the routine called.

This implementation is built on :mod:`cobra`, expressing the formulations
through cobrapy's own solver interface rather than raw matrices. A MATLAB
implementation built on the COBRA Toolbox is also available; see the
documentation for the differences between the two.

The usual order of operations is:

1. :func:`spectra_cc` to drop the blocked reactions of a universal model, or
   skip it and use :func:`spectra_ccme`, which folds the check into the
   extraction.
2. :func:`spectra_me` to extract a context-specific model around a set of core
   reactions, choosing a network inference formulation with `problem_type`.

Reference
---------
S, P. K., Sridhar, S., Alsmadi, N., Mahadevan, R., & Bhatt, N. P. (2026).
Generalist method to reconstruct metabolic networks from multi-omics data at
large-scale. bioRxiv. https://doi.org/10.64898/2026.04.02.716249
"""

from ._orientation import STOICHIOMETRY, TOPOLOGY
from .consistency import (
    blocked_reaction_ids,
    consistent_reaction_ids,
    spectra_cc,
)
from .exceptions import (
    SpectraError,
    SpectraInfeasibleCoreError,
    SpectraSolverError,
)
from .extraction import (
    CORE_DIRECTION,
    GROWTH_OPTIM,
    MIN_NET_LP,
    MIN_NET_MILP,
    PATHWAY_EXCLUSION,
    PROBLEM_TYPES,
    TRADE_OFF,
    spectra_ccme,
    spectra_me,
)
from .flux_reducer import flux_reducer
from .formulations import growth_optim, min_net_lp, min_net_milp, trade_off
from .verification import ExtractionReport, check_extraction

__version__ = "0.1.0.dev0"

__all__ = (
    # Consistency.
    "spectra_cc",
    "consistent_reaction_ids",
    "blocked_reaction_ids",
    # Extraction.
    "spectra_me",
    "spectra_ccme",
    # Formulations, for use on their own.
    "min_net_lp",
    "min_net_milp",
    "trade_off",
    "growth_optim",
    "flux_reducer",
    # Verification.
    "check_extraction",
    "ExtractionReport",
    # Option names.
    "STOICHIOMETRY",
    "TOPOLOGY",
    "MIN_NET_LP",
    "MIN_NET_MILP",
    "TRADE_OFF",
    "GROWTH_OPTIM",
    "PROBLEM_TYPES",
    "CORE_DIRECTION",
    "PATHWAY_EXCLUSION",
    # Exceptions.
    "SpectraError",
    "SpectraSolverError",
    "SpectraInfeasibleCoreError",
    "__version__",
)
