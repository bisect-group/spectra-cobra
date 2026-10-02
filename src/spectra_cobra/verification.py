"""Check that an extracted model is actually fit for purpose.

An extraction promises two things: that every core reaction is in the result,
and that the result is flux consistent. The first always holds. The second
does **not** always hold for the LP formulations, because of how they decide
which reactions to keep.

``minNetLP`` and ``growthOptim`` read their answer off a flux vector, keeping
any reaction whose flux clears a cutoff near the solver's own tolerance. An L1
objective has many optimal solutions and spreads tiny fluxes over thousands of
reactions, so discarding everything below that cutoff throws away flux that
was balancing a hub metabolite. The kept flux vector is then only
approximately mass balanced, and a reaction whose own flux was near ``tol``
can turn out unable to carry it — a core reaction included.

The mixed-integer formulations do not have this problem: a binary at one
forces its reaction to carry at least ``tol`` and a binary at zero forces
exactly zero, so the kept flux vector is exactly balanced.

:func:`check_extraction` reports on a result so the difference is visible
rather than silent.
"""

from logging import getLogger
from typing import TYPE_CHECKING, Iterable, List, Optional

from ._orientation import STOICHIOMETRY
from .consistency import consistent_reaction_ids

if TYPE_CHECKING:
    from cobra.core import Model


logger = getLogger(__name__)


class ExtractionReport:
    """What an extracted model does and does not deliver.

    Attributes
    ----------
    n_reactions : int
        How many reactions the extracted model has.
    missing_core : list of str
        Core reactions absent from the extracted model. Should always be
        empty; an extraction that drops a core reaction is a bug.
    blocked : list of str
        Reactions in the extracted model that cannot carry flux in it.
    blocked_core : list of str
        The core reactions among `blocked`. These are the damaging ones: the
        reaction is present, so a glance at the model suggests the extraction
        worked, but it cannot carry flux and so cannot play the role it was
        chosen for.
    is_consistent : bool
        Whether the extracted model is flux consistent.
    is_valid : bool
        Whether the extracted model delivers both promises, i.e. no missing
        core reactions and nothing blocked.

    """

    def __init__(
        self,
        n_reactions: int,
        missing_core: List[str],
        blocked: List[str],
        blocked_core: List[str],
    ) -> None:
        """Record the findings.

        Parameters
        ----------
        n_reactions : int
            How many reactions the extracted model has.
        missing_core : list of str
            Core reactions absent from the extracted model.
        blocked : list of str
            Reactions that cannot carry flux in the extracted model.
        blocked_core : list of str
            The core reactions among `blocked`.

        """
        self.n_reactions = n_reactions
        self.missing_core = missing_core
        self.blocked = blocked
        self.blocked_core = blocked_core

    @property
    def is_consistent(self) -> bool:
        """Return whether the extracted model is flux consistent.

        Returns
        -------
        bool
            True when no reaction in it is blocked.

        """
        return not self.blocked

    @property
    def is_valid(self) -> bool:
        """Return whether the extraction delivered on both promises.

        Returns
        -------
        bool
            True when every core reaction is present and able to carry flux,
            and no other reaction is blocked either.

        """
        return not self.missing_core and not self.blocked

    def __repr__(self) -> str:
        """Return a one-line summary.

        Returns
        -------
        str
            The summary.

        """
        return (
            f"<ExtractionReport {self.n_reactions} reactions, "
            f"{len(self.missing_core)} core missing, "
            f"{len(self.blocked)} blocked "
            f"({len(self.blocked_core)} of them core), "
            f"valid={self.is_valid}>"
        )

    def summary(self) -> str:
        """Return a readable multi-line summary.

        Returns
        -------
        str
            The summary, suitable for printing or logging.

        """
        lines = [f"{self.n_reactions} reactions in the extracted model"]
        if self.missing_core:
            lines.append(
                f"  {len(self.missing_core)} core reaction(s) MISSING: "
                f"{self.missing_core[:5]}"
            )
        else:
            lines.append("  every core reaction is present")
        if self.blocked:
            lines.append(f"  {len(self.blocked)} reaction(s) cannot carry flux")
            if self.blocked_core:
                lines.append(
                    f"  {len(self.blocked_core)} of those are CORE reactions: "
                    f"{self.blocked_core[:5]}"
                )
        else:
            lines.append("  flux consistent: every reaction can carry flux")
        return "\n".join(lines)


def check_extraction(
    extracted: "Model",
    core_reactions: Iterable,
    tol: float = 1e-4,
    consistency_type: str = STOICHIOMETRY,
    detection_cutoff: Optional[float] = None,
    seed: Optional[int] = None,
) -> ExtractionReport:
    """Check an extracted model against what the extraction promised.

    Parameters
    ----------
    extracted : cobra.Model
        The model returned by :func:`~spectra_cobra.spectra_me` or
        :func:`~spectra_cobra.spectra_ccme`.
    core_reactions : iterable
        The core reactions that were asked for, as identifiers or
        :class:`cobra.Reaction` objects. For
        :func:`~spectra_cobra.spectra_ccme`, leave out the ones it reported as
        blocked, since it was never going to include those.
    tol : float, optional
        The flux threshold the extraction was run with (default 1e-4). Use the
        same value, or the answer will not describe the model you have.
    consistency_type : {"stoichiometry", "topology"}, optional
        The consistency type the extraction was run with
        (default "stoichiometry").
    detection_cutoff : float, optional
        The absolute flux at which a reaction counts as carrying flux
        (default ``0.99 * tol``).
    seed : int, optional
        A seed for the consistency check's random coefficients (default None).

    Returns
    -------
    ExtractionReport
        What the model does and does not deliver.

    Examples
    --------
    >>> report = check_extraction(extracted, core, tol=1e-4)
    >>> report.is_valid
    True

    Notes
    -----
    This solves a handful of LPs, so it is cheap next to the extraction
    itself, and worth running whenever the result matters. See
    :mod:`spectra_cobra.verification` for why an LP-based extraction can
    return a model with blocked reactions in the first place.

    """
    core_ids = [c if isinstance(c, str) else c.id for c in core_reactions]
    present = {rxn.id for rxn in extracted.reactions}

    missing_core = [c for c in core_ids if c not in present]

    consistent, _ = consistent_reaction_ids(
        extracted, tol, consistency_type, detection_cutoff, seed
    )
    blocked = [rxn.id for rxn in extracted.reactions if rxn.id not in consistent]
    blocked_core = [c for c in core_ids if c in set(blocked)]

    report = ExtractionReport(
        n_reactions=len(extracted.reactions),
        missing_core=missing_core,
        blocked=blocked,
        blocked_core=blocked_core,
    )
    if not report.is_valid:
        logger.warning(
            "The extracted model is not valid: %d core reaction(s) missing, "
            "%d reaction(s) blocked, %d of those core.",
            len(missing_core),
            len(blocked),
            len(blocked_core),
        )
    return report
