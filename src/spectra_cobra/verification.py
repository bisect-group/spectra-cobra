"""Check that an extracted model is actually fit for purpose.

An extraction promises two things: that every core reaction is in the result,
and that the result is flux consistent. The first always holds. The second
does **not** always hold for the LP formulations, because of how they decide
which reactions to keep.

``minNetLP`` and ``growthOptim`` read their answer off a flux vector, keeping
any reaction whose flux clears a cutoff near the solver's own tolerance. A few
reactions end up carrying flux just under that cutoff while doing load-bearing
balancing work, and dropping them leaves a mass-balance residual of the same
order. That residual is above the solver's feasibility tolerance, so the kept
flux vector is not actually feasible in the extracted model, and the chains
that relied on those sub-cutoff reactions can then carry no flux at all -- a
core reaction included. On Recon3D this is measurable rather than theoretical;
see the package documentation.

The mixed-integer formulations should not have this problem, since a binary at
zero forces its reaction's flux to exactly zero and so discards nothing.

:func:`check_extraction` reports on a result so the difference is visible
rather than silent.
"""

from logging import getLogger
from typing import TYPE_CHECKING, Iterable, List, Optional

from cobra.flux_analysis import find_blocked_reactions

from ._orientation import STOICHIOMETRY
from .consistency import consistent_reaction_ids
from .exceptions import SpectraError

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

    See Also
    --------
    is_valid : Whether both promises were kept.
    is_consistent : Whether nothing in the model is blocked.

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
    method: str = "fva",
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
        (default "stoichiometry"). Only used by ``method="spectra"``.
    detection_cutoff : float, optional
        The absolute flux at which a reaction counts as carrying flux
        (default ``0.99 * tol`` for ``method="spectra"``, ``model.tolerance``
        for ``method="fva"``).
    seed : int, optional
        A seed for the consistency check's random coefficients (default None).
    method : {"fva", "spectra"}, optional
        How to decide whether a reaction is blocked (default "fva"). Only
        ``"spectra"`` supports ``consistency_type="topology"``.
        ``"fva"`` maximises and minimises each reaction in turn with
        :func:`cobra.flux_analysis.find_blocked_reactions`, which is slow but
        authoritative. ``"spectra"`` reuses this package's own consistency
        check, which is far quicker but **can report a dead reaction as
        consistent** on a model holding numerically marginal reactions —
        exactly the models an LP extraction produces. Use it only when the
        speed matters more than the answer.

    Returns
    -------
    ExtractionReport
        What the model does and does not deliver.

    Raises
    ------
    SpectraError
        If `method` is not one of the two accepted values.

    Examples
    --------
    >>> report = check_extraction(extracted, core, tol=1e-4)
    >>> report.is_valid
    True

    Notes
    -----
    The default is the slow, trustworthy option on purpose. Measured on
    Recon3D extractions, the ``"spectra"`` method reported no blocked core
    reactions where FVA found four, so a fast check here would hand back
    false assurance about the very thing being checked.

    """
    if method not in ("fva", "spectra"):
        raise SpectraError(f'method must be "fva" or "spectra", not {method!r}.')
    if method == "fva" and consistency_type != STOICHIOMETRY:
        # find_blocked_reactions always assumes a steady state, so it cannot
        # answer the question the accumulation condition asks. Saying so beats
        # quietly checking something other than what was requested.
        raise SpectraError(
            f'method="fva" cannot check consistency_type={consistency_type!r}, '
            f"because cobrapy's find_blocked_reactions always assumes a steady "
            f'state. Pass method="spectra" to check the accumulation condition.'
        )

    core_ids = [c if isinstance(c, str) else c.id for c in core_reactions]
    present = {rxn.id for rxn in extracted.reactions}

    missing_core = [c for c in core_ids if c not in present]

    if method == "fva":
        blocked = find_blocked_reactions(extracted, zero_cutoff=detection_cutoff)
    else:
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
