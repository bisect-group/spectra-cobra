"""Copying a model when the solver will not come with it.

:meth:`cobra.Model.copy` deep-copies the attached solver, and optlang
implements that for some solvers by writing an LP file and reading it
back. That round trip is not lossless. An identifier that happens to be a
keyword of the LP format -- ``St``, which begins the constraint section,
and which the CarveMe universal reconstruction uses for sulfur diffusion
-- produces a file the solver then refuses to parse, and the copy fails
with an error that says nothing about why.

The network itself is fine; only the solver's own serialisation is not.
So on failure the model is rebuilt from its dictionary form, which
carries reactions, metabolites, genes, bounds, GPRs and annotations but
no solver, and the solver is attached afterwards from the interface the
original was using.
"""

from logging import getLogger
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from cobra.core import Model


logger = getLogger(__name__)


def copy_model(model: "Model") -> "Model":
    """Return a copy of a model, with or without the solver's help.

    Parameters
    ----------
    model : cobra.Model
        The model to copy.

    Returns
    -------
    cobra.Model
        An independent copy, using the same solver interface.

    Notes
    -----
    Tries cobrapy's own copy first, since it is faster and carries
    everything. The fallback is only reached when the solver cannot
    serialise itself, and it is reported when it happens so that the
    cause is visible rather than merely survived.

    """
    try:
        return model.copy()
    except Exception as error:  # noqa: BLE001 - whatever the solver raises
        from cobra.io.dict import model_from_dict, model_to_dict

        logger.warning(
            "%s could not copy its solver (%s: %s), so the model was rebuilt "
            "without it. The usual cause is an identifier that is a keyword "
            "of the LP file format, such as a reaction named 'St'.",
            type(model.solver).__module__,
            type(error).__name__,
            error,
        )
        copied = model_from_dict(model_to_dict(model))
        copied.solver = model.problem
        copied.tolerance = model.tolerance
        return copied
