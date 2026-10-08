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
from typing import TYPE_CHECKING, Set

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


def subset_model(
    model: "Model",
    keep_ids: Set[str],
    remove_genes: bool = True,
    keep_orphan_metabolites: bool = False,
) -> "Model":
    """Return a model holding only the given reactions.

    Parameters
    ----------
    model : cobra.Model
        The model to take the subset of. It is left unchanged.
    keep_ids : set of str
        The identifiers of the reactions to keep.
    remove_genes : bool, optional
        Whether to drop the genes left without a reaction (default True).
    keep_orphan_metabolites : bool, optional
        Whether to keep the metabolites left without a reaction (default
        False). They cost nothing, and a metabolic task that names one cannot
        be set up against a model that has dropped it.

    Returns
    -------
    cobra.Model
        A new model carrying the kept reactions and everything attached to
        them, in the order the source had them.

    Notes
    -----
    The model is *built* from the kept reactions rather than copied whole and
    then cut down. Deleting is quadratic: every variable removal makes
    ``optlang.Container._reindex`` rebuild the whole name-to-index map, so
    discarding 9227 of Recon3D's 10600 reactions costs 22.6 s, of which 15.2 s
    is reindexing alone. Building instead is linear in what is kept -- 0.16 s
    for the same result, and the gap widens with model size: on Harvey (81094
    reactions) keeping 10% takes 1.8 s against more than 600 s. The two agree
    byte for byte, down to element ordering.

    ``deepcopy`` carries every attribute, so nothing has to be enumerated and
    nothing is silently dropped. It needs ``_model`` detached first --
    ``Reaction.__getstate__`` does not null it, unlike ``Object.__getstate__``,
    so a plain deepcopy would drag the entire source model along. The
    originals are restored in a ``finally``.

    Orphaned metabolites always go; orphaned genes only when asked, and since
    a built model never has any they are put back when they are wanted.

    """
    from copy import deepcopy

    from cobra import Model
    from cobra.core.dictlist import DictList
    from cobra.core.group import Group

    keep = set(keep_ids)
    reactions = [rxn for rxn in model.reactions if rxn.id in keep]
    metabolites = {met for rxn in reactions for met in rxn._metabolites}
    genes = {gene for rxn in reactions for gene in rxn._genes}

    detached = [(obj, obj._model) for obj in (*reactions, *metabolites, *genes)]
    for obj, _ in detached:
        obj._model = None
    try:
        copies = deepcopy(reactions)
    finally:
        for obj, owner in detached:
            obj._model = owner

    subset = Model(model.id)
    subset.name = model.name
    subset.compartments = dict(model.compartments)
    for attribute in ("annotation", "notes"):
        try:
            setattr(subset, attribute, dict(getattr(model, attribute)))
        except (AttributeError, TypeError, ValueError):
            pass
    subset.add_reactions(copies)

    if keep_orphan_metabolites:
        absent = [
            met for met in model.metabolites if met.id not in subset.metabolites
        ]
        if absent:
            subset.add_metabolites(deepcopy(absent))

    # ``Reaction.__getstate__`` serialises the GPR to a string, so the rule is
    # re-parsed when the reactions are added and the genes come back as fresh
    # objects carrying nothing but an identifier. Their metadata has to be put
    # back from the source, or a model whose genes are annotated -- iJO1366,
    # for one -- loses every annotation it had.
    for gene in subset.genes:
        try:
            original = model.genes.get_by_id(gene.id)
        except KeyError:
            continue
        gene.name = original.name
        try:
            gene.annotation = dict(original.annotation)
            gene.notes = dict(original.notes)
        except (AttributeError, TypeError, ValueError):
            pass

    if not remove_genes:
        # The kept reactions only reach the genes they mention, so the rest
        # have to be put back by hand to match ``remove_genes=False``.
        missing = [gene for gene in model.genes if gene.id not in subset.genes]
        if missing:
            for gene in deepcopy(missing):
                gene._model = subset
                subset.genes.append(gene)

    # ``add_reactions`` interns metabolites and genes in the order it first
    # meets them, which is not the order they had in the source model. Restore
    # it, so that an extraction keeping everything is indistinguishable from
    # the original rather than merely equivalent to it.
    for attribute in ("metabolites", "genes"):
        position = {obj.id: i for i, obj in enumerate(getattr(model, attribute))}
        current = list(getattr(subset, attribute))
        if all(obj.id in position for obj in current):
            setattr(
                subset,
                attribute,
                DictList(sorted(current, key=lambda obj: position[obj.id])),
            )

    # Groups survive the loss of all their members, as they would have done
    # when the reactions were deleted out from under them.
    groups = []
    for group in getattr(model, "groups", []):
        members = []
        for member in group.members:
            for container in (
                subset.reactions,
                subset.metabolites,
                subset.genes,
            ):
                if member.id in container:
                    members.append(container.get_by_id(member.id))
                    break
        groups.append(
            Group(group.id, name=group.name, members=members, kind=group.kind)
        )
    if groups:
        subset.add_groups(groups)

    # Nothing above reaches model-level state, so it is carried over by hand.
    try:
        subset.solver = model.solver.interface
    except (AttributeError, TypeError, ValueError):
        pass
    try:
        subset.tolerance = model.tolerance
    except (AttributeError, TypeError, ValueError):
        pass
    coefficients = {
        subset.reactions.get_by_id(rxn.id): rxn.objective_coefficient
        for rxn in reactions
        if rxn.objective_coefficient
    }
    if coefficients:
        subset.objective = coefficients
    # Unconditional: with every objective reaction gone the coefficients are
    # empty but the direction still has to follow the original.
    try:
        subset.objective.direction = model.objective.direction
    except (AttributeError, TypeError, ValueError):
        pass
    return subset
