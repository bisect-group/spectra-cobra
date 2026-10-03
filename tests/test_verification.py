"""Test the extraction verification helper."""

from cobra import Model

from spectra_cobra import ExtractionReport, check_extraction, spectra_me


def test_reports_a_good_extraction(toy_model: Model) -> None:
    """A sound extraction is reported as valid."""
    extracted = spectra_me(toy_model, ["R3"], seed=0)
    report = check_extraction(extracted, ["R3"])

    assert report.is_valid
    assert report.is_consistent
    assert report.missing_core == []
    assert report.blocked == []
    assert report.blocked_core == []
    assert report.n_reactions == len(extracted.reactions)


def test_detects_a_missing_core_reaction(toy_model: Model) -> None:
    """A core reaction absent from the model is reported."""
    extracted = spectra_me(toy_model, ["R3"], seed=0)
    # R7 is on a different pathway and is not in this extraction.
    report = check_extraction(extracted, ["R3", "R7"])

    assert not report.is_valid
    assert report.missing_core == ["R7"]
    # Nothing is blocked, though: the model itself is fine.
    assert report.is_consistent


def test_detects_blocked_reactions(blocked_model: Model) -> None:
    """A model carrying a blocked reaction is reported as inconsistent.

    The extraction routines do not produce one on a model this small, so the
    check is given a model that is inconsistent by construction.
    """
    report = check_extraction(blocked_model, ["R1"])

    assert not report.is_valid
    assert not report.is_consistent
    assert set(report.blocked) == {"R2", "R3"}
    # R1 is unblocked, so the core reaction itself is fine here.
    assert report.blocked_core == []


def test_detects_a_blocked_core_reaction(blocked_model: Model) -> None:
    """A core reaction that is present but cannot carry flux is called out.

    This is the damaging case: the reaction is in the model, so its presence
    alone suggests the extraction worked.
    """
    report = check_extraction(blocked_model, ["R1", "R3"])

    assert not report.is_valid
    assert report.missing_core == [], "R3 is present..."
    assert "R3" in report.blocked_core, "...but it cannot carry flux"


def test_accepts_reaction_objects(toy_model: Model) -> None:
    """Core reactions may be given as objects rather than identifiers."""
    extracted = spectra_me(toy_model, ["R3"], seed=0)
    by_id = check_extraction(extracted, ["R3"])
    by_object = check_extraction(extracted, [extracted.reactions.R3])

    assert by_id.missing_core == by_object.missing_core
    assert by_id.is_valid == by_object.is_valid


def test_summary_mentions_the_problem(blocked_model: Model) -> None:
    """The readable summary names what went wrong."""
    report = check_extraction(blocked_model, ["R3"])
    summary = report.summary()

    assert "cannot carry flux" in summary
    assert "CORE" in summary


def test_repr_is_a_one_liner(toy_model: Model) -> None:
    """The repr summarises without newlines."""
    extracted = spectra_me(toy_model, ["R3"], seed=0)
    text = repr(check_extraction(extracted, ["R3"]))

    assert "\n" not in text
    assert "ExtractionReport" in text
    assert "valid=True" in text


def test_topology_mode_is_honoured(blocked_model: Model) -> None:
    """Checking under the accumulation condition unblocks the dead end."""
    steady = check_extraction(
        blocked_model, ["R1"], consistency_type="stoichiometry", method="spectra"
    )
    accumulating = check_extraction(
        blocked_model, ["R1"], consistency_type="topology", method="spectra"
    )

    assert not steady.is_consistent
    assert accumulating.is_consistent


def test_fva_refuses_topology(blocked_model: Model) -> None:
    """FVA cannot answer the accumulation question, and says so.

    cobrapy's find_blocked_reactions always assumes a steady state, so
    silently accepting topology would check something other than what was
    asked for.
    """
    import pytest

    from spectra_cobra import SpectraError

    with pytest.raises(SpectraError, match="topology"):
        check_extraction(blocked_model, ["R1"], consistency_type="topology")


def test_report_is_constructible_directly() -> None:
    """The report can be built without running a check."""
    report = ExtractionReport(
        n_reactions=5, missing_core=[], blocked=["a"], blocked_core=["a"]
    )
    assert not report.is_valid
    assert not report.is_consistent
    assert report.blocked_core == ["a"]


def test_rejects_an_unknown_method(toy_model: Model) -> None:
    """An unrecognised check method is an error, not a silent fallback."""
    import pytest

    from spectra_cobra import SpectraError

    extracted = spectra_me(toy_model, ["R3"], seed=0)
    with pytest.raises(SpectraError, match="method"):
        check_extraction(extracted, ["R3"], method="guess")


def test_both_methods_agree_on_a_clean_model(toy_model: Model) -> None:
    """On a well-behaved model FVA and the fast check give the same verdict."""
    extracted = spectra_me(toy_model, ["R3"], seed=0)
    by_fva = check_extraction(extracted, ["R3"], method="fva")
    by_spectra = check_extraction(extracted, ["R3"], method="spectra")

    assert by_fva.is_valid == by_spectra.is_valid
    assert set(by_fva.blocked) == set(by_spectra.blocked)


def test_fva_is_the_default(blocked_model: Model) -> None:
    """The default method is the authoritative one."""
    default = check_extraction(blocked_model, ["R1"])
    by_fva = check_extraction(blocked_model, ["R1"], method="fva")

    assert set(default.blocked) == set(by_fva.blocked)
