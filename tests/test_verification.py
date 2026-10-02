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
    steady = check_extraction(blocked_model, ["R1"], consistency_type="stoichiometry")
    accumulating = check_extraction(blocked_model, ["R1"], consistency_type="topology")

    assert not steady.is_consistent
    assert accumulating.is_consistent


def test_report_is_constructible_directly() -> None:
    """The report can be built without running a check."""
    report = ExtractionReport(
        n_reactions=5, missing_core=[], blocked=["a"], blocked_core=["a"]
    )
    assert not report.is_valid
    assert not report.is_consistent
    assert report.blocked_core == ["a"]
