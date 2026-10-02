"""Exceptions raised by SPECTRA."""


class SpectraError(Exception):
    """Base class for every error raised by SPECTRA."""


class SpectraSolverError(SpectraError):
    """Raised when an LP or MILP does not reach an acceptable status.

    The MATLAB implementation warns and returns an empty or ``NaN`` solution
    in this situation, which the callers then misread: ``abs(NaN) >= tol`` is
    false, so a loop waiting for reactions to be explained never terminates.
    Raising instead surfaces the failure where it happens.
    """

    def __init__(self, message: str, status: str = None) -> None:
        """Record the message and the offending solver status.

        Parameters
        ----------
        message : str
            What was being solved and why the status is a problem.
        status : str, optional
            The status the solver reported (default None).

        """
        super().__init__(message)
        self.status = status


class SpectraInfeasibleCoreError(SpectraError):
    """Raised when the core reactions cannot all carry flux at once.

    `spectra_me` requires every core reaction to be able to carry at least
    `tol` of flux simultaneously. If the core set is mutually inconsistent, or
    if any core reaction is blocked in the input model, no such flux
    distribution exists. Use `spectra_ccme` for a model that is not known to
    be flux consistent, as it drops the blocked core reactions and reports
    them instead of failing.
    """
