"""The network inference formulations SPECTRA extracts models with.

Each function here takes a model, a per-reaction direction read off the
Stage-1 core direction discovery phase, and per-reaction (or per-group)
weights, and returns the reactions that belong in the extracted model. They
differ in what they optimise:

==================  ==========================================  ==============
Formulation         Objective                                   Weights
==================  ==========================================  ==============
``min_net_lp``      minimise weighted total absolute flux       non-negative
``min_net_milp``    minimise weighted reaction count            non-negative
``trade_off``       maximise net weight of included reactions   any real
``optim_biomass``   maximise objective minus weighted flux      non-negative
==================  ==========================================  ==============

A reaction's *direction* is 1 if it must carry positive oriented flux, -1 if
it must carry negative oriented flux, and 0 if it is free to be dropped. The
directed reactions are the core set; the rest are what the formulations choose
between.
"""

import math
from dataclasses import dataclass
from logging import getLogger
from typing import (
    TYPE_CHECKING,
    Dict,
    Iterable,
    List,
    Mapping,
    Optional,
    Sequence,
    Set,
    Tuple,
    Union,
)

import numpy as np
import scipy.sparse as sp
from cobra.util.solver import linear_reaction_coefficients
from optlang.interface import FEASIBLE, OPTIMAL, TIME_LIMIT
from optlang.symbolics import Zero

from ._copy import subset_model
from ._orientation import (
    STOICHIOMETRY,
    TOPOLOGY,
    apply_direction_bounds,
    is_reversible,
    oriented_bounds,
    reaction_signs,
    relaxed_mass_balance,
)
from .exceptions import SpectraError, SpectraSolverError, SpectraValueError

if TYPE_CHECKING:
    from cobra.core import Model


logger = getLogger(__name__)

#: Statuses a MILP may stop at and still carry a usable solution: optimal,
#: or a feasible solution found before the time limit.
ACCEPTABLE_MILP_STATUSES = (OPTIMAL, FEASIBLE, TIME_LIMIT)

#: The fraction of `tol` below which a flux counts as zero when reading a
#: mixed-integer solution's continuous flux support.
INCLUSION_CUTOFF_FACTOR = 1e-7

#: The integrality tolerance the mixed-integer formulations ask for.
MILP_INTEGRALITY_TOLERANCE = 1e-9

#: The smallest tolerance the solvers accept.
SOLVER_TOLERANCE_FLOOR = 1e-9

#: The solver tolerance the top-level entry points ask for.
DEFAULT_SOLVER_TOLERANCE = 1e-9

#: Continuation schedule for Split-Variable minNetLP IRLS refinement.
IRLS_DELTAS = (1e-3, 1e-5, 1e-6, 1e-8)


@dataclass(frozen=True)
class _LpToleranceConfig:
    """Derived numerical tolerances for LP formulations and support checks."""

    tol: float
    eff_solver_tol: float
    supp_tol: float
    dead_end_zero_tol: float
    rescue_tol: float
    primal_res_tol: float
    cert_unblocked_tol: float
    weight_cap: float
    lock_tol: float
    user_set_cutoff: bool


def _validate_positive_finite(value: float, name: str) -> float:
    """Validate that `value` is a positive finite real scalar."""
    if isinstance(value, bool) or not isinstance(value, (int, float, np.number)):
        raise SpectraValueError(f"{name} must be a positive finite scalar, got {value!r}.")
    val = float(value)
    if not math.isfinite(val) or val <= 0.0:
        raise SpectraValueError(f"{name} must be a positive finite scalar, got {value!r}.")
    return val


def _validate_weights_map(
    weights: Mapping[str, float],
    required_keys: Iterable[str],
    allow_negative: bool,
    formulation_name: str,
    check_all_entries: bool = True,
) -> None:
    """Validate weight finiteness, required keys, and sign constraints."""
    req = list(required_keys)
    missing = sorted(key for key in req if key not in weights)
    if missing:
        raise SpectraValueError(
            f"{len(missing)} indicator group(s) or reaction(s) have no weight "
            f"for {formulation_name}, for example {missing[:5]}."
        )
    keys_to_check = list(weights.keys()) if check_all_entries else req
    for key in keys_to_check:
        raw = weights[key]
        if isinstance(raw, bool) or not isinstance(raw, (int, float, np.number)):
            raise SpectraValueError(
                f"Weight for {key!r} in {formulation_name} must be a finite numeric value."
            )
        val = float(raw)
        if not math.isfinite(val):
            raise SpectraValueError(
                f"Weight for {key!r} in {formulation_name} must be finite, got {raw!r}."
            )
        if not allow_negative and val < 0.0:
            raise SpectraValueError(
                f"weights must be non-negative for {formulation_name!r} "
                f"(got {val:g} on {key!r}); use 'tradeOff' for signed weights."
            )


def _resolve_lp_tolerances(
    model: "Model",
    tol: float,
    feasibility_tolerance: Optional[float] = None,
    inclusion_cutoff: Optional[float] = None,
) -> _LpToleranceConfig:
    """Derive scale-aware tolerances from `(tol, feasibility_tolerance, inclusion_cutoff)`."""
    tol_val = _validate_positive_finite(tol, "tol")

    if feasibility_tolerance is None:
        model_tol = getattr(model, "tolerance", None)
        if (
            model_tol is not None
            and not isinstance(model_tol, bool)
            and isinstance(model_tol, (int, float, np.number))
            and math.isfinite(float(model_tol))
            and float(model_tol) > 0.0
        ):
            eff_solver_tol = max(SOLVER_TOLERANCE_FLOOR, float(model_tol))
        else:
            eff_solver_tol = SOLVER_TOLERANCE_FLOOR
    else:
        ft_val = _validate_positive_finite(feasibility_tolerance, "feasibility_tolerance")
        eff_solver_tol = max(SOLVER_TOLERANCE_FLOOR, ft_val)

    if inclusion_cutoff is not None:
        supp_tol = _validate_positive_finite(inclusion_cutoff, "inclusion_cutoff")
        user_set_cutoff = True
    else:
        supp_tol = min(0.1 * tol_val, max(eff_solver_tol, 1e-5 * tol_val))
        user_set_cutoff = False

    dead_end_zero_tol = min(eff_solver_tol, supp_tol)
    rescue_tol = 1e-2 * min(eff_solver_tol, supp_tol)
    primal_res_tol = min(0.1 * tol_val, 1e2 * eff_solver_tol)
    if user_set_cutoff:
        cert_unblocked_tol = min(tol_val, max(supp_tol, eff_solver_tol))
    else:
        cert_unblocked_tol = min(tol_val, max(1e2 * eff_solver_tol, supp_tol))
    weight_cap = max(1e4, 1.0 / max(eff_solver_tol, 1e-4 * tol_val))
    lock_tol = max(0.99 * tol_val, tol_val - 10.0 * eff_solver_tol)

    return _LpToleranceConfig(
        tol=tol_val,
        eff_solver_tol=eff_solver_tol,
        supp_tol=supp_tol,
        dead_end_zero_tol=dead_end_zero_tol,
        rescue_tol=rescue_tol,
        primal_res_tol=primal_res_tol,
        cert_unblocked_tol=cert_unblocked_tol,
        weight_cap=weight_cap,
        lock_tol=lock_tol,
        user_set_cutoff=user_set_cutoff,
    )


def _build_network_arrays(
    model: "Model",
) -> Tuple[sp.csc_matrix, np.ndarray, np.ndarray, List[str], bool, List["object"]]:
    """Extract sparse stoichiometry, bounds, and extra constraints in O(nnz)."""
    met_index = {met: idx for idx, met in enumerate(model.metabolites)}
    met_ids = {met.id for met in model.metabolites}
    n_mets = len(met_index)
    n_rxns = len(model.reactions)

    rows: List[int] = []
    cols: List[int] = []
    vals: List[float] = []
    lb = np.empty(n_rxns, dtype=float)
    ub = np.empty(n_rxns, dtype=float)
    rxn_ids: List[str] = []

    for j, rxn in enumerate(model.reactions):
        rxn_ids.append(rxn.id)
        lb[j] = float(rxn.lower_bound)
        ub[j] = float(rxn.upper_bound)
        for met, coeff in rxn._metabolites.items():
            if coeff != 0.0:
                rows.append(met_index[met])
                cols.append(j)
                vals.append(float(coeff))

    s_csc = sp.csc_matrix((vals, (rows, cols)), shape=(n_mets, n_rxns), dtype=float)

    extra_cons = []
    has_inhomog_met = False
    for cons in model.constraints:
        if cons.name in met_ids:
            if (cons.lb is not None and cons.lb != 0.0) or (
                cons.ub is not None and cons.ub != 0.0
            ):
                has_inhomog_met = True
        else:
            extra_cons.append(cons)

    has_coupling = bool(extra_cons) or has_inhomog_met
    return s_csc, lb, ub, rxn_ids, has_coupling, extra_cons


def _zero_flux_is_feasible(model: "Model", steady_state: bool = True) -> bool:
    """Return whether the zero flux vector `v = 0` satisfies all bounds and constraints."""
    for rxn in model.reactions:
        if rxn.lower_bound > 0.0 or rxn.upper_bound < 0.0:
            return False

    met_ids = {met.id for met in model.metabolites}
    for cons in model.constraints:
        if cons.lb is not None and cons.lb > 0.0:
            return False
        if cons.name in met_ids and not steady_state:
            continue
        if cons.ub is not None and cons.ub < 0.0:
            return False
    return True


def _read_flux_vector(model: "Model") -> np.ndarray:
    """Read the current reaction fluxes in model order as a NumPy vector."""
    n_rxns = len(model.reactions)
    if n_rxns == 0:
        return np.zeros(0, dtype=float)
    try:
        primals = model.solver.primal_values
        return np.array(
            [
                float(primals[rxn.forward_variable.name])
                - float(primals[rxn.reverse_variable.name])
                for rxn in model.reactions
            ],
            dtype=float,
        )
    except Exception:  # noqa: BLE001 - fallback for custom/stub solver interfaces
        return np.array([float(rxn.flux) for rxn in model.reactions], dtype=float)


def _check_row_residual(s_mat: sp.spmatrix, v: np.ndarray, steady_state: bool) -> float:
    """Compute the maximum stoichiometric or topological row residual."""
    if s_mat.shape[0] == 0 or s_mat.shape[1] == 0 or v.size == 0:
        return 0.0
    sv = s_mat.dot(v)
    if sv.size == 0:
        return 0.0
    if steady_state:
        return float(np.max(np.abs(sv)))
    return float(np.max(np.maximum(0.0, -sv)))


def _degree2_chain_multipliers(
    s_csc: sp.csc_matrix, open_mask: np.ndarray, paid_non_core: np.ndarray
) -> np.ndarray:
    """Count paid non-core reactions per degree-2 linear chain in 0 LPs."""
    n_mets, n_rxns = s_csc.shape
    if n_rxns == 0 or n_mets == 0 or not np.any(paid_non_core):
        return np.ones(n_rxns, dtype=float)

    open_indices = np.flatnonzero(open_mask)
    if open_indices.size == 0:
        return np.ones(n_rxns, dtype=float)

    s_csr = s_csc[:, open_indices].tocsr()
    s_csr.eliminate_zeros()
    row_nnz = np.diff(s_csr.indptr)
    deg2_rows = np.flatnonzero(row_nnz == 2)
    if deg2_rows.size == 0:
        return np.ones(n_rxns, dtype=float)

    parent = np.arange(n_rxns, dtype=np.int64)
    indptr = s_csr.indptr
    indices = s_csr.indices
    for r in deg2_rows:
        start = indptr[r]
        if indptr[r + 1] - start != 2:
            continue
        ri = int(open_indices[indices[start]])
        rj = int(open_indices[indices[start + 1]])
        while parent[ri] != ri:
            parent[ri] = parent[parent[ri]]
            ri = int(parent[ri])
        while parent[rj] != rj:
            parent[rj] = parent[parent[rj]]
            rj = int(parent[rj])
        if ri != rj:
            parent[ri] = rj

    roots = np.empty(n_rxns, dtype=np.int64)
    for i in range(n_rxns):
        r = i
        while parent[r] != r:
            r = int(parent[r])
        curr = i
        while parent[curr] != r:
            nxt = int(parent[curr])
            parent[curr] = r
            curr = nxt
        roots[i] = r

    comp_paid_counts = np.bincount(roots[paid_non_core], minlength=n_rxns)
    multipliers = np.ones(n_rxns, dtype=float)
    multipliers[paid_non_core] = np.maximum(
        1.0, comp_paid_counts[roots[paid_non_core]].astype(float)
    )
    return multipliers


def _prepare_merged_weights(
    s_csc: sp.csc_matrix,
    lb: np.ndarray,
    ub: np.ndarray,
    core_pinned_mask: np.ndarray,
    raw_weights: np.ndarray,
    steady_state: bool,
) -> np.ndarray:
    """Compute 'With-Chain' median-normalized non-negative weights in 0 LPs."""
    n = raw_weights.size
    if n == 0:
        return np.zeros(0, dtype=float)

    open_mask = (np.abs(lb) > 0.0) | (np.abs(ub) > 0.0)
    user_w = np.asarray(raw_weights, dtype=float)
    paid_non_core = open_mask & (~core_pinned_mask) & (user_w > 0.0)
    u_pos = user_w[paid_non_core]
    u_med = float(np.median(u_pos)) if u_pos.size > 0 else 1.0

    u_norm = np.zeros(n, dtype=float)
    pos_all = open_mask & (user_w > 0.0)
    if np.any(pos_all):
        u_norm[pos_all] = np.clip(user_w[pos_all] / u_med, 1e-4, 1e4)

    m_chain = (
        _degree2_chain_multipliers(s_csc, open_mask, paid_non_core)
        if steady_state
        else np.ones(n, dtype=float)
    )
    return u_norm * m_chain


def _strip_dead_ends(
    s_csc: sp.csc_matrix,
    lb: np.ndarray,
    ub: np.ndarray,
    keep_mask: np.ndarray,
    core_mask: np.ndarray,
    steady_state: bool,
    dead_end_zero_tol: float,
) -> np.ndarray:
    """Iteratively prune structural dead-end non-core reactions in 0 LPs."""
    active = np.array(keep_mask, dtype=bool, copy=True)
    core = np.asarray(core_mask, dtype=bool)
    n_rows = s_csc.shape[0]
    if n_rows == 0:
        return active

    while True:
        cols = np.flatnonzero(active)
        if cols.size == 0:
            break
        s_sub = s_csc[:, cols].tocoo()
        if s_sub.nnz == 0:
            break
        nz = s_sub.data != 0.0
        if not np.all(nz):
            r_idx = s_sub.row[nz]
            c_idx = s_sub.col[nz]
            vals = s_sub.data[nz]
        else:
            r_idx = s_sub.row
            c_idx = s_sub.col
            vals = s_sub.data
        if vals.size == 0:
            break

        row_nnz = np.bincount(r_idx, minlength=n_rows)
        lbc = lb[cols]
        ubc = ub[cols]

        can_pos = ((vals > 0.0) & (ubc[c_idx] >= dead_end_zero_tol)) | (
            (vals < 0.0) & (lbc[c_idx] <= -dead_end_zero_tol)
        )
        can_neg = ((vals > 0.0) & (lbc[c_idx] <= -dead_end_zero_tol)) | (
            (vals < 0.0) & (ubc[c_idx] >= dead_end_zero_tol)
        )

        pos_cnt = np.bincount(r_idx, weights=can_pos.astype(np.int64), minlength=n_rows)
        neg_cnt = np.bincount(r_idx, weights=can_neg.astype(np.int64), minlength=n_rows)

        act_rows = row_nnz > 0
        if steady_state:
            dead_rows = act_rows & ((row_nnz == 1) | (pos_cnt == 0) | (neg_cnt == 0))
        else:
            dead_rows = act_rows & (pos_cnt == 0) & (neg_cnt > 0)

        if not np.any(dead_rows):
            break

        dead_local_cols = np.unique(c_idx[dead_rows[r_idx]])
        dead_global = cols[dead_local_cols]
        dead_non_core = dead_global[~core[dead_global]]
        if dead_non_core.size == 0:
            break
        active[dead_non_core] = False

    return active


def _extra_constraints_satisfied(
    model: "Model",
    extra_cons: Sequence["object"],
    cols: np.ndarray,
    scaled_u: np.ndarray,
    primal_res_tol: float,
) -> bool:
    """Check whether `scaled_u` on `cols` (0 elsewhere) satisfies extra constraints."""
    if not extra_cons:
        return True

    var_values: Dict["object", float] = {}
    for rxn in model.reactions:
        var_values[rxn.forward_variable] = 0.0
        var_values[rxn.reverse_variable] = 0.0
    for col_idx, flux_val in zip(cols, scaled_u):
        rxn = model.reactions[int(col_idx)]
        val = float(flux_val)
        var_values[rxn.forward_variable] = max(val, 0.0)
        var_values[rxn.reverse_variable] = max(-val, 0.0)

    for cons in extra_cons:
        coefs = cons.get_linear_coefficients(cons.variables)
        cu = 0.0
        for var, coef in coefs.items():
            if coef == 0.0:
                continue
            if var not in var_values:
                return False
            cu += float(coef) * var_values[var]
        if cons.ub is not None and cu - float(cons.ub) > primal_res_tol:
            return False
        if cons.lb is not None and float(cons.lb) - cu > primal_res_tol:
            return False
    return True


def _project_extra_constraints(
    sub_mod: "Model", extra_cons: Sequence["object"]
) -> None:
    """Project extra linear constraints onto `sub_mod` (with dropped reactions at 0)."""
    if not extra_cons:
        return
    sub_vars = sub_mod.variables
    prob = sub_mod.problem
    new_cons = []
    new_coefs = []
    for cons in extra_cons:
        raw_coefs = cons.get_linear_coefficients(cons.variables)
        mapped = {
            sub_vars[var.name]: float(coef)
            for var, coef in raw_coefs.items()
            if coef != 0.0 and var.name in sub_vars
        }
        c_copy = prob.Constraint(Zero, name=cons.name, lb=cons.lb, ub=cons.ub)
        new_cons.append(c_copy)
        if mapped:
            new_coefs.append((c_copy, mapped))
    if new_cons:
        sub_mod.add_cons_vars(new_cons)
        sub_mod.solver.update()
        for c_copy, mapped in new_coefs:
            c_copy.set_linear_coefficients(mapped)


def _certify_support(
    model: "Model",
    s_csc: sp.csc_matrix,
    lb: np.ndarray,
    ub: np.ndarray,
    rxn_ids: Sequence[str],
    v_sol: np.ndarray,
    core_mask: np.ndarray,
    steady_state: bool,
    has_coupling: bool,
    extra_cons: Sequence["object"],
    feasibility_tolerance: Optional[float],
    optimality_tolerance: Optional[float],
    tol_cfg: _LpToleranceConfig,
) -> Tuple[np.ndarray, int]:
    """Certify that the LP support is free of solver-noise dead ends."""
    extra_lps = 0
    supp_tol = tol_cfg.supp_tol
    rescue_tol = tol_cfg.rescue_tol
    res_limit = 0.5 * tol_cfg.primal_res_tol
    unblock_tol = tol_cfg.cert_unblocked_tol
    zero_tol = tol_cfg.dead_end_zero_tol
    tol = tol_cfg.tol

    keep_mask = core_mask | (np.abs(v_sol) >= supp_tol)
    if not has_coupling:
        keep_mask = _strip_dead_ends(
            s_csc, lb, ub, keep_mask, core_mask, steady_state, zero_tol
        )
    cols = np.flatnonzero(keep_mask)
    if cols.size == 0:
        return keep_mask, extra_lps

    u = v_sol[cols]
    res_inf = _check_row_residual(s_csc[:, cols], u, steady_state)
    if not tol_cfg.user_set_cutoff and res_inf > res_limit:
        keep_mask = core_mask | (np.abs(v_sol) >= rescue_tol)
        if not has_coupling:
            keep_mask = _strip_dead_ends(
                s_csc, lb, ub, keep_mask, core_mask, steady_state, zero_tol
            )
        cols = np.flatnonzero(keep_mask)
        if cols.size == 0:
            return keep_mask, extra_lps
        u = v_sol[cols]
        res_inf = _check_row_residual(s_csc[:, cols], u, steady_state)

    # Zero-LP Conic Ray-Scaling Certificate (valid when homogeneous mass balance holds)
    has_inhomog_met = has_coupling and not extra_cons
    if not has_inhomog_met and res_inf <= res_limit:
        pos = u > 0.0
        neg = u < 0.0
        eps = float(np.finfo(float).eps)
        max_scale = max(50.0, tol / max(supp_tol, eps))
        alpha_ub = float(np.min(ub[cols[pos]] / u[pos])) if np.any(pos) else max_scale
        alpha_lb = float(np.min(lb[cols[neg]] / u[neg])) if np.any(neg) else max_scale
        alpha_res = (res_limit / res_inf) if res_inf > 0.0 else max_scale
        alpha = max(1.0, min(max_scale, alpha_ub, alpha_lb, alpha_res))

        cpl_ok = _extra_constraints_satisfied(
            model, extra_cons, cols, alpha * u, tol_cfg.primal_res_tol
        )
        if cpl_ok:
            u_scaled = np.abs(u) * alpha
            core_loc = core_mask[cols]
            if np.all(u_scaled[core_loc] >= 0.99 * tol) and np.all(
                u_scaled[~core_loc] >= unblock_tol
            ):
                return keep_mask, extra_lps

    from .consistency import consistent_reaction_ids

    sub_ids = [rxn_ids[int(i)] for i in cols]
    sub_mod = subset_model(model, set(sub_ids), remove_genes=True)
    if extra_cons:
        _project_extra_constraints(sub_mod, extra_cons)

    c_type = STOICHIOMETRY if steady_state else TOPOLOGY
    cons_ids, cc_lps = consistent_reaction_ids(
        sub_mod,
        tol=unblock_tol,
        consistency_type=c_type,
        feasibility_tolerance=feasibility_tolerance,
        optimality_tolerance=optimality_tolerance,
    )
    extra_lps += cc_lps
    is_cons_sub = np.array([rxn_ids[int(i)] in cons_ids for i in cols], dtype=bool)
    dead_global = cols[(~is_cons_sub) & (~core_mask[cols])]
    keep_mask[dead_global] = False
    return keep_mask, extra_lps


def _certify_milp_support(
    model: "Model",
    s_csc: sp.csc_matrix,
    lb: np.ndarray,
    ub: np.ndarray,
    rxn_ids: Sequence[str],
    raw_included: Set[str],
    v_sol: Optional[np.ndarray],
    core_pinned_mask: np.ndarray,
    has_groups: bool,
    steady_state: bool,
    has_coupling: bool,
    extra_cons: Sequence["object"],
    feasibility_tolerance: Optional[float],
    optimality_tolerance: Optional[float],
    tol_cfg: _LpToleranceConfig,
) -> Tuple[Set[str], int]:
    """Certify MILP support consistency on inconsistent models (spectra_ccme)."""
    extra_lps = 0
    n = len(rxn_ids)
    keep_mask = np.array([rid in raw_included for rid in rxn_ids], dtype=bool)
    has_inhomog_met = has_coupling and not extra_cons

    if not has_groups and v_sol is not None and v_sol.size >= n:
        keep_cand = core_pinned_mask | (np.abs(v_sol[:n]) >= tol_cfg.supp_tol)
        if not has_coupling:
            keep_cand = _strip_dead_ends(
                s_csc,
                lb,
                ub,
                keep_cand,
                core_pinned_mask,
                steady_state,
                tol_cfg.dead_end_zero_tol,
            )
        cols = np.flatnonzero(keep_cand)
        if cols.size == 0:
            return set(), extra_lps

        u = v_sol[cols]
        res_inf = _check_row_residual(s_csc[:, cols], u, steady_state)
        res_limit = 0.5 * tol_cfg.primal_res_tol
        if not has_inhomog_met and res_inf <= res_limit:
            pos = u > 0.0
            neg = u < 0.0
            eps = float(np.finfo(float).eps)
            max_scale = max(50.0, tol_cfg.tol / max(tol_cfg.supp_tol, eps))
            alpha_ub = float(np.min(ub[cols[pos]] / u[pos])) if np.any(pos) else max_scale
            alpha_lb = float(np.min(lb[cols[neg]] / u[neg])) if np.any(neg) else max_scale
            alpha_res = (res_limit / res_inf) if res_inf > 0.0 else max_scale
            alpha = max(1.0, min(max_scale, alpha_ub, alpha_lb, alpha_res))
            cpl_ok = _extra_constraints_satisfied(
                model, extra_cons, cols, alpha * u, tol_cfg.primal_res_tol
            )
            if cpl_ok:
                u_scaled = np.abs(u) * alpha
                core_loc = core_pinned_mask[cols]
                if np.all(u_scaled[core_loc] >= 0.99 * tol_cfg.tol) and np.all(
                    u_scaled[~core_loc] >= tol_cfg.cert_unblocked_tol
                ):
                    return {rxn_ids[int(i)] for i in cols}, extra_lps
        keep_mask = keep_cand

    cols = np.flatnonzero(keep_mask)
    if cols.size == 0:
        return set(), extra_lps

    from .consistency import consistent_reaction_ids

    sub_ids = [rxn_ids[int(i)] for i in cols]
    sub_mod = subset_model(model, set(sub_ids), remove_genes=True)
    if extra_cons:
        _project_extra_constraints(sub_mod, extra_cons)

    c_type = STOICHIOMETRY if steady_state else TOPOLOGY
    cons_ids, cc_lps = consistent_reaction_ids(
        sub_mod,
        tol=tol_cfg.cert_unblocked_tol,
        consistency_type=c_type,
        feasibility_tolerance=feasibility_tolerance,
        optimality_tolerance=optimality_tolerance,
    )
    extra_lps += cc_lps
    is_cons_sub = np.array([rxn_ids[int(i)] in cons_ids for i in cols], dtype=bool)
    keep_mask[cols[(~is_cons_sub) & (~core_pinned_mask[cols])]] = False
    return {rxn_ids[int(i)] for i in np.flatnonzero(keep_mask)}, extra_lps


def _check_status(model: "Model", description: str, milp: bool = False) -> str:
    """Raise unless the solver reached a status carrying a usable solution."""
    status = model.solver.status
    acceptable = ACCEPTABLE_MILP_STATUSES if milp else (OPTIMAL,)
    if status not in acceptable:
        raise SpectraSolverError(
            f"The {description} terminated with status {status!r} instead of "
            f"reaching {'an acceptable status' if milp else 'optimality'}.",
            status=status,
        )
    if status != OPTIMAL:
        if not _has_solution(model):
            raise SpectraSolverError(
                f"The {description} terminated with status {status!r} before "
                f"finding any feasible solution. Allow it more time, or use a "
                f"linear formulation such as minNetLP.",
                status=status,
            )
        logger.warning(
            "The %s stopped at status %r rather than optimality; the returned "
            "model is feasible but may not be minimal.",
            description,
            status,
        )
    return status


def _has_solution(model: "Model") -> bool:
    """Return whether the solver is holding a solution that can be read."""
    try:
        value = model.solver.objective.value
    except Exception:  # noqa: BLE001
        return False
    if value is None or math.isnan(value):
        return False
    try:
        variables = model.solver.variables
        if not len(variables):
            return True
        variables[0].primal
    except AttributeError as error:
        if "primal" not in str(error).lower() and "'X'" not in str(error):
            return True
        return False
    except Exception:  # noqa: BLE001
        return False
    return True


def _included_reactions(
    model: "Model", tol: float, inclusion_cutoff: Optional[float] = None
) -> Set[str]:
    """Return the reactions carrying flux in the current solution.

    Parameters
    ----------
    model : cobra.Model
        A model whose solver holds a solution.
    tol : float
        The flux threshold the formulation was given.
    inclusion_cutoff : float, optional
        Minimum absolute flux to retain a non-core reaction
        (default ``tol * 1e-7``).

    Returns
    -------
    set of str
        The identifiers of the reactions to keep.

    """
    if inclusion_cutoff is None:
        cutoff = tol * INCLUSION_CUTOFF_FACTOR
    else:
        cutoff = _validate_positive_finite(inclusion_cutoff, "inclusion_cutoff")
    included = {rxn.id for rxn in model.reactions if abs(rxn.flux) >= cutoff}
    _warn_on_noisy_support(model, included, tol)
    return included


def _warn_on_noisy_support(model: "Model", included: Set[str], tol: float) -> None:
    """Warn when the solution's support is partly below solver precision."""
    precision = model.tolerance
    if not precision:
        return
    noisy = sum(
        1 for rxn in model.reactions if rxn.id in included and abs(rxn.flux) < precision
    )
    if not noisy:
        return
    logger.warning(
        "%d of %d reactions in this solution carry less flux than the "
        "solver can resolve (model.tolerance=%.3g), so whether they belong "
        "in the model is not something the solution settles. The extracted "
        "model may be unable to carry flux at all, core reactions included. ",
        noisy,
        len(included),
        precision,
    )


def _included_from_indicators(
    directions: Dict[str, int],
    indicators: Dict[str, "object"],
    groups: Mapping[str, Sequence[str]],
) -> Set[str]:
    """Return the reactions a mixed-integer solution keeps from its binaries."""
    kept = {rxn_id for rxn_id, direction in directions.items() if direction != 0}
    for key, var in indicators.items():
        if var.primal > 0.5:
            kept.update(groups[key])
    return kept


@dataclass(frozen=True)
class MilpSolution:
    """What a mixed-integer formulation found.

    Parameters
    ----------
    included : set of str
        The reactions the extracted model should contain.
    selected : set of str
        The **indicator groups** the solver switched on, named by their keys.

    """

    included: Set[str]
    selected: Set[str]


def _milp_solution(
    model: "Model",
    directions: Dict[str, int],
    indicators: Dict[str, "object"],
    groups: Mapping[str, Sequence[str]],
    tol: float,
    inclusion_cutoff: Optional[float] = None,
) -> MilpSolution:
    """Read a mixed-integer solution off both its indicators and its flux."""
    included = _included_from_indicators(directions, indicators, groups)
    carrying = _included_reactions(model, tol, inclusion_cutoff)
    selected = {key for key, var in indicators.items() if var.primal > 0.5}
    return MilpSolution(included=included | carrying, selected=selected)


def _free_reaction_ids(model: "Model", directions: Dict[str, int]) -> List[str]:
    """Return the identifiers of the reactions free to be dropped."""
    return [rxn.id for rxn in model.reactions if directions.get(rxn.id, 0) == 0]


def _resolve_indicator_groups(
    free_ids: Sequence[str],
    indicator_reactions: Optional[Iterable[str]],
    indicator_groups: Optional[Mapping[str, Iterable[str]]],
) -> Tuple[Dict[str, List[str]], Set[str]]:
    """Work out which reactions each binary governs."""
    if indicator_reactions is not None and indicator_groups is not None:
        raise SpectraValueError(
            "Pass indicator_reactions or indicator_groups, not both: one "
            "gives each named reaction its own binary and the other shares "
            "one binary between several, so giving both leaves it ambiguous "
            "which applies."
        )

    order = {rxn_id: index for index, rxn_id in enumerate(free_ids)}
    if indicator_groups is not None:
        groups: Dict[str, List[str]] = {}
        seen: Dict[str, str] = {}
        for key, members in indicator_groups.items():
            wanted = list(dict.fromkeys(members))
            if not wanted:
                raise SpectraValueError(
                    f"Indicator group {key!r} is empty, so its binary would "
                    f"govern nothing while still being counted."
                )
            unknown = [r for r in wanted if r not in order]
            if unknown:
                raise SpectraValueError(
                    f"indicator_groups must name free reactions of the model, "
                    f"but group {key!r} names {len(unknown)} that are not, for "
                    f"example {sorted(unknown)[:5]}. A reaction given a "
                    f"direction is already forced in and cannot also be "
                    f"selected over."
                )
            for rxn_id in wanted:
                if rxn_id in seen:
                    raise SpectraValueError(
                        f"{rxn_id!r} is in both group {seen[rxn_id]!r} and "
                        f"group {key!r}, so which binary governs it is "
                        f"ambiguous. Every reaction belongs to at most one."
                    )
                seen[rxn_id] = key
            groups[key] = sorted(wanted, key=order.__getitem__)
        return groups, set(free_ids) - set(seen)

    if indicator_reactions is None:
        return {rxn_id: [rxn_id] for rxn_id in free_ids}, set()

    wanted_set = set(indicator_reactions)
    unknown_set = wanted_set - set(free_ids)
    if unknown_set:
        raise SpectraValueError(
            f"indicator_reactions must be free reactions of the model, "
            f"but {len(unknown_set)} are not, for example "
            f"{sorted(unknown_set)[:5]}. A reaction given a direction is "
            f"already forced in and cannot also be selected over."
        )
    chosen = [rxn_id for rxn_id in free_ids if rxn_id in wanted_set]
    return {rxn_id: [rxn_id] for rxn_id in chosen}, set(free_ids) - wanted_set


def _absolute_flux_coefficients(
    model: "Model", rxn_ids: Iterable[str], weights: Mapping[str, float]
) -> Dict["object", float]:
    """Objective coefficients summing the weighted absolute flux."""
    coefficients: Dict["object", float] = {}
    for rxn_id in rxn_ids:
        reaction = model.reactions.get_by_id(rxn_id)
        weight = float(weights[rxn_id])
        coefficients[reaction.forward_variable] = weight
        coefficients[reaction.reverse_variable] = weight
    return coefficients


def _add_exclusion_constraints(
    model: "Model",
    indicators: Dict[str, "object"],
    previous_solutions: Optional[List[Set[str]]],
) -> None:
    """Forbid each previously found indicator selection from recurring."""
    if not previous_solutions:
        return

    prob = model.problem
    constraints = []
    coefficients = []
    for index, solution in enumerate(previous_solutions):
        in_solution = [indicators[r] for r in solution if r in indicators]
        if not in_solution:
            continue
        constraint = prob.Constraint(
            Zero,
            name=f"spectra_exclude_{index}",
            ub=len(in_solution) - 1,
        )
        constraints.append(constraint)
        coefficients.append((constraint, {var: 1.0 for var in in_solution}))

    model.add_cons_vars(constraints)
    model.solver.update()
    for constraint, terms in coefficients:
        constraint.set_linear_coefficients(terms)
    constraints.reverse()


def _min_net_lp_core(
    model: "Model",
    directions: Dict[str, int],
    weights: Mapping[str, float],
    tol: float,
    steady_state: bool = True,
    inclusion_cutoff: Optional[float] = None,
    feasibility_tolerance: Optional[float] = None,
    optimality_tolerance: Optional[float] = None,
    run_irls: bool = True,
    pass0_flux_ready: Optional[Union[Mapping[str, float], np.ndarray]] = None,
) -> Tuple[Set[str], Dict[str, float], int]:
    """Run Split-Variable minNetLP (with optional IRLS + support certification)."""
    rxn_list = list(model.reactions)
    rxn_ids = [rxn.id for rxn in rxn_list]
    n = len(rxn_ids)

    _validate_weights_map(
        weights, rxn_ids, allow_negative=False, formulation_name="minNetLP"
    )
    tol_cfg = _resolve_lp_tolerances(
        model, tol, feasibility_tolerance, inclusion_cutoff
    )
    _apply_tolerances(model, feasibility_tolerance, optimality_tolerance)

    core_pinned_mask = np.array(
        [directions.get(rid, 0) != 0 for rid in rxn_ids], dtype=bool
    )

    if run_irls and not np.any(core_pinned_mask) and _zero_flux_is_feasible(model, steady_state):
        return set(), {rid: 0.0 for rid in rxn_ids}, 0

    s_csc, orig_lb, orig_ub, _, has_coupling, extra_cons = _build_network_arrays(model)
    raw_w = np.array([float(weights[rid]) for rid in rxn_ids], dtype=float)
    eff_weights = _prepare_merged_weights(
        s_csc, orig_lb, orig_ub, core_pinned_mask, raw_w, steady_state
    )

    signs = reaction_signs(model)
    free_indices = np.flatnonzero(~core_pinned_mask)
    deltas = IRLS_DELTAS if run_irls else ()
    n_lps = 0

    with model, relaxed_mass_balance(model, steady_state):
        apply_direction_bounds(model, directions, signs, tol_cfg.tol)

        col_lb = np.array([float(r.lower_bound) for r in rxn_list], dtype=float)
        col_ub = np.array([float(r.upper_bound) for r in rxn_list], dtype=float)
        open_free = (~core_pinned_mask) & ((col_lb < 0.0) | (col_ub > 0.0))
        paid_mask = open_free & (eff_weights > 0.0)

        model.objective = model.problem.Objective(Zero, direction="min")

        v_ori: Optional[np.ndarray] = None
        if pass0_flux_ready is not None:
            if isinstance(pass0_flux_ready, Mapping):
                cand_arr = np.array(
                    [float(pass0_flux_ready.get(rid, float("nan"))) for rid in rxn_ids],
                    dtype=float,
                )
            else:
                cand_arr = np.asarray(pass0_flux_ready, dtype=float).reshape(-1)
            if cand_arr.size >= n and not np.any(np.isnan(cand_arr[:n])):
                v_ori = cand_arr[:n].copy()

        if v_ori is None:
            init_coeffs = {}
            for idx in free_indices:
                rxn = rxn_list[int(idx)]
                w_val = float(eff_weights[int(idx)])
                init_coeffs[rxn.forward_variable] = w_val
                init_coeffs[rxn.reverse_variable] = w_val
            if init_coeffs:
                model.objective.set_linear_coefficients(init_coeffs)

            model.slim_optimize()
            n_lps += 1
            _check_status(model, "minNetLP problem")
            v_ori = _read_flux_vector(model)

        if deltas and np.any(paid_mask):
            m_cap = tol_cfg.weight_cap
            supp_tol = tol_cfg.supp_tol
            res_tol = tol_cfg.primal_res_tol
            prev_support = int(np.sum(np.abs(v_ori[paid_mask]) >= supp_tol))

            for delta in deltas:
                c_step = (
                    tol_cfg.tol
                    * eff_weights
                    * np.minimum(m_cap, 1.0 / (np.abs(v_ori) + delta))
                )
                step_coeffs = {}
                for idx in free_indices:
                    rxn = rxn_list[int(idx)]
                    w_val = float(c_step[int(idx)])
                    step_coeffs[rxn.forward_variable] = w_val
                    step_coeffs[rxn.reverse_variable] = w_val
                if step_coeffs:
                    model.objective.set_linear_coefficients(step_coeffs)

                model.slim_optimize()
                n_lps += 1
                if model.solver.status == OPTIMAL:
                    v_cand = _read_flux_vector(model)
                    if _check_row_residual(s_csc, v_cand, steady_state) <= res_tol:
                        v_ori = v_cand
                    else:
                        break
                else:
                    break

                cur_support = int(np.sum(np.abs(v_ori[paid_mask]) >= supp_tol))
                if cur_support == prev_support and delta <= 0.1 * tol_cfg.tol:
                    break
                prev_support = cur_support

    flux_map = {rxn_ids[i]: float(v_ori[i]) for i in range(n)}
    if not run_irls:
        keep_mask = core_pinned_mask | (np.abs(v_ori) >= tol_cfg.supp_tol)
        kept = {rxn_ids[int(i)] for i in np.flatnonzero(keep_mask)}
        return kept, flux_map, n_lps

    keep_mask, cert_lps = _certify_support(
        model=model,
        s_csc=s_csc,
        lb=orig_lb,
        ub=orig_ub,
        rxn_ids=rxn_ids,
        v_sol=v_ori,
        core_mask=core_pinned_mask,
        steady_state=steady_state,
        has_coupling=has_coupling,
        extra_cons=extra_cons,
        feasibility_tolerance=feasibility_tolerance,
        optimality_tolerance=optimality_tolerance,
        tol_cfg=tol_cfg,
    )
    n_lps += cert_lps
    kept = {rxn_ids[int(i)] for i in np.flatnonzero(keep_mask)}
    return kept, flux_map, n_lps


def min_net_lp(
    model: "Model",
    directions: Dict[str, int],
    weights: Dict[str, float],
    tol: float,
    steady_state: bool = True,
    inclusion_cutoff: Optional[float] = None,
    feasibility_tolerance: Optional[float] = None,
    optimality_tolerance: Optional[float] = None,
) -> Set[str]:
    """Extract a model by minimising the weighted total absolute flux.

    Parameters
    ----------
    model : cobra.Model
        The model to extract from.
    directions : dict of {str: int}
        The oriented direction each reaction must carry flux in.
    weights : dict of {str: float}
        Non-negative reaction cost weights. Higher weights penalise inclusion
        more strongly; zero is free; negative weights raise :class:`ValueError`.
    tol : float
        The minimum absolute flux a directed reaction has to carry.
    steady_state : bool, optional
        Whether to enforce ``S v = 0`` rather than ``S v >= 0`` (default True).
    inclusion_cutoff : float, optional
        Minimum absolute flux to retain a non-core reaction (default derived
        from ``(tol, feasibility_tolerance)``).
    feasibility_tolerance : float, optional
        The solver's primal feasibility tolerance, or None to leave its own
        in place (the default).
    optimality_tolerance : float, optional
        The solver's dual feasibility tolerance, or None to leave its own in
        place (the default).

    Returns
    -------
    set of str
        The identifiers of the reactions in the extracted model.

    Raises
    ------
    SpectraValueError
        If `weights` contains negative or non-finite values, or if `tol` or
        `inclusion_cutoff` is not a positive finite scalar.

    Notes
    -----
    Computes "With-Chain" median-normalised non-negative weights (multiplying
    by degree-2 paid non-core chain length in ``steady_state=True`` mode and
    ``1`` in ``steady_state=False`` mode), solves a warm-started Split-Variable
    iteratively reweighted L1 sequence (``deltas = [1e-3, 1e-5, 1e-6, 1e-8]``),
    and runs tolerance-derived support certification.

    """
    kept, _, _ = _min_net_lp_core(
        model=model,
        directions=directions,
        weights=weights,
        tol=tol,
        steady_state=steady_state,
        inclusion_cutoff=inclusion_cutoff,
        feasibility_tolerance=feasibility_tolerance,
        optimality_tolerance=optimality_tolerance,
        run_irls=True,
    )
    return kept


def _min_net_milp_core(
    model: "Model",
    directions: Dict[str, int],
    weights: Mapping[str, float],
    tol: float,
    steady_state: bool = True,
    time_limit: Optional[float] = 7200.0,
    previous_solutions: Optional[List[Set[str]]] = None,
    indicator_reactions: Optional[Iterable[str]] = None,
    inclusion_cutoff: Optional[float] = None,
    indicator_groups: Optional[Mapping[str, Iterable[str]]] = None,
    feasibility_tolerance: Optional[float] = None,
    optimality_tolerance: Optional[float] = None,
    integrality_tolerance: Optional[float] = MILP_INTEGRALITY_TOLERANCE,
    warm_start: Optional[Mapping[str, float]] = None,
) -> Tuple[MilpSolution, np.ndarray]:
    """Run minNetMILP and return both `MilpSolution` and continuous flux vector."""
    tol_val = _validate_positive_finite(tol, "tol")
    if inclusion_cutoff is not None:
        _validate_positive_finite(inclusion_cutoff, "inclusion_cutoff")

    free_ids = _free_reaction_ids(model, directions)
    groups, always_keep = _resolve_indicator_groups(
        free_ids, indicator_reactions, indicator_groups
    )
    _validate_weights_map(
        weights, groups.keys(), allow_negative=False, formulation_name="minNetMILP"
    )
    _apply_tolerances(
        model,
        feasibility_tolerance,
        optimality_tolerance,
        integrality_tolerance,
    )

    has_pinned = any(directions.get(rxn.id, 0) != 0 for rxn in model.reactions)
    if (
        not has_pinned
        and not previous_solutions
        and indicator_reactions is None
        and indicator_groups is None
        and all(float(weights[key]) > 0.0 for key in groups)
        and _zero_flux_is_feasible(model, steady_state)
    ):
        return (
            MilpSolution(included=set(), selected=set()),
            np.zeros(len(model.reactions), dtype=float),
        )

    signs = reaction_signs(model)
    prob = model.problem

    with model, relaxed_mass_balance(model, steady_state):
        apply_direction_bounds(model, directions, signs, tol_val)

        indicators = {}
        to_add = []
        coefficients = []
        for key, members in groups.items():
            var = prob.Variable(f"spectra_z_{key}", type="binary")
            indicators[key] = var
            to_add.append(var)
            for rxn_id in members:
                reaction = model.reactions.get_by_id(rxn_id)
                lower = prob.Constraint(
                    Zero, name=f"spectra_z_lower_{rxn_id}", lb=0.0
                )
                upper = prob.Constraint(
                    Zero, name=f"spectra_z_upper_{rxn_id}", ub=0.0
                )
                to_add.extend([lower, upper])
                forward = reaction.forward_variable
                backward = reaction.reverse_variable
                coefficients.append(
                    (
                        lower,
                        {
                            forward: 1.0,
                            backward: -1.0,
                            var: -float(reaction.lower_bound),
                        },
                    )
                )
                coefficients.append(
                    (
                        upper,
                        {
                            forward: 1.0,
                            backward: -1.0,
                            var: -float(reaction.upper_bound),
                        },
                    )
                )
        model.add_cons_vars(to_add)
        model.solver.update()
        for constraint, terms in coefficients:
            constraint.set_linear_coefficients(terms)
        to_add.reverse()
        _add_exclusion_constraints(model, indicators, previous_solutions)

        model.objective = prob.Objective(Zero, direction="min")
        model.objective.set_linear_coefficients(
            {indicators[key]: float(weights[key]) for key in groups}
        )
        _set_time_limit(model, time_limit)
        _apply_tolerances(
            model,
            feasibility_tolerance,
            optimality_tolerance,
            integrality_tolerance,
        )
        _apply_warm_start(model, warm_start, indicators, groups, tol_val)

        model.slim_optimize()
        _check_status(model, "minNetMILP problem", milp=True)
        v_sol = _read_flux_vector(model)
        solution = _milp_solution(
            model, directions, indicators, groups, tol_val, inclusion_cutoff
        )
        return (
            MilpSolution(
                included=solution.included | always_keep, selected=solution.selected
            ),
            v_sol,
        )


def min_net_milp(
    model: "Model",
    directions: Dict[str, int],
    weights: Dict[str, float],
    tol: float,
    steady_state: bool = True,
    time_limit: Optional[float] = 7200.0,
    previous_solutions: Optional[List[Set[str]]] = None,
    indicator_reactions: Optional[Iterable[str]] = None,
    inclusion_cutoff: Optional[float] = None,
    indicator_groups: Optional[Mapping[str, Iterable[str]]] = None,
    feasibility_tolerance: Optional[float] = None,
    optimality_tolerance: Optional[float] = None,
    integrality_tolerance: Optional[float] = MILP_INTEGRALITY_TOLERANCE,
    warm_start: Optional[Mapping[str, float]] = None,
) -> MilpSolution:
    """Extract a model by minimising the weighted reaction count.

    Parameters
    ----------
    model : cobra.Model
        The model to extract from.
    directions : dict of {str: int}
        The oriented direction each reaction must carry flux in.
    weights : dict of {str: float}
        Non-negative reaction (or group) cost weights. Higher weights make a
        reaction less likely to be kept; negative weights raise
        :class:`ValueError`.
    tol : float
        The minimum absolute flux a directed reaction has to carry.
    steady_state : bool, optional
        Whether to enforce ``S v = 0`` rather than ``S v >= 0`` (default True).
    time_limit : float, optional
        The maximum time to spend solving, in seconds (default 7200). A
        solution found before the limit is used even if not proven optimal.
    previous_solutions : list of set of str, optional
        Indicator sets to exclude, so that a different model is returned
        (default None).
    indicator_reactions : iterable of str, optional
        Restrict the binaries to these reactions; the rest stay continuous
        and are kept whatever they do (default None).
    inclusion_cutoff : float, optional
        Minimum absolute flux to retain a non-core reaction
        (default ``tol * 1e-7``).
    indicator_groups : dict of {str: iterable of str}, optional
        Give the reactions of each group **one shared binary** instead of
        one each. Mutually exclusive with `indicator_reactions`.
    feasibility_tolerance : float, optional
        The solver's primal feasibility tolerance, or None to leave its own
        in place (the default).
    optimality_tolerance : float, optional
        The solver's dual feasibility tolerance, or None to leave its own in
        place (the default).
    integrality_tolerance : float, optional
        How far a binary may sit from a whole number and still count as one
        (default :data:`MILP_INTEGRALITY_TOLERANCE`, 1e-9).
    warm_start : mapping of {str: float}, optional
        Initial feasible flux distribution for MIP warm-starting.

    Returns
    -------
    MilpSolution
        The reactions to keep, and the indicator groups the solver
        switched on.

    """
    sol, _ = _min_net_milp_core(
        model=model,
        directions=directions,
        weights=weights,
        tol=tol,
        steady_state=steady_state,
        time_limit=time_limit,
        previous_solutions=previous_solutions,
        indicator_reactions=indicator_reactions,
        inclusion_cutoff=inclusion_cutoff,
        indicator_groups=indicator_groups,
        feasibility_tolerance=feasibility_tolerance,
        optimality_tolerance=optimality_tolerance,
        integrality_tolerance=integrality_tolerance,
        warm_start=warm_start,
    )
    return sol


def _trade_off_core(
    model: "Model",
    directions: Dict[str, int],
    weights: Mapping[str, float],
    tol: float,
    steady_state: bool = True,
    time_limit: Optional[float] = 7200.0,
    previous_solutions: Optional[List[Set[str]]] = None,
    indicator_reactions: Optional[Iterable[str]] = None,
    inclusion_cutoff: Optional[float] = None,
    indicator_groups: Optional[Mapping[str, Iterable[str]]] = None,
    feasibility_tolerance: Optional[float] = None,
    optimality_tolerance: Optional[float] = None,
    integrality_tolerance: Optional[float] = MILP_INTEGRALITY_TOLERANCE,
    warm_start: Optional[Mapping[str, float]] = None,
) -> Tuple[MilpSolution, np.ndarray]:
    """Run tradeOff and return both `MilpSolution` and continuous flux vector."""
    tol_val = _validate_positive_finite(tol, "tol")
    if inclusion_cutoff is not None:
        _validate_positive_finite(inclusion_cutoff, "inclusion_cutoff")

    free_ids = _free_reaction_ids(model, directions)
    groups, always_keep = _resolve_indicator_groups(
        free_ids, indicator_reactions, indicator_groups
    )
    _validate_weights_map(
        weights, groups.keys(), allow_negative=True, formulation_name="tradeOff"
    )
    _apply_tolerances(
        model,
        feasibility_tolerance,
        optimality_tolerance,
        integrality_tolerance,
    )

    has_pinned = any(directions.get(rxn.id, 0) != 0 for rxn in model.reactions)
    if (
        not has_pinned
        and not previous_solutions
        and indicator_reactions is None
        and indicator_groups is None
        and all(float(weights[key]) < 0.0 for key in groups)
        and _zero_flux_is_feasible(model, steady_state)
    ):
        return (
            MilpSolution(included=set(), selected=set()),
            np.zeros(len(model.reactions), dtype=float),
        )

    signs = reaction_signs(model)
    prob = model.problem

    with model, relaxed_mass_balance(model, steady_state):
        apply_direction_bounds(model, directions, signs, tol_val)

        indicators = {}
        forward_bins = {}
        reverse_bins = {}
        to_add = []
        coefficients = []
        for key, members in groups.items():
            included = prob.Variable(f"spectra_z_{key}", type="binary")
            indicators[key] = included
            to_add.append(included)

            for rxn_id in members:
                reaction = model.reactions.get_by_id(rxn_id)
                sign = signs[rxn_id]
                lower, upper = oriented_bounds(reaction, sign)
                forward = reaction.forward_variable
                backward = reaction.reverse_variable
                oriented = {forward: float(sign), backward: -float(sign)}

                if not is_reversible(reaction, sign):
                    irr_lower = prob.Constraint(
                        Zero, name=f"spectra_irr_lower_{rxn_id}", lb=0.0
                    )
                    irr_upper = prob.Constraint(
                        Zero, name=f"spectra_irr_upper_{rxn_id}", ub=0.0
                    )
                    to_add.extend([irr_lower, irr_upper])
                    coefficients.append(
                        (irr_lower, {**oriented, included: -float(tol_val)})
                    )
                    coefficients.append(
                        (irr_upper, {**oriented, included: -float(upper)})
                    )
                    continue

                forward_on = prob.Variable(f"spectra_a_{rxn_id}", type="binary")
                reverse_on = prob.Variable(f"spectra_b_{rxn_id}", type="binary")
                forward_bins[rxn_id] = forward_on
                reverse_bins[rxn_id] = reverse_on
                split = prob.Constraint(
                    Zero, name=f"spectra_split_{rxn_id}", lb=0.0, ub=0.0
                )
                rev_lower = prob.Constraint(
                    Zero, name=f"spectra_rev_lower_{rxn_id}", lb=0.0
                )
                rev_upper = prob.Constraint(
                    Zero, name=f"spectra_rev_upper_{rxn_id}", ub=0.0
                )
                to_add.extend([forward_on, reverse_on, split, rev_lower, rev_upper])
                coefficients.append(
                    (
                        split,
                        {forward_on: 1.0, reverse_on: 1.0, included: -1.0},
                    )
                )
                coefficients.append(
                    (
                        rev_lower,
                        {
                            **oriented,
                            forward_on: -float(tol_val),
                            reverse_on: -float(lower),
                        },
                    )
                )
                coefficients.append(
                    (
                        rev_upper,
                        {
                            **oriented,
                            forward_on: -float(upper),
                            reverse_on: float(tol_val),
                        },
                    )
                )

        model.add_cons_vars(to_add)
        model.solver.update()
        for constraint, terms in coefficients:
            constraint.set_linear_coefficients(terms)
        to_add.reverse()
        _add_exclusion_constraints(model, indicators, previous_solutions)

        model.objective = prob.Objective(Zero, direction="max")
        model.objective.set_linear_coefficients(
            {indicators[key]: float(weights[key]) for key in groups}
        )
        _set_time_limit(model, time_limit)
        _apply_tolerances(
            model,
            feasibility_tolerance,
            optimality_tolerance,
            integrality_tolerance,
        )
        if indicator_groups is None:
            _apply_trade_off_warm_start(
                model,
                warm_start,
                indicators,
                groups,
                signs,
                forward_bins,
                reverse_bins,
                tol_val,
            )

        model.slim_optimize()
        _check_status(model, "tradeOff problem", milp=True)
        v_sol = _read_flux_vector(model)
        solution = _milp_solution(
            model, directions, indicators, groups, tol_val, inclusion_cutoff
        )
        return (
            MilpSolution(
                included=solution.included | always_keep, selected=solution.selected
            ),
            v_sol,
        )


def trade_off(
    model: "Model",
    directions: Dict[str, int],
    weights: Dict[str, float],
    tol: float,
    steady_state: bool = True,
    time_limit: Optional[float] = 7200.0,
    previous_solutions: Optional[List[Set[str]]] = None,
    indicator_reactions: Optional[Iterable[str]] = None,
    inclusion_cutoff: Optional[float] = None,
    indicator_groups: Optional[Mapping[str, Iterable[str]]] = None,
    feasibility_tolerance: Optional[float] = None,
    optimality_tolerance: Optional[float] = None,
    integrality_tolerance: Optional[float] = MILP_INTEGRALITY_TOLERANCE,
    warm_start: Optional[Mapping[str, float]] = None,
) -> MilpSolution:
    """Extract a model by maximising the net weight of included reactions.

    Parameters
    ----------
    model : cobra.Model
        The model to extract from.
    directions : dict of {str: int}
        The oriented direction each reaction must carry flux in.
    weights : dict of {str: float}
        Signed reaction (or group) evidence weights: **any finite real number**.
        Positive weights reward inclusion, negative weights penalise inclusion,
        and zero is neutral.
    tol : float
        The minimum absolute flux a kept reaction has to carry.
    steady_state : bool, optional
        Whether to enforce ``S v = 0`` rather than ``S v >= 0`` (default True).
    time_limit : float, optional
        The maximum time to spend solving, in seconds (default 7200).
    previous_solutions : list of set of str, optional
        Indicator sets to exclude, so that a different model is returned
        (default None).
    indicator_reactions : iterable of str, optional
        Restrict the binaries to these reactions; the rest stay continuous
        and are kept whatever they do (default None).
    inclusion_cutoff : float, optional
        Minimum absolute flux to retain a non-core reaction
        (default ``tol * 1e-7``).
    indicator_groups : dict of {str: iterable of str}, optional
        Give the reactions of each group **one shared binary**. When provided,
        `weights` looks up each group's score by its group key, and switching a
        group on requires all of its members to carry at least `tol` of flux.
        Mutually exclusive with `indicator_reactions`.
    feasibility_tolerance : float, optional
        The solver's primal feasibility tolerance, or None to leave its own
        in place (the default).
    optimality_tolerance : float, optional
        The solver's dual feasibility tolerance, or None to leave its own in
        place (the default).
    integrality_tolerance : float, optional
        How far a binary may sit from a whole number and still count as one
        (default :data:`MILP_INTEGRALITY_TOLERANCE`, 1e-9).
    warm_start : mapping of {str: float}, optional
        Initial feasible flux distribution for MIP warm-starting.

    Returns
    -------
    MilpSolution
        The reactions to keep, and the indicator groups the solver selected.

    """
    sol, _ = _trade_off_core(
        model=model,
        directions=directions,
        weights=weights,
        tol=tol,
        steady_state=steady_state,
        time_limit=time_limit,
        previous_solutions=previous_solutions,
        indicator_reactions=indicator_reactions,
        inclusion_cutoff=inclusion_cutoff,
        indicator_groups=indicator_groups,
        feasibility_tolerance=feasibility_tolerance,
        optimality_tolerance=optimality_tolerance,
        integrality_tolerance=integrality_tolerance,
        warm_start=warm_start,
    )
    return sol


def _optim_biomass_core(
    model: "Model",
    directions: Dict[str, int],
    weights: Mapping[str, float],
    tol: float,
    steady_state: bool = True,
    inclusion_cutoff: Optional[float] = None,
    feasibility_tolerance: Optional[float] = None,
    optimality_tolerance: Optional[float] = None,
) -> Tuple[Set[str], Dict[str, float], int]:
    """Run optimBiomass and return `(kept_ids, flux_map, n_lps)`."""
    rxn_ids = [rxn.id for rxn in model.reactions]
    n = len(rxn_ids)

    _validate_weights_map(
        weights, rxn_ids, allow_negative=False, formulation_name="optimBiomass"
    )
    tol_cfg = _resolve_lp_tolerances(
        model, tol, feasibility_tolerance, inclusion_cutoff
    )
    _apply_tolerances(model, feasibility_tolerance, optimality_tolerance)

    core_pinned_mask = np.array(
        [directions.get(rid, 0) != 0 for rid in rxn_ids], dtype=bool
    )
    objective_ids = {
        rxn.id for rxn, coef in linear_reaction_coefficients(model).items() if coef
    }

    if (
        not np.any(core_pinned_mask)
        and _zero_flux_is_feasible(model, steady_state)
        and (
            not objective_ids
            or all(float(weights[rid]) == 0.0 for rid in objective_ids)
        )
    ):
        return set(), {rid: 0.0 for rid in rxn_ids}, 0

    s_csc, orig_lb, orig_ub, _, has_coupling, extra_cons = _build_network_arrays(model)
    penalised_ids = [rid for rid in rxn_ids if rid not in objective_ids]
    signs = reaction_signs(model)

    with model, relaxed_mass_balance(model, steady_state):
        apply_direction_bounds(model, directions, signs, tol_cfg.tol)
        model.objective = model.problem.Objective(Zero, direction="max")
        coefficients = {
            variable: -weight
            for variable, weight in _absolute_flux_coefficients(
                model, penalised_ids, weights
            ).items()
        }
        for rxn_id in objective_ids:
            reaction = model.reactions.get_by_id(rxn_id)
            w_obj = float(weights[rxn_id])
            coefficients[reaction.forward_variable] = w_obj
            coefficients[reaction.reverse_variable] = -w_obj
        if coefficients:
            model.objective.set_linear_coefficients(coefficients)

        model.slim_optimize()
        n_lps = 1
        _check_status(model, "optimBiomass problem")
        v_ori = _read_flux_vector(model)

    flux_map = {rxn_ids[i]: float(v_ori[i]) for i in range(n)}
    keep_mask, cert_lps = _certify_support(
        model=model,
        s_csc=s_csc,
        lb=orig_lb,
        ub=orig_ub,
        rxn_ids=rxn_ids,
        v_sol=v_ori,
        core_mask=core_pinned_mask,
        steady_state=steady_state,
        has_coupling=has_coupling,
        extra_cons=extra_cons,
        feasibility_tolerance=feasibility_tolerance,
        optimality_tolerance=optimality_tolerance,
        tol_cfg=tol_cfg,
    )
    n_lps += cert_lps
    kept = {rxn_ids[int(i)] for i in np.flatnonzero(keep_mask)}
    return kept, flux_map, n_lps


def optim_biomass(
    model: "Model",
    directions: Dict[str, int],
    weights: Dict[str, float],
    tol: float,
    steady_state: bool = True,
    inclusion_cutoff: Optional[float] = None,
    feasibility_tolerance: Optional[float] = None,
    optimality_tolerance: Optional[float] = None,
) -> Set[str]:
    """Extract a model by maximising objective flux while penalising total flux.

    Parameters
    ----------
    model : cobra.Model
        The model to extract from. Nonzero linear coefficients in
        ``model.objective`` identify the biomass/objective reactions; if none
        exist, ``optim_biomass`` reduces cleanly to weighted L1 flux
        minimisation.
    directions : dict of {str: int}
        The oriented direction each reaction must carry flux in.
    weights : dict of {str: float}
        Non-negative reaction weights. For objective reactions the weight
        rewards flux; for all other reactions it penalises absolute flux.
        Negative weights raise :class:`ValueError`.
    tol : float
        The minimum absolute flux a directed reaction has to carry.
    steady_state : bool, optional
        Whether to enforce ``S v = 0`` rather than ``S v >= 0`` (default True).
    inclusion_cutoff : float, optional
        Minimum absolute flux to retain a non-core reaction (default derived
        from ``(tol, feasibility_tolerance)``).
    feasibility_tolerance : float, optional
        The solver's primal feasibility tolerance, or None to leave its own
        in place (the default).
    optimality_tolerance : float, optional
        The solver's dual feasibility tolerance, or None to leave its own in
        place (the default).

    Returns
    -------
    set of str
        The identifiers of the reactions in the extracted model.

    Raises
    ------
    SpectraValueError
        If `weights` contains negative or non-finite values, or if `tol` or
        `inclusion_cutoff` is not a positive finite scalar.

    Notes
    -----
    Maximises
    :math:`\\sum_{b \\in c} w_b v_b - \\sum_{i \\notin c} w_i |v_i|`, where
    :math:`c` are the reactions with nonzero linear objective coefficients,
    followed by tolerance-derived support certification.

    """
    kept, _, _ = _optim_biomass_core(
        model=model,
        directions=directions,
        weights=weights,
        tol=tol,
        steady_state=steady_state,
        inclusion_cutoff=inclusion_cutoff,
        feasibility_tolerance=feasibility_tolerance,
        optimality_tolerance=optimality_tolerance,
    )
    return kept


def _apply_tolerances(
    model: "Model",
    feasibility: Optional[float] = None,
    optimality: Optional[float] = None,
    integrality: Optional[float] = None,
) -> None:
    """Set the solver tolerances a formulation asked for."""
    wanted = (
        ("feasibility", feasibility),
        ("optimality", optimality),
        ("integrality", integrality),
    )
    for name, value in wanted:
        if value is None:
            continue
        val = _validate_positive_finite(value, f"{name}_tolerance")
        if val < SOLVER_TOLERANCE_FLOOR:
            logger.warning(
                "The %s tolerance cannot go below %.0e, so %.0e has been "
                "raised to that. Asking for less would leave the solver at "
                "its own tolerance instead, which is looser than either.",
                name,
                SOLVER_TOLERANCE_FLOOR,
                val,
            )
            val = SOLVER_TOLERANCE_FLOOR
        try:
            setattr(model.solver.configuration.tolerances, name, val)
        except Exception as error:  # noqa: BLE001
            logger.warning(
                "Could not set the %s tolerance to %.0e on the %s solver: %s",
                name,
                val,
                model.problem.__name__,
                error,
            )
            continue
        if name == "feasibility":
            model._tolerance = val


def _apply_warm_start(
    model: "Model",
    fluxes: Optional[Mapping[str, float]],
    indicators: Mapping[str, "object"],
    groups: Mapping[str, Iterable[str]],
    tol: float,
) -> bool:
    """Hand the solver a clean feasible starting point for minNetMILP."""
    if not fluxes:
        logger.warning(
            "no warm start: the solver must find its own first feasible "
            "solution, which it may not manage at this size"
        )
        return False

    threshold = min(1e-12, tol * INCLUSION_CUTOFF_FACTOR)
    group_members = {key: list(members) for key, members in groups.items()}
    on = {
        key: any(abs(fluxes.get(rxn_id, 0.0)) > threshold for rxn_id in members)
        for key, members in group_members.items()
    }
    off_rxns = {
        rxn_id
        for key, members in group_members.items()
        if not on[key]
        for rxn_id in members
    }
    try:
        for key, variable in indicators.items():
            variable._internal_variable.setAttr("Start", 1.0 if on[key] else 0.0)
        for reaction in model.reactions:
            flux = fluxes.get(reaction.id)
            if flux is None:
                continue
            cleaned = 0.0 if reaction.id in off_rxns else float(flux)
            reaction.forward_variable._internal_variable.setAttr(
                "Start", max(cleaned, 0.0)
            )
            reaction.reverse_variable._internal_variable.setAttr(
                "Start", max(-cleaned, 0.0)
            )
        model.solver.problem.update()
    except (AttributeError, TypeError, ValueError) as error:
        logger.warning(
            "the %s solver took no warm start (%s), so it must find its own "
            "first feasible solution",
            model.problem.__name__,
            error,
        )
        return False
    logger.warning(
        "warm start: %d of %d binaries begin switched on, from %d fluxes",
        sum(on.values()),
        len(on),
        len(fluxes),
    )
    return True


def _apply_trade_off_warm_start(
    model: "Model",
    fluxes: Optional[Mapping[str, float]],
    indicators: Mapping[str, "object"],
    groups: Mapping[str, Iterable[str]],
    signs: Mapping[str, float],
    forward_bins: Mapping[str, "object"],
    reverse_bins: Mapping[str, "object"],
    tol: float,
) -> bool:
    """Hand the solver a clean feasible starting point (`|x0| >= 0.99 * tol`) for tradeOff."""
    if not fluxes:
        logger.warning(
            "no warm start: the solver must find its own first feasible "
            "solution, which it may not manage at this size"
        )
        return False

    threshold = 0.99 * tol
    group_members = {key: list(members) for key, members in groups.items()}
    on = {
        key: all(abs(fluxes.get(rxn_id, 0.0)) >= threshold for rxn_id in members)
        for key, members in group_members.items()
    }
    off_rxns = {
        rxn_id
        for key, members in group_members.items()
        if not on[key]
        for rxn_id in members
    }
    try:
        for key, variable in indicators.items():
            variable._internal_variable.setAttr("Start", 1.0 if on[key] else 0.0)
        for rxn_id, a_var in forward_bins.items():
            ori = signs[rxn_id] * float(fluxes.get(rxn_id, 0.0))
            is_on = rxn_id not in off_rxns
            a_var._internal_variable.setAttr(
                "Start", 1.0 if (is_on and ori >= threshold) else 0.0
            )
        for rxn_id, b_var in reverse_bins.items():
            ori = signs[rxn_id] * float(fluxes.get(rxn_id, 0.0))
            is_on = rxn_id not in off_rxns
            b_var._internal_variable.setAttr(
                "Start", 1.0 if (is_on and ori <= -threshold) else 0.0
            )
        for reaction in model.reactions:
            flux = fluxes.get(reaction.id)
            if flux is None:
                continue
            cleaned = 0.0 if reaction.id in off_rxns else float(flux)
            reaction.forward_variable._internal_variable.setAttr(
                "Start", max(cleaned, 0.0)
            )
            reaction.reverse_variable._internal_variable.setAttr(
                "Start", max(-cleaned, 0.0)
            )
        model.solver.problem.update()
    except (AttributeError, TypeError, ValueError) as error:
        logger.warning(
            "the %s solver took no warm start (%s), so it must find its own "
            "first feasible solution",
            model.problem.__name__,
            error,
        )
        return False
    logger.warning(
        "warm start: %d of %d binaries begin switched on, from %d fluxes",
        sum(on.values()),
        len(on),
        len(fluxes),
    )
    return True


def _set_time_limit(model: "Model", time_limit: Optional[float]) -> None:
    """Ask the solver to stop after the given time, if it supports that."""
    if time_limit is None:
        return
    try:
        model.solver.configuration.timeout = int(math.ceil(time_limit))
    except (AttributeError, TypeError, ValueError) as error:
        logger.warning(
            "Could not set a time limit of %.4gs on the %s solver: %s",
            time_limit,
            model.problem.__name__,
            error,
        )
