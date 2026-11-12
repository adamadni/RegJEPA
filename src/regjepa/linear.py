"""Frozen-feature linear readout; this does not train an encoder or infer product ratios."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .evaluation import CandidateTable, choice_loss


def linear_scores(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Fixed per-row reductions keep identical inputs tied across batch positions."""
    features = np.ascontiguousarray(values, dtype=np.float64)
    coefficients = np.asarray(weights, dtype=np.float64)
    if features.ndim != 2 or coefficients.shape != (features.shape[1],):
        raise ValueError("Linear score dimensions do not match")
    # A BLAS GEMV may use different summation kernels for different row positions.
    # The resulting ~1e-16 differences can incorrectly split exact feature ties.
    # Do not merge outcome rows or introduce a label-dependent score tolerance.
    return np.einsum("ij,j->i", features, coefficients, optimize=False)


def centered_features(values: np.ndarray, table: CandidateTable) -> np.ndarray:
    data = np.asarray(values, dtype=np.float64)
    if data.ndim != 2 or data.shape[0] != table.n_rows or not np.isfinite(data).all():
        raise ValueError("Features must be finite and aligned to complete candidate sets")
    centered = data.copy()
    for indices in table.rows:
        centered[indices] -= centered[indices].mean(axis=0)
    return centered


@dataclass(frozen=True)
class LinearTransform:
    """Train-only within-molecule RMS scales; centering preserves linear choices exactly."""

    keep: np.ndarray
    scale: np.ndarray
    standardize: bool

    @classmethod
    def fit(cls, values: np.ndarray, table: CandidateTable, *, standardize: bool,
            minimum_rms: float = 1e-12) -> LinearTransform:
        centered = centered_features(values, table)
        second_moment = np.mean([np.mean(centered[indices] ** 2, axis=0) for indices in table.rows], axis=0)
        rms = np.sqrt(second_moment)
        keep = rms > minimum_rms
        scale = np.where(keep, rms, 1.0) if standardize else np.ones(centered.shape[1])
        return cls(keep, scale, standardize)

    def apply(self, values: np.ndarray, table: CandidateTable) -> np.ndarray:
        centered = centered_features(values, table)
        if centered.shape[1] != len(self.keep):
            raise ValueError("Feature width differs from the training transform")
        return np.ascontiguousarray(centered[:, self.keep] / self.scale[self.keep])


def linear_objective(weights: np.ndarray, values: np.ndarray, table: CandidateTable,
                     penalties: np.ndarray) -> tuple[float, np.ndarray]:
    if weights.shape != (values.shape[1],) or penalties.shape != weights.shape or np.any(penalties <= 0):
        raise ValueError("One positive L2 penalty is required per fitted coefficient")
    nll, derivative = choice_loss(linear_scores(values, weights), table)
    return nll + float(0.5 * np.dot(penalties, weights ** 2)), values.T @ derivative + penalties * weights


@dataclass(frozen=True)
class LinearScorer:
    transform: LinearTransform
    weights: np.ndarray
    diagnostics: dict

    def predict(self, values: np.ndarray, table: CandidateTable) -> np.ndarray:
        return linear_scores(self.transform.apply(values, table), self.weights)


def fit_linear(values: np.ndarray, table: CandidateTable, penalties: float | np.ndarray,
               *, standardize: bool, max_iterations: int = 1000,
               gradient_tolerance: float = 1e-7, accepted_gradient: float = 1e-5) -> LinearScorer:
    """Convex conditional-choice fit with analytic gradients; fail on nonconvergence."""
    from scipy.optimize import minimize

    transform = LinearTransform.fit(values, table, standardize=standardize)
    features = transform.apply(values, table)
    raw_penalties = np.broadcast_to(np.asarray(penalties, dtype=np.float64), (values.shape[1],))
    if not np.isfinite(raw_penalties).all() or np.any(raw_penalties <= 0):
        raise ValueError("Penalties must be finite and positive")
    selected_penalties = raw_penalties[transform.keep]
    if not features.shape[1]:
        nll, _ = choice_loss(np.zeros(table.n_rows), table)
        return LinearScorer(transform, np.zeros(0), {"converged": True, "iterations": 0,
                            "active_coefficients": 0, "objective": nll, "gradient_inf_norm": 0.0})
    result = minimize(linear_objective, np.zeros(features.shape[1]), args=(features, table, selected_penalties),
                      method="L-BFGS-B", jac=True,
                      options={"maxiter": max_iterations, "gtol": gradient_tolerance,
                               "ftol": 1e-14, "maxls": 50, "maxcor": 10})
    gradient = float(np.max(np.abs(result.jac)))
    converged = bool(result.success and np.isfinite(result.fun) and np.isfinite(result.x).all()
                     and gradient <= accepted_gradient)
    diagnostics = {"converged": converged, "iterations": int(result.nit),
                   "active_coefficients": features.shape[1], "objective": float(result.fun),
                   "gradient_inf_norm": gradient, "optimizer_message": str(result.message)}
    if not converged:
        raise RuntimeError(f"Linear scorer did not converge: {diagnostics}")
    return LinearScorer(transform, result.x.copy(), diagnostics)
