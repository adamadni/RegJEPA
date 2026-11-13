"""Linear and small nonlinear listwise probes; neither trains a molecular encoder."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from regjepa.evaluation import CandidateTable, choice_loss
from regjepa.linear import linear_objective, linear_scores


@dataclass
class Transform:
    head: str
    keep: np.ndarray
    mean: np.ndarray
    scale: np.ndarray

    @classmethod
    def fit(cls, values: np.ndarray, table: CandidateTable, head: str, morgan: np.ndarray) -> Transform:
        if values.shape[0] != table.n_rows or not np.isfinite(values).all():
            raise ValueError("Invalid feature matrix")
        if head not in {"linear", "mlp"}:
            raise ValueError("Unknown head")
        # Shared preprocessing permits a head comparison without discarding the
        # molecular context that a nonlinear head can use through interactions.
        mean = np.mean([values[rows].mean(axis=0) for rows in table.rows], axis=0)
        data = values - mean
        rms = np.sqrt(np.mean([np.mean(data[rows] ** 2, axis=0) for rows in table.rows], axis=0))
        keep = rms > 1e-12
        scale = np.where(keep, rms, 1.0)
        if head == "linear":
            scale[morgan] = 1.0
        return cls(head, keep, mean, scale)

    def apply(self, values: np.ndarray, table: CandidateTable) -> np.ndarray:
        if values.shape != (table.n_rows, len(self.keep)) or not np.isfinite(values).all():
            raise ValueError("Unaligned/nonfinite transform inputs")
        data = values - self.mean
        return np.ascontiguousarray(data[:, self.keep] / self.scale[self.keep])


def unpack(parameters: np.ndarray, width: int, hidden: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    end = width * hidden
    return parameters[:end].reshape(width, hidden), parameters[end:end + hidden], parameters[end + hidden:]


def nonlinear_scores(values: np.ndarray, parameters: np.ndarray, hidden: int) -> tuple[np.ndarray, np.ndarray]:
    weights, bias, output = unpack(parameters, values.shape[1], hidden)
    # Fixed reductions preserve identical rows across batch positions.
    activation = np.tanh(np.einsum("ij,jk->ik", values, weights, optimize=False) + bias)
    return linear_scores(activation, output), activation


def nonlinear_objective(parameters: np.ndarray, values: np.ndarray, table: CandidateTable,
                        penalty: float, hidden: int) -> tuple[float, np.ndarray]:
    scores, activation = nonlinear_scores(values, parameters, hidden)
    nll, derivative = choice_loss(scores, table)
    weights, _, output = unpack(parameters, values.shape[1], hidden)
    delta = derivative[:, None] * output[None, :] * (1.0 - activation ** 2)
    gradient = np.concatenate([(values.T @ delta).ravel(), delta.sum(axis=0), activation.T @ derivative])
    return nll + 0.5 * penalty * float(parameters @ parameters), gradient + penalty * parameters


@dataclass
class Scorer:
    transform: Transform
    parameters: np.ndarray
    hidden: int
    diagnostics: dict

    def predict(self, values: np.ndarray, table: CandidateTable) -> np.ndarray:
        data = self.transform.apply(values, table)
        if self.transform.head == "linear":
            return linear_scores(data, self.parameters)
        return nonlinear_scores(data, self.parameters, self.hidden)[0]


def fit(values: np.ndarray, table: CandidateTable, morgan: np.ndarray, head: str,
        penalties: np.ndarray, config: dict) -> Scorer:
    transform = Transform.fit(values, table, head, morgan)
    data = transform.apply(values, table)
    selected = penalties[transform.keep]
    if not np.isfinite(penalties).all() or np.any(penalties <= 0):
        raise ValueError("All penalties must be finite and positive")
    hidden = config["mlp_hidden"] if head == "mlp" else 0
    seeds = config["mlp_restart_seeds"] if head == "mlp" else [None]
    attempts, solutions = [], []
    for seed in seeds:
        if head == "linear":
            initial = np.zeros(data.shape[1])
            objective = lambda p: linear_objective(p, data, table, selected)
        else:
            if not np.all(penalties == penalties[0]):
                raise ValueError("MLP uses one shared parameter penalty, not blockwise penalties")
            rng = np.random.default_rng(seed)
            initial = np.concatenate([rng.normal(0, 1 / np.sqrt(max(data.shape[1], 1)), (data.shape[1], hidden)).ravel(),
                                      np.zeros(hidden), rng.normal(0, 1 / np.sqrt(hidden), hidden)])
            objective = lambda p: nonlinear_objective(p, data, table, float(penalties[0]), hidden)
        if not len(initial):
            value, gradient = objective(initial)
            parameters = initial
            diagnostic = {"converged": True, "iterations": 0, "objective": value,
                          "gradient_inf_norm": 0.0, "optimizer_message": "no usable linear features", "seed": seed}
        else:
            result = minimize(objective, initial, jac=True, method="L-BFGS-B",
                              options={"maxiter": config["max_iterations"], "gtol": config["gradient_tolerance"],
                                       "ftol": 1e-14, "maxls": 50, "maxcor": 10,
                                       "maxfun": config["max_iterations"] * 50})
            gradient = float(np.max(np.abs(result.jac)))
            converged = bool(result.success and np.isfinite(result.fun) and np.isfinite(result.x).all()
                             and gradient <= config["accepted_gradient"])
            diagnostic = {"converged": converged, "iterations": int(result.nit), "objective": float(result.fun),
                          "gradient_inf_norm": gradient, "optimizer_message": str(result.message), "seed": seed}
            if not converged:
                raise RuntimeError(f"Follow-up scorer failed convergence: {diagnostic}")
            parameters = result.x.copy()
        attempts.append(diagnostic)
        solutions.append(parameters)
    best = min(range(len(attempts)), key=lambda i: attempts[i]["objective"])
    return Scorer(transform, solutions[best], hidden,
                  {"selected_restart": best, "attempts": attempts, "selection": "minimum regularized training objective",
                   "active_features": int(transform.keep.sum()), "parameter_count": len(solutions[best])})
