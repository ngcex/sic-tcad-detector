"""Automatic, adaptive device-design optimisation on the fast engines (issue #10).

* ``run_design`` (also ``etna.run_design``): the facade, a SimResult with the Pareto
  front (or the best design), alternative-physics re-evaluation and breakdown check.
* ``etna.design.space``: variables, lattices, unit-cube encoding.
* ``etna.design.problem``: objectives, constraints, the vectorised Problem.
* ``etna.design.bo``: GP Bayesian optimisation (constrained EI, ParEGO, kriging-believer
  batches), cache, checkpoint / resume.
* ``etna.design.gp``: Matern-5/2 ARD Gaussian process (numpy + scipy).
* ``etna.design.pareto``: non-dominated sort, exact 2D / 3D hypervolume.
* ``etna.design.device``: batched engine evaluation of DeviceConfigs (no devsim).
* ``etna.design.dosimeter``: bounds from the engine box and tolerances, robust mode.
* ``etna.design.testproblems``: analytic problems with known optima.
"""

from etna.design.api import ALT_PHYSICS, METRIC_SENSE, clear_cache, run_design
from etna.design.bo import DesignHistory, EvaluationCache, optimize
from etna.design.device import (
    ALL_METRICS,
    CHEAP_METRICS,
    DesignSettings,
    evaluate_configs,
)
from etna.design.dosimeter import (
    DEFAULT_VARIABLES,
    TOLERANCES,
    VARIABLES,
    RobustSpec,
    area_bounds,
    design_space,
    make_problem,
)
from etna.design.problem import Constraint, Objective, Problem
from etna.design.space import DesignSpace, Variable

__all__ = [
    "run_design", "clear_cache", "ALT_PHYSICS", "METRIC_SENSE", "optimize",
    "DesignHistory", "EvaluationCache", "DesignSettings", "evaluate_configs",
    "ALL_METRICS", "CHEAP_METRICS", "DEFAULT_VARIABLES", "TOLERANCES", "VARIABLES",
    "RobustSpec", "area_bounds", "design_space", "make_problem", "Objective",
    "Constraint", "Problem", "DesignSpace", "Variable",
]  # fmt: skip
