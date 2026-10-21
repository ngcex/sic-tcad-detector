"""Uncertainty quantification and global sensitivity analysis over the ETNA chain.

- ``etna.uq.params``: the registry of uncertain inputs (distribution, bounds, provenance,
  label) and ``ParameterSpace`` (correlations through a Gaussian copula; correlated
  blocks are single Sobol' factors).
- ``etna.uq.sampling``: Sobol' / Latin hypercube / MC designs, correlated normal draws,
  the extendable Saltelli design.
- ``etna.uq.sensitivity``: Sobol' first / total indices with bootstrap CIs; Ishigami.
- ``etna.uq.observables``: batched evaluators on the numba engines (no devsim).
- ``etna.uq.run``: ``run_uq`` (bands + indices, adaptive to a CI target).
- ``etna.uq.plots``: ``uq_figure`` in the etna.viz style (import it explicitly).

See docs/USAGE.md §20.
"""

from etna.uq.observables import OBSERVABLES
from etna.uq.params import (
    DEFAULT_PARAMS,
    ParameterSpace,
    UncertainParameter,
    default_parameters,
    tied_removal,
)
from etna.uq.run import indices_table, run_uq
from etna.uq.sampling import SaltelliDesign, correlated_normal, sample
from etna.uq.sensitivity import ishigami, ishigami_indices, sobol_indices

__all__ = [
    "run_uq",
    "indices_table",
    "OBSERVABLES",
    "DEFAULT_PARAMS",
    "ParameterSpace",
    "UncertainParameter",
    "default_parameters",
    "tied_removal",
    "SaltelliDesign",
    "sample",
    "correlated_normal",
    "sobol_indices",
    "ishigami",
    "ishigami_indices",
]
